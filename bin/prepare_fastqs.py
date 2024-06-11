#!/usr/bin/env python3

import argparse
import itertools
import logging
import os
import pathlib
import re
import shutil
import subprocess
import sys

from collections import Counter


logging.basicConfig(
	level=logging.DEBUG,
	format='[%(asctime)s] %(message)s'
)

logger = logging.getLogger(__name__)


def check_pairwise(r1, r2):
	""" Checks if two sets of read files contain the same prefixes.

	Input:
	 - r1/r2: lists of tuples (prefix, filename) of R1/R2 paired-end read files

	Raises error if one prefix is not found in either r1 or r2.

	"""
	r1_r2 = tuple(prefix[:-1] for prefix, _ in itertools.chain(r1, r2))
	for prefix in set(r1_r2):
		if r1_r2.count(prefix) != 2:
			raise ValueError(f"Missing mates for prefix {prefix}.")


def transfer_file(source, dest, remote_input=False):
	""" Transfers a file depending on its location and compressed state.

	Input:
	 - path to source file
	 - path to destination file
	 - whether source file is considered to be located on a remote file system

	"""
	resolved_src = pathlib.Path(source).resolve()

	if os.path.splitext(source)[1][1:] in ("gz", "bz2"):
		if remote_input:
			# if file is on remote file system, copy it to destination
			logging.debug('transfer_file: source=%s, dest=%s, remote_input=%s, action=copy', source, dest, remote_input)
			shutil.copyfile(resolved_src, dest)
		else:
			# if file is compressed and on local fs, just symlink it
			logging.debug('transfer_file: source=%s, dest=%s, remote_input=%s, action=symlink', source, dest, remote_input)
			pathlib.Path(dest).symlink_to(resolved_src)
	else:
		# if file is not compressed, gzip it to destination
		logging.debug('transfer_file: source=%s, dest=%s, remote_input=%s, action=gzip', source, dest, remote_input)
		with open(dest, "wt") as _out:
			subprocess.run(("gzip", "-c", resolved_src), stdout=_out)


def transfer_multifiles(files, dest, remote_input=False, compression=None):
	""" Transfers a set of files depending on their location and compressed state.

	Input:
	 - list of source file paths
	 - path to destination file
	 - whether source files are considered to be located on a remote file system
	 - the compression type of the files (supported: None, gz, bz2)
	"""
	
	if len(files) > 1:
		src_files = tuple(os.path.abspath(f) for f in files)  # tuple(f.resolve() for f in files)
		cat_cmd = ("cat", ) + src_files

		if compression in ("gz", "bz2"):
			# multiple compressed files can just be concatenated
			logging.debug('transfer_multifiles: compression=%s, remote_input=%s, action=concatenate', compression, remote_input)
			with open(dest, "wt") as _out:
				subprocess.run(cat_cmd, stdout=_out)
		else:
			# multiple uncompressed files will be cat | gzipped
			logging.debug('transfer_multifiles: compression=%s, remote_input=%s, action=concatenate+gzip', compression, remote_input)
			cat_pr = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
			with open(dest, "wt") as _out:
				subprocess.run(("gzip", "-c", "-"), stdin=cat_pr.stdout, stdout=_out)
			
	else:
		logging.debug('transfer_multifiles: single file, source=%s, dest=%s, remote_input=%s, action=defer->transfer_file', files[0], dest, remote_input)
		transfer_file(files[0], dest, remote_input=remote_input)


def process_sample(
	sample, fastqs, output_dir,
	fastq_suffix_pattern,
	remove_suffix=None, remote_input=False,
	add_suffix=None,
):
	""" Checks if a set of fastq files in a directory is a valid collection
	and transfers files to a destination dir upon success.

	Input:
	 - sample_id
	 - list of fastq files
	 - path to output directory
	 - suffix to strip off from filenames (e.g. _001)
	 - whether fastq files are located on remote file system
	"""

	if len(fastqs) == 1:
		# remove potential "single(s)" string from single fastq file name prefix
		sample_sub = re.sub(r"[._]singles?", "", sample)
		# 20221018: and attach it at the end of the sample name
		# - this might be a temporary fix, but @93a73d0
		# single-end samples without .singles-suffix cause problems 
		# with fastqc results in the collate step
		if add_suffix:
			sample_sub += f".{add_suffix}"
		sample = sample_sub + ".singles"

		sample_dir = os.path.join(output_dir, sample)
		pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)

		dest_compression = fastqs[0][fastqs[0].rfind(".") + 1:]

		dest = os.path.join(sample_dir, f"{sample}_R1.fastq.{dest_compression}")
		transfer_file(fastqs[0], dest, remote_input=remote_input)

		yield sample, False

	elif fastqs:

		# check if all fastq files are compressed the same way
		# suffixes-counter will have True/False counts for compressed/uncompressed
		# if True only - all files are compressed
		# if False only - all files are uncompressed
		suffixes = Counter(
			f[f.rfind("."):] in (".gz", ".bz2") for f in fastqs
		)

		# mix -> error out
		if len(suffixes) > 1:
			raise ValueError(f"sample: {sample} has mixed compressed and uncompressed input files. Please check.")

		all_compressed, _ = suffixes.most_common()[0]
		
		if all_compressed:
			suffixes = Counter(
				f[f.rfind(".") + 1:] for f in fastqs
			)
			if len(suffixes) > 1:
				raise ValueError(f"sample: {sample} has mixed gzip and bzip2 files. Please check.")
			dest_compression = compression = suffixes.most_common()[0][0]
		else:
			dest_compression, compression = "gz", None

		# extract the file name prefixes
		prefixes = [re.sub(fastq_suffix_pattern, "", os.path.basename(f)) for f in fastqs]
		if remove_suffix:
			# remove suffix pattern if requested
			prefixes = [re.sub(remove_suffix + r"$", "", p) for p in prefixes]

		print("PRE", prefixes, file=sys.stderr)

		# partition fastqs into R1, R2, and 'other' sets
		r1 = [(p, f) for p, f in zip(prefixes, fastqs) if re.search(r"[._R]1$", p)]
		r2 = [(p, f) for p, f in zip(prefixes, fastqs) if re.search(r"[._R]2$", p)]
		others = sorted(list(set(fastqs).difference({f for _, f in r1}).difference({f for _, f in r2})))

		# check if R1/R2 sets have equal sizes or are empty
		# R1 empty: potential scRNAseq (or any protocol with barcode reads in R1)
		# R2 empty: typical single end reads with (R?)1 suffix
		assert len(r2) == 0 or len(r1) == 0 or (r1 and len(r1) == len(r2)), "R1/R2 sets are not of the same length"

		# if R1 and R2 are of equal size, check if the prefixes match
		if len(r1) == len(r2) and r1:
			check_pairwise(r1, r2)

		# sort R1/R2 for concatenation, get rid off prefixes
		r1 = sorted(f for _, f in r1)
		r2 = sorted(f for _, f in r2)

		print("R1", r1, file=sys.stderr)
		print("R2", r2, file=sys.stderr)
		print("others", others, file=sys.stderr, flush=True)

		if add_suffix:
			sample += f".{add_suffix}"

		sample_dir = os.path.join(output_dir, sample)

		if r1 or r2:

			pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)

			if r1:
				# if R1 is not empty, transfer R1-files
				dest = os.path.join(sample_dir, f"{sample}_R1.fastq.{dest_compression}")
				transfer_multifiles(r1, dest, remote_input=remote_input, compression=compression)
			if r2:
				# if R2 is not empty, transfer R2-files,
				# if R1 is empty, rename R2 to R1 so that files can be processed as normal single-end
				target_r = "R2" if r1 else "R1"
				dest = os.path.join(sample_dir, f"{sample}_{target_r}.fastq.{dest_compression}")
				transfer_multifiles(r2, dest, remote_input=remote_input, compression=compression)
		
			yield sample, bool(r1 and r2)

		if others:
			# if single-end reads exist,
			# transfer them to <sample>.singles
			# these will be processed independently and merged with the paired-end reads
			# at a later stage
			sample_dir = sample_dir + ".singles"
			sample = sample + ".singles"
			pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
			dest = os.path.join(sample_dir, f"{sample}_R1.fastq.{dest_compression}")
			transfer_multifiles(others, dest, remote_input=remote_input, compression=compression)

			yield sample, bool(r1 or r2)
		

def is_fastq(f, valid_fastq_suffixes, valid_compression_suffixes):
	""" Checks if a file is a fastq file (compressed or uncompressed.)

	Input:
	 - filename

	Output:
	 - true if file is fastq else false

	"""
	file_extensions = re.split(r"[._]", os.path.basename(f))[-2:]
	# try:
	# 	fastq_suffix, *compression_suffix = filename_tokens[-2:]
	# except ValueError:
	# 	return False
	valid_compression = len(file_extensions) == 2 and file_extensions[-1] in valid_compression_suffixes
	valid_fastq = any(
        (
            valid_compression and file_extensions[0] in valid_fastq_suffixes,
            file_extensions and file_extensions[-1] in valid_fastq_suffixes,
        )
    )

	# valid_compression = not compression_suffix or compression_suffix[0] in valid_compression_suffixes

	logger.info('OBJECT: %s FASTQ: %s COMPRESSION: %s ISFILE: %s' % (f, valid_fastq, valid_compression, os.path.isfile(f)))

	return os.path.isfile(f) and valid_fastq

	# if not compression_suffix:
	# 	return fq_suffix in valid_fastq_suffixes
	# else:
	# 	return compression_suffix[0] in valid_compression_suffixes and fq_suffix in valid_fastq_suffixes


	# prefix, suffix = os.path.splitext(f)
	# if suffix in valid_fastq_suffixes:
	# 	return True
	# if suffix in valid_compression_suffixes:
	# 	_, suffix = os.path.splitext(prefix)
	# 	return suffix in valid_fastq_suffixes
	# return False


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input_dir", type=str, default=".")
	ap.add_argument("-o", "--output_dir", type=str, default="prepared_samples")
	ap.add_argument("-p", "--prefix", type=str, required=True)
	ap.add_argument("--remote-input", action="store_true")
	ap.add_argument("--remove-suffix", type=str, default=None)
	ap.add_argument("--valid-fastq-suffixes", type=str, default="fastq,fq")
	ap.add_argument("--valid-compression-suffixes", type=str, default="gz,bz2")
	ap.add_argument("--add_sample_suffix", type=str)

	args = ap.parse_args()

	valid_fastq_suffixes = tuple(f"{suffix}" for suffix in args.valid_fastq_suffixes.split(","))
	print(valid_fastq_suffixes)
	valid_compression_suffixes = tuple(f"{suffix}" for suffix in args.valid_compression_suffixes.split(","))
	print(valid_compression_suffixes)

	fastq_file_suffix_pattern = r"[._](" + \
		args.valid_fastq_suffixes.replace(",", "|") + \
		")([._](" + \
		args.valid_compression_suffixes.replace(",", "|") + \
		"))?$"

	fastq_mate_pattern = r"([._R][12])$"

	def collect_fastq_files(input_dir, valid_fastq_suffixes, valid_compression_suffixes):
		return sorted(
				os.path.join(input_dir, f)
				for f in os.listdir(input_dir)
				if is_fastq(os.path.join(input_dir, f), valid_fastq_suffixes, valid_compression_suffixes)
			)

	try:
		pwd, dirs, _ = next(os.walk(args.input_dir))
	except StopIteration:
		raise ValueError(f"Could not find input directory {args.input_dir} ({os.path.abspath(args.input_dir)})")

	samples = {}
	pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

	for sample_dir in dirs:
		sample, sample_dir = sample_dir, os.path.join(pwd, sample_dir)

		samples.setdefault(sample, []).extend(
			collect_fastq_files(sample_dir, valid_fastq_suffixes, valid_compression_suffixes)			
		)

	root_fastqs = collect_fastq_files(args.input_dir, valid_fastq_suffixes, valid_compression_suffixes)

	if samples and root_fastqs:
		raise ValueError(f"Found {len(root_fastqs)} fastq files in input directory together with {len(samples)} sample directories. Please check input data.")
	elif root_fastqs:
		for f in root_fastqs:
			sample = re.sub(fastq_file_suffix_pattern, "", os.path.basename(f))
			sample = re.sub(fastq_mate_pattern, "", sample)
			samples.setdefault(sample, []).append(f)

	# check and transfer the files
	with open("sample_library_info.txt", "wt") as lib_out:
		for sample, fastqs in samples.items():
			try:
				renamed = process_sample(
					sample, fastqs, args.output_dir,
					fastq_file_suffix_pattern,
					remove_suffix=args.remove_suffix, remote_input=args.remote_input,
					add_suffix=args.add_sample_suffix,
				)
			except Exception as e:
				raise ValueError(f"Encountered problems processing sample '{sample}': {e}.\nPlease check your file names.")
			else:
				for sample, is_paired in renamed:
					print(sample, int(is_paired), sep="\t", file=lib_out)

if __name__ == "__main__":
	main()
