#!/usr/bin/env python3

import argparse
import gzip
import os
import pathlib
import shutil
import subprocess
import sys


def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("--input_dir", "-i", type=str, default=".")
	ap.add_argument("--output_dir", "-o", type=str, default="fastq")
	args = ap.parse_args()

	pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

	for f in sorted(os.listdir(args.input_dir)):
		full_f = os.path.join(args.input_dir, f)
		if pathlib.Path(full_f).is_symlink() and f != os.path.basename(__file__):
			link_target = pathlib.Path(full_f).resolve()
			sample = os.path.basename(os.path.dirname(link_target))
			print(f, sample)
			if not sample:
				raise NotImplementedError("Flat-directories not implemented.")
			sample_dir = os.path.join(args.output_dir, sample)
			print(sample_dir)
			pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
			dest = os.path.join(sample_dir, f)
			if f.endswith(".gz"):
				shutil.copyfile(link_target, dest)
			else:
				with open(dest + ".gz", "wt") as _out:
					subprocess.run(("gzip", "-c", full_f), stdout=_out)






	...


if __name__ == "__main__":
	main()