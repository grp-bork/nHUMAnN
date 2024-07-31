# nHUMAnN workflow
<table>
  <tr width="100%">
    <td width="150px">
      <a href="https://www.bork.embl.de/"><img src="https://www.bork.embl.de/assets/img/normal_version.png" alt="Bork Group Logo" width="150px" height="auto"></a>
    </td>
    <td width="425px" align="center">
      <b>Developed by the <a href="https://www.bork.embl.de/">Bork Group</a></b><br>
      Raise an <a href="https://github.com/grp-bork/nHUMAnN/issues">issue</a> or <a href="mailto:N4M@embl.de">contact us</a><br><br>
      See our <a href="https://www.bork.embl.de/services.html">other Software & Services</a>
    </td>
    <td width="500px">
      Contributors:<br>
      <ul>
        <li>
          <a href="https://github.com/cschu/">Christian Schudoma</a> <a href="https://orcid.org/0000-0003-1157-1354"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
        <li>
          <a href="https://github.com/danielpodlesny/">Daniel Podlesny</a> <a href="https://orcid.org/0000-0002-5685-0915"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td colspan="3" align="center">The development of this workflow was supported by <a href="https://www.nfdi4microbiota.de/">NFDI4Microbiota <img src="https://github.com/user-attachments/assets/1e78f65e-9828-46c0-834c-0ed12ca9d5ed" alt="NFDI4Microbiota icon" width="20px" height="20px"></a> 
</td>
  </tr>
</table>

---
#### Description
The `nHUMAnN workflow` is a nextflow workflow for running `HUMAnN3` based on `Metaphlan4` profiles via joint index generation. The workflow includes optional read preprocessing and host/human decontamination steps provided by the [nevermore](https://github.com/cschu/nevermore) workflow library.

Due to compatibility issues between current `CHOCOPhlAn` databases and recent versions of `HUMAnN3`, `nHUMAnN` makes use of a patched `HUMAnN3` version obtainable as a [Docker container](registry.git.embl.de/schudoma/humann3-docker:latest).

#### Citation
This workflow: badge tbd.

Also cite:
```
Beghini F, McIver LJ, Blanco-Míguez A, et al. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. Elife. 2021;10:e65088. Published 2021 May 4. doi:10.7554/eLife.65088
```
---
# Overview

![nevermore_workflow](https://raw.githubusercontent.com/grp-bork/nHUMAnN/main/docs/nevermore.svg)
![HUMAnN3_subworkflow](https://raw.githubusercontent.com/grp-bork/nHUMAnN/main/docs/humann3.svg)

---
# Requirements
The easiest way to handle dependencies is via Singularity/Docker containers. Alternatively, conda environments, software module systems or native installations can be used.

## Preprocessing

Preprocessing and QA is done with `bbmap`, `fastqc`, and `multiqc`.

### Decontamination/Host removal

Decontamination is done with `kraken2` and additionally requires `seqtk`. 

#### Kraken2 database

Host removal requires a `kraken2` host database.

### Metaphlan Profiling

The default supported `MetaPhlAn` version is 4.

#### CHOCOPhlAn database for Metaphlan4

Get the `mpa_vOct22_CHOCOPhlAnSGB_202212` database from [here](http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.tar), unpack the tarball, and point the `--mp4_db` parameter to the database's root directory. 

In `params.yml`:

```
mp4_db: "/path/to/mpa_vOct22_CHOCOPhlAnSGB_202212/"
```

On the command line:

```
--mp4_db "/path/to/mpa_vOct22_CHOCOPhlAnSGB_202212/"
```

#### HUMAnN3 databases

Get the annotated `CHOCOPhlAn` db from [here](http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz) and the annotated uniref db from [here](http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref50_annotated_v201901b_full.tar.gz), unpack the tarballs and set the respective parameters.

  In `params.yml`:

```
humann_nuc_db: "/path/to/full_chocophlan_db/"
humann_prot_db: "/path/to/uniref90_annotated_v201901b_full/"
```

On the command line:

```
--humann_nuc_db "/path/to/full_chocophlan_db/"
--humann_prot_db "/path/to/uniref90_annotated_v201901b_full/"
```


---
# Usage
## Cloud-based Workflow Manager (CloWM)
This workflow will be available on the `CloWM` platform (coming soon).

## Command-Line Interface (CLI)
The workflow run is controlled by environment-specific parameters (see [run.config](https://raw.githubusercontent.com/grp-bork/nHUMAnN/main/config/run.config)) and study-specific parameters (see [params.yml](https://raw.githubusercontent.com/grp-bork/nHUMAnN/main/config/params.yml)). The parameters in the `params.yml` can be specified on the command line as well.

You can either clone this repository from GitHub and run it as follows
```
git clone https://github.com/grp-bork/nHUMAnN.git
nextflow run /path/to/nhumann [-resume] -c /path/to/run.config -params-file /path/to/params.yml
```

Or, you can have nextflow pull it from github and run it from the `$HOME/.nextflow` directory.
```
nextflow run cschu/nHUMAnN [-resume] -c /path/to/run.config -params-file /path/to/params.yml
```

## Input files
Fastq files are supported and can be either uncompressed (but shouldn't be!) or compressed with `gzip` or `bzip2`. Sample data must be arranged in one directory per sample.

### Per-sample input directories
All files in a sample directory will be associated with the name of the sample folder. Paired-end mate files need to have matching prefixes. Mates 1 and 2 can be specified with suffixes `_[12]`, `_R[12]`, `.[12]`, `.R[12]`. Lane IDs or other read id modifiers have to precede the mate identifier. Files with names not containing either of those patterns will be assigned to be single-ended. Samples consisting of both single and paired end files are assumed to be paired end with all single end files being orphans (quality control survivors). 
