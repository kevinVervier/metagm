# metagm: the metagenomic platform from Team162

In this page, we provide examples illustrating the different options offered by the platform.

## Before starting
To use all the scripts described in this page, users need to add the following to their `~/.profile` file:

```
# add the metagm library to your path
export PATH=/nfs/team162/kv4/github/metagm/metagm/wrapper:$PATH
#add Bracken local install (needed for now as the Pathinf team has not installed it yet)
export PATH=/nfs/team162/kv4/bin/Bracken/:$PATH
```

## metagm_build module

This section describes the process of building databases for different softwares (Kraken, Bracken, Mash, ...) using a list of reference genomes.

![metagm_build_pipeline](img/metagm_build_pipeline.PNG)

### Inputs (mandatory)
There is two mandatory positional arguments for this function:
* a `genomes` list (text file):
  * __mandatory__ first column contains the absolute paths to genome assemblies (`.fa` or `.fna`)
  * second column contains the genome names (if not provided, the file name will be used)
  * third column contains the taxids (if not provided, a [taxonomic assignment](https://github.com/kevinVervier/metagm/blob/master/README.md#taxonomic-assignment) step is performed)
* an `output` folder to store all the files produced by the script

### Options
The script `metagm_build.py` also offers options to add flexibility in the building database process:
* 

### Outputs

The `output` folder contains a directory for each task that is performed:
* `output/merge_final` contains quality control (QC) results for all the genomes
  * `output/merge_final/ValidatedGenomes.txt` is the list of all genomes that pass QC
  * `output/merge_final/FilteredGenomes.txt` is the list of all genomes that fail QC
  * `output/merge_final/log.txt` provides details on why a genome failed QC
* `output/Kraken` contains the [Kraken2 database](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#kraken-2-databases), as well as the [Bracken files](https://github.com/jenniferlu717/Bracken#step-1-generate-the-bracken-database-file-databasexmerskmer_distrib) (if requested)
* `output/Mash` contains the [Mash sketch](https://mash.readthedocs.io/en/latest/sketches.html) file


### Comments

* You can access all the information in the command help `metagm_build.py -h `
* Building a Kraken database requires a large memory amount (>100Gb). The function will automatically submit the final `kraken_build` job on `teramem` queue which means users need to run it on farm3 (or farm4/5 when available).
* If neither `--QC`, `--Kraken`, `--Bracken`, or `--Mash` is provided, the `metagm_build.py` is not going to do anything.

## metagm_classify module

This section describes the process of classifying metagenomics sequencing reads to get both taxonomic and functional profiles.

![metagm_classify_pipeline](img/metagm_classify_pipeline.PNG)

### Inputs (mandatory)

### Options

### Outputs

### Comments

## Classes

# Taxonomy

## Taxonomic assignment
In this section, we describe the process used to assign a given genome to the current taxonomy.
As mentioned in the section 'metagm_build module', users have the option to rely on either the [gtdb](https://gtdb.ecogenomic.org/) taxonomy done with `gtdb-tk classify_wf` [function](https://github.com/Ecogenomics/GtdbTk).

## How to build a tree in NCBI format using gtdb metadata

Kraken software relies on taxonomic information presented in the 'NCBI format' (a pair of `nodes.dmp` and `names.dmp`). 
Unfortunately, gtdb does not provide (yet?) its taxonomy in such format. Therefore, we downloaded all the gtdb [archeal metadata](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar122_metadata.tsv) and [bacterial metadata](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv) containing all taxonomic paths for every genome found in the database.

# Statistical analysis

This section provides some R snippets to do analysis using the metagm_classify output files. 

More details are given in the R vignette '.Rmd', also found in this repository.

# Miscalleneous

## Other utility functions/scripts

# TODO
