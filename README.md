# metagm: the metagenomic platform from Team162

In this page, we provide examples illustrating the different options offered by the platform.

## metagm_build module

This section describes the process of building databases for different softwares (Kraken, Bracken, Mash, ...) using a list of reference genomes.

![metagm_build_pipeline](img/metagm_build_pipeline.PNG)

### Inputs (mandatory)
There is two mandatory positional arguments for this function:
* a `genomes` list
* an `output` folder

### Options
The `output` folder contains a directory for each task that is performed:
* `output/Kraken` contains the [Kraken2 database](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#kraken-2-databases), as well as the [Bracken files](https://github.com/jenniferlu717/Bracken#step-1-generate-the-bracken-database-file-databasexmerskmer_distrib) (if requested)
* `output/Mash` contains the [Mash sketch](https://mash.readthedocs.io/en/latest/sketches.html) file

### Outputs

### Comments

* Building a Kraken database requires a large memory amount (>100Gb). The function will automatically submit the final `kraken_build` job on `teramem` queue which means users need to run it on farm3 (or farm4/5 when available).

## metagm_classify module

This section describes the process of classifying metagenomics sequencing reads to get both taxonomic and functional profiles.

![metagm_classify_pipeline](img/metagm_classify_pipeline.PNG)

### Inputs (mandatory)

### Options

### Outputs

### Comments

## Classes

# Taxonomy

In this section, we describe the process used to assign a given genome to the current taxonomy.
As mentioned in the section 'metagm_build module', users have the option to rely on either the [gtdb](https://gtdb.ecogenomic.org/) taxonomy done with `gtdb-tk classify_wf` [function](https://github.com/Ecogenomics/GtdbTk).

### How to build a tree in NCBI format using gtdb metadata

Kraken software relies on taxonomic information presented in the 'NCBI format' (a pair of `nodes.dmp` and `names.dmp`). 
Unfortunately, gtdb does not provide (yet?) its taxonomy in such format. Therefore, we downloaded all the gtdb [archeal metadata](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar122_metadata.tsv) and [bacterial metadata](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv) containing all taxonomic paths for every genome found in the database.

# Statistical analysis

This section provides some R snippets to do analysis using the metagm_classify output files. 

More details are given in the R vignette '.Rmd', also found in this repository.

# Miscalleneous

## Other utility functions/scripts

# TODO
