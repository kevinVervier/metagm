# metagm: the metagenomic platform from Team162

In this page, we provide examples illustrating the different options offered by the platform.

## metagm_build module

This section describes the process of building databases for different softwares (Kraken, Bracken, Mash, ...) using a list of reference genomes.

![metagm_build_pipeline](img/metagm_build_pipeline.PNG)

### Inputs (mandatory)

### Options

### Outputs


## metagm_classify module

This section describes the process of classifying metagenomics sequencing reads to get both taxonomic and functional profiles.

![metagm_classify_pipeline](img/metagm_classify_pipeline.PNG)

### Inputs (mandatory)

### Options

### Outputs

## Taxonomy

In this section, we describe the process used to assign a given genome to the current taxonomy.
As mentioned in the section 'metagm_build module', users have the option to rely on either the [gtdb](https://gtdb.ecogenomic.org/) taxonomy done with `gtdb-tk classify_wf` [function](https://github.com/Ecogenomics/GtdbTk).

## statistical analysis

This section provides some R snippets to do analysis using the metagm_classify output files. 

More details are given in the R vignette '.Rmd', also found in this repository.
