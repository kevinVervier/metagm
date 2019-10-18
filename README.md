# metagm: the metagenomic platform from Team162

In this page, we provide examples illustrating the different options offered by the platform.

## Before starting
__To use all the scripts described in this page, users need to add__ the following to their `~/.profile` file:

```
# add the metagm library to your path
export PATH=/nfs/team162/kv4/github/metagm/metagm/wrapper:$PATH
#add Bracken local install (needed for now as the Pathinf team has not installed it yet)
export PATH=/nfs/team162/kv4/bin/Bracken/:$PATH
```
and then run the command `source ~/.profile`

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
The script `metagm_build.py` also offers options in the building database process:
* `-h`: help display
* `-v`: verbose mode for additional details on each step
* `--QC`: run [quality control](https://github.com/kevinVervier/metagm/blob/master/README.md#quality-control) on the list of genomes before building any database. If not done, the scirpt assumes that all the genomes have already been cheked.
* `--taxoAssign`: run [taxonomic assignment](https://github.com/kevinVervier/metagm/blob/master/README.md#taxonomic-assignment) step on all the genomes. It will automatically be done if Kraken/Bracken databases are built.
* `--KrakenDB`: build a Kraken2 database using the list of genomes
* `--BrackenDB`: build a Bracken database using the list of genomes. Requires a Kraken database to exist, and will therefore automatically creates one.
* `--MashDB`: create a [Mash sketch](https://mash.readthedocs.io/en/latest/sketches.html) of all the genomes.
* `--ncbi`: rely on NCBI taxonomy instead of gtdb (_default: false_)
* `-b`: define the number of genomes to be analyzed in each batch (_default: 10_)
* `-t`: define the number of threads used in each job (_default: 2_)
* `-q`: define to which queue the jobs are submitted (_default: long_)
* `-m`: define the amount of memory requested for each job (_default: 64_)
* `--maxcontamination`: define the maximal value on checkm contamination to filter out a genome during QC (_default: 5_)
* `--mincompleteness`: define the minimal value on checkm completeness to filter out a genome during QC (_default: 90_)

### Outputs

The `output` folder contains a directory for each task that is performed:
* `output/merge_final` contains quality control (QC) results for all the genomes
  * `output/merge_final/ValidatedGenomes.txt` is the list of all genomes that pass QC
  * `output/merge_final/FilteredGenomes.txt` is the list of all genomes that fail QC
  * `output/merge_final/log.txt` provides details on why a genome failed QC
* `output/tmp$i` folders contain the [checkm lineage_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow) output for the batch `i` of genomes (if parallelized)
* `output/Kraken` contains the [Kraken2 database](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#kraken-2-databases), as well as the [Bracken files](https://github.com/jenniferlu717/Bracken#step-1-generate-the-bracken-database-file-databasexmerskmer_distrib) (if requested)
* `output/Mash` contains the [Mash sketch](https://mash.readthedocs.io/en/latest/sketches.html) file
* `output/genome_with_[gtdb|ncbi]_taxid.txt.txt` contains the list of taxids assigned to the genomes list if a [taxonomic assignment](https://github.com/kevinVervier/metagm/blob/master/README.md#taxonomic-assignment) step is performed

### Examples

The following examples illustrate various featuresfrom the `metagm_build.py` script. Depending how busy _farm_ is, these examples can take some time to run.

#### Quality control + taxonomic assignment on a list of genomes

```
metagm_build.py /nfs/team162/kv4/bin/list_example_pipeline.txt ./ --QC --taxoAssign
```

The command applies:
1. [Quality control](https://github.com/kevinVervier/metagm/blob/master/README.md#quality-control) on 10 genomes.
 * according to `./merge_final/ValidatedGenomes.txt`, there is XYZ genomes that passed QC
 * according to `./merge_final/FilteredGenomes.txt`, there is XYZ genomes that failed QC
 * according to `./merge_final/log.txt`, the genomes were filtered because of XYZ
2. [taxonomic assignment](https://github.com/kevinVervier/metagm/blob/master/README.md#taxonomic-assignment) on validated genomes only, using GTDB taxonomy (default).
 * the taxonomic assignment can be found in `./genome_with_gtdb_taxid.txt`

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

Current taxonomic tree build from GTDB metadata is stored here: `/nfs/pathogen005/team162/gtdb_taxonomy` 

## How to build a tree in NCBI format using gtdb metadata

Kraken software relies on taxonomic information presented in the 'NCBI format' (a pair of `nodes.dmp` and `names.dmp`). 
Unfortunately, gtdb does not provide (yet?) its taxonomy in such format. Therefore, we need to build a GTDB-based taxonomy that can be saved in the NCBI format.
1. download all the gtdb [archeal metadata](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar122_metadata.tsv) and [bacterial metadata](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv) containing all taxonomic paths for every genome found in the database.
2. extract unique paths in the metadata files to make the following steps faster
```
# merge bacterial and archeal metadata
cat ar122_metadata.tsv bac120_metadata.tsv > bac_ar_metadata.tsv
# get unique taxonomic paths
cut -f17 bac_ar_metadata.tsv | sort -u > tmp.txt
# remove last line that contains header
head -n -1 tmp.txt > taxo_path_gtdb_bac120_ar122_unique.txt
#remove temporary file
rm tmp.txt
```
3. create a Taxonomy tree in Python

```
#import library
import pickle # save binary object, like trees
import sys  # system command library
sys.path.append('/nfs/team162/kv4/github/metagm') # add the metagm library to your session
from metagm.phylogeny.TaxonomyTree import TaxonomyTree # class to generate trees
# this class only needs a list of taxonomic paths and will generate a taxonomy in NCBI format
tree = TaxonomyTree('taxo_path_gtdb_bac120_ar122_unique.txt')
# save the tree in NCBI format at the given location
tree.saveNCBIFormat('taxonomy_bac_ar')
# also save the tree in a Pickle format (Python binary) for later
tree.savePickleFormat('/nfs/team162/kv4/bin/gtdb_metadata/taxonomy_bac_ar/taxo.pyc')
```

## How to add nodes to an existing taxonomic tree

If a tree needs to be updated by adding new nodes, it is not necessary to re-build it from scratch.
User can provide a Pickle tree and new taxonomic paths:

```
#import library
import pickle  # load binary object, like trees
import csv # read csv files
import sys  # system command library
sys.path.append('/nfs/team162/kv4/github/metagm') # add the metagm library to your session
from metagm.phylogeny.TaxonomyTree import TaxonomyTree # class to generate trees

#load existing tree to add nodes
with open('taxonomy_bac_ar/taxo.pyc', "rb") as input_file:
    tree = pickle.load(input_file)

#Here we want to add non bacterial and non archeal nodes (fungi, virus and other eukaryotes)
extraNodes = '/nfs/team162/kv4/Kraken1019/nonBacterial_taxopath.txt'
with open(extraNodes) as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for i, row in enumerate(reader):
        print(row[0])
        # add a node if it is new
        tree.add_node(row[0])

# save the updated tree in NCBI format at the given location
tree.saveNCBIFormat('taxonomy')
# also save the tree in a Pickle format (Python binary) for later
tree.savePickleFormat('taxonomy/taxo.pyc')
```

# Statistical analysis

This section provides some R snippets to do analysis using the metagm_classify output files. 

More details are given in the R vignette '.Rmd', also found in this repository.

# MySQL Knowledge database

This section describes how to log in the lab mySQL database.
For obvious security reasons, the password is not provided here and needs to be requested to Nick/Hilary.

# Miscalleneous

## Other utility functions/scripts

### Quality control

# TODO
