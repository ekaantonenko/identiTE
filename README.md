# identiTE
Workflow for precise location and pangenome-identification of Transposable Elements


## Installation
### Prerequisites

You need to install the following packages (here is an example using conda environment):

```bash
conda create -n bioenv python=3.11
conda activate bioenv
conda install -c conda-forge -c bioconda pybedtools biopython=1.79 pandas
conda install -c bioconda blast
conda install yq
conda install -c conda-forge jq
```

## Usage

You need to upload the following data and specify the paths in the config file for the ***reference*** genome:
  - The reference sequence (.fasta format)
  - The gene annotation
  - The pericentromeres coordinates (.bed format)

You need to upload the following data and specify the paths in the config file for the ***query*** genomes:
  - The library file (either GraffiTE .vcf output or a .bed files are accepted)
  - The list of genomes to be analyzed
  - The sequences (.fasta format)
  - The gene annotations (liftoff, ensemble, bed formats accepted)

Use the flag --config to specify the config file, config_default.yaml is used otherwise.

Use the flags --step1, --step2, --step3, --step4 to specify the workflow steps to run; all steps are run if not specified.

Example of running on a cluster using SLURM job manager:

```bash
sbatch annotate_TIPs.slurm --config config_test.yaml --step1 --step2
```

### Step 1

Blasting of all TEs and their flanks in all genomes. 

Output:
  - blastedSummary/blastedSummary.csv file: summary of all TEs in all genomes, including start and end coordinates, quality of alignment of identified TEs in the query genomes.
  - blastedSummary/binarySummary.csv file: binary table of 0 (confirmed absence), 1 (confirmed presence), NA. 

### Step 2

Splitting the results of Step 1 into separate .bed files for each genome. 

### Steps 3 and 4

Annotation for the methylation status.



