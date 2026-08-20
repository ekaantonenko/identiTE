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

### Running on a cluster using SLURM job manager
```bash
sbatch annotate_TIPs.slurm --config config_poly.yaml --step1 --step2
```
