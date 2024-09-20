
# Perturb-Seq Analysis Pipeline

This is a very simple pipeline processes Perturb-Seq data by running the 10x Genomics **Cell Ranger Multi** pipeline followed by the **Seurat** pipeline for downstream analysis. It is implemented using **Nextflow DSL2** to handle the workflow orchestration.

## Features

-   Runs **Cell Ranger Multi** to process multiplexed gene expression data.
-   Uses **Seurat** for further analysis on the resulting matrix.
-   Modular design: individual components (Cell Ranger and Seurat) are organized into separate Nextflow modules.
-   Includes a help message feature to guide the user on how to use the pipeline.

## Requirements

### Software

-   **Nextflow** v20.10.0 or higher: Required for workflow orchestration.
-   **Docker** or **Singularity**: Used for containerization (alternatively, you can configure a conda or local environment if not using containers).
-   **Cell Ranger** (8.0.1 or compatible version): Installed within a container to process gene expression data.
-   **Seurat**: Installed within a container for downstream analysis.
-   **Git**: For cloning the repository and version control.

_**All Docker files are provided in the `containers` folder of the pipeline.**_

### Hardware

#### For **Cell Ranger**:

-   **CPU**: Cell Ranger is designed to run on a multi-core server. A minimum of **8 cores** is recommended, but **16 or more cores** will significantly improve performance.
-   **Memory**: Cell Ranger requires a significant amount of memory, depending on the size of the dataset. For processing 1k cells, at least **16GB of RAM** is recommended (it WILL break under 16GB). For larger datasets (e.g., 10k+ cells), **64 GB or more** is required.
-   **Disk Space**: Depending on the size of the input FASTQ files and reference genome, you may need **20 GB to 500 GB** of free disk space. Ensure that you have enough storage to handle the temporary files and output data generated by Cell Ranger.

#### For **Seurat**:

-   **CPU**: Seurat can efficiently use **multiple cores**, but it is more memory-bound. At least **4 cores** are recommended for medium-sized datasets.
-   **Memory**: Seurat requires at least **16 GB of RAM** for small datasets (e.g., 1k cells). For larger datasets, **32 GB or more** may be required to handle the memory load during the analysis.
-   **Disk Space**: Seurat typically requires **1 GB to 50 GB** of free disk space for input/output files, depending on the dataset size.

### Recommended Hardware for Full Pipeline:

For smooth execution of the full pipeline (Cell Ranger + Seurat):

-   **CPU**: 16 cores or more.
-   **Memory**: 16GB or more.
-   **Disk Space**: At least 30 GB of free space to accommodate both intermediate and output files from Cell Ranger and Seurat.

## Installation

1.  Install Nextflow:
    
    
    `curl -s https://get.nextflow.io | bash` 
    
2.  Ensure Docker or Singularity is installed on your system to run the containers.
    
3.  Clone this repository:
        
    `git clone https://github.com/srx3001/Perturb-seq` 
    

## Input Files

This pipeline expects the following input files:

-   **Analysis ID**: A unique identifier for the analysis, simultaneously Cell Ranger output name.
-   **CSV Configuration File**: A `.csv` file that contains the configuration for Cell Ranger multi.
-   **Run Directory**: The directory that contains input FASTQ files or other data required by Cell Ranger.
-   **Reference Directory**: The directory containing the reference genome for alignment and counting.
-   **Seurat Results Output Directory**: Directory where the Seurat output will be stored.

### Example `params` Configuration

`--analysis_id my_analysis
--csv_file_name config.csv
--rundir /path/to/fastq/
--reference /path/to/reference/
--seurat_results out_folder_name` 

## Usage

Run the pipeline with Nextflow using the following command:

`nextflow run main.nf --analysis_id my_analysis --csv_file_name config.csv --rundir /path/to/fastq --reference /path/to/reference --seurat_results /path/to/seurat_results` 

### Help Message

You can view the help message by running:

`nextflow run main.nf --help` 

### Workflow Description

1.  **Cell Ranger Multi**:
    
    -   The pipeline begins by running the **Cell Ranger Multi** module to process multiplexed gene expression data. This module takes in the provided CSV configuration, analysis ID, run directory, and reference genome to generate a matrix of filtered feature-barcode data.
2.  **Seurat**:
    
    -   The output matrix from Cell Ranger is then passed to the **Seurat** module for downstream analysis. The results from Seurat are saved in the specified output directory.



    
## Outputs

-   **Cell Ranger Multi Outputs**:
    -  Cell Ranger matrix file
    -  Cell Ranger web summary file
    -  Cell Ranger metrics file
-   **Seurat Outputs**:
    - Seurat plots and top differentially expressed genes for Guide and cluster comparisons


## Modular Components

### `./help`

Contains helper functions:

-   **`printHeader`**: Displays the pipeline header and information.
-   **`helpMessage`**: Displays the help message for the pipeline.

### `./modules/cell_ranger/main.nf`

This module wraps **Cell Ranger Multi** for processing gene expression data.

### `./modules/seurat/main.nf`

This module wraps **Seurat** for performing downstream analysis.

## Sample Run

To run the pipeline with sample data, follow these steps:

### 1. Download the Required Files

#### Download the Reference Data

Download the reference genome from 10x Genomics:


`wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` 

Once downloaded, extract the reference genome:

`tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz` 

This will extract the reference data into a folder called `refdata-gex-GRCh38-2024-A`.

#### Download the FASTQ Files

Download the sample FASTQ files from 10x Genomics:

`wget https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_fastqs.tar` 

Extract the FASTQ files:

`tar -xvf 1k_CRISPR_5p_gemx_Multiplex_fastqs.tar` 

This will create a folder called `1k_CRISPR_5p_gemx_Multiplex_fastqs`.

### 2. Update Configuration File

Edit the configuration CSV file (`1k_CRISPR_5p_gemx_Multiplex_config.csv`) founds in the repositories test folder to contain the correct absolute paths for your reference and FASTQ files.

Here’s an example of how to edit the file:

`[gene-expression]
reference,/absolute/path/to/refdata-gex-GRCh38-2024-A
fastqs,/absolute/path/to/1k_CRISPR_5p_gemx_Multiplex_fastqs
sample_id,1k_CRISPR_5p_gemx_Multiplex
create-bam,true` 

-   Replace `/absolute/path/to/refdata-gex-GRCh38-2024-A` with the actual path where you extracted the reference genome.
-   Replace `/absolute/path/to/1k_CRISPR_5p_gemx_Multiplex_fastqs` with the actual path where you extracted the FASTQ files.

### 3. Run the Pipeline

After updating the configuration file with the correct paths, you can run the pipeline by executing the following command:


`nextflow run main.nf \
    --analysis_id '1k_CRISPR_5p_analysis' \
    --csv_file_name '1k_CRISPR_5p_gemx_Multiplex_config.csv' \
    --rundir '/absolute/path/to/1k_CRISPR_5p_gemx_Multiplex_fastqs' \
    --reference '/absolute/path/to/refdata-gex-GRCh38-2024-A' \
    --seurat_results '/absolute/path/to/seurat_results'` 

Make sure to replace `/absolute/path/to/` with the actual paths on your system.

### 4. Review Output

-   The pipeline will output processed Cell Ranger Multi data, including the filtered feature-barcode matrix and analysis metrics.
-   The Seurat analysis results will be saved in the `seurat_results` directory.


## Tolls used are made by:

-   **Cell Ranger**: 10x Genomics
-   **Seurat**: Satija Lab