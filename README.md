# PROSS - Protein Repair One Stop Shop (In-House Adapted Version)

This repository contains an in-house adapted version of PROSS (Protein Repair One Stop Shop), a computational method for improving the stability and expression levels of proteins through sequence design without compromising their functionality.

## Overview

This adapted version enhances the original PROSS pipeline by enabling parallel processing via GNU Parallel and integrating additional features tailored for in-house use. It automates the protein stabilization process by orchestrating several bioinformatics tools and databases.

## Dependencies

Before using this software, ensure that the following dependencies are installed and properly configured:

- **NCBI BLAST 2.2.27+**: Used for sequence alignment and generating PSSMs.
- **UniRef90 Database**: Required for running BLAST searches. Must be formatted using `makeblastdb`.
- **Rosetta**: Any installation of Rosetta software is required for protein modeling.
- **Python Packages**:
  - **rosetta_finder**: Helps in locating the Rosetta binary.
  - **Biopython**: Used for processing PDB files and sequences.
- **GNU Parallel**: Enables parallel processing; must be installed and available in your system's PATH.

## Installation

### Install NCBI BLAST 2.2.27+

Download and install BLAST 2.2.27+ from the [NCBI website](https://www.ncbi.nlm.nih.gov/books/NBK52640/).

### Install and Format UniRef90 Database

Download the UniRef90 database from the [UniProt website](https://www.uniprot.org/downloads). Format the database using `makeblastdb`:

```bash
makeblastdb -in uniref90.fasta -dbtype prot -out /path/to/uniref90_db
```

### Install Rosetta

Obtain and install Rosetta from the [Rosetta Commons](https://www.rosettacommons.org/software/license-and-download). Follow the official installation instructions.

### Install Python Packages

Install the required Python packages using pip:

```bash
pip install rosetta_finder biopython
```

### Install GNU Parallel

Install GNU Parallel via your package manager or from source:

- **Ubuntu/Debian**:

  ```bash
  sudo apt-get install parallel
  ```

- **macOS (using Homebrew)**:

  ```bash
  brew install parallel
  ```

- **From Source**:

  Download from the [official website](https://www.gnu.org/software/parallel/) and follow the installation instructions.

## Usage

The main script for running the PROSS pipeline is `run_pross_parallel.sh`.

The script can be executed with the following options:

```bash
Usage: ./run_pross_parallel.sh <OPTIONS>

Required Parameters:
      -s <init_pdb>          Input PDB file.

Optional Parameters:
      -f <res_to_fix>        Comma-separated list of residue IDs to fix.
      -r <res_to_restrict>   Comma-separated list of residue IDs to restrict mutations.
      -p <pssm_fn>           PSSM filename (skip PSSM generation if provided).
      -l <ligand_fas>        FASTA file of ligands.
      -j <nproc>             Number of processors for parallel execution (default: maximum available).
      -c <chain_id>          Chain ID to process (default: A).
      -d true|false          Delay execution (default: false).
      -o <output_dir>        Output directory (default: derived from job name and instance).
      -u <uniref90_db>       Path/prefix to UniRef90 database (must be formatted with `makeblastdb`).
      -B <blast_bin>         Path to NCBI BLAST binaries (default: directory containing `psiblast`).
      -h                     Print help message and exit.
```

### Parameters Description

- `-s <init_pdb>`: **(Required)** Path to the initial PDB file of the protein to be stabilized.
- `-f <res_to_fix>`: Comma-separated residue IDs (e.g., `10,20,30`) that should remain fixed during the design process.
- `-r <res_to_restrict>`: Comma-separated residue IDs where mutations should be restricted.
- `-p <pssm_fn>`: Filename of a precomputed Position-Specific Scoring Matrix (PSSM) to use. If provided, the PSSM generation step will be skipped.
- `-l <ligand_fas>`: FASTA file containing ligand sequences, if applicable.
- `-j <nproc>`: Number of processors to use for parallel execution (default is the maximum available as reported by `nproc`).
- `-c <chain_id>`: Chain identifier from the PDB file to process (default is `A`).
- `-d true|false`: Whether to delay execution (useful for scheduling; default is `false`).
- `-o <output_dir>`: Directory where output files will be saved (default is derived from the job name and instance).
- `-u <uniref90_db>`: Path or prefix to the formatted UniRef90 database (must be formatted using `makeblastdb`).
- `-B <blast_bin>`: Path to the directory containing NCBI BLAST binaries. Defaults to the directory containing `psiblast` in your PATH.
- `-h`: Display the help message and exit.

### Example Commands

#### Basic Execution

```bash
./run_pross_parallel.sh -s input_structure.pdb
```

#### Specify Residues to Fix and Restrict

```bash
./run_pross_parallel.sh -s input_structure.pdb -f 10,20,30 -r 15,25,35
```

#### Use Multiple Processors and Specify Output Directory

```bash
./run_pross_parallel.sh -s input_structure.pdb -j 8 -o /path/to/output_directory
```

#### Provide Precomputed PSSM

```bash
./run_pross_parallel.sh -s input_structure.pdb -p precomputed_pssm.txt
```

### Running the Test Case

A test case can be executed to verify the installation and functionality of the PROSS pipeline.

Run the test script:

```bash
bash tests/scripts/run_3afp_hf3_A.sh
```

This script uses the PDB file `3afp_hf3_A` and runs the PROSS pipeline.

#### PSSM Files Location

The PSSM files for the test case are located at:

```
tests/outputs/3afp_hf3_A/pssm_msa
```

#### Report File

After the test run completes, a report file will be generated at:

```
./tests/outputs/3afp_hf3_A/report/3afp_hf3_A_PROSS_report.txt
```

This report file contains the PDB file paths of the PROSS designs generated during the test run.

## Pipeline Steps

The script automates the following steps:

1. **Preparation of PSSM Profile**:
   - If a PSSM file is not provided via the `-p` option, the script generates one using `psiblast`.
   - Requires the UniRef90 database formatted with `makeblastdb`.

2. **Constraint File Generation**:
   - Generates a constraints file (`*.cst`) from the initial PDB structure using an AWK script.

3. **Structure Refinement**:
   - Refines the input structure using Rosetta scripts to prepare for mutation scanning.

4. **Filter Scan**:
   - Performs a filter scan to evaluate the impact of mutations at each residue position.
   - Utilizes parallel processing to speed up computation.

5. **Resfile Generation**:
   - Merges individual resfiles generated during the filter scan into comprehensive resfiles for design.

6. **Design Generation**:
   - Produces designed protein sequences with enhanced stability based on the resfiles.
   - Outputs models and scores for each design level.

7. **Report Generation**:
   - Compiles a report summarizing the designs and their respective paths.

## Output

The script generates the following outputs in the specified output directory:

- **Designed Sequences**: Protein sequences with predicted enhanced stability.
- **Model Structures**: PDB files of the designed protein models.
- **Log Files**: Detailed logs of the execution and processing steps.
- **Reports**: Summary reports of the design results and statistics.

For the test case, after running `bash tests/scripts/run_3afp_hf3_A.sh`, the outputs will be located under `tests/outputs/3afp_hf3_A/`.

The report file `3afp_hf3_A_PROSS_report.txt` in `tests/outputs/3afp_hf3_A/report/` contains the paths to the generated PROSS designs.

## Notes

- Ensure all dependencies are correctly installed and accessible in your system's PATH.
- The script uses GNU Parallel to speed up computations. Adjust the `-j` parameter according to your system's capabilities.
- If processing large proteins or datasets, sufficient computational resources are required.
- The default number of processors (`nproc`) is set to the maximum available, as determined by the `nproc` command.
- The `res_to_fix` and `res_to_restrict` options accept comma-separated residue IDs directly, not files.
- For testing purposes, PSSM files are located in `tests/outputs/3afp_hf3_A/pssm_msa`.

## Troubleshooting

- **BLAST Errors**: Ensure that BLAST 2.2.27+ is installed and that the UniRef90 database is correctly formatted with `makeblastdb`.
- **Rosetta Issues**: Verify that Rosetta is installed and that `rosetta_finder` can locate the binaries.
- **Python Package Errors**: Confirm that `rosetta_finder` and `biopython` are installed in your Python environment.
- **GNU Parallel Not Found**: Make sure GNU Parallel is installed and in your PATH.
- **PSSM Generation Issues**: If you experience issues during PSSM generation, consider providing a precomputed PSSM using the `-p` option.

## License

This in-house adapted version of PROSS is subject to the original PROSS licensing terms. Please refer to the original authors for licensing information.

## Acknowledgments

- **PROSS Developers**: The original PROSS method was developed by the Fleishman Lab at the Weizmann Institute of Science.
- **Rosetta Commons**: For providing the Rosetta software suite.
- **GNU Parallel**: For enabling efficient parallel processing.

## References

- **PROSS Web Server**: [https://pross.weizmann.ac.il/](https://pross.weizmann.ac.il/)
- **Rosetta Software**: [https://www.rosettacommons.org/](https://www.rosettacommons.org/)
- **NCBI BLAST**: [https://blast.ncbi.nlm.nih.gov/Blast.cgi](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- **UniProt UniRef Databases**: [https://www.uniprot.org/downloads](https://www.uniprot.org/downloads)
- **GNU Parallel**: [https://www.gnu.org/software/parallel/](https://www.gnu.org/software/parallel/)

---

**Please ensure all software and databases are used in accordance with their respective licenses and terms of use.**