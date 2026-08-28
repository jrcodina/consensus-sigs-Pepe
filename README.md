## ConsensusSigs

In this repository you will find the instructions to run ConsensusSigs. ConsensusSigs is a tool that combines 4 signature extraction algorithms - [SigProfiler](https://github.com/AlexandrovLab/SigProfilerAssignment), [SigMiner](https://github.com/ShixiangWang/sigminer), [deconstructSigs](https://github.com/raerose01/deconstructSigs) and [MESiCA](https://github.com/Adarya/MESiCA) - and performs an agreement metric.

The input mutational catalog must be in MAF format. For more information about this format visit the [NIH documentation page.](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)

Variant data is processed to generate single base substitution (SBS) context matrices in the SBS96 format for assignments. The reference genome used for analysis is GRCh37. Signature deconvolution isbconducted using the [COSMIC SBS v3.4 reference signature set](https://cancer.sanger.ac.uk/signatures/sbs/).

Individual SBS signature assignments is further grouped into eight clinically relevant mutational processes, using the following mapping:

-   **APOBEC:** SBS2, SBS13.
-   **Ultraviolet (UV) exposure:** SBS7a, SBS7b, SBS7c, SBS7d, SBS38.
-   **Tobacco smoke:** SBS4.
-   **Homologous recombination deficiency (HRD):** SBS3.
-   **Mismatch repair (MMR) deficiency**: SBS6, SBS14, SBS15, SBS20, SBS21, SBS26, SBS44.
-   **POLE-related mutation:** SBS10a, SBS10b.
-   **Clock_SBS1:** SBS1.
-   **Clock_SBS5:** SBS5.

For each sample, a signature group is considered present if at least one of its associated SBS signatures is detected (i.e., non-zero contribution) by the respective algorithm. Unlike the other three algorithms, MESiCA directly provides signature group assignments, eliminating the need for post-processing or manual mapping.

To facilitate comparison across methods, the results from each tool are dichotomized to indicate the presence or absence of a given signature in each sample. A signature group is assigned a value of 1 (present) if its reported contribution exceeded zero, or 0 (absent) otherwise. Using this approach, a consensus matrix is generated, where each cell represents the number of tools that assigned a given signature group to a sample.

The consensus is then used to determine the most relevant signature for each sample. This consolidation was performed by retaining only those signatures for which at least three of the four algorithms showed agreement. If no signature group reached this threshold, no single dominant signature was assigned.

## Running with Docker

First, pull the container image:

```sh
docker pull jip2013/consensus-sigs:v0.2
```

You can then run the consensus-sigs pipeline using the provided container. Mount your data directory and specify the input and output paths as follows:

```sh
docker run --rm -v ./data:/data consensus-sigs -m /data/data_mutations_extended.txt -o /data
```

- `-v ./data:/data` mounts your local `data` directory to `/data` inside the container.
- `-m` specifies the input MAF file path (inside the container).
- `-o` specifies the output directory (inside the container).

## Running Algorithms Individually (Without Docker)

If you prefer to run the algorithms outside the container, you can use the provided environment lock files to set up the necessary dependencies for R or Python:

### R Environment (using `renv.lock`)

1. Install [renv](https://rstudio.github.io/renv/) in R if you don’t have it:
   ```r
   install.packages("renv")
   ```
2. In the project directory, restore the environment:
   ```r
   renv::restore()
   ```
3. You can now run any of the R scripts (e.g., `scripts/deconstructsigs_assignment.R`, `scripts/mesica_assignment.R`) with all required packages installed.

### Python Environment (using `poetry.lock`)

1. Install [Poetry](https://python-poetry.org/docs/#installation) if you don’t have it.
2. In the project directory, install dependencies:
   ```sh
   poetry install
   ```
3. Activate the Poetry environment (for zsh/bash):
   ```sh
   $(poetry env activate)
   ```
4. You can now run Python scripts (e.g., for SigProfilerAssignment) with all required packages.
