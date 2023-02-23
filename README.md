## Colocalization analysis pipeline

### 1. Lift coordinates

Execute `./lift.sh` to lift GWAS variant coordinates from GRCh36 to GRCh38.

### 2. Prepare input files

Run `prep_input_coloc.R` to create input files formatted for eQTL Catalogue queries and coloc.

### 3. Query eQTL catalogue and perform colocalization

Run the script `run_coloc.R` to obtain data from eQTL catalogue and perform colocalization.

This step was performed in an interactive job with 8 CPUs and 64GB of RAM.
