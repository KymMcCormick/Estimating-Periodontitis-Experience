# Applying this update to the GitHub repository

1. Open the local cloned repository in RStudio using the `.Rproj` file.
2. Copy the folders/files from this update bundle into the repository root, allowing replacements of `README.md`, `NEWS.md`, `Makefile`, `DESCRIPTION`, and `docs/output-data-dictionary.md`.
3. Confirm that the new scripts are visible in `scripts/`:
   - `03_diagnose_incident_attribution_model.R`
   - `04_summarise_epe_distribution.R`
   - `05_assess_hba1c_associational_coherence.R`
   - `06_assess_diabetes_discrimination_auc.R`
4. In RStudio, run `source("scripts/00_setup.R")`.
5. If the EPE analytic dataset already exists locally, test one downstream script, for example `source("scripts/04_summarise_epe_distribution.R")`.
6. Use the RStudio Git pane to stage the changed code/documentation files only.
7. Commit with a message such as `Add EPE validation and reporting analyses`.
8. Push to GitHub.

Do not commit raw NHANES files or generated processed datasets.
