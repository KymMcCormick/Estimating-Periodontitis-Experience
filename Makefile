setup:
	Rscript scripts/00_setup.R

prepare:
	Rscript scripts/01_prepare_nhanes_core.R

build-epe:
	Rscript scripts/02_build_epe_iteration01.R

incident-diagnostics:
	Rscript scripts/03_diagnose_incident_attribution_model.R

distribution:
	Rscript scripts/04_summarise_epe_distribution.R

hba1c-coherence:
	Rscript scripts/05_assess_hba1c_associational_coherence.R

diabetes-auc:
	Rscript scripts/06_assess_diabetes_discrimination_auc.R

validation: incident-diagnostics distribution hba1c-coherence diabetes-auc

all: prepare build-epe validation
