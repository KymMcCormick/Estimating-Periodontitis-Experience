prepare:
	Rscript scripts/01_prepare_nhanes_core.R

build-epe:
	Rscript scripts/02_build_epe_iteration01.R

all: prepare build-epe
