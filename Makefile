all:
	@echo "qc"
	@echo "qc-avg"

qc:
	-rm -rf tmp/{work,complete}/BeadArray_qc
	bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/BeadArray_qc.R

qc-avg:
	-rm -rf tmp/{work,complete}/BeadArray_qc_average
	bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/BeadArray_qc_average.R

test:
	-rm -rf tmp/{work,complete}/BeadArray_test
	bin/util/get-mapinfo -u chip_barcode | sort -R | head -7 | bin/slurm-submit code/BeadArray_test.R
