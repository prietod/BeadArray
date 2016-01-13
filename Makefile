all:
	@echo "cancel-all"
	@echo "test"
	@echo "qc"
	@echo "qc-clean-all"
	@echo "qc-avg"
	@echo "qc-avg-clean-all"

cancel-all:
	scancel -phii02 --signal=9 --full --user $$LOGNAME

qc:
	bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/BeadArray_qc.R

qc-avg:
	bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/BeadArray_qc_average.R

test:
	scancel -phii02 --signal=9 --full --user $$LOGNAME --name BeadArray_test
	-rm -rf tmp/{work,log,complete}/BeadArray_test &
	sleep 10
	bin/util/get-mapinfo -u chip_barcode | sort -R | head -7 | SUBMIT_RANDOM_SECONDS=10 bin/slurm-submit code/BeadArray_test.R
	sleep 10; find tmp/log/BeadArray_test/current/ -type f -name '*.log' | xargs tail -f

qc-clean-all:
	-rm -rf tmp/{work,complete}/BeadArray_qc

qc-avg-clean-all:
	-rm -rf tmp/{work,complete}/BeadArray_qc_average
