all:
	@echo "clean"
	@echo "combine"
	@echo "combine-qc"
	@echo "combine-qc-average"
	@echo "cancel-all"
	@echo "test"
	@echo "qc"
	@echo "qc-clean-all"
	@echo "qc-average"
	@echo "qc-average-clean-all"
	@echo "test-fix-data"

clean:
	-rm -rf tmp/work tmp/log tmp/slurm-array

cancel-all:
	scancel -phii02 --signal=9 --full --user $$LOGNAME

combine: combine-qc combine-qc-average

test-fix-data:
	mkdir -p tmp/test/fix-data
	bin/fix-data test/fix-data/input/{Avg_Signal.txt,BEAD_STDERR.txt,Avg_NBEADS.txt,Detection_Pval.txt} \
		> tmp/test/fix-data/result.txt
	diff tmp/test/fix-data/result.txt test/fix-data/expected/result.txt

test-fix-data-7-cols:
	mkdir -p tmp/test/fix-data
	bin/fix-data test/fix-data/input/{Avg_Signal-7-cols.txt,BEAD_STDERR-7-cols.txt,Avg_NBEADS-7-cols.txt,Detection_Pval-7-cols.txt} \
		> tmp/test/fix-data/result-7-cols.txt
	diff tmp/test/fix-data/result-7-cols.txt test/fix-data/expected/result-7-cols.txt

combine-qc:
	mkdir -p tmp/work/BeadArray_qc/combined
	cat static/qc_header.txt > tmp/work/BeadArray_qc/combined/raw_qc_details.txt
	find tmp/work/BeadArray_qc/chips -name '*_raw_qc_details.txt' \
		| xargs -n1 -I{} cat {} >> tmp/work/BeadArray_qc/combined/raw_qc_details.txt
	bin/run-sample-filter tmp/work/BeadArray_qc/combined/raw_qc_details.txt \
		tmp/work/BeadArray_qc/combined/qc_exclude_sample_list.txt

combine-qc-average:
	mkdir -p tmp/work/BeadArray_qc_average/combined
	cat static/qc_header.txt > tmp/work/BeadArray_qc_average/combined/raw_qc_average_details.txt
	find tmp/work/BeadArray_qc_average/chips -name '*_raw_qc_average_details.txt' \
		| xargs -n1 -I{} cat {} >> tmp/work/BeadArray_qc_average/combined/raw_qc_average_details.txt
	bin/run-sample-filter tmp/work/BeadArray_qc_average/combined/raw_qc_average_details.txt \
		tmp/work/BeadArray_qc/combined/qc_average_exclude_sample_list.txt

qc:
	bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/BeadArray_qc.R

qc-average:
	bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/BeadArray_qc_average.R

test:
	scancel -phii02 --signal=9 --full --user $$LOGNAME --name BeadArray_test
	-rm -rf tmp/{work,log,complete}/BeadArray_test &
	sleep 10
	bin/util/get-mapinfo -u chip_barcode | sort -R | head -7 | SUBMIT_RANDOM_SECONDS=10 bin/slurm-submit code/BeadArray_test.R
	sleep 10; find tmp/log/BeadArray_test/current/ -type f -name '*.log' | xargs tail -f

qc-clean-all:
	-rm -rf tmp/{work,complete}/BeadArray_qc

qc-average-clean-all:
	-rm -rf tmp/{work,complete}/BeadArray_qc_average
