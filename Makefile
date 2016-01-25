.PHONY: test

slurm_partition = hii02
jinfiniti_data_dir = /hiidata/teddy/data/jinfiniti/gene_expression
raw_data_dir = tmp/data/raw
filtered_data_dir = $$(pwd)/tmp/data/filtered
chip_list = tmp/work/all/chip-list.txt
raw_data_exclude_patterns = static/exclude-data-raw.txt
filtered_data_exclude_patterns = tmp/work/all/exclude-data-filtered.txt

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

stage-1: setup-dirs generate-chip-list copy-raw-data
	cat $(chip_list) | bin/slurm-submit-array bin/run-R-chips-qc code/BeadArray_qc.R
	cat $(chip_list) | bin/slurm-submit-array bin/run-R-chips-qc code/BeadArray_qc_average.R

stage-2: combine-qc combine-qc-average run-sample-filter copy-filtered-data

stage-3: run-approach-a-step-1 run-method-n-step-1

setup-dirs:
	[[ -d tmp/work/all ]] || mkdir -p tmp/work/all
	[[ -d $(raw_data_dir) ]] || mkdir -p $(raw_data_dir)
	[[ -d $(filtered_data_dir) ]] || mkdir -p $(filtered_data_dir)
	[[ -d tmp/work/BeadArray_qc.R/combined ]] || mkdir -p tmp/work/BeadArray_qc.R/combined
	[[ -d tmp/work/BeadArray_qc_average.R/combined ]] || mkdir -p tmp/work/BeadArray_qc_average.R/combined

generate-chip-list:
	bin/util/generate-chip-list $(jinfiniti_data_dir)/map_info.txt > $(chip_list)

copy-raw-data:
	bin/util/copy-data $(chip_list) $(raw_data_exclude_patterns) $(jinfiniti_data_dir) $(raw_data_dir)

combine-qc:
	cat static/qc_header.txt > tmp/work/BeadArray_qc.R/combined/raw-qc-details.txt
	find tmp/work/BeadArray_qc.R/chip -name '*_raw_qc_details.txt' \
		| xargs -n1 -I{} cat {} >> tmp/work/BeadArray_qc.R/combined/raw-qc-details.txt

combine-qc-average:
	cat static/qc_header.txt > tmp/work/BeadArray_qc_average.R/combined/raw-qc-details.txt
	find tmp/work/BeadArray_qc_average.R/chip -name '*_raw_qc_details.txt' \
		| xargs -n1 -I{} cat {} >> tmp/work/BeadArray_qc_average.R/combined/raw-qc-details.txt

run-sample-filter:
	bin/run-sample-filter $$(pwd)/tmp/work/BeadArray_qc.R/combined raw-qc-details.txt
	bin/run-sample-filter $$(pwd)/tmp/work/BeadArray_qc_average.R/combined raw-qc-details.txt
	diff $$(pwd)/tmp/work/BeadArray_qc{,_average}.R/combined/raw-qc-details.txt
	diff $$(pwd)/tmp/work/BeadArray_qc{,_average}.R/combined/Exclude_sample_list.txt

copy-filtered-data:
	awk '{print $$1}' $$(pwd)/tmp/work/BeadArray_qc.R/combined/Exclude_sample_list.txt \
		| sed -n '2,$$p' > $(filtered_data_exclude_patterns)
	bin/util/copy-data $(chip_list) $(filtered_data_exclude_patterns) $(raw_data_dir) $(filtered_data_dir)

run-approach-a-step-1:
	ls $(filtered_data_dir) | \
		bin/slurm-submit-array \
			bin/run-R-chips-qc-results \
			code/BeadArray_approach_a_step_1.R

run-method-n-step-1:
	for n in {1..5}; do \
	 ls $(filtered_data_dir) | \
		bin/slurm-submit-array \
			bin/run-R-chips-results \
			code/BeadArray_method_$${n}_step_1.R; \
			sleep 10; done

clean:
	-rm -rf tmp/work tmp/log tmp/slurm-array

cancel-all:
	scancel -p$(slurm_partition) --signal=9 --full --user $$LOGNAME

test:
	mkdir -p tmp/test/fix-data
	bin/fix-data test/fix-data/input/{Avg_Signal.txt,BEAD_STDERR.txt,Avg_NBEADS.txt,Detection_Pval.txt} \
		> tmp/test/fix-data/result.txt
	diff tmp/test/fix-data/result.txt test/fix-data/expected/result.txt
	mkdir -p tmp/test/fix-data
	bin/fix-data test/fix-data/input/{Avg_Signal-7-cols.txt,BEAD_STDERR-7-cols.txt,Avg_NBEADS-7-cols.txt,Detection_Pval-7-cols.txt} \
		> tmp/test/fix-data/result-7-cols.txt
	diff tmp/test/fix-data/result-7-cols.txt test/fix-data/expected/result-7-cols.txt

q:
	squeue -p$(slurm_partition) --user $$LOGNAME

test-array-job:
	printf "foo\nbar\n" | SUBMIT_RANDOM_SECONDS="1" bin/slurm-submit-array bin/test-array-job


# combine: combine-qc combine-qc-average
#
# test-fix-data:
# 	mkdir -p tmp/test/fix-data
# 	bin/fix-data test/fix-data/input/{Avg_Signal.txt,BEAD_STDERR.txt,Avg_NBEADS.txt,Detection_Pval.txt} \
# 		> tmp/test/fix-data/result.txt
# 	diff tmp/test/fix-data/result.txt test/fix-data/expected/result.txt
#
# test-fix-data-7-cols:
# 	mkdir -p tmp/test/fix-data
# 	bin/fix-data test/fix-data/input/{Avg_Signal-7-cols.txt,BEAD_STDERR-7-cols.txt,Avg_NBEADS-7-cols.txt,Detection_Pval-7-cols.txt} \
#
# qc:
# 	bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/BeadArray_qc.R
#
# qc-average:
# 	bin/util/get-mapinfo -u chip_barcode | bin/slurm-submit code/BeadArray_qc_average.R
#
# test:
# 	scancel -phii02 --signal=9 --full --user $$LOGNAME --name BeadArray_test
# 	-rm -rf tmp/{work,log,complete}/BeadArray_test &
# 	sleep 10
# 	bin/util/get-mapinfo -u chip_barcode | sort -R | head -7 | SUBMIT_RANDOM_SECONDS=10 bin/slurm-submit code/BeadArray_test.R
# 	sleep 10; find tmp/log/BeadArray_test/current/ -type f -name '*.log' | xargs tail -f
#
# qc-clean-all:
# 	-rm -rf tmp/{work,complete}/BeadArray_qc
#
# qc-average-clean-all:
# 	-rm -rf tmp/{work,complete}/BeadArray_qc_average
