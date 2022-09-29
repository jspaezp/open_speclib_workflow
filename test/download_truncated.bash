
files="20181210_15cm_3_HelaDilution_30min_0ng_28Hz_R1.raw  20181210_15cm_3_HelaDilution_30min_0ng_41Hz_R1.raw   20181210_15cm_3_HelaDilution_30min_1000ng_28Hz_R1.raw  20181210_15cm_3_HelaDilution_30min_1000ng_41Hz_R1.raw"

for curr_file in $files ; do
    gsutil cp gs://truncated_speclib_workflow/${curr_file%.raw}.mzML .
done