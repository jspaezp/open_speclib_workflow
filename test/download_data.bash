#!/bin/bash

set -x
set -e

files="20181210_15cm_3_HelaDilution_30min_0ng_28Hz_R1.raw  20181210_15cm_3_HelaDilution_30min_0ng_41Hz_R1.raw   20181210_15cm_3_HelaDilution_30min_1000ng_28Hz_R1.raw  20181210_15cm_3_HelaDilution_30min_1000ng_41Hz_R1.raw"

for curr_file in $files ; do
    if ! [ -f ${curr_file} ]; then
        curl ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/12/PXD015087/${curr_file} --output ${curr_file}
    fi
    gsutil cp ${curr_file}  gs://open_speclib_workflow
done