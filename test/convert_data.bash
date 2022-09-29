#!/bin/bash

set -x
set -e

files="20181210_15cm_3_HelaDilution_30min_0ng_28Hz_R1.raw  20181210_15cm_3_HelaDilution_30min_0ng_41Hz_R1.raw   20181210_15cm_3_HelaDilution_30min_1000ng_28Hz_R1.raw  20181210_15cm_3_HelaDilution_30min_1000ng_41Hz_R1.raw"

for curr_file in $files ; do
    mkdir -p data
    gsutil cp gs://open_speclib_workflow/${curr_file} ./data/.

    sudo docker run -it --rm -v $PWD/data:/data proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses \
        wine msconvert \
        /data/${curr_file} \
        -o /data  \
        --filter "peakPicking true 1-" \
        --zlib \
        --32 \
        --filter "scanNumber [10000,12000]" \
        --filter "zeroSamples removeExtra" \
        --filter "threshold count 3000 most-intense"
    # TODO change it so it extracts based on time and not scan number

    gsutil cp ./data/*.mzML gs://truncated_speclib_workflow
    rm -rf data
done