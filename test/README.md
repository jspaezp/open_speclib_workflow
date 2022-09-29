
# Data provenance

Data from the kuster lab, downloaded from pride

`test/download_data.bash`

Then data was processed in a google cloud instance to make it a .mzML, compress it and extract a subsection (2000 spectra) (`test/convert_data.bash`).

Then data was re-downloaded:
`test/download_truncated.bash`

