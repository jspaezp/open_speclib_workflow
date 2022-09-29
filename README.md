
# Workflow to run a mass spec search and generate spectral libraries

1. Yaml validation
2. Discover system comet and bibliospec, or add config option
3. moularize
  - https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html


## Configuration

The config is split into 3 main configurations:
1. The default configuration, which is applied to all the files if it is not
        overridden by the other configurations.
2. The 'experiment' configuration, which is applied to all the files in the
        experiment.
        - The 'experiment' bundles files that share search parameters and database.
3. The 'combine' configuration, which is applied to all the files in the
        combine.
        - The 'combine' bundles experiments that should be in the same spectral
        library.

The section

```yaml
title: Supercoolrun
default:
  fdr_level: protein # currently not used
  comet:
    peptide_mass_tolerance: 20
    num_threads: 10
  fasta:
    type: uniprot
    value: UP000005640
  BiblioSpec:
    cutoff: 0.99
combine:
  - name: SuperCoolCombination
    experiments:
      - Firstexp
      - Secondexp
  - name: Only first exp
    experiments:
      - Firstexp
experiments:
  - name: Firstexp
    files:
      - /path/to/raw/file1.mzML
      - /path/to/raw/otherfile.mzML
    fasta:
      type: uniprot
      value: UP000005640
    comet:
      peptide_mass_tolerance: 20
  - name: Secondexp
    files:
      - /path/to/others.mzML
      - /path/to/others2.mzML
    fasta:
      type: uniprot
      value: UP000000589
    comet:
      peptide_mass_tolerance: 50
    BiblioSpec:
      cutoff: 0.95

```

## Running

You need to provide a yaml configuration with the `--configfile` option

You can run this to check the workflow with that dummy configuration:

```shell
poetry run snakemake --verbose --cores 1 --directory $PWD -s snakefile.smk --configfile ./reference_files/run.yml --dry-run
```

## Project directory structure

```
.
├── README.md ####################### This file #########################################
├── build_tools
│   └── docker_build_bibliospec ##### Has a container that builds BiblioSpec ############
│       ├── DOCKERFILE
│       ├── build.bash
│       └── pwiz #################### proteowizard ######################################
│           ├── .....
├── extras
│   └── msconvert_gsutil.sh # Utility script to convert files inside google cloud services
├── open_speclib_workflow # placeholder directory
├── poetry.lock
├── pyproject.toml
├── raw_files
│   ├── 20181221_N15cm_3_Hela_2ug_10min_R1.mzML
│   ├── 20181221_N15cm_3_Hela_2ug_10min_R1.raw
│   ├── HUMAN.fasta
│   ├── README.md
│   └── comet.params
├── reference_files ################# Files to use as a reference for the workflow ######
│   ├── comet.params.new
│   ├── fakefiles
│   │   ├── file1.mzML
│   │   ├── foo.mzML
│   │   ├── otherfile.mzML
│   │   └── others2.mzML
│   └── run.yml
├── requirements.txt
├── snakefile.smk ################### <<<<<<< main snakemake file >>>>>>>> ###############
└── test
    ├── 20181210_15cm_3_HelaDilution_30min_0ng_28Hz_R1.mzML
    ├── 20181210_15cm_3_HelaDilution_30min_1000ng_28Hz_R1.mzML
    ├── 20181210_15cm_3_HelaDilution_30min_1000ng_41Hz_R1.mzML
    ├── DOCKERFILE
    ├── README.md
    ├── convert_data.bash
    ├── dockerized_run.bash
    ├── download_data.bash
    ├── download_truncated.bash
    └── run.yml

```

## Testing

### Full workflow run

Run the dockerized version in the `tests` directory

### Check management

`tox` runs a dry run on the dummy data inside `reference_files`
