
**This in still in progress, it will not run**

2. Yaml validation


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

```shell
poetry run snakemake --dry-run -s snakefile.smk --configfile ./reference_files/run.yml
```