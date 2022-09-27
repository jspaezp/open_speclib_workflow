from loguru import logger
import requests

DEFAULT_CONFIG = config["default"]
combine_configs = {x["name"]: DEFAULT_CONFIG.copy() for x in config["combine"]}

for x in config["combine"]:
    combine_configs[x["name"]].update(x)

logger.info(combine_configs)

experiment_configs = {x["name"]: DEFAULT_CONFIG.copy() for x in config["experiments"]}
for x in config["experiments"]:
    experiment_configs[x["name"]].update(x)

logger.info(experiment_configs)


rule all:
    input:
        # expand("{sample}.{param}.output.pdf", sample=config["samples"], param=config["yourparam"])
        "foo",
        psms = "./experiment_assets/{experiment}/mokapot/mokapot.psms.txt",
        pin="./experiment_assets/{experiment}/comet/{raw_file}.comet.pin",

for exp_name in experiment_configs:
    shell(
        f"mkdir -p ./experiment_assets/{exp_name}/{{raw,mzml,comet,fasta,mokapot,bibliospec}}"
    )


rule convert_raw:
    input:
        "./experiment_assets/{experiment}/raw/{raw_file}",
    output:
        "./experiment_assets/{experiment}/mzml/{raw_file}.mzML"
    run:
        raise NotImplementedError


rule get_fasta:
    output:
        "./experiment_assets/{experiment}/fasta/{fasta_file}.fasta",
    run:
        fasta_name = experiment_configs[experiment]["fasta"]["name"]
        fasta_type = experiment_configs[experiment]["fasta"]["type"]
        if fasta_type == "url":
            shell(f"wget {fasta_name} -O {output}")
        elif fasta_type == "file":
            shell(f"cp {fasta_name} {output}")
        elif fasta_type == "unirpot":

            url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3A{PROTEOME}%29%20AND%20%28reviewed%3Atrue%29"
            url = url.format(PROTEOME=fasta_name)

            all_fastas = requests.get(url).text
            with open(output, "w") as f:
                f.write(all_fastas)
        else:
            raise Exception(f"Unknown fasta type {fasta_type}")



def update_comet_config(infile, outfile, config_dict: dict):
    with open(infile, "r") as f:
        with open(output, "w+", encoding="utf-8") as of:
            for line in f:
                for key in config_dict:
                    if line.startswith(key):
                        line = f"{key} = {config_dict[key]}\n"

                of.write(line)


rule generate_comet_config:
    output:
        param_file = "./experiment_assets/{experiment}/{experiment}.comet.params",
    run:
        comet_config = experiment_configs[experiment]["comet"]
        update_comet_config("reference_files/comet.params.new", output.param_file, comet_config)


# rule decoy_fasta:
#   """
#   Generates a concatenated decoy database fasta file
#   """
#   input: "./experiment_assets/{experiment}/fasta/{fasta_file}.fasta"
#   output: "./experiment_assets/{experiment}/fasta/{fasta_file}.decoy.fasta"
#   run:
#     pyteomics.fasta.write_decoy_db(
#       source=None,
#       output=None,
#       mode='reverse',
#       prefix="DECOY_",
#       decoy_only=False)

"""
Future versions of comet might support this better.
Pre-making the index
rule comet_index:
  input:
    fasta = "./experiment_assets/{experiment}/fasta/{fasta_file}.fasta",
    comet_config = "./comet_params/{experiment}.comet.params"
  output: "./experiment_assets/{experiment}/comet/{fasta_file}.fasta.idx"
  run:
    shell(f"comet -P{input.comet_config} -D{input.fasta}")
"""


rule comet_exe:
    """
    Downloads comet from github
    """
    output:
        "./bin/comet.exe",
    params:
        source="https://github.com/UWPR/Comet/releases/download/v2022.01.2/comet.linux.exe",
    run:
        shell("mkdir -p ./bin")
        shell("wget {params.source} -O {output}")
        shell("chmod +x {output}")

rule build_bibliospec:
    """
    Builds bibliospec and stores the binaries

    Maybe we shoudl contact the Brendan to have them make a release for linux ...
    """
    output:
        "./bin/BlibBuild",
        "./bin/BlibFilter",
        "./bin/BlibSearch",
        "./bin/BlibToMs2",
    threads: 8
    run:
        build_cmd = (
            "bash ./quickbuild.sh "
            f"-j{threads} --abbreviate-paths "
            "--hash optimization=space "
            "address-model=64 "
            "pwiz pwiz_tools/BiblioSpec "
            "executables gcc"
        )

        shell("mkdir -p ./bin")
        shell("cd ./bin && git clone --depth 1 https://github.com/ProteoWizard/pwiz.git")
        shell(f"cd ./bin/pwiz && {build_cmd}")
        shell("cp pwiz/build-*/BiblioSpec/*.tar.bz2 ./bin/.")
        shell("tar -x --bzip2 -f ./bin/bibliospec-bin-*-*-release-*.tar.bz2")
        shell("cd ./pwiz && bash clean.sh")



rule comet_search:
    """Uses comet to search a single mzml file

    Every run takes ~1.5 CPU hours, so in a 20 CPU machine, every file takes 5 minutes.
    Usually 16gb of memory is more than enough for any file.
    """
    input:
        fasta="./experiment_assets/{experiment}/fasta/{fasta_file}.fasta",
        mzml="./experiment_assets/{experiment}/mzml/{raw_file}.mzml",
        comet_config = "./experiment_assets/{experiment}/{experiment}.comet.params",
        comet_executable="./bin/comet.exe",
    output:
        mzid_decoy="./experiment_assets/{experiment}/comet/{raw_file}.comet.decoy.mzid",
        mzid="./experiment_assets/{experiment}/comet/{raw_file}.comet.mzid",
        pepxml_decoy="./experiment_assets/{experiment}/comet/{raw_file}.comet.decoy.pep.xml",
        pepxml="./experiment_assets/{experiment}/comet/{raw_file}.comet.pep.xml",
        sqt_decoy="./experiment_assets/{experiment}/comet/{raw_file}.comet.decoy.sqt",
        sqt="./experiment_assets/{experiment}/comet/{raw_file}.comet.sqt",
        txt_decoy="./experiment_assets/{experiment}/comet/{raw_file}.comet.decoy.txt",
        txt="./experiment_assets/{experiment}/comet/{raw_file}.comet.txt",
        pin="./experiment_assets/{experiment}/comet/{raw_file}.comet.pin",
    params:
        base_name="./experiment_assets/{experiment}/comet/{raw_file}.comet",
    threads: 20
    run:
        update_dict = {
            "num_threads": threads,
            "output_sqtfile": 1,
            "output_txtfile": 1,
            "output_pepxmlfile": 1,
            "output_mzidentmlfile": 1,
            "output_percolatorfile": 1,
            "decoy_search": 2,
        }
        with tmpdir.tempfile() as tmp_config:
            update_comet_config(input.comet_config, tmp_config, update_dict)
            shell_cmd = f"{comet_executable} -P{tmp_config} -D{input.fasta} -N{params.base_name} {input.mzml}"
            shell(shell_cmd)

rule mokapot:
    input:
        pin_files = lambda x: xml
    output:
        decoy_peptides = "./experiment_assets/{experiment}/mokapot/mokapot.decoy.peptides.txt",
        decoy_proteins = "./experiment_assets/{experiment}/mokapot/mokapot.decoy.proteins.txt",
        decoy_psms = "./experiment_assets/{experiment}/mokapot/mokapot.decoy.psms.txt",
        peptides = "./experiment_assets/{experiment}/mokapot/mokapot.peptides.txt",
        proteins = "./experiment_assets/{experiment}/mokapot/mokapot.proteins.txt",
        psms = "./experiment_assets/{experiment}/mokapot/mokapot.psms.txt",
    params:
        out_dir = "./experiment_assets/{experiment}/mokapot/",
    run:
        shell_cmd = (
            'mokapot --keep_decoys '
            '--enzyme "[KR]" '
            '--proteins HUMAN.fasta '
            '--decoy_prefix "DECOY_" '
            '--aggregate '
            f'--dest_dir {out_dir} '
            " ".join(input.pin_files)
        )
        shell(shell_cmd)

