#RUNID1 = ["1901449_nanopore_20220517", "1901449_nanopore_20220523"]
#RUNID2 = ["20220517_1140_MN38629_FAR64540_5c4e90ff", "20220523_1017_MN38629_FAR64621_15a8c2b1"]
import pandas as pd

configfile: "config.yaml"

samples = pd.read_csv(config["samples"], sep = "\t").set_index("sample", drop=False)

def get_path(wildcards):
    return samples.loc[wildcards.sample, 'path']


rule all:
	input:
			expand("/storage/scratch01/groups/bu/myeloma_altum/agostoli/nanoplot/output/{sample}",zip, sample=samples['sample'])
rule guppy:
	input:
			get_path
	output: 
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/guppy/output/{sample}/sequencing_summary.txt"
	shell:
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/ont-guppy/bin/guppy_basecaller --model_file template_r9.4.1_450bps_hac.jsn --save_path {output} --input_path {input} --device cuda:all:100% --disable_pings -c dna_r9.4.1_450bps_fast.cfg --min_qscore 7 --recursive -x 'cuda:all' --num_callers 4 --gpu_runners_per_device 4 --chunks_per_runner 1024 --chunk_size 1000 --compress_fastq"
rule nanoplot:
	input: 
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/guppy/output/{sample}/sequencing_summary.txt" 
	output:
			directory("/storage/scratch01/groups/bu/myeloma_altum/agostoli/nanoplot/output/{sample}")
	shell: 
			"NanoPlot --summary {input} --loglength -o {output}"
