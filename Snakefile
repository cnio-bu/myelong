import pandas as pd

configfile: "/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/myelong/envs/config.yaml"

samples = pd.read_csv(config["samples"], sep = "\t").set_index("sample", drop=False)

def get_path(wildcards):
    return samples.loc[wildcards.sample, 'path']


rule all:
	input:
			expand("/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/nanoplot/output/{sample}",zip, sample=samples['sample'])
rule guppy:
	input:
			get_path
	output: 
			directory("/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/guppy/output/{sample}")
	benchmark: 
			"benchmarks/{sample}.tsv"
	threads:
			8
	resources: 
			mem_mb= 32000,
			runtime= 1440
	conda:
			"guppy.yaml"   
	shell:
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/ont-guppy/bin/guppy_basecaller --model_file template_r9.4.1_450bps_hac.jsn --save_path {output} --input_path {input} --device cuda:all:100% --disable_pings -c dna_r9.4.1_450bps_fast.cfg --min_qscore 7 --recursive -x 'cuda:all' --num_callers 4 --gpu_runners_per_device 4 --chunks_per_runner 1024 --chunk_size 1000 --compress_fastq"

rule nanoplot:
	input: 
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/guppy/{sample}"
	output:
			directory("/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/nanoplot/{sample}")
	benchmark: "benchmarks/{sample}.tsv"
	threads: 
			8
	resources: 
			mem_mb= 32000,
			runtime= 1440
	conda:
			"nanoplot.yaml"   
	shell: 
			"NanoPlot --summary {input}/sequencing_summary.txt --loglength -o {output}"

rule concat:
	input: 
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/guppy/{sample}"
	output:
			"{sample}_fastq.gz"
	benchmark: "benchmarks/{sample}.tsv"
	threads: 
			8
	resources: 
			mem_mb= 32000,
			runtime= 1440     
	shell:
			"cat {input}/pass/*.fastq.gz > {output}"
			
rule minimap2:
	input: 
			fastq="{sample}_fastq.gz",
			genome="/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/GRCh38.primary_assembly.genome.fa.gz"
	output:
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/minimap/{sample}.sam"
	params:
			extra="-a -x map-ont"
	benchmark: "benchmarks/{sample}.tsv"
	threads: 
			8
	resources: 
			mem_mb= 32000,
			runtime= 1440
	conda:
			"minimap2.yaml"     
	shell: 
			"minimap2 {input.genome} {input.fastq} > {output}"

rule samtools_conversion_to_bam:
	input:
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/minimap/{sample}.sam"
	output:
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/samtools/{sample}.bam"
	conda:
			"envs/samtools.yaml"
	benchmark: "benchmarks/{sample}.tsv"
	threads: 
			8
	resources: 
			mem_mb= 32000,
			runtime= 1440
	conda:
			"samtools.yaml"      
	shell:
			"samtools view -bo {input} {output}"

rule samtools_sort:
    input:
        	"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/samtools/{sample}.bam"
    output:
        	"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/samtools/{sample}.sorted.bam"
    benchmark: "benchmarks/{sample}.tsv"
	threads: 
			8
	resources: 
			mem_mb= 32000,
			runtime= 1440
	conda:
			"samtools.yaml"
	shell:
			"samtools sort -o {output} {input}"

rule samtools_index:
    input:
        	"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/samtools/{sample}.sorted.bam"
    output:
        	"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/samtools/{sample}.sorted.bam.bai"
    benchmark: "benchmarks/{sample}.tsv"
	threads: 
			8
	resources: 
			mem_mb= 32000,
			runtime= 1440   
	conda:
			"samtools.yaml"
	shell:
        	"samtools index {input}"
rule sniffles:
	input:
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/samtools/{sample}.sorted.bam",
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/samtools/{sample}.sorted.bam.bai"
	output:
			"/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/sniffles/{sample}.vcf"
	conda:
			"envs/sniffles.yaml"
	benchmark: "benchmarks/{sample}.tsv"
	threads: 
			8
	resources: 
			mem_mb= 32000,
			runtime= 1440   
	shell:
			"sniffles --input {input} --vcf {output}"