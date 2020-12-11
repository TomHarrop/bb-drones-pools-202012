#!/usr/bin/env python3

import pandas

honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.12')
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'
phase_honeybee_vcf = 'shub://TomHarrop/phase-honeybee-vcf:phase_honeybee_vcf_v0.0.2'

ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
fai = f'{ref}.fai'

sample_info = 'config/bb_gt_table.csv'
cnv_map = 'config/cnv_map.txt'

sample_df = pandas.read_csv(sample_info,
                            index_col='sample')
all_samples = sorted(set(sample_df.index))


rule target:
    input:
        'output/030_phasing/phased.vcf.gz'

rule phase:
    input:
        ref = ref,
        vcf = 'output/020_filtered-genotypes/filtered.vcf.gz',
        bam = 'output/010_genotypes/merged.bam',
        csv = 'config/phasing.csv'
    output:
        'output/030_phasing/phased.vcf.gz'
    log:
        'output/logs/phase.log'
    params:
        outdir = 'output/030_phasing'
    threads:
        workflow.cores
    container:
        phase_honeybee_vcf
    shell:
        'phase_honeybee_vcf '
        '--threads {threads} '
        '--ref {input.ref} '
        '--vcf {input.vcf} '
        '--bam {input.bam} '
        '--samples_csv {input.csv} '
        '--outdir {params.outdir} '
        '2> {log}'

rule filter:
    input:
        vcf = 'output/010_genotypes/calls.vcf.gz'
    output:
        temp('output/020_filtered-genotypes/filtered.vcf')
    params:
        min_maf = 0.05,
        f_missing = 0.2
    log:
        'output/logs/filter.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '-m2 -M2 -v snps '  # biallelic snps only
        '--min-af {params.min_maf}:nonmajor '
        '--exclude "F_MISSING>{params.f_missing}" '
        '{input.vcf} '
        '> {output} '
        '2> {log}'


# genotype
checkpoint genotype:
    input:
        csv = sample_info,
        ref = ref,
        cnv_map = cnv_map,
        reads = expand('data/reads/{sample}_r{r}.fastq.gz',
                       sample=all_samples,
                       r=[1, 2])
    output:
        cutoffs = 'output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/calls.vcf.gz',
        bam = 'output/010_genotypes/merged.bam',
        ref = 'output/010_genotypes/015_ref/ref.fasta',
        fai = 'output/010_genotypes/015_ref/ref.fasta.fai'
    params:
        wd = 'output/010_genotypes',
    log:
        'output/logs/genotype.log'
    threads:
        workflow.cores
    singularity:
        honeybee_genotype_pipeline
    shell:
        'honeybee_genotype_pipeline '
        '--ref {input.ref} '
        '--samples_csv {input.csv} '
        '--outdir {params.wd} '
        '--cnv_map {input.cnv_map} '
        '--threads {threads} '
        '--csd '
        '--restart_times 1 '
        '&>> {log}'


# generic index rule
rule index_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        gz = 'output/{folder}/{file}.vcf.gz',
        tbi = 'output/{folder}/{file}.vcf.gz.tbi'
    log:
        'output/logs/{folder}/{file}_index-vcf.log'
    # wildcard_constraints:
    #     file = 'filtered|phased_drones|phased_pools'
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} 2> {log} '
        '; '
        'tabix -p vcf {output.gz} 2>> {log}'
