import pandas as pd
import pathlib

# USAGE: snakemake --snakefile varcallpe.smk --configfile config.yaml --cores=5 -p


###############################################################################
#                                  CONSTANTS                                  #
###############################################################################

# WARNING!!! Global declarations are repeatedly evaluated/executed for each
# rule. Therefore, do not include anything in this section that is dynamic. As
# these are constants within the pipeline, we will use all capitals in the
# naming convention.

# Adapter Sequences ###########################################################

ADAPTERS = str()
for adapter in config['trimming']['trim adapter']:
    ADAPTERS += '-a {} '.format(adapter)
ADAPTERS = ADAPTERS.strip()

# Sample and Reference Associations ###########################################

OUTDIR = config['pipeline']['smk dir']
REFERENCE_IDS = set()
SAMPLE_IDS = set()
SAMPLES_TO_REFERENCES = dict()
REFERENCES_TO_SAMPLES = dict()
COMBO_PREFIXES = set()

for reference in config['reference genomes']:
    REFERENCE_IDS.add(reference)

for sample in config['sample key']:
    SAMPLE_IDS.add(sample)
    SAMPLES_TO_REFERENCES.update({sample: config['sample key'][sample]['reference id']})
    for reference in config['sample key'][sample]['reference id']:
        COMBO_PREFIXES.add('{}__{}'.format(sample, reference))
        if reference not in REFERENCES_TO_SAMPLES.keys():
            REFERENCES_TO_SAMPLES.update({reference: {sample}})
        else:
            REFERENCES_TO_SAMPLES[reference].add(sample)

# Output Directories ##########################################################

MAP_SET_INSTRUCTIONS_DIR = pathlib.Path(OUTDIR, '01_map_set_instructions')
LINK_REFERENCE_GENOMES_DIR = pathlib.Path(OUTDIR, '02_link_reference_genomes')
LINK_FASTQ_READS_DIR = pathlib.Path(OUTDIR, '03_link_fastq_files')
TRIM_READS_DIR =  pathlib.Path(OUTDIR, '04_trim_adapter_sequences')
FLC_READS_DIR = pathlib.Path(OUTDIR, '05_filter_low_complexity')
REPAIR_FASTQ_DIR = pathlib.Path(OUTDIR, '06_repair_fastq')
MAP_READS_DIR = pathlib.Path(OUTDIR, '07_map_reads')
ONLY_MAPPED_READS_DIR = pathlib.Path(OUTDIR, '08_write_mapped_only_sam')
PIDV_READS_DIR = pathlib.Path(OUTDIR, '09_percent_identity_variance_filter_sam')
BAM_DIR = pathlib.Path(OUTDIR, '10_sam_to_bam_mapped_reads')
NAME_SORT_BAM_DIR = pathlib.Path(OUTDIR, '11_name_sort_mapped_bam')
FIXMATE_BAM_DIR = pathlib.Path(OUTDIR, '12_fixmate_bam')
COORDINATE_SORT_BAM_DIR = pathlib.Path(OUTDIR, '13_coordinate_sort_mapped_bam')
MARK_DUPLICATES_BAM_DIR = pathlib.Path(OUTDIR, '14_mark_duplicates_bam')
BAM_TO_MPILEUP_DIR = pathlib.Path(OUTDIR, '15_bam_to_mpileup')
CALL_SNPS_DIR = pathlib.Path(OUTDIR, '16_call_snps')
CALL_INDELS_DIR = pathlib.Path(OUTDIR, '17_call_indels')
READ_COUNT_DIR = pathlib.Path(OUTDIR, '18_write_read_counts')


###############################################################################
#                                    RULES                                    #
###############################################################################

rule all:
    input:
        expand(pathlib.Path(MAP_SET_INSTRUCTIONS_DIR, '{instructions}.map'), instructions=COMBO_PREFIXES),
        expand(pathlib.Path(LINK_REFERENCE_GENOMES_DIR, '{references}.fna'), references=REFERENCE_IDS),
        expand(pathlib.Path(LINK_FASTQ_READS_DIR, '{samples}.R{pairs}.fastq'), samples=SAMPLE_IDS, pairs=[1, 2]),
        expand(pathlib.Path(CALL_SNPS_DIR, '{instructions}.snps.vcf'), instructions=COMBO_PREFIXES),
        expand(pathlib.Path(CALL_INDELS_DIR, '{instructions}.indels.vcf'), instructions=COMBO_PREFIXES),
        expand(pathlib.Path(READ_COUNT_DIR, '{instructions}.counts'), instructions=COMBO_PREFIXES)


# 01. Mapping Instructions ####################################################

# We expand the outputs here as the files are generated without wildcards and
# no input files are required. We will require the output files in All.

rule map_set_instructions:
    output:
       expand(pathlib.Path(MAP_SET_INSTRUCTIONS_DIR, '{instructions}.map'), instructions=COMBO_PREFIXES),
       expand(pathlib.Path(MAP_SET_INSTRUCTIONS_DIR, '{instructions}.map.cmd'), instructions=COMBO_PREFIXES)
    run:
        for sample in SAMPLES_TO_REFERENCES.keys():
            for ref in SAMPLES_TO_REFERENCES[sample]:
                instructions = '{}__{}.map'.format(sample, ref)
                cmd_file = pathlib.Path(MAP_SET_INSTRUCTIONS_DIR, instructions + '.cmd')
                cmd = 'touch {}'.format(pathlib.Path(MAP_SET_INSTRUCTIONS_DIR, instructions))
                shell('echo "' + cmd + '" > {}'.format(cmd_file))
                shell(cmd)


# 02. Linking Reference genomes ###############################################

# We expand the outputs here as the files are generated without wildcards and
# no input files are required. We will require the output files in All. We will
# link in all index files associated with the FASTA reference files.

rule link_reference_genomes:
    output:
        expand(pathlib.Path(LINK_REFERENCE_GENOMES_DIR, '{references}.fna'), references=REFERENCE_IDS),
        expand(pathlib.Path(LINK_REFERENCE_GENOMES_DIR, '{references}_linking.cmd'), references=REFERENCE_IDS)
    run:
        for reference in config['reference genomes']:
            cmd_file = pathlib.Path(LINK_REFERENCE_GENOMES_DIR, reference + '_linking.cmd')
            cmd = ' '.join([
                'ln -s',
                config['reference genomes'][reference] + '*',
                pathlib.Path(LINK_REFERENCE_GENOMES_DIR).resolve().as_posix()
            ])
            shell('echo "' + cmd + '" > {}'.format(cmd_file))
            shell(cmd)


# 03. Linking FASTQ Files #####################################################

# We expand the outputs here as the files are generated without wildcards and
# no input files are required. We will require the output files in All.

rule link_fastq_files:
    output:
        expand(pathlib.Path(LINK_FASTQ_READS_DIR, '{samples}.R{pairs}.fastq'), samples=SAMPLE_IDS, pairs=[1, 2]),
        expand(pathlib.Path(LINK_FASTQ_READS_DIR, '{samples}.R{pairs}.fastq.cmd'), samples=SAMPLE_IDS, pairs=[1, 2])
    run:
        for sample in config['sample key']:
            cmd_file = pathlib.Path(LINK_FASTQ_READS_DIR, '{}.R1.fastq.cmd'.format(sample))
            cmd = ' '.join([
                'ln -s',
                config['sample key'][sample]['fastq r1'],
                pathlib.Path(LINK_FASTQ_READS_DIR, '{}.R1.fastq'.format(sample)).resolve().as_posix()
            ])
            shell('echo "' + cmd + '" > {}'.format(cmd_file))
            shell(cmd)
            cmd_file = pathlib.Path(LINK_FASTQ_READS_DIR, '{}.R2.fastq.cmd'.format(sample))
            cmd = ' '.join([
                'ln -s',
                config['sample key'][sample]['fastq r2'],
                pathlib.Path(LINK_FASTQ_READS_DIR, '{}.R2.fastq'.format(sample)).resolve().as_posix()
            ])
            shell('echo "' + cmd + '" > {}'.format(cmd_file))
            shell(cmd)


# 04. Trimming Adapter Sequences ##############################################

# We use wildcards to process the linked-in FASTQ files, as they have
# previously been placed in the name space. Adapter trimming may be skipped if
# the --no-trim flag is passed. The cutadapt command has some output to STDERR
# for its processes. This section will accept fastq.gz files.

if config['pipeline']['trim'] is True:

    rule trim_adapter_sequences:
        input:
            fastq = pathlib.Path(LINK_FASTQ_READS_DIR, '{samples}.R{pairs}.fastq')
        output:
            trim = temp(pathlib.Path(TRIM_READS_DIR, '{samples}.R{pairs}.trim.fastq')),
            cmd = pathlib.Path(TRIM_READS_DIR, '{samples}.R{pairs}.trim.fastq.cmd')
        log:
            pathlib.Path(TRIM_READS_DIR, '{samples}.R{pairs}.trim_adapter_sequences.log')
        run:
            cmd = ' '.join([
                'cutadapt',
                ADAPTERS,
                '-o {output.trim}',
                '--quality-base {}'.format(config['trimming']['trim quality base']),
                '--quality-cutoff {}'.format(config['trimming']['trim quality cutoff']),
                '--minimum-length {}'.format(config['trimming']['trim minimum length']),
                '--max-n {}'.format(config['trimming']['trim max n']),
                '{input.fastq}',
                '>',
                '{log}',
                '2>&1'

            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

elif config['pipeline']['trim'] is False:

    rule trim_adapter_sequences:
        input:
            fastq = pathlib.Path(LINK_FASTQ_READS_DIR, '{samples}.R{pairs}.fastq')
        output:
            trim = temp(pathlib.Path(TRIM_READS_DIR, '{samples}.R{pairs}.trim.fastq')),
            cmd = pathlib.Path(TRIM_READS_DIR, '{samples}.R{pairs}.trim.fastq.cmd')
        run:
            cmd = ' '.join([
                'cp',
                '{input.fastq}',
                '{output.trim}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)


# 05. Low Complexity Filtering ################################################

# We use wildcards to process the FASTQ files. Low complexity filtering may be
# skipped if the --no-lcf flag is passed. This section will accept fastq.gz
# files.

if config['pipeline']['lcf'] is True:

    rule filter_low_complexity:
        input:
            fastq = pathlib.Path(TRIM_READS_DIR, '{samples}.R{pairs}.trim.fastq')
        output:
            flc = temp(pathlib.Path(FLC_READS_DIR, '{samples}.R{pairs}.flc.fastq')),
            cmd = pathlib.Path(FLC_READS_DIR, '{samples}.R{pairs}.flc.fastq.cmd')
        log:
            pathlib.Path(FLC_READS_DIR, '{samples}.R{pairs}.filter_low_complexity.log')
        run:
            cmd = ' '.join([
                'vsearch',
                '--fastx_mask {input.fastq}',
                '--min_unmasked_pct {}'.format(config['low complexity filter']['lcf min unmasked pct']),
                '--fastqout {output.flc}',
                '--hardmask',
                '2> {log}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

elif config['pipeline']['lcf'] is False:

    rule filter_low_complexity:
        input:
            fastq = pathlib.Path(TRIM_READS_DIR, '{samples}.R{pairs}.trim.fastq')
        output:
            flc = temp(pathlib.Path(FLC_READS_DIR, '{samples}.R{pairs}.flc.fastq')),
            cmd = pathlib.Path(FLC_READS_DIR, '{samples}.R{pairs}.flc.fastq.cmd')
        run:
            cmd = ' '.join([
                'cp',
                '{input.fastq}',
                '{output.flc}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)


# 06. Repair FASTQ Files ######################################################

# We use wildcards to process the FASTQ files. Repairing FASTQ files may be
# skipped if both --no-trim and --no-lcf flags were passed---i.e. no alteration
# of the original FASTQ files; else we will need repair the pair order of FASTQ
# reads. If the FASTQ files were altered by either filtering script then the
# resultant FASTQ files will be uncompressed (repair.sh does not like fastq.gz
# files).

if (config['pipeline']['trim'] is False) and (config['pipeline']['lcf'] is False):

    rule repair_fastq:
        input:
            fastq = pathlib.Path(FLC_READS_DIR, '{samples}.R{pairs}.flc.fastq')
        output:
            repair = temp(pathlib.Path(REPAIR_FASTQ_DIR, '{samples}.R{pairs}.repair.fastq')),
            cmd = pathlib.Path(REPAIR_FASTQ_DIR, '{samples}.R{pairs}.repair.fastq.cmd')
        run:
            cmd = ' '.join([
                'cp',
                '{input.fastq}',
                '{output.repair}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

else:

    rule repair_fastq:
        input:
            r1 = pathlib.Path(FLC_READS_DIR, '{samples}.R1.flc.fastq'),
            r2 = pathlib.Path(FLC_READS_DIR, '{samples}.R2.flc.fastq')
        output:
            r1 = temp(pathlib.Path(REPAIR_FASTQ_DIR, '{samples}.R1.repair.fastq')),
            r2 = temp(pathlib.Path(REPAIR_FASTQ_DIR, '{samples}.R2.repair.fastq')),
            s = temp(pathlib.Path(REPAIR_FASTQ_DIR, '{samples}.singleton.repair.fastq')),
            cmd = pathlib.Path(REPAIR_FASTQ_DIR, '{samples}.repair.fastq.cmd')
        log:
            pathlib.Path(REPAIR_FASTQ_DIR, '{samples}.repair_fastq.log')
        run:
            cmd = ' '.join([
                'repair.sh',
                'in={input.r1}',
                'in2={input.r2}',
                'out={output.r1}',
                'out2={output.r2}',
                'outs={output.s}',
                '2> {log}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)


# 07. Map Reads to Reference Genome ###########################################

# We use wildcards to process the FASTQ files, namely we require a sample to
# have post-repair R1/R2 pairs for paired-end mapping.

rule map_reads:
    input:
        map = pathlib.Path(MAP_SET_INSTRUCTIONS_DIR, '{sample}__{reference}.map'),
        r1 = pathlib.Path(REPAIR_FASTQ_DIR, '{sample}.R1.repair.fastq'),
        r2 = pathlib.Path(REPAIR_FASTQ_DIR, '{sample}.R2.repair.fastq'),
        reference = pathlib.Path(LINK_REFERENCE_GENOMES_DIR, '{reference}.fna')
    output:
        bam = temp(pathlib.Path(MAP_READS_DIR, '{sample}__{reference}.bam')),
        cmd = pathlib.Path(MAP_READS_DIR, '{sample}__{reference}.bam.cmd')
    log:
        pathlib.Path(MAP_READS_DIR, '{sample}__{reference}.map_reads.log')
    run:
        cmd = ' '.join([
            'bwa mem',
            '{input.reference}',
            '{input.r1} {input.r2}',
            '2> {log} | samtools view -bS',
            '>',
            output.bam
        ])
        shell('echo "' + cmd + '" > {output.cmd}')
        shell(cmd)


# 08. Write Mapped Only SAM ###################################################

# Write a version of the SAM  file that contains only mapped reads.

rule write_mapped_only_sam:
    input:
        pathlib.Path(MAP_READS_DIR, '{instructions}.bam')
    output:
        sam = temp(pathlib.Path(ONLY_MAPPED_READS_DIR, '{instructions}.mapped.sam')),
        cmd = pathlib.Path(ONLY_MAPPED_READS_DIR, '{instructions}.mapped.sam.cmd'),
    run:
        cmd = ' '.join([
            'samtools view',
            '-h',
            '-O SAM',
            '-F 0x4',
            '{input}',
            '>',
            '{output.sam}'
        ])
        shell('echo "' + cmd + '" > {output.cmd}')
        shell(cmd)


# 09. Percent Identity Variance Filter ########################################

# The PIDV filter step will produce a SAM that can have a mixture of PE and
# singleton reads of mapped and unmapped reads. Therefore, output from the PIDV
# filter step---or the skipped version---will be a SAM file ready for samtools
# view filtering based on SAM/BAM flags. This step is skipped if the --no-pidv
# flag is passed.

if config['pipeline']['pidv'] is True:

    rule percent_identity_variance_filter_sam:
        input:
            pathlib.Path(ONLY_MAPPED_READS_DIR, '{instructions}.mapped.sam')
        output:
            sam = temp(pathlib.Path(PIDV_READS_DIR, '{instructions}.pidv.sam')),
            log = pathlib.Path(PIDV_READS_DIR, '{instructions}.pidv.log'),
            cmd = pathlib.Path(PIDV_READS_DIR, '{instructions}.pidv.sam.cmd')
        run:
            cmd = ' '.join([
                'filter_pidv_pe.py',
                '--sam {input}',
                '--out {output.sam}',
                '--log {output.log}',
                '--maxscore {}'.format(config['pidv filter']['max pidv pct'])
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

elif config['pipeline']['pidv'] is False:

    rule percent_identity_variance_filter_sam:
        input:
            pathlib.Path(ONLY_MAPPED_READS_DIR, '{instructions}.mapped.sam')
        output:
            sam = temp(pathlib.Path(PIDV_READS_DIR, '{instructions}.pidv.sam')),
            cmd = pathlib.Path(PIDV_READS_DIR, '{instructions}.pidv.sam.cmd')
        run:
            cmd = ' '.join([
                'cp',
                '{input}',
                '{output.sam}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)


# 10. SAM to BAM Conversion ###################################################

# We now convert the SAM file to BAM for input to resolving the mpileup file
# for variant calling purposes. The resultant BAM file will be determined by
# the samtools view arguments provided in the configuration file. For example,
# if --only-properly-paired is passed, then we will only output properly paired
# reads in the BAM file.

flags = {'-F 0x4'}  # only mapped reads

if config['mapping']['no secondary'] is True:
    flags.add('-F 0x100')

if config['mapping']['no supplemental'] is True:
    flags.add('-F 0x800')

if config['mapping']['no singletons'] is True:
    flags.add('-F 0x8')

if config['mapping']['only properly paired'] is True:
    flags.add('-f 0x2')

flags_str = ' '.join(flags)

rule sam_to_bam_mapped_reads:
    input:
        pathlib.Path(PIDV_READS_DIR, '{instructions}.pidv.sam')
    output:
        bam = temp(pathlib.Path(BAM_DIR, '{instructions}.bam')),
        cmd = pathlib.Path(BAM_DIR, '{instructions}.bam.cmd')
    run:
        cmd = ' '.join([
            'samtools view',
            '-h',
            '-O BAM',
            flags_str,
            '{input}',
            '>',
            '{output.bam}'
        ])
        shell('echo "' + cmd + '" > {output.cmd}')
        shell(cmd)


# 11. Name Sort BAM ###########################################################
# 12. BAM Fixmate  ############################################################
# 13. Coordinate Sort BAM #####################################################
# 14. BAM Mark Duplicates #####################################################

# This entire section runs conceptually as a block; all steps are required for
# marking and removing mapped duplicates. Deduplication may be skipped entirely
# if --no-dedup is passed. In order to mark duplicates, we will fix the mate
# pair score information and associated fields in the BAM using samtools
# fixmate. Prior to running fixmate, the BAM must be name sorted; however, the
# BAM must be sorted on mapped read coordinates prior to running samtools
# markdup. The final BAM will have had duplicates removed.

if config['pipeline']['dedup'] is True:

    rule name_sort_mapped_bam:
        input:
            pathlib.Path(BAM_DIR, '{instructions}.bam')
        output:
            bam = temp(pathlib.Path(NAME_SORT_BAM_DIR, '{instructions}.nsort.bam')),
            cmd = pathlib.Path(NAME_SORT_BAM_DIR, '{instructions}.nsort.bam.cmd')
        run:
            cmd = ' '.join([
                'samtools sort',
                '-n',
                '-O BAM',
                '-o {output.bam}',
                '{input}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

    rule fixmate_bam:
        input:
            pathlib.Path(NAME_SORT_BAM_DIR, '{instructions}.nsort.bam')
        output:
            bam = temp(pathlib.Path(FIXMATE_BAM_DIR, '{instructions}.fixmate.bam')),
            cmd = pathlib.Path(FIXMATE_BAM_DIR, '{instructions}.fixmate.bam.cmd')
        run:
            cmd = ' '.join([
                'samtools fixmate',
                '-m',
                '-O BAM',
                '{input}',
                '{output.bam}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

    rule coordinate_sort_mapped_bam:
        input:
            pathlib.Path(FIXMATE_BAM_DIR, '{instructions}.fixmate.bam')
        output:
            bam = temp(pathlib.Path(COORDINATE_SORT_BAM_DIR, '{instructions}.sort.bam')),
            cmd = pathlib.Path(COORDINATE_SORT_BAM_DIR, '{instructions}.sort.bam.cmd')
        run:
            cmd = ' '.join([
                'samtools sort',
                '-O BAM',
                '-o {output.bam}',
                '{input}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

    rule mark_duplicates_bam:
        input:
            pathlib.Path(COORDINATE_SORT_BAM_DIR, '{instructions}.sort.bam')
        output:
            bam = temp(pathlib.Path(MARK_DUPLICATES_BAM_DIR, '{instructions}.dedup.bam')),
            cmd = pathlib.Path(MARK_DUPLICATES_BAM_DIR, '{instructions}.dedup.bam.cmd')
        log:
            pathlib.Path(MARK_DUPLICATES_BAM_DIR, '{instructions}.mark_duplicates_bam.log')
        run:
            cmd = ' '.join([
                'samtools markdup',
                '-O BAM',
                '-s',  # report stats
                '-r',  # remove duplicates
                '{input}',
                '{output.bam}',
                '2>',
                '{log}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

elif config['pipeline']['dedup'] is False:

    rule name_sort_mapped_bam:
        input:
            pathlib.Path(BAM_DIR, '{instructions}.bam')
        output:
            bam = temp(pathlib.Path(NAME_SORT_BAM_DIR, '{instructions}.nsort.bam')),
            cmd = pathlib.Path(NAME_SORT_BAM_DIR, '{instructions}.nsort.bam.cmd')
        run:
            cmd = ' '.join([
                'cp',
                '{input}',
                '{output.bam}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

    rule fixmate_bam:
        input:
            pathlib.Path(NAME_SORT_BAM_DIR, '{instructions}.nsort.bam')
        output:
            bam = temp(pathlib.Path(FIXMATE_BAM_DIR, '{instructions}.fixmate.bam')),
            cmd = pathlib.Path(FIXMATE_BAM_DIR, '{instructions}.fixmate.bam.cmd')
        run:
            cmd = ' '.join([
                'cp',
                '{input}',
                '{output.bam}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

    rule coordinate_sort_mapped_bam:
        input:
            pathlib.Path(FIXMATE_BAM_DIR, '{instructions}.fixmate.bam')
        output:
            bam = temp(pathlib.Path(COORDINATE_SORT_BAM_DIR, '{instructions}.sort.bam')),
            cmd = pathlib.Path(COORDINATE_SORT_BAM_DIR, '{instructions}.sort.bam.cmd')
        run:
            cmd = ' '.join([
                'samtools sort',
                '-O BAM',
                '-o {output.bam}',
                '{input}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)

    rule mark_duplicates_bam:
        input:
            pathlib.Path(COORDINATE_SORT_BAM_DIR, '{instructions}.sort.bam')
        output:
            bam = temp(pathlib.Path(MARK_DUPLICATES_BAM_DIR, '{instructions}.dedup.bam')),
            cmd = pathlib.Path(MARK_DUPLICATES_BAM_DIR, '{instructions}.dedup.bam.cmd')
        run:
            cmd = ' '.join([
                'cp',
                '{input}',
                '{output.bam}'
            ])
            shell('echo "' + cmd + '" > {output.cmd}')
            shell(cmd)


# 15. BAM to Mpileup ##########################################################

# In order to call variants, we must produce an input mpileup file for each BAM
# file. We match the BAM file with the associated reference genome the reads
# were aligned against. By default the allowed depth is not constrained, but
# this value may be altered as --mpileup-max-depth in the configuration file.

rule bam_to_mpileup:
    input:
        bam = pathlib.Path(MARK_DUPLICATES_BAM_DIR, '{sample}__{reference}.dedup.bam'),
        reference = pathlib.Path(LINK_REFERENCE_GENOMES_DIR, '{reference}.fna')
    output:
        mpileup = temp(pathlib.Path(BAM_TO_MPILEUP_DIR, '{sample}__{reference}.mpileup')),
        cmd = pathlib.Path(BAM_TO_MPILEUP_DIR, '{sample}__{reference}.mpileup.cmd'),
    log:
        pathlib.Path(BAM_TO_MPILEUP_DIR, '{sample}__{reference}.bam_to_mpileup.log'),
    run:
        cmd = ' '.join([
            'samtools mpileup',
            '-A',  # do not discard anomalous read pairs
            '-d {}'.format(config['variants']['mpileup max depth']),
            '-f {input.reference}',
            '{input.bam}',
            '>',
            '{output.mpileup}',
            '2>',
            '{log}'
        ])
        shell('echo "' + cmd + '" > {output.cmd}')
        shell(cmd)


# 16. Variant Call SNPs #######################################################

# We now call single nucleotide variants on the mpileup file using VarScan2.

rule call_snps:
    input:
        pathlib.Path(BAM_TO_MPILEUP_DIR, '{instructions}.mpileup')
    output:
        vcf = pathlib.Path(CALL_SNPS_DIR, '{instructions}.snps.vcf'),
        cmd = pathlib.Path(CALL_SNPS_DIR, '{instructions}.snps.vcf.cmd')
    log:
        log = pathlib.Path(CALL_SNPS_DIR, '{instructions}.call_snps.log')
    run:
        cmd = ' '.join([
            'varscan mpileup2snp',
            '{input}',
            '--min-coverage {}'.format(config['variants']['snp min coverage']),
            '--min-reads2 {}'.format(config['variants']['snp min reads2']),
            '--min-avg-qual {}'.format(config['variants']['snp min avg qual']),
            '--min-var-freq {}'.format(config['variants']['snp min var freq']),
            '--min-freq-for-hom {}'.format(config['variants']['snp min freq for hom']),
            '--p-value {}'.format(config['variants']['snp p value']),
            '--strand-filter {}'.format(config['variants']['snp strand filter']),
            '--output-vcf 1',
            '>',
            '{output.vcf}',
            '2>',
            '{log}'
        ])
        shell('echo "' + cmd + '" > {output.cmd}')
        shell(cmd)


# 17. Variant Call INDELs #####################################################

# We now call short insertion/deletion variants on the mpileup file using
# VarScan2.

rule call_indels:
    input:
        pathlib.Path(BAM_TO_MPILEUP_DIR, '{instructions}.mpileup')
    output:
        vcf = pathlib.Path(CALL_INDELS_DIR, '{instructions}.indels.vcf'),
        cmd = pathlib.Path(CALL_INDELS_DIR, '{instructions}.indels.vcf.cmd')
    log:
        log = pathlib.Path(CALL_INDELS_DIR, '{instructions}.call_indels.log')
    run:
        cmd = ' '.join([
            'varscan mpileup2indel',
            '{input}',
            '--min-coverage {}'.format(config['variants']['indel min coverage']),
            '--min-reads2 {}'.format(config['variants']['indel min reads2']),
            '--min-avg-qual {}'.format(config['variants']['indel min avg qual']),
            '--min-var-freq {}'.format(config['variants']['indel min var freq']),
            '--min-freq-for-hom {}'.format(config['variants']['indel min freq for hom']),
            '--p-value {}'.format(config['variants']['indel p value']),
            '--strand-filter {}'.format(config['variants']['indel strand filter']),
            '--output-vcf 1',
            '>',
            '{output.vcf}',
            '2>',
            '{log}'
        ])
        shell('echo "' + cmd + '" > {output.cmd}')
        shell(cmd)


# 18. Write Read Counts #######################################################

# We can write read coverage metrics for the variant calls based on the input
# mpileup file.

rule write_read_counts:
    input:
        pathlib.Path(BAM_TO_MPILEUP_DIR, '{instructions}.mpileup')
    output:
        counts = pathlib.Path(READ_COUNT_DIR, '{instructions}.counts'),
        cmd = pathlib.Path(READ_COUNT_DIR, '{instructions}.counts.cmd')
    log:
        pathlib.Path(READ_COUNT_DIR, '{instructions}.write_read_counts.log')
    run:
        cmd = ' '.join([
            'varscan readcounts',
            '{input}',
            '--min-coverage {}'.format(config['variants']['counts_min_coverage']),
            '--min-base-qual {}'.format(config['variants']['counts min base qual']),            
            '--output-file {output.counts}',
            '2>',
            '{log}'
        ])
        shell('echo "' + cmd + '" > {output.cmd}')
        shell(cmd)


# __END__

# T.N. Wylie  <twylie@wustl.edu>
# Last Update: Mon Jun 27 13:06:22 CDT 2022
