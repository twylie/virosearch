#! /usr/bin/python3.7

import argparse
import pathlib
from virosearch.references.get_ref import GetRef
from virosearch.references.index_ref import IndexRef
from virosearch.support.stdout_message import StdoutMessage
from virosearch.variants.var_call_pe import VarCallPE
from virosearch.variants.var_call_dir import VarCallDir
from virosearch.variants.var_call_washu import VarCallWashu

version = 'ViroSearch (Revision 18)'


###############################################################################
#                                  FUNCTIONS                                  #
###############################################################################

def collect_cli_arguments():

    description = """
    -- References
       getref         Download reference files from Entrez by acc id.
       indexref       Index reference files in a ViroSearch DB directory.

    -- Mapping and Variant Calling
       varcallpe      Paired end mapping and variant calling.
       washupp        WashU Only: Add parallel processing setup."""

    # Top-level argument parser.

    parser = argparse.ArgumentParser(
        description='*** VIROSEARCH ***',
        prog='virosearch',
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Optional arguments.

    parser.add_argument(
        '-h',
        '--help',
        action='help',
        help='Display the extended usage statement.'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=version,
        help='Display the software version number.'
    )

    # Sub-parsers.

    subparsers = parser.add_subparsers(
        title='Subcommands:',
        description=description,
        required=True,
        dest='command'
    )

    # getref ##################################################################

    subparser_getref = subparsers.add_parser(
        'getref',
        description='*** Download reference files from Entrez. ***'
    )

    getref_required_group = subparser_getref.add_argument_group('required arguments')

    getref_required_group.add_argument(
        '--email',
        metavar='STR',
        type=str,
        action='store',
        help='The email address to be used when downloading genome sequences and associated metadata. Entrez requires an email address for anonymous connection.',
        required=True,
    )

    getref_required_group.add_argument(
        '--virosearch-db',
        metavar='PATH',
        type=pathlib.Path,
        action='store',
        help='The path to write the resultant ViroSearch DB reference genome directory. May be a new or preexisting directory. An existing directory will trigger evaluation of underlying files and assumes directory extension is desired.',
        required=True,
    )

    getref_required_group.add_argument(
        '--acc-id',
        metavar='STR',
        type=str,
        action='store',
        help='Supply a list (space-delimited) of accession ids for all of the reference genome sequences you wish to include in the ViroSearch DB directory. Already existing references will be skipped. Running the getref command multiple times with new accession ids will extend/add references to the ViroSearch DB directory.',
        required=True,
        nargs='+'
    )

    # indexref ################################################################

    subparser_indexref = subparsers.add_parser(
        'indexref',
        description='*** Prepare/index files in a ViroSearch DB directory. ***'
    )

    indexref_required_group = subparser_indexref.add_argument_group('required arguments')

    indexref_required_group.add_argument(
        '--virosearch-db',
        metavar='PATH',
        type=pathlib.Path,
        action='store',
        help='The path to an existing ViroSearch DB reference genome directory. All underlying reference genomes in this directory will be evaluated for indexes required for alignment; missing indexes will be created as needed.',
        required=True,
    )

    # washupp #################################################################

    subparser_washupp = subparsers.add_parser(
        'washupp',
        description='*** Add parallel processing support for WashU users. ***'
    )

    washupp_required_group = subparser_washupp.add_argument_group('required arguments')

    washupp_required_group.add_argument(
        '--results-dir',
        metavar='PATH',
        type=pathlib.Path,
        action='store',
        help='Path to existing ViroSearch varcallpe --results-dir directory. This directory prepares ViroSearch to run the variant calling pipeline; washupp will extend this directory for parallel processing.',
        required=True
    )

    washupp_required_group.add_argument(
        '--docker-volumes',
        metavar='PATH',
        type=str,
        nargs='+',
        action='store',
        help='A space-delimited list of volume paths to map under Docker. Any volume used for processing data under ViroSearch should be included here.',
        required=True
    )

    washupp_required_group.add_argument(
        '--lsf-compute-group',
        metavar='STR',
        type=str,
        action='store',
        help='Provide a valid LSF compute group name for processing such as compute-twylie.',
        required=True
    )

    washupp_optional_group = subparser_washupp.add_argument_group('optional parameters')

    washupp_optional_group.add_argument(
        '--docker-image',
        metavar='STR',
        type=str,
        action='store',
        help='Name of the ViroSearch Docker image to pull from DockerHub. [twylie/virosearch:latest]',
        required=False,
        default='twylie/virosearch:latest'
    )

    washupp_optional_group.add_argument(
        '--lsf-memory',
        metavar='STR',
        type=str,
        action='store',
        help='Amount of memory to request from LSF per child process. [16G]',
        required=False,
        default='16G'
    )

    washupp_optional_group.add_argument(
        '--lsf-queue',
        metavar='STR',
        type=str,
        action='store',
        help='Provide a valid LSF compute queue name for processing. [general]',
        required=False,
        default='general'
    )

    washupp_optional_group.add_argument(
        '--lsf-latency-wait',
        metavar='INT',
        type=int,
        action='store',
        help='Wait given seconds if an output file of a job is not present after the job finished. [15]',
        required=False,
        default='15'
    )

    washupp_optional_group.add_argument(
        '--lsf-restart-times',
        metavar='INT',
        type=int,
        action='store',
        help='Number of times to restart failing jobs. [3]',
        required=False,
        default='3'
    )

    washupp_optional_group.add_argument(
        '--lsf-log-dir',
        metavar='PATH',
        type=pathlib.Path,
        action='store',
        help='Path to write LSF log files from child processes. [lsf_logs]',
        required=False,
        default='lsf_logs'
    )

    # varcallpe ###############################################################

    subparser_varcallpe = subparsers.add_parser(
        'varcallpe',
        description='*** Mapping and variant calling of samples. ***'
    )

    varcallpe_required_group = subparser_varcallpe.add_argument_group('required arguments')

    varcallpe_required_group.add_argument(
        '--virosearch-db',
        metavar='PATH',
        type=pathlib.Path,
        action='store',
        help='The path to the existing ViroSearch DB reference genome directory. Used for mapping reads against.',
        required=True,
    )

    varcallpe_required_group.add_argument(
        '--results-dir',
        metavar='PATH',
        type=pathlib.Path,
        action='store',
        help='The path to the write the processing/results directory. This directory will prepare ViroSearch to run the variant calling pipeline.',
        required=True,
    )

    varcallpe_required_group.add_argument(
        '--sample-key',
        metavar='PATH',
        type=pathlib.Path,
        action='store',
        help='Path to an existing, formatted sample key (tab-delimited) file. This file associated samples with reference genomes for mapping and variant calling.',
        required=True,
    )

    varcallpe_optional_group = subparser_varcallpe.add_argument_group('optional parameters')

    varcallpe_optional_group.add_argument(
        '--smk-dir',
        metavar='PATH',
        type=pathlib.Path,
        action='store',
        help='Directory path for Snakemake to write pipeline output files. By default, the directory is within the --results-dir directory. [smk]',
        required=False
    )

    varcallpe_optional_group.add_argument(
        '--trim',
        action='store_true',
        help='Runs the adapter trimming portion of the pipeline. By default this is True and does not require specifying. [True]'
    )

    varcallpe_optional_group.add_argument(
        '--no-trim',
        dest='trim',
        action='store_false',
        help='Skips adapter trimming step. [False]'
    )

    varcallpe_optional_group.add_argument(
        '--lcf',
        action='store_true',
        help='Runs the low complexity filtering portion of the pipeline. By default this is True and does not require specifying. [True]'
    )

    varcallpe_optional_group.add_argument(
        '--no-lcf',
        dest='lcf',
        action='store_false',
        help='Skips the low complexity filter portion of the pipeline. [False]'
    )

    varcallpe_optional_group.add_argument(
        '--pidv',
        action='store_true',
        help='Runs the percent identity variation filtering portion of the pipeline. By default this is True and does not require specifying. [True]'
    )

    varcallpe_optional_group.add_argument(
        '--no-pidv',
        dest='pidv',
        action='store_false',
        help='Skips percent identity variation filtering step. [False]'
    )

    varcallpe_optional_group.add_argument(
        '--dedup',
        action='store_true',
        help='Runs the mapped read deduplication portion of the pipeline. By default this is True and does not require specifying. [True]'
    )

    varcallpe_optional_group.add_argument(
        '--no-dedup',
        dest='dedup',
        action='store_false',
        help='Skips mapped read deduplication step. [False]'
    )

    varcallpe_optional_group.set_defaults(
        trim=True,
        lcf=True,
        pidv=True,
        dedup=True
    )

    varcallpe_optional_group.add_argument(
        '--no-secondary',
        action='store_true',
        help='Remove reads flagged as secondary alignments prior to variant calling. [False]'
    )

    varcallpe_optional_group.add_argument(
        '--no-supplemental',
        action='store_true',
        help='Remove reads flagged as supplemental alignments prior to variant calling. [False]'
    )

    varcallpe_optional_group.add_argument(
        '--no-singletons',
        action='store_true',
        help='Remove reads flagged as singleton alignments prior to variant calling. [False]'
    )

    varcallpe_optional_group.add_argument(
        '--only-properly-paired',
        action='store_true',
        help='Enforce properly paired read alignments when calling variants. [False]'
    )

    varcallpe_optional_group.add_argument(
        '--mpileup-max-depth',
        metavar='INT',
        action='store',
        default=0,
        type=int,
        help='Max per-file depth as used by mpileup to avoid excessive memory usage. Default of 0 does not define a maximum depth. [0]',
    )

    varcallpe_optional_group.add_argument(
        '--trim-adapter',
        metavar='STR',
        action='store',
        help='Space delimited list of adapter sequences to trim. The default sequences are Illumina TruSeq single index and TruSeq CD index sequences. [AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]',
        type=str,
        default=['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
        nargs='+'
    )

    varcallpe_optional_group.add_argument(
        '--trim-quality-base',
        metavar='INT',
        action='store',
        default=33,
        type=int,
        help='Sets FASTQ quality encoding of 33 or (older) 64 for trimming purposes. [33]',
    )

    varcallpe_optional_group.add_argument(
        '--trim-quality-cutoff',
        metavar='INT',
        action='store',
        default=10,
        type=int,
        help='Trim 3\'-end when quality drops below given value; may alter sequence length. [10]',
    )

    varcallpe_optional_group.add_argument(
        '--trim-minimum-length',
        metavar='INT',
        action='store',
        default=50,
        type=int,
        help='Minimum required read length after trimming, or will be removed from downstream processing. [50]'
    )

    varcallpe_optional_group.add_argument(
        '--trim-max-n',
        metavar='FLOAT',
        action='store',
        default=0.5,
        type=int,
        help='Maximum percent of Ns allowed per read post-trimming; may remove sequences from downstream processing. [0.5]'
    )

    varcallpe_optional_group.add_argument(
        '--lcf-min-unmasked-pct',
        metavar='INT',
        action='store',
        default=50,
        type=int,
        help='Minimum required unmask percent for read filtering; may remove sequences from downstream processing. [50]'
    )

    varcallpe_optional_group.add_argument(
        '--max-pidv-pct',
        metavar='FLOAT',
        action='store',
        default=0.15,
        type=float,
        help='Maximum percent identity variance (compared to reference) allowed for nucleotide alignments; may remove sequences from downstream processing. [0.15]'
    )

    varcallpe_optional_group.add_argument(
        '--snp-min-coverage',
        metavar='INT',
        action='store',
        default=8,
        type=int,
        help='Minimum read depth at a position to make a SNP call. [8]'
    )

    varcallpe_optional_group.add_argument(
        '--snp-min-reads2',
        metavar='INT',
        action='store',
        default=2,
        type=int,
        help='Minimum supporting reads at a position to call SNP variants. [2]'
    )

    varcallpe_optional_group.add_argument(
        '--snp-min-avg-qual',
        metavar='INT',
        action='store',
        default=15,
        type=int,
        help='Minimum base quality at a position to count a read. [15]'
    )

    varcallpe_optional_group.add_argument(
        '--snp-min-var-freq',
        metavar='FLOAT',
        action='store',
        default=0.01,
        type=float,
        help='Minimum variant allele frequency threshold for SNP calls. [0.01]'
    )

    varcallpe_optional_group.add_argument(
        '--snp-min-freq-for-hom',
        metavar='FLOAT',
        action='store',
        default=0.75,
        type=float,
        help='Minimum frequency to call homozygote. [0.75]'
    )

    varcallpe_optional_group.add_argument(
        '--snp-p-value',
        metavar='FLOAT',
        action='store',
        default=0.99,
        type=float,
        help='Default p-value threshold for calling SNP variants. [0.99]'
    )

    varcallpe_optional_group.add_argument(
        '--snp-strand-filter',
        metavar='INT',
        action='store',
        default=1,
        type=int,
        help=r'Ignore variants with more than 90 percent support on one strand for SNPs. (1=True; 0=False) [1]'
    )

    varcallpe_optional_group.add_argument(
        '--indel-min-coverage',
        metavar='INT',
        action='store',
        default=8,
        type=int,
        help='Minimum read depth at a position to make a INDEL call. [8]'
    )

    varcallpe_optional_group.add_argument(
        '--indel-min-reads2',
        metavar='INT',
        action='store',
        default=2,
        type=int,
        help='Minimum supporting reads at a position to call INDEL variants. [2]'
    )

    varcallpe_optional_group.add_argument(
        '--indel-min-avg-qual',
        metavar='INT',
        action='store',
        default=15,
        type=int,
        help='Minimum base quality at a position to count a read. [15]'
    )

    varcallpe_optional_group.add_argument(
        '--indel-min-var-freq',
        metavar='FLOAT',
        action='store',
        default=0.01,
        type=float,
        help='Minimum variant allele frequency threshold for INDEL calls. [0.01]'
    )

    varcallpe_optional_group.add_argument(
        '--indel-min-freq-for-hom',
        metavar='FLOAT',
        action='store',
        default=0.75,
        type=float,
        help='Minimum frequency to call homozygote. [0.75]'
    )

    varcallpe_optional_group.add_argument(
        '--indel-p-value',
        metavar='FLOAT',
        action='store',
        default=0.99,
        type=float,
        help='Default p-value threshold for calling INDEL variants. [0.99]'
    )

    varcallpe_optional_group.add_argument(
        '--indel-strand-filter',
        metavar='INT',
        action='store',
        default=1,
        type=int,
        help='Ignore variants with more than 90 percent support on one strand for INDELs. (1=True; 0=False) [1]'
    )

    varcallpe_optional_group.add_argument(
        '--counts-min-coverage',
        metavar='INT',
        action='store',
        default=1,
        type=int,
        help='Minimum read depth at a position to make a call and include in count coverage. [1]'
    )

    varcallpe_optional_group.add_argument(
        '--counts-min-base-qual',
        metavar='INT',
        action='store',
        default=20,
        type=int,
        help='Minimum base quality at a position to count a read and include in count coverage. [20]'
    )

    return parser.parse_args()


def execute_get_ref_cmd(arguments):

    # Created a ViroSearch DB directory and downloads required sequence files
    # from Entrez based on a list of accession IDs.

    StdoutMessage(
        command='getref',
        function='execute_get_ref_cmd',
        message='BEGIN'
    )

    get_ref = GetRef(arguments)
    get_ref.make_db_dir()
    get_ref.download_entrez_references()

    StdoutMessage(
        command='getref',
        function='execute_get_ref_cmd',
        message='END'
    )

    return


def execute_index_ref_cmd(arguments):

    # When pointed to a preexisting ViroSearch DB directory, this function will
    # index (e.g. bwa index) the underlying sequence files as needed for
    # downstream processing.

    StdoutMessage(
        command='indexref',
        function='execute_index_ref_cmd',
        message='BEGIN'
    )

    index_ref = IndexRef(arguments)
    index_ref.index_references()

    StdoutMessage(
        command='indexref',
        function='execute_index_ref_cmd',
        message='END'
    )

    return


def execute_var_call_pe_cmd(arguments):

    # This function runs a Snakemake pipeline that maps and calls variants
    # given a sample key of FASTQ files and an associated ViroSearch formatted
    # reference directory.

    StdoutMessage(
        command='varcallpe',
        function='execute_var_call_pe_cmd',
        message='BEGIN'
    )

    var_call_pe = VarCallPE(arguments)
    var_call_pe.eval_dir_db()
    var_call_pe.eval_reference_indexes()
    var_call_pe.make_processing_outdir()
    var_call_pe.make_smk_dir()
    var_call_pe.eval_samplekey()
    var_call_pe.write_config_file()
    var_call_pe.copy_varcallpe_recipe()
    var_call_pe.write_example_smk_cmd()

    StdoutMessage(
        command='varcallpe',
        function='execute_var_call_pe_cmd',
        message='END'
    )

    return


def execute_washupp_cmd(arguments):

    # This function takes a pre-existing results directory and adds parallel
    # processing support for people at Washington University School of Medicine
    # using RIS LSF compute.

    StdoutMessage(
        command='washupp',
        function='execute_washupp_cmd',
        message='BEGIN'
    )

    var_call_dir = VarCallDir(arguments.results_dir)
    var_call_dir.load_directory()

    var_call_washu = VarCallWashu(arguments)
    var_call_washu.make_lsf_logs_directory()
    var_call_washu.copy_submit_lsf_script()
    var_call_washu.write_yaml_config_file()
    var_call_washu.write_example_smk_cmd(var_call_dir)

    StdoutMessage(
        command='washupp',
        function='execute_washupp_cmd',
        message='END'
    )

    return


def eval_varcallpe_arguments(arguments):

    if arguments.trim_quality_base:
        if (arguments.trim_quality_base != 33) and (arguments.trim_quality_base != 64):
            StdoutMessage(
                command='varcallpe',
                function='eval_varcallpe_arguments',
                message='Option --trim-quality-base must be 33 or 64 in value.',
                fatal=True
            )

    return


###############################################################################
#                                     MAIN                                    #
###############################################################################

if __name__ == '__main__':

    arguments = collect_cli_arguments()

    # Determine which sub-command is being used, and execute that block of
    # code. Selections are as follows:

    if arguments.command == 'getref':
        execute_get_ref_cmd(arguments)
    elif arguments.command == 'indexref':
        execute_index_ref_cmd(arguments)
    elif arguments.command == 'washupp':
        execute_washupp_cmd(arguments)
    elif arguments.command == 'varcallpe':
        eval_varcallpe_arguments(arguments)
        execute_var_call_pe_cmd(arguments)

# __END__

# T.N. Wylie  <twylie@wustl.edu>
# Last Update: Mon Jun 27 11:29:51 CDT 2022
