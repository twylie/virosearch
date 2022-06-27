from ..support.stdout_message import StdoutMessage
import yaml
import pathlib


###############################################################################
#                                VarCallConfig                                #
###############################################################################

# This class will take CLI arguments and format and write a configuration file
# (YAML) suitable for importing into Snakmake for pipeline executation.


class VarCallConfig:

    def __init__(self, verbose=True):

        self.verbose = verbose
        self.config_file = None
        self.config_file_yaml = None

        return

    def populate_arguments(self, arguments, dir_db, sample_key):

        # EXAMPLE:
        #
        # command='varcallpe',
        # counts_min_base_qual=20,
        # counts_min_coverage=1,
        # dedup=True,
        # indel_min_avg_qual=15,
        # indel_min_coverage=8,
        # indel_min_freq_for_hom=0.75,
        # indel_min_reads2=2,
        # indel_min_var_freq=0.01,
        # indel_p_value=0.99,
        # indel_strand_filter=1,
        # lcf=True,
        # lcf_min_unmasked_pct=50,
        # max_pidv_pct=0.15,
        # mpileup_max_depth=0,
        # no_secondary=False,
        # no_singletons=False,
        # no_supplemental=False,
        # only_properly_paired=False,
        # pidv=True,
        # results_dir=PosixPath('/tmp/floop'),
        # sample_key=PosixPath('pe_samplekey.tsv'),
        # snp_min_avg_qual=15,
        # snp_min_coverage=8,
        # snp_min_freq_for_hom=0.75,
        # snp_min_reads2=2,
        # snp_min_var_freq=0.01,
        # snp_p_value=0.99,
        # snp_strand_filter=1,
        # trim=True,
        # trim_adapter=['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
        # 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
        # trim_max_n=0.5,
        # trim_minimum_length=50,
        # trim_quality_base=33,
        # trim_quality_cutoff=10,
        # virosearch_db=PosixPath('ref_EBV')
        # smk_dir=PosixPath('/tmp/floop/smk')

        self.dir_db = dir_db
        self.smk_dir = arguments.smk_dir
        self.sample_key = sample_key
        self.sample_key_path = arguments.sample_key
        self.command = arguments.command
        self.counts_min_base_qual = arguments.counts_min_base_qual
        self.counts_min_coverage = arguments.counts_min_coverage
        self.dedup = arguments.dedup
        self.indel_min_avg_qual = arguments.indel_min_avg_qual
        self.indel_min_coverage = arguments.indel_min_coverage
        self.indel_min_freq_for_hom = arguments.indel_min_freq_for_hom
        self.indel_min_reads2 = arguments.indel_min_reads2
        self.indel_min_var_freq = arguments.indel_min_var_freq
        self.indel_p_value = arguments.indel_p_value
        self.indel_strand_filter = arguments.indel_strand_filter
        self.lcf = arguments.lcf
        self.lcf_min_unmasked_pct = arguments.lcf_min_unmasked_pct
        self.max_pidv_pct = arguments.max_pidv_pct
        self.mpileup_max_depth = arguments.mpileup_max_depth
        self.no_secondary = arguments.no_secondary
        self.no_singletons = arguments.no_singletons
        self.no_supplemental = arguments.no_supplemental
        self.only_properly_paired = arguments.only_properly_paired
        self.pidv = arguments.pidv
        self.results_dir = arguments.results_dir
        self.snp_min_avg_qual = arguments.snp_min_avg_qual
        self.snp_min_coverage = arguments.snp_min_coverage
        self.snp_min_freq_for_hom = arguments.snp_min_freq_for_hom
        self.snp_min_reads2 = arguments.snp_min_reads2
        self.snp_min_var_freq = arguments.snp_min_var_freq
        self.snp_p_value = arguments.snp_p_value
        self.snp_strand_filter = arguments.snp_strand_filter
        self.trim = arguments.trim
        self.trim_adapter = arguments.trim_adapter
        self.trim_max_n = arguments.trim_max_n
        self.trim_minimum_length = arguments.trim_minimum_length
        self.trim_quality_base = arguments.trim_quality_base
        self.trim_quality_cutoff = arguments.trim_quality_cutoff
        self.virosearch_db = arguments.virosearch_db

        StdoutMessage(
            command='VarCallConfig',
            function='populate_arguments',
            message='Config populated with CLI arguments.'
        )

        return

    def write_config_file(self, config_file):

        # We will format and write the config file as simple YAML. Snakemake
        # will load the configuration file prior to running the pipeline and
        # this is how it gets the arguments passed from the CLI.

        if config_file.is_file():
            StdoutMessage(
                command='VarCallConfig',
                function='write_config_file',
                message='Configuration file already exists.',
                fatal=True
            )

        with open(config_file, 'w') as fho:

            # Required arguments.

            fho.write('required:\n')
            fho.write('  sample key: \'{}\'\n'.format(self.sample_key_path.resolve()))
            fho.write('  virosearch db: \'{}\'\n'.format(self.virosearch_db.resolve()))
            fho.write('  results dir: \'{}\'\n'.format(self.results_dir.resolve()))

            # Pipeline specific arguments.

            fho.write('pipeline:\n')
            fho.write('  smk dir: \'{}\'\n'.format(self.smk_dir))
            fho.write('  command: \'{}\'\n'.format(self.command))
            fho.write('  trim: {}\n'.format(self.trim))
            fho.write('  lcf: {}\n'.format(self.lcf))
            fho.write('  pidv: {}\n'.format(self.pidv))
            fho.write('  dedup: {}\n'.format(self.dedup))

            # Mapping and variant calling arguments.

            fho.write('mapping:\n')
            fho.write('  no secondary: {}\n'.format(self.no_secondary))
            fho.write('  no supplemental: {}\n'.format(self.no_supplemental))
            fho.write('  no singletons: {}\n'.format(self.no_singletons))
            fho.write('  only properly paired: {}\n'.format(self.only_properly_paired))

            # Trimming parameters.

            fho.write('trimming:\n')
            fho.write('  trim adapter:\n')
            for adapter in self.trim_adapter:
                fho.write('    - \'{}\'\n'.format(adapter))
            fho.write('  trim max n: {}\n'.format(self.trim_max_n))
            fho.write('  trim minimum length: {}\n'.format(self.trim_minimum_length))
            fho.write('  trim quality base: {}\n'.format(self.trim_quality_base))
            fho.write('  trim quality cutoff: {}\n'.format(self.trim_quality_cutoff))

            # Low complexity filter.

            fho.write('low complexity filter:\n')
            fho.write('  lcf min unmasked pct: {}\n'.format(self.lcf_min_unmasked_pct))

            # Variant read counts.

            fho.write('variants:\n')
            fho.write('  counts min base qual: {}\n'.format(self.counts_min_base_qual))
            fho.write('  counts_min_coverage: {}\n'.format(self.counts_min_coverage))
            fho.write('  indel min avg qual: {}\n'.format(self.indel_min_avg_qual))
            fho.write('  indel min coverage: {}\n'.format(self.indel_min_coverage))
            fho.write('  indel min freq for hom: {}\n'.format(self.indel_min_freq_for_hom))
            fho.write('  indel min reads2: {}\n'.format(self.indel_min_reads2))
            fho.write('  indel min var freq: {}\n'.format(self.indel_min_var_freq))
            fho.write('  indel p value: {}\n'.format(self.indel_p_value))
            fho.write('  indel strand filter: {}\n'.format(self.indel_strand_filter))
            fho.write('  snp min avg qual: {}\n'.format(self.snp_min_avg_qual))
            fho.write('  snp min coverage: {}\n'.format(self.snp_min_coverage))
            fho.write('  snp min freq for hom: {}\n'.format(self.snp_min_freq_for_hom))
            fho.write('  snp min reads2: {}\n'.format(self.snp_min_reads2))
            fho.write('  snp min var freq: {}\n'.format(self.snp_min_var_freq))
            fho.write('  snp p value: {}\n'.format(self.snp_p_value))
            fho.write('  snp strand filter: {}\n'.format(self.snp_strand_filter))
            fho.write('  mpileup max depth: {}\n'.format(self.mpileup_max_depth))

            # Percent identity variance filtering.

            fho.write('pidv filter:\n')
            fho.write('  max pidv pct: {}\n'.format(self.max_pidv_pct))

            # Reference genomes used for mapping.

            fho.write('reference genomes:\n')
            for acc_id in self.dir_db.for_accids():
                fho.write('  {}: \'{}\'\n'.format(acc_id, self.dir_db.fna_file_path(acc_id)))

            # Sample key information.

            fho.write('sample key:\n')
            for i in self.sample_key.df.index:
                r1 = pathlib.Path(self.sample_key.df.loc[i]['fastq r1']).resolve().as_posix()
                r2 = pathlib.Path(self.sample_key.df.loc[i]['fastq r2']).resolve().as_posix()
                fho.write('  {}:\n'.format(self.sample_key.df.loc[i]['sample id']))
                fho.write('    fastq r1: \'{}\'\n'.format(r1))
                fho.write('    fastq r2: \'{}\'\n'.format(r2))
                fho.write('    reference id:\n')
                for ref in self.sample_key.df.loc[i]['reference id'].split():
                    fho.write('      - \'{}\'\n'.format(ref))

        if config_file.is_file():
            StdoutMessage(
                command='VarCallConfig',
                function='write_config_file',
                message='Wrote YAML configuration: {}'.format(config_file.resolve())
            )
        else:
            StdoutMessage(
                command='VarCallConfig',
                function='write_config_file',
                message='Could not write YAML configuration file.',
                fatal=True
            )

        self.config_file = config_file

        return

    def load_yaml_config(self, config_file):

        if config_file.is_file():
            StdoutMessage(
                command='VarCallConfig',
                function='load_yaml_config',
                message='Config YAML file exists.'
            )
        else:
            StdoutMessage(
                command='VarCallConfig',
                function='load_yaml_config',
                message='Could not find YAML configuration: {}'.format(config_file.resolve()),
                fatal=True
            )

        # Load the pre-existing YAML configuration file.

        with open(config_file.resolve(), 'r') as fhi:
            config_yaml = yaml.safe_load(fhi)

        StdoutMessage(
            command='VarCallConfig',
            function='load_yaml_config',
            message='Loaded the YAML configuration file.'
        )

        self.config_file = config_file
        self.config_file_yaml = config_yaml

        return


# __END__
