from ..references.dir_db import DirDB
from ..support.stdout_message import StdoutMessage
from .var_call_config import VarCallConfig
from .var_call_dir import VarCallDir
from .var_call_smk_dir_pe import VarCallSmkDirPE
from .var_call_sample_key import VarCallSampleKey
import pathlib
import shutil


###############################################################################
#                                  VarCallPE                                  #
###############################################################################

# This class will perform mapping and variant calling using Snakemake. We will
# copy in and execute the appropriate Snakefile recipe provided all of the
# required information from the CLI. We will point the instance at a
# preexisting ViroSearch DB directory which should have had appropriate
# reference genomes formatted and indexed.


class VarCallPE:

    def __init__(self, arguments, verbose=True):

        self.arguments = arguments
        self.verbose = verbose
        self.var_call_dir = None
        self.dir_db = DirDB(self.arguments.virosearch_db)
        self.samplekey = None
        self.config = None
        self.snakefile_path = None
        self.snakemake_cmd_path = None

        return

    def eval_dir_db(self):

        # Does the top-level ViroSearch DB directory exist?

        StdoutMessage(
            command='VarCallPE',
            function='eval_dir_db',
            message='Evaluating the ViroSearch DB directory.'
        )

        if self.dir_db.is_dir_path():
            StdoutMessage(
                command='VarCallPE',
                function='eval_dir_db',
                message='Parent ViroSearch DB directory present.'
            )
        else:
            StdoutMessage(
                command='VarCallPE',
                function='eval_dir_db',
                message='Parent ViroSearch DB directory does not exist.',
                fatal=True
            )

        # Minimally, the 'referencesDB' sub-directory should be present.

        StdoutMessage(
            command='VarCallPE',
            function='eval_dir_db',
            message=f'Evaluating --virosearch-db: {self.arguments.virosearch_db}'
        )

        if self.dir_db.is_references_db():
            StdoutMessage(
                command='VarCallPE',
                function='eval_dir_db',
                message='ViroSearch referencesDB directory present.'
            )
        else:
            StdoutMessage(
                command='VarCallPE',
                function='eval_dir_db',
                message='ViroSearch referencesDB directory does not exist.',
                fatal=True
            )

        return

    def eval_reference_indexes(self):

        # Scan the preexisting ViroSearch DB directory. Assess if the
        # references are indexed or not.

        StdoutMessage(
            command='VarCallPE',
            function='eval_reference_indexes',
            message=f'Review BWA MEM indexing of references.'
        )

        if len(self.dir_db.for_accids()) > 0:

            for accid in self.dir_db.for_accids():
                if self.dir_db.is_bwa_index(accid):
                    StdoutMessage(
                        command='VarCallPE',
                        function='eval_reference_indexes',
                        message=f'BWA MEM indexes exist: {accid}'
                    )
                else:
                    StdoutMessage(
                        command='VarCallPE',
                        function='eval_reference_indexes',
                        message=f'BWA MEM indexes missing: {accid}',
                        fatal=True
                    )

        else:

            StdoutMessage(
                command='VarCallPE',
                function='eval_reference_indexes',
                message=f'No reference genome files.',
                fatal=True
            )

        return

    def make_processing_outdir(self):

        StdoutMessage(
            command='VarCallPE',
            function='make_processing_outdir',
            message='Making the processing outdir directory.'
        )

        self.var_call_dir = VarCallDir(self.arguments.results_dir)
        self.var_call_dir.make_directory()

        return

    def make_smk_dir(self):

        StdoutMessage(
            command='VarCallPE',
            function='make_smk_dir',
            message='Making the Snakemale processing directory.'
        )

        self.smk_dir = VarCallSmkDirPE(self.arguments.smk_dir, self.arguments.results_dir)
        self.smk_dir.make_directory()
        self.arguments.smk_dir = self.smk_dir.smk_dir.resolve().as_posix()

        return

    def write_config_file(self):

        StdoutMessage(
            command='VarCallPE',
            function='write_config_file',
            message='Writing the YAML configuration file.'
        )

        self.config = VarCallConfig()
        self.config.populate_arguments(self.arguments, self.dir_db, self.samplekey)
        config_file = pathlib.Path(self.arguments.results_dir, 'config.yaml')
        self.config.write_config_file(config_file)

        return

    def copy_varcallpe_recipe(self):

        StdoutMessage(
            command='VarCallPE',
            function='copy_varcallpe_recipe',
            message='Copying the Snakefile recipe for pipeline execution.'
        )

        smk_recipe = '/usr/lib/python3.7/virosearch/recipes/varcallpe.smk'
        shutil.copy(smk_recipe, self.arguments.results_dir)
        self.snakefile_path = pathlib.Path(self.arguments.results_dir, 'varcallpe.smk')

        return

    def eval_samplekey(self):

        StdoutMessage(
            command='VarCallPE',
            function='eval_samplekey',
            message='Evaluating the provided samplekey.'
        )

        if self.arguments.sample_key.is_file():

            StdoutMessage(
                command='VarCallPE',
                function='eval_samplekey',
                message='Samplekey file is present.'
            )

            self.samplekey = VarCallSampleKey(self.arguments.sample_key)
            self.samplekey.populate_samples()
            self.samplekey.eval_reference_ids(self.dir_db)
            self.samplekey.eval_pe_samplekey_format()
            tsv_file = pathlib.Path(self.arguments.results_dir, 'sample_key.tsv')
            self.samplekey.write_samplekey_tsv_file(tsv_file)

        else:

            StdoutMessage(
                command='VarCallPE',
                function='eval_samplekey',
                message='Samplekey does not exist: {}'.format(self.arguments.sample_key.resolve().as_posix()),
                fatal=True
            )

        return

    def write_example_smk_cmd(self):

        # We can provide an example Snakemake execution shell command based on
        # the configuration settings. This version of the snakemake command is
        # for non-parallel processing.

        # EXAMPLE:

        # snakemake \
        # --snakefile  varcallpe.smk \
        # --configfile config.yaml \
        # --cores=5 \
        # -p

        StdoutMessage(
            command='VarCallPE',
            function='write_example_smk_cmd',
            message='Write an example Snakemake pipeline execution command. [cmd.sh]'
        )

        cmd = ' '.join([
            'snakemake',
            '--snakefile {}'.format(self.snakefile_path.resolve().as_posix()),
            '--configfile {}'.format(self.config.config_file.resolve().as_posix()),
            '--cores=1',
            '-p'
        ])

        cmd = '#! /bin/sh\n' + cmd

        self.snakemake_cmd_path = pathlib.Path(self.arguments.results_dir, 'cmd.sh')

        with open(self.snakemake_cmd_path, 'w') as fho:
            fho.write(cmd + '\n')

        StdoutMessage(
            command='VarCallPE',
            function='write_example_smk_cmd',
            message='Making the example command executable.'
        )

        self.snakemake_cmd_path.chmod(0o775)

        return


# __END__
