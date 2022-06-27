import pathlib
import shutil
from ..support.stdout_message import StdoutMessage
from virosearch.variants.var_call_washu_config import VarCallWashuConfig


###############################################################################
#                                 VarCallWashu                                #
###############################################################################

# --- Only for Washington University School of Medicine users. ---

# This class is responsible for helping extend a ViroSearch results directory
# with LSF parallel processing support for WashU users.


class VarCallWashu:

    def __init__(self, arguments, verbose=True):

        self.verbose = verbose
        self.arguments = arguments
        self.snakemake_cmd_path = None
        self.submit_lsf_path = None
        self.lsf_log_dir = None

        return

    def make_lsf_logs_directory(self):

        StdoutMessage(
            command='VarCallWashu',
            function='make_lsf_logs_directory',
            message='Evaluating LSF log directory settings.'
        )

        if self.arguments.lsf_log_dir.as_posix() == 'lsf_logs':

            StdoutMessage(
                command='VarCallWashu',
                function='make_lsf_logs_directory',
                message='Attempting default LSF log directory within results directory.'
            )

            lsf_log_dir = pathlib.Path(self.arguments.results_dir, 'lsf_logs')

            if lsf_log_dir.is_dir():
                StdoutMessage(
                    command='VarCallWashu',
                    function='make_lsf_logs_directory',
                    message='The default lsf_logs directory already exists.',
                    fatal=True
                )
            else:
                StdoutMessage(
                    command='VarCallWashu',
                    function='make_lsf_logs_directory',
                    message='Writing default LSF log directory. [lsf_logs]'
                )

                lsf_log_dir.mkdir()

                self.lsf_log_dir = lsf_log_dir

        else:

            if self.arguments.lsf_log_dir.is_dir():
                StdoutMessage(
                    command='VarCallWashu',
                    function='make_lsf_logs_directory',
                    message='The external LSF logs directory already exists.',
                    fatal=True
                )
            else:
                StdoutMessage(
                    command='VarCallWashu',
                    function='make_lsf_logs_directory',
                    message='Writing external LSF log directory: {}'.format(self.arguments.lsf_log_dir.resolve().as_posix())
                )

                self.arguments.lsf_log_dir.mkdir()

                self.lsf_log_dir = self.arguments.lsf_log_dir

        return

    def copy_submit_lsf_script(self):

        StdoutMessage(
            command='VarCallWashu',
            function='copy_submit_lsf_script',
            message='Copying LSF submission script to results directory.'
        )

        shutil.copy('/usr/lib/python3.7/virosearch/washu/submit_lsf.py', self.arguments.results_dir)

        self.submit_lsf_path = pathlib.Path(self.arguments.results_dir, 'submit_lsf.py')

        return

    def write_yaml_config_file(self):

        StdoutMessage(
            command='VarCallWashu',
            function='write_yaml_config_file',
            message='Writing WashU YAML configuration file.'
        )

        config = pathlib.Path(self.arguments.results_dir, 'washu.yaml')
        var_call_washu_config = VarCallWashuConfig()
        var_call_washu_config.populate_arguments(self)
        var_call_washu_config.write_config_file(config)

        return

    def write_example_smk_cmd(self, var_call_dir):

        # We can provide an example Snakemake execution shell command based on
        # the configuration settings. This version of the snakemake command is
        # for parallel processing at Washington University's RIS. Settings will
        # be for LSF commands.

        # EXAMPLE:

        # snakemake \
        # --snakefile varcallpe.smk \
        # --configfile config.yaml \
        # --cluster submit_lsf.py \
        # --cores 100 \
        # --local-cores 1 \
        # --restart-times 3 \
        # --latency-wait 10 \
        # -p

        StdoutMessage(
            command='VarCallWashu',
            function='write_example_smk_cmd',
            message='Write an example WashU Snakemake pipeline execution command. [washu.sh]'
        )

        cmd = ' '.join([
            'snakemake',
            '--snakefile {}'.format(var_call_dir.snakefile_path.resolve().as_posix()),
            '--configfile {}'.format(var_call_dir.config.resolve().as_posix()),
            '--cluster {}'.format(self.submit_lsf_path.resolve().as_posix()),
            '--jobs 100',
            '--local-cores 1',
            '--restart-times 3',
            '--latency-wait 10',
            '-p'
        ])

        cmd = '#! /bin/sh\n' + cmd

        self.snakemake_cmd_path = pathlib.Path(self.arguments.results_dir, 'washu.sh')

        with open(self.snakemake_cmd_path, 'w') as fho:
            fho.write(cmd + '\n')

        StdoutMessage(
            command='VarCallWashu',
            function='write_example_smk_cmd',
            message='Making the WashU example command executable.'
        )

        self.snakemake_cmd_path.chmod(0o775)

        return


# __END__
