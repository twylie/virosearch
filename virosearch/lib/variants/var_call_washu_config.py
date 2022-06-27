from ..support.stdout_message import StdoutMessage
import yaml


###############################################################################
#                              VarCallWashuConfig                             #
###############################################################################

# --- Only for Washington University School of Medicine users. ---

# This class will provide support for formatting and writing the WashU
# configuration file used for LSF parallel processing at WashU.

class VarCallWashuConfig:

    def __init__(self, verbose=True):

        self.verbose = verbose
        self.config_file = None
        self.config_file_yaml = None

        return

    def populate_arguments(self, var_call_washu):

        # EXAMPLE:
        #
        # 'command': 'washupp',
        # 'results_dir': PosixPath('results'),
        # 'docker_volumes': ['foo', 'foo'],
        # 'lsf_compute_group': 'twylie-group',
        # 'docker_image': 'twylie/virosearch:latest',
        # 'lsf_memory': '16G',
        # 'lsf_queue': 'general',
        # 'lsf_latency_wait': 15,
        # 'lsf_restart_times': 3,
        # 'lsf_log_dir': PosixPath('lsf_logs')},

        self.results_dir = var_call_washu.arguments.results_dir
        self.lsf_log_dir = var_call_washu.lsf_log_dir
        self.lsf_latency_wait = var_call_washu.arguments.lsf_latency_wait
        self.lsf_restart_times = var_call_washu.arguments.lsf_restart_times
        self.lsf_queue = var_call_washu.arguments.lsf_queue
        self.lsf_memory = var_call_washu.arguments.lsf_memory
        self.docker_image = var_call_washu.arguments.docker_image
        self.lsf_compute_group = var_call_washu.arguments.lsf_compute_group
        self.command = var_call_washu.arguments.command
        self.docker_volumes = var_call_washu.arguments.docker_volumes

        StdoutMessage(
            command='VarCallWashuConfig',
            function='populate_arguments',
            message='WashU config populated with CLI arguments.'
        )

        return

    def write_config_file(self, config_file):

        # We will format and write the config file as simple YAML. Snakemake
        # will load the configuration each time a child process is created by
        # the submit_lsf.py script.

        if config_file.is_file():
            StdoutMessage(
                command='VarCallWashuConfig',
                function='write_config_file',
                message='Washu configuration file already exists.',
                fatal=True
            )

        with open(config_file, 'w') as fho:

            fho.write('docker:\n')
            fho.write('  image: \'{}\'\n'.format(self.docker_image))
            fho.write('  volumes:\n')
            for volume in self.docker_volumes:
                fho.write('    - \'{}\'\n'.format(volume))
            fho.write('lsf:\n')
            fho.write('  memory: \'{}\'\n'.format(self.lsf_memory))
            fho.write('  results dir: \'{}\'\n'.format(self.results_dir.resolve().as_posix()))
            fho.write('  compute group: \'{}\'\n'.format(self.lsf_compute_group))
            fho.write('  queue: \'{}\'\n'.format(self.lsf_queue))
            fho.write('  latency wait: {}\n'.format(self.lsf_latency_wait))
            fho.write('  restart times: {}\n'.format(self.lsf_restart_times))
            fho.write('  lsf log dir: \'{}\'\n'.format(self.lsf_log_dir.resolve().as_posix()))

            if config_file.is_file():
                StdoutMessage(
                    command='VarCallWashuConfig',
                    function='write_config_file',
                    message='Wrote WashU YAML configuration: {}'.format(config_file.resolve())
                )
            else:
                StdoutMessage(
                    command='VarCallWashuConfig',
                    function='write_config_file',
                    message='Could not write WashU YAML configuration file.',
                    fatal=True
                )

            self.config_file = config_file

        return

    def load_yaml_config(self, config_file):

        if config_file.is_file():
            StdoutMessage(
                command='VarCallWashuConfig',
                function='load_yaml_config',
                message='Config WashU YAML file exists.'
            )
        else:
            StdoutMessage(
                command='VarCallWashuConfig',
                function='load_yaml_config',
                message='Could not find WashU YAML configuration: {}'.format(config_file.resolve()),
                fatal=True
            )

        # Load the pre-existing YAML configuration file.

        with open(config_file.resolve(), 'r') as fhi:
            config_yaml = yaml.safe_load(fhi)

        StdoutMessage(
            command='VarCallWashuConfig',
            function='load_yaml_config',
            message='Loaded the WashU YAML configuration file.'
        )

        self.config_file = config_file
        self.config_file_yaml = config_yaml

        return


# __END__
