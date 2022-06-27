import os
import pathlib
from ..support.stdout_message import StdoutMessage


###############################################################################
#                                VarCallSmkDirPE                              #
###############################################################################

# This class is responsible for setting up a structured Snakemake
# processing/output directory for paired-end variant calling. By default, the
# output directory is called 'smk' within the parent results directory;
# however, a user may pass a different path for writing the Snakemake
# processing output files.


class VarCallSmkDirPE:

    def __init__(self, smk_dir, results_dir, verbose=True):

        self.smk_dir = smk_dir
        self.results_dir = results_dir
        self.verbose = verbose
        self.dirs = set()
        self.files = set()

        return

    def make_directory(self):

        # By default, the Snakemake processing directory will be under the
        # parent results directory, unless there is an argument passed to
        # change this behavior.

        if not self.smk_dir:
            self.smk_dir = pathlib.Path(self.results_dir, 'smk')

        # We're making a brand new instance of a Snakemake processing directory
        # structure.

        if self.smk_dir.is_dir():
            StdoutMessage(
                command='VarCallSmkDirPE',
                function='make_directory',
                message='The Snakemake directory already exists: {}'.format(self.smk_dir.resolve().as_posix()),
                fatal=True
            )
        else:
            StdoutMessage(
                command='VarCallSmkDirPE',
                function='make_directory',
                message='Creating new Snakemake processing directory: {}'.format(self.smk_dir.resolve().as_posix())
            )
            self.smk_dir.mkdir()

        return

    def smk_dir(self):
        return self.smk_dir

    def load_directory(self):

        # Recursively traverse the given Snakemake output directory and collect
        # file paths of interest. As the directory structure is standardized, we
        # may use this class to evaluate and navigate the file structure of the
        # directory.

        StdoutMessage(
            command='VarCallSmkDirPE',
            function='load_directory',
            message='Loading Snakemake processing directory.'
        )

        smk = pathlib.Path(self.smk_dir).resolve()

        for dir_path, dir_names, file_names in os.walk(smk):

            if dir_names:
                for dir_ in dir_names:
                    self.dirs.add(dir_)

            if file_names:

                for file_ in file_names:
                    file_path = pathlib.Path(dir_path, file_)
                    file_name = file_path.name
                    self.files.add(file_name)

        return


# __END__
