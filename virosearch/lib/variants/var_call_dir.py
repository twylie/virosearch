import os
import pathlib
from ..support.stdout_message import StdoutMessage
import re


###############################################################################
#                                  VarCallDir                                 #
###############################################################################

# This class is responsible for setting up a structured ViroSearch variant
# calling directory. This structure is used by the various Snakemake pipelines
# to read and write variant calling information.


class VarCallDir:

    def __init__(self, outdir, verbose=True):

        self.outdir = outdir
        self.verbose = verbose
        self.dirs = set()
        self.files = set()
        self.file_paths = set()
        self.config = None
        self.snakefile_path = None

        return

    def make_directory(self):

        # We're making a brand new instance of a variant calling processing
        # directory structure.

        if self.outdir.is_dir():
            StdoutMessage(
                command='VarCallDir',
                function='make_directory',
                message='The outdir already exists: {}'.format(self.outdir.resolve().as_posix()),
                fatal=True
            )
        else:
            StdoutMessage(
                command='VarCallDir',
                function='make_directory',
                message='Creating new outdir processing directory.'
            )
            self.outdir.mkdir()

        return

    def load_directory(self):

        # Recursively traverse the given ViroSearch output directory and collect
        # file paths of interest. As the directory structure is standardized, we
        # may use this class to evaluate and navigate the file structure of the
        # directory. A generic example of this tree follows:
        #
        # outdir
        # ├── config.yaml
        # ├── sample_key.tsv
        # └── varcallpe.smk

        StdoutMessage(
            command='VarCallDir',
            function='load_directory',
            message='Loading variant calling results directory.'
        )

        varcall_dir = pathlib.Path(self.outdir).resolve()

        for dir_path, dir_names, file_names in os.walk(varcall_dir):

            if dir_names:
                for dir_ in dir_names:
                    self.dirs.add(dir_)

            if file_names:

                for file_ in file_names:
                    file_path = pathlib.Path(dir_path, file_)
                    file_name = file_path.name
                    self.files.add(file_name)
                    self.file_paths.add(file_path)

        # We do some minimal checks to see if the loaded directory looks like a
        # variant calling results directory.

        if ('sample_key.tsv' in self.files) and ('config.yaml' in self.files):
            StdoutMessage(
                command='VarCallDir',
                function='load_directory',
                message='Passed results directory evaluation step.'
            )
        else:
            StdoutMessage(
                command='VarCallDir',
                function='load_directory',
                message='Not a results directory; missing sample_key.tsv and config.yaml files.',
                fatal=True
            )

        # We can load the paths to the existing files of interest if they
        # are present.

        for path_ in self.file_paths:
            file_path = path_.resolve().as_posix()
            if re.search('config.yaml$', file_path):
                self.config = path_
            elif re.search('varcallpe.smk$', file_path):
                self.snakefile_path = path_

        return


# __END__
