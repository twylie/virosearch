from .dir_db import DirDB
from subprocess import call
from ..support.stdout_message import StdoutMessage


###############################################################################
#                                   IndexRef                                  #
###############################################################################

# This class is responsible for indexing underlying sequence files for
# downstream processing. We will point the instance at a preexisting ViroSearch
# DB directory. If all of the required index files are present, we will just
# pass without doing anything. If the index files are missing, we will generate
# all of the required index files.


class IndexRef:

    def __init__(self, arguments, verbose=True):

        self.arguments = arguments
        self.verbose = verbose
        self.dir_db = DirDB(self.arguments.virosearch_db.resolve())

        self.__eval_dir_db()

        return

    def __eval_dir_db(self):

        # Does the top-level ViroSearch DB directory exist?

        StdoutMessage(
            command='IndexRef',
            function='eval_dir_db',
            message=f'Evaluating the ViroSearch DB directory.'
        )

        if self.dir_db.is_dir_path():
            StdoutMessage(
                command='IndexRef',
                function='eval_dir_db',
                message=f'Parent ViroSearch directory present.'
            )
        else:
            StdoutMessage(
                command='IndexRef',
                function='eval_dir_db',
                message=f'Parent ViroSearch directory does not exist.',
                fatal=True
            )

        # Minimally, the 'referencesDB' sub-directory should be present.

        StdoutMessage(
            command='IndexRef',
            function='eval_dir_db',
            message=f'Evaluating --virosearch-db: {self.arguments.virosearch_db}'
        )

        if self.dir_db.is_references_db():
            StdoutMessage(
                command='IndexRef',
                function='eval_dir_db',
                message=f'ViroSearch referencesDB directory present.'
            )
        else:
            StdoutMessage(
                command='IndexRef',
                function='eval_dir_db',
                message=f'ViroSearch referencesDB directory does not exist.',
                fatal=True
            )

        return

    def index_references(self):

        # Scan the preexisting ViroSearch DB directory. Assess if we need to
        # index the sequence files or not.

        # --- DEPENDENCIES ---
        #
        # bwa

        StdoutMessage(
            command='IndexRef',
            function='index_references',
            message=f'BWA MEM indexing of references.'
        )

        if len(self.dir_db.for_accids()) > 0:

            for accid in self.dir_db.for_accids():
                if self.dir_db.is_bwa_index(accid):
                    StdoutMessage(
                        command='IndexRef',
                        function='index_references',
                        message=f'BWA MEM indexes exist: skipping {accid}'
                    )
                else:
                    StdoutMessage(
                        command='IndexRef',
                        function='index_references',
                        message=f'BWA MEM indexes missing: {accid}'
                    )
                    refseq_fna = self.dir_db.fna_file_path(accid)
                    call(['bwa', 'index', refseq_fna])

        else:

            StdoutMessage(
                command='IndexRef',
                function='index_references',
                message=f'No reference genome files to index.',
                fatal=True
            )

        return


# __END__
