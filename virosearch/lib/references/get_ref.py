from Bio import Entrez
import pathlib
import time
from ..support.stdout_message import StdoutMessage
from .dir_db import DirDB


###############################################################################
#                                    GetRef                                   #
###############################################################################

# This class is responsible for setting up a structured ViroSearch database
# directory and automatically downloading a list of reference sequence files,
# given a list of accession IDs. If the ViroSearch directory already exists, we
# will be appending new reference files to the directory structure.


class GetRef:

    def __init__(self, arguments, verbose=True):

        self.arguments = arguments
        self.verbose = verbose
        self.directory_mode = None
        self.dir_db = DirDB(self.arguments.virosearch_db.resolve())

        return

    def __download_genbank(self, accid, sleep=5):

        # Make the accession ID specific directory.

        accid_dir = pathlib.Path(self.db_dir.resolve(), f'{accid}')

        if not accid_dir.is_dir():
            StdoutMessage(
                command='GetRef',
                function='download_genbank',
                message=f'Making directory: {accid_dir}'
            )
            accid_dir.mkdir()

        # Download the indicated accession ID's GenBank file.

        genbank_file = pathlib.Path(self.db_dir.resolve(), f'{accid}/{accid}.gbk')

        if not genbank_file.is_file():

            StdoutMessage(
                command='GetRef',
                function='download_genbank',
                message=f'Downloading GenBank file: {accid}'
            )

            net_handle = Entrez.efetch(db="nucleotide", id=accid, rettype="gb", retmode="text")

            with open(genbank_file, 'w') as fho:
                fho.write(net_handle.read())

            net_handle.close()
            time.sleep(sleep)

            StdoutMessage(
                command='GetRef',
                function='download_genbank',
                message=f'Wrote GenBank file: {genbank_file}'
            )

        else:

            StdoutMessage(
                command='GetRef',
                function='download_genbank',
                message=f'GenBank file exists: skipping {accid}'
            )

        return

    def __download_cds_amino_acid(self, accid, sleep=5):

        # Make the accession ID specific directory.

        accid_dir = pathlib.Path(self.db_dir.resolve(), f'{accid}')

        if not accid_dir.is_dir():
            StdoutMessage(
                command='GetRef',
                function='download_cds_amino_acid',
                message=f'Making directory: {accid_dir}'
            )
            accid_dir.mkdir()

        # Download the indicated accession ID's CDS amino acid file.

        cds_faa_file = pathlib.Path(self.db_dir.resolve(), f'{accid}/{accid}.cds.faa')

        if not cds_faa_file.is_file():

            StdoutMessage(
                command='GetRef',
                function='download_cds_amino_acid',
                message=f'Downloading CDS amino acid file: {accid}'
            )

            net_handle = Entrez.efetch(db="nuccore", id=accid, rettype="fasta_cds_aa", retmode="text")

            with open(cds_faa_file, 'w') as fho:
                fho.write(net_handle.read())

            net_handle.close()
            time.sleep(sleep)

            StdoutMessage(
                command='GetRef',
                function='download_cds_amino_acid',
                message=f'Wrote cds.faa file: {cds_faa_file}'
            )

        else:

            StdoutMessage(
                command='GetRef',
                function='download_cds_amino_acid',
                message=f'The cds.faa exists: skipping {accid}'
            )

        return

    def __download_cds_nucleotide(self, accid, sleep=5):

        # Make the accession ID specific directory.

        accid_dir = pathlib.Path(self.db_dir.resolve(), f'{accid}')

        if not accid_dir.is_dir():
            StdoutMessage(
                command='GetRef',
                function='download_cds_nucleotide',
                message=f'Making directory: {accid_dir}'
            )
            accid_dir.mkdir()

        # Download the indicated accession ID's CDS nucleotide file.

        cds_fna_file = pathlib.Path(self.db_dir.resolve(), f'{accid}/{accid}.cds.fna')

        if not cds_fna_file.is_file():

            StdoutMessage(
                command='GetRef',
                function='download_cds_nucleotide',
                message=f'Downloading CDS nucleotide file: {accid}'
            )

            net_handle = Entrez.efetch(db="nuccore", id=accid, rettype="fasta_cds_na", retmode="text")

            with open(cds_fna_file, 'w') as fho:
                fho.write(net_handle.read())

            net_handle.close()
            time.sleep(sleep)

            StdoutMessage(
                command='GetRef',
                function='download_cds_nucleotide',
                message=f'Wrote cds.fna file: {cds_fna_file}'
            )

        else:

            StdoutMessage(
                command='GetRef',
                function='download_cds_nucleotide',
                message=f'The cds.fna exists: skipping {accid}'
            )

        return

    def __download_refseq(self, accid, sleep=5):

        # Make the accession ID specific directory.

        accid_dir = pathlib.Path(self.db_dir.resolve(), f'{accid}')

        if not accid_dir.is_dir():
            StdoutMessage(
                command='GetRef',
                function='download_refseq',
                message=f'Making directory: {accid}'
            )
            accid_dir.mkdir()

        # Download the indicated accession ID's RefSeq nucleotide file.

        refseq_fna_file = pathlib.Path(self.db_dir.resolve(), f'{accid}/{accid}.fna')

        if not refseq_fna_file.is_file():

            StdoutMessage(
                command='GetRef',
                function='download_refseq',
                message=f'Downloading RefSeq nucleotide file: {accid}'
            )

            net_handle = Entrez.efetch(db="nuccore", id=accid, rettype="fasta", retmode="text")

            with open(refseq_fna_file, 'w') as fho:
                fho.write(net_handle.read())
            net_handle.close()
            time.sleep(sleep)

            StdoutMessage(
                command='GetRef',
                function='download_refseq',
                message=f'Wrote RefSeq file: {refseq_fna_file}'
            )

        else:

            StdoutMessage(
                command='GetRef',
                function='download_refseq',
                message=f'The RefSeq fna exists: skipping {accid}'
            )

        return

    def make_db_dir(self):

        # Evaluate if the working directory already exists. If so, we will be
        # appending; if not, we will create the ViroSearch directory from
        # scratch. Current directory structure:
        #
        # path/
        # path/referencesDB/

        StdoutMessage(
            command='GetRef',
            function='make_db_dir',
            message='Assessing database directory.'
        )

        if not self.dir_db.dir_path().is_dir():

            self.directory_mode = 'new'
            self.arguments.virosearch_db.mkdir()

            StdoutMessage(
                command='GetRef',
                function='make_db_dir',
                message=f'Making directory: {self.arguments.virosearch_db.resolve()}'
            )

            db_dir = pathlib.Path(self.arguments.virosearch_db.resolve(), 'referencesDB')
            db_dir.mkdir()
            self.db_dir = db_dir

            StdoutMessage(
                command='GetRef',
                function='make_db_dir',
                message=f'Making directory: {db_dir}'
            )

        else:

            self.directory_mode = 'appending'

            StdoutMessage(
                command='GetRef',
                function='make_db_dir',
                message=f'Evaluating: {self.arguments.virosearch_db.resolve()}'
            )

            db_dir = pathlib.Path(self.arguments.virosearch_db.resolve(), 'referencesDB')
            self.db_dir = db_dir

            StdoutMessage(
                command='GetRef',
                function='make_db_dir',
                message=f'Appending: {db_dir}'
            )

        return

    def download_entrez_references(self):

        # Based on the list of accession IDs passed from the CLI, we will
        # automatically download the required reference sequence files from
        # NCBI GenBank. It is important to put a small wait-time in the
        # download, so NCBI doesn't get upset about overloading their servers.
        # If we are in "appending" mode, the new databases will be downloaded
        # as needed; else, we simply skip over preexisting files.

        # Currently, the desired files per accession ID are:
        #
        # accid.gbk (GenBank Entry)
        # accid.fna (RefSeq nucleotide sequence)
        # accid.cds.fna (CDS nucleotide sequence)
        # accid.cds.faa (CDS amino acid sequence)

        # Set the Entrez email for the handshake.

        Entrez.email = self.arguments.email

        # Download the reference files.

        StdoutMessage(
            command='GetRef',
            function='download_entrez_references',
            message=f'Downloading data from Entrez.'
        )

        for accid in self.arguments.acc_id:
            self.__download_genbank(accid=accid)
            self.__download_cds_amino_acid(accid=accid)
            self.__download_cds_nucleotide(accid=accid)
            self.__download_refseq(accid=accid)

        return


# __END__
