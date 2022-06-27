from ..support.stdout_message import StdoutMessage
import pandas as pd


###############################################################################
#                               VarCallSampleKey                              #
###############################################################################

# This class will populate an object that maps samplekey information used for
# running the variant calling pipelines. This class handles both single-end
# and paired-end samplekey inputs.


class VarCallSampleKey:

    def __init__(self, samplekey, verbose=True):

        self.samplekey = samplekey
        self.verbose = verbose
        self.df = None
        self.samplekey_tsv = None

        return

    def populate_samples(self):

        # NOTE: Entries are tab delimited, while the references field is space
        # delimited. The sample id field is free-form; however, the references
        # field must match the accession id format as provided by ViroSearch
        # the reference directory. Depending on the end mode (SE or PE), we may
        # expect one or two associated FASTQ files.
        #
        # TestSampleOne  HPV-6 HPV-18  test.r1.fastq test.r2.fastq
        # TestSampleTwo  HPV-18  test.r1.fastq test.r2.fastq

        self.df = pd.read_csv(
            self.samplekey,
            names=['sample id', 'reference id', 'fastq r1', 'fastq r2'],
            sep='\t'
        )

        if len(self.df['sample id']) != len(self.df['sample id'].unique()):

            StdoutMessage(
                command='VarCallSampleKey',
                function='populate_samples',
                message='Provided sample ids are not unique',
                fatal=True
            )

        else:

            StdoutMessage(
                command='VarCallSampleKey',
                function='populate_samples',
                message='Provided sample ids are uniqued.'
            )

        StdoutMessage(
            command='VarCallSampleKey',
            function='populate_samples',
            message='Samplekey information populated.'
        )

        return

    def eval_reference_ids(self, dir_db):

        # Reference ids in the samplekey must match those in the provided
        # ViroSearch referencesDB.

        references = set()
        for reference in self.df['reference id']:
            for i in reference.split(' '):
                references.add(i)

        samplekey_ref_set = set(dir_db.references.keys())

        # TODO: Maybe don't jump directly into this class---use method?

        for samplekey_ref in references:
            if samplekey_ref not in samplekey_ref_set:
                StdoutMessage(
                    command='VarCallSampleKey',
                    function='eval_reference_ids',
                    message='Reference id {} not in {}'.format(
                        samplekey_ref,
                        dir_db.dir_.resolve().as_posix()
                    ),
                    fatal=True
                )
            else:
                StdoutMessage(
                    command='VarCallSampleKey',
                    function='eval_reference_ids',
                    message='Reference id {} is a match.'.format(samplekey_ref)
                )

        return

    def eval_pe_samplekey_format(self):

        # Paired-end format expects 4 columns with 2 columns reserved for FASTQ
        # R1 and R2 file paths. Fail anything that doesn't match this for PE
        # processing.

        column_count = len(self.df.columns)

        if column_count != 4:
            StdoutMessage(
                command='VarCallSampleKey',
                function='eval_pe_samplekey_format',
                message='Expected 4 columns but found {}.'.format(column_count),
                fatal=True
            )
        elif self.df['sample id'].isna().any():
            StdoutMessage(
                command='VarCallSampleKey',
                function='eval_pe_samplekey_format',
                message='Column 1 (sample id) has null values present.',
                fatal=True
            )
        elif self.df['reference id'].isna().any():
            StdoutMessage(
                command='VarCallSampleKey',
                function='eval_pe_samplekey_format',
                message='Column 2 (reference id) has null values present.',
                fatal=True
            )
        elif self.df['fastq r1'].isna().any():
            StdoutMessage(
                command='VarCallSampleKey',
                function='eval_pe_samplekey_format',
                message='Column 3 (fastq r1) has null values present.',
                fatal=True
            )
        elif self.df['fastq r2'].isna().any():
            StdoutMessage(
                command='VarCallSampleKey',
                function='eval_pe_samplekey_format',
                message='Column 4 (fastq r2) has null values present.',
                fatal=True
            )
        else:
            StdoutMessage(
                command='VarCallSampleKey',
                function='eval_pe_samplekey_format',
                message='Passed samplekey format evaluation.'
            )

        # Neither the sample id nor the reference id names may contain a dunder
        # (double-underscore, '__'), as we use that delimiter within the
        # pipeline for sample-to-reference mapping.

        if (self.df['sample id'].str.contains('__').any()) or (self.df['reference id'].str.contains('__').any()):
            StdoutMessage(
                command='VarCallSampleKey',
                function='eval_pe_samplekey_format',
                message='An illegal double-underscore is present in the sample or reference ids.',
                fatal=True
            )

        return

    def write_samplekey_tsv_file(self, tsv_file):

        self.df.to_csv(tsv_file.resolve().as_posix(), sep='\t')
        self.samplekey_tsv = tsv_file

        StdoutMessage(
            command='VarCallSampleKey',
            function='write_samplekey_tsv_file',
            message='Wrote samplekey to processing directory.'
        )

        return


# __END__
