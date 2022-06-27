#! /usr/bin/python3.7

import argparse
import pysam
from viromatch.lib.hits import CalcPercentVariationSam

version = '1.1'


def eval_cli_arguments():

    parser = argparse.ArgumentParser(
        description='Given mapped reads in a SAM file, filter (paired) hits based on percent identity variance.',
        prog='filter_pidv_pe.py',
        add_help=False
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

    # Required arguments.

    required_group = parser.add_argument_group('required')

    required_group.add_argument(
        '--sam',
        metavar='FILE',
        action='store',
        help='Path to input SAM mapped read file to be filtered.',
        required=True
    )

    required_group.add_argument(
        '--log',
        metavar='FILE',
        action='store',
        help='Path to log file listing filtered read ids and PID info.',
        required=True
    )

    required_group.add_argument(
        '--out',
        metavar='FILE',
        action='store',
        help='Path to write filtered SAM file.',
        required=True
    )

    required_group.add_argument(
        '--maxscore',
        metavar='FLOAT',
        action='store',
        default=0.15,
        help='Maximum variation score for filtering. [0.15]'
    )

    return parser.parse_args()


def filter_percent_identity(arguments):

    # First pass through the SAM file determines which read IDs have at least
    # one failure based on pidv. The second pass re-evaluates reads given that
    # one of the read pairs has a known failure; SAM/BAM flags will be updated
    # to reflect the changes and a VS:f tag will be added for pidv values.

    failed_read_ids = set()

    # First Pass ##############################################################

    sam = pysam.AlignmentFile(arguments.sam, 'r')

    for line in sam:

        if line.is_unmapped is False:

            read_id = line.query_name
            variation_score = CalcPercentVariationSam(line.to_string()).calc_percent_variation()

            if variation_score > float(arguments.maxscore):
                failed_read_ids.add(read_id)

    # Second Pass #############################################################

    fh_log = open(arguments.log, 'w')
    fh_out = open(arguments.out, 'w')

    sam = pysam.AlignmentFile(arguments.sam, 'r')

    fh_out.write(str(sam.header))

    fh_log.write('\t'.join(['read id', 'pe suffix', 'pid variance', 'alteration']) + '\n')

    for line in sam:

        read_id = line.query_name
        variation_score = CalcPercentVariationSam(line.to_string()).calc_percent_variation()

        if read_id in failed_read_ids:

            if line.is_read1 is True:
                pair_end_suffix = 'R1'
            elif line.is_read2 is True:
                pair_end_suffix = 'R2'
            else:
                pair_end_suffix = 'unset'

            if variation_score > float(arguments.maxscore):

                if line.mate_is_unmapped is False:
                    line.flag = line.flag + 8  # set as singleton
                line.flag = line.flag + 4  # set as unmapped
                line.set_tag('VS', variation_score)
                fh_out.write(line.to_string() + '\n')
                fh_log.write(f'{read_id}\t{pair_end_suffix}\t{variation_score}\tfail pidv\n')

            else:

                if line.mate_is_unmapped is False:
                    line.flag = line.flag + 8  # set as singleton
                line.set_tag('VS', variation_score)
                fh_out.write(line.to_string() + '\n')
                fh_log.write(f'{read_id}\t{pair_end_suffix}\t{variation_score}\tpass pidv\n')

        else:

            line.set_tag('VS', variation_score)
            fh_out.write(line.to_string() + '\n')

    fh_log.close()
    fh_out.close()

    return


if __name__ == '__main__':

    # Given an input SAM file of mapped reads, we will be calculating percent
    # identity vrianc (PIDV) score---based on CIGAR string and not a strict
    # Levenshtein distance calculation---and filtering based on a given
    # maxscore threshold. Output will include a log file of filtered read ids
    # and associated PID scores and a filtered SAM file with altered SAM flags
    # to indicate failed/unmapped reads. Reads that fail the PIDV evaluation
    # will have a VS:f tag indicating the pidv value for the read. Of note,
    # this only makes sense for MAPPED SAM hits, as the NM (and other mapped)
    # fields are required to calculate a score between a query hit and the
    # reference subject sequence.

    arguments = eval_cli_arguments()
    filter_percent_identity(arguments)


# __END__

# T.N. Wylie  <twylie@wustl.edu>
# Last Update: Mon Jun  6 15:44:46 CDT 2022
