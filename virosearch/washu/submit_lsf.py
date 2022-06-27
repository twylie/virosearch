#! /usr/bin/python3.7

# This script may be used by Snakemake to submit sub-processing to the WashU
# LSF system in parallel. We will be using Docker containers to run these
# processes. The master Snakemake job will poll to see when the sub-processes
# are finished and also keep track of jobs based on the original DAG.

###############################################################################
#                                DO NOT EDIT!!!                               #
###############################################################################

import os
import sys
import yaml
from snakemake.utils import read_job_properties

with open('washu.yaml', 'r') as fh:
    config = yaml.load(fh, Loader=yaml.FullLoader)

jobscript = sys.argv[-1]
props = read_job_properties(jobscript)
job_id = props['jobid']

volumes = list()
for volume in config['docker']['volumes']:
    volume += ':' + volume
    volumes.append(volume)
volumes = str(' '.join(volumes))

lsf_error = os.path.join(config['lsf']['lsf log dir'], 'LSF.err.' + str(job_id))
lsf_out = os.path.join(config['lsf']['lsf log dir'], 'LSF.out.' + str(job_id))
lsf_cmd = os.path.join(config['lsf']['lsf log dir'], 'LSF.cmd.' + str(job_id))

cmd = ' '.join([
    'LSF_DOCKER_VOLUMES="{}"'.format(volumes),
    'bsub',
    '-M {}'.format(config['lsf']['memory']),
    '-R "select[mem>{}] rusage[mem={}]"'.format(config['lsf']['memory'], config['lsf']['memory']),
    '-G {}'.format(config['lsf']['compute group']),
    '-q {}'.format(config['lsf']['queue']),
    '-e {}'.format(lsf_error),
    '-o {}'.format(lsf_out),
    '-a \'docker({})\''.format(config['docker']['image']),
    jobscript
])

with open(lsf_cmd, 'w') as fho:
    fho.write(cmd + '\n')

os.system(cmd)

# __END__

# T.N. Wylie  <twylie@wustl.edu>
# Last Update: Tue Jun 21 11:08:51 CDT 2022
