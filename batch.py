"""
This script runs the entire suite of seasonal flu builds via AWS batch
One batch run is created per lineage
"""

import subprocess
import argparse

def get_cpus(jobs):
    count = 1
    if jobs >= 72:
        count = 72
    elif jobs >= 36:
        count = 36
    elif jobs >= 16:
        count = 16
    elif jobs >= 8:
        count = 8
    elif jobs >= 4:
        count = 4
    elif jobs >= 2:
        count = 2
    return count

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run flu builds on aws', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--system', type = str, default = 'local', help='where to run, local or batch')
    parser.add_argument('-l', '--lineages', nargs='+', type = str,  help ="flu lineages to include", default=['h3n2', 'h1n1pdm', 'vic', 'yam'])
    parser.add_argument('-r', '--resolutions', nargs='+', type = str,  help ="flu resolutions to include", default=['1y'])
    parser.add_argument('-s', '--segments', nargs='+', type = str, help ="flu segments to include", default=['genome', 'ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns'])
    params = parser.parse_args()

    targets = []
    for lineage in params.lineages:
        for resolution in params.resolutions:
            for segment in params.segments:
                targets.append('auspice/seattleflu_flu_seasonal_%s_%s_%s.json'%(lineage, segment, resolution))

    cpus = get_cpus(len(targets))
    memory = 1800 * cpus
    if params.system == 'local':
        call = ['nextstrain', 'build', '.', '--jobs', '1']
    elif params.system == 'batch':
        call = ['nextstrain', 'build', '--aws-batch', '--detach', '--aws-batch-cpus', str(cpus), '--aws-batch-memory', str(memory), '.', '--jobs', str(cpus)]
    call.extend(targets)

    if targets:
        print(' '.join(call))
        pro = subprocess.call(call)
