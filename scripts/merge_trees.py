'''
Merge multiple Newick trees into a single Newick tree leaving a polytomy at the root
'''

import argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Merge multiple Newick trees into a single Newick tree leaving a polytomy at the root",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--trees', nargs='+', type=str, required=True, help="list of Newick tree files")
    parser.add_argument('--outgroup', type=str, required=True, help="remove outgroup with name")
    parser.add_argument('--output', required=True, help="file to write merged Newick tree")
    args = parser.parse_args()

    newick_strings = []
    for filename in args.trees:
        with open(filename, 'r') as file:
            newick_string = file.read().replace('\n', '')
            newick_string = re.sub(r'NODE_\d+', '', newick_string)
            newick_string = re.sub(r':\d+\.\d+', '', newick_string)
            match = re.search(r'^\(' + re.escape(args.outgroup) + ',(.+)\);', newick_string)
            if match:
                newick_string = match.groups()[0]
            match = re.search(r'^\((.+),' + re.escape(args.outgroup) + '\);', newick_string)
            if match:
                newick_string = match.groups()[0]
            newick_strings.append(newick_string)

    joined_newick_string = "(" + ",".join(newick_strings) + ");"

    with open (args.output, "w") as file:
        file.write(joined_newick_string + '\n')
