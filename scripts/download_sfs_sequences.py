import requests
import argparse
import json
import os
from urllib.parse import urljoin


def get_full_lineage(short_lineage):
    """
    Given a *short_lineage*, return the full lineage required to find exact
    lineage match within ID3C.
    """
    lineage_map = {
        'h1n1pdm': 'Influenza.A.H1N1',
        'h3n2': 'Influenza.A.H3N2',
        'vic': 'Influenza.B.Vic',
        'yam': 'Influenza.B.Yam'
    }

    return lineage_map[short_lineage]


def generate_full_url(base_url, lineage, segment):
    """
    Generate the full URL for the API endpoint to get sequences of a specific
    *lineage* and *segment*
    """
    params = "/".join([lineage, segment])
    return urljoin(base_url, params)


def get_sequences_from_id3c(url, username, password, lineage, segment, output):
    """
    GET sequences from ID3C server with provided *lineage* and *segment*
    """
    r = requests.get(url, auth=(username,password), stream=True)

    with open(output, 'w+') as fasta_file:
        for line in r.iter_lines():
            if line:
                sequence = json.loads(line)
                fasta_file.write("".join([">", sequence['sample'], "\n", sequence['seq'].lower(), "\n"]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Downloads SFS sequences from ID3C of given lineage and segment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--output", type = str,
                        help="The output file for sequences, expected to be FASTA file")
    parser.add_argument("--lineage", type = str,
                        help="The lineage wildcard from Snakefile")
    parser.add_argument("--segment", type = str,
                        help="The segment wildcard from Snakefile")

    args = parser.parse_args()

    id3c_url = urljoin(os.environ["ID3C_URL"], "v1/shipping/genomic-data/")
    id3c_username = os.environ["ID3C_USERNAME"]
    id3c_password = os.environ["ID3C_PASSWORD"]

    lineage = get_full_lineage(args.lineage)
    url = generate_full_url(id3c_url, lineage, args.segment)

    get_sequences_from_id3c(url, id3c_username, id3c_password,
                            lineage, args.segment, args.output)
