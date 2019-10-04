import os
import argparse
import json
import csv
from urllib.parse import urljoin
import requests

def get_metadata_from_id3c(id3c_url, id3c_username, id3c_password, output):
    r = requests.get(id3c_url, auth=(id3c_username, id3c_password), stream=True)
    stream = r.iter_lines()

    with open(output, 'w+') as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        for i, record in enumerate(map(json.loads, stream)):
            # Write the TSV header
            if i == 0:
                tsv_writer.writerow(record.keys())
            tsv_writer.writerow(record.values())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Downloads SFS metadata from ID3C",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--output", type = str, default="data/seattle_metadata.tsv",
        help="The output file for metadata, expected to be TSV file")
    args = parser.parse_args()

    id3c_url = urljoin(os.environ["ID3C_URL"], "v1/shipping/augur-build-metadata")
    id3c_username = os.environ["ID3C_USERNAME"]
    id3c_password = os.environ["ID3C_PASSWORD"]

    get_metadata_from_id3c(id3c_url, id3c_username, id3c_password, args.output)
