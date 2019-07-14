"""
    Given auspice-ready JSONs (i.e. from `augur export`) and a node--data JSON
    which defines which nodes should be hidden, produce an auspice-compatable
    JSON which hides these nodes (exploded-tree-like).
    NOTE: as of augur v6 this capability will be build into `augur export`
    NOTE: designed to work with augur v5 & auspice v1
    original author: James Hadfield                                 July 2019
"""
import argparse
import json


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract sample sequences by name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, metavar="JSON", help="Tree JSON from augur export")
    parser.add_argument("--output", required=True, metavar="JSON", help="Tree JSON for auspice")
    return parser.parse_args()


def post_order_traversal_iterative(tree, fn):
    s1 = [tree] 
    s2 = [] # Items added in postorder

    # Run while first stack is not empty 
    while s1:
        node = s1.pop()
        s2.append(node)
        if "children" in node:
            for child in node["children"]:
                s1.append(child)
    while s2:
        node = s2.pop()
        fn(node)


def get_cluster(node):
    try:
        return node["attr"]["cluster"]
    except KeyError:
        return False

def mark_clusters_as_hidden(node):
    if not "children" in node:
        return # ignore terminal nodes
    child_clusters = {get_cluster(child) for child in node["children"]}
    if len(child_clusters) == 1:
        node["attr"]["cluster"] = child_clusters.pop()
    else:
        node["hidden"] = "always"


if __name__ == '__main__':
    args = parse_args()

    with open(args.input, "rU") as fh:
        tree = json.load(fh)

    post_order_traversal_iterative(tree, mark_clusters_as_hidden)

    with open(args.output, "w") as fh:
        json.dump(tree, fh, indent=2, sort_keys=True)