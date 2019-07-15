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

def make_stub(node):
    
    if "cluster" not in node["attr"] or node["attr"]["cluster"] == False:
        return node

    stub = {
        "strain": "stub_{}".format(get_cluster(node)),
        "attr": {
            "num_date": node["attr"]["num_date"],
            "div": node["attr"]["div"]
        },
        "hidden": "always",
        "children": [node],
        "clade": int("100000{}".format(get_cluster(node))) # required in auspice v1
    }

    return stub


def mark_clusters_as_hidden(node):
    if not "children" in node:
        return # ignore terminal nodes
    child_clusters = {get_cluster(child) for child in node["children"]}
    if len(child_clusters) == 1:
        node["attr"]["cluster"] = child_clusters.pop()
    else:
        node["hidden"] = "always"
        # the children here will each lead to a different cluster
        # and we want to make them "stubs", i.e. not show the long
        # branch length leading to the start of the cluster
        node["children"] = [make_stub(child) for child in node["children"]]


def flatten_tree(tree):
    flat = []
    stack = [tree]
    while(len(stack) > 0):  
        node = stack.pop()
        flat.append(node)
        if "children" in node:
            for child in node["children"]:
                stack.append(child)
    return flat


def shift_hidden_div_and_time(flat_tree):
    """
        Shift the num_date and div attributes of the hidden nodes
        so that the non-hidden nodes take up the entire screen in
        auspice. (Auspice calculates the view based on all the nodes,
        whether they are hidden or not, so a deep MRCA -- even if
        hidden -- doesn't look great.)

        This uses the trick that both flat_tree & tree contain references
        to the same object, so modifications to one can affect the other
    """
    min_date = min([n["attr"]["num_date"] for n in flat_tree if get_cluster(n)])
    min_div = min([n["attr"]["div"] for n in flat_tree if get_cluster(n)]) # cumulative!
    for n in flat_tree:
        if n["attr"]["num_date"] < min_date:
            n["attr"]["num_date"] = min_date
            if "num_date_confidence" in n["attr"]:
                del n["attr"]["num_date_confidence"]
        if n["attr"]["div"] < min_div:
            n["attr"]["div"] = 0
        else:
            n["attr"]["div"] -= min_div


if __name__ == '__main__':
    args = parse_args()

    with open(args.input, "rU") as fh:
        tree = json.load(fh)

    post_order_traversal_iterative(tree, mark_clusters_as_hidden)
    flat_tree = flatten_tree(tree)
    shift_hidden_div_and_time(flat_tree)

    with open(args.output, "w") as fh:
        json.dump(tree, fh, indent=2, sort_keys=True)