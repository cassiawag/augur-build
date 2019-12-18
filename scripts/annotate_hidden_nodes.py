"""
    Given auspice-ready JSONs (i.e. from `augur export`) produce an
    auspice-compatable JSON which hides nodes that are basal to monophyletic
    clusters (as defined by the 'cluster' node attribute).
    original author: James Hadfield                                 July 2019
"""
import argparse
import json

## Magics / hardcoded parameters
STUB_DIV_LENGTH = 0.0002
STUB_TIME_LENGTH = 0.05 # ~18 days

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract sample sequences by name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, metavar="JSON", help="Unified JSON from augur export v2")
    parser.add_argument("--output", required=True, metavar="JSON", help="Unified JSON for auspice")
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
        return node["node_attrs"]["cluster"]["value"]
    except KeyError:
        return False

def make_stub(node):
    cluster = get_cluster(node)
    # if the node doesn't have a cluster then it doesn't need to be a stub!
    if not cluster:
        return node

    stub = {
        "name": "stub_{}".format(cluster),
        "node_attrs": {
            "num_date": {
                "value": node["node_attrs"]["num_date"]["value"] - STUB_TIME_LENGTH
            },
            "div": node["node_attrs"]["div"] - STUB_DIV_LENGTH,
            "cluster": {
                "value": cluster
            },
            "hidden": "always"
        },
        "children": [node]
    }

    return stub


def mark_clusters_as_hidden(node):
    if not "children" in node:
        return # ignore terminal nodes
    child_clusters = {get_cluster(child) for child in node["children"]}
    # if (node["name"] == "NODE_0000000"): import pdb; pdb.set_trace()
    if child_clusters == {False}:
        # nodes where all children don't have a cluster should be hidden
        node["node_attrs"]["hidden"] = "always"
    elif len(child_clusters) == 1:
        # clusters are not annotated on internal nodes, but this sets them
        # on internal nodes if all children (a) have a cluster set and
        # (b) the cluster is the same
        node["node_attrs"]["cluster"] = {
            "value": child_clusters.pop()
        }
    else:
        # nodes with different clusters amoung the children are those we want to
        # "chop off"
        node["node_attrs"]["hidden"] = "always"
        # the children here will each lead to a different cluster
        # and we want to make them "stubs", i.e. not show the long
        # branch length leading to the start of the cluster
        node["children"] = [make_stub(child) for child in node["children"]]


def flatten_tree(tree):
    flat = []
    stack = [tree]
    while len(stack) > 0:
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

        Note: plenty of optimisation potential here
    """

    # find the date of the earliest node which has a cluster and modify all
    # nodes earlier than that to have that date
    min_date = min([n["node_attrs"]["num_date"]["value"] for n in flat_tree if get_cluster(n)])
    for n in flat_tree:
        if n["node_attrs"]["num_date"]["value"] < min_date:
            # setting the `num_date` here will also overwrite any confidences if they are set
            n["node_attrs"]["num_date"] = {"value": min_date}

    # find the min divergence (cumulative!) for each cluster
    min_div_per_cluster = {}
    for n in flat_tree:
        cluster = get_cluster(n)
        if cluster:
            if cluster not in min_div_per_cluster:
                min_div_per_cluster[cluster] = n["node_attrs"]["div"]
            elif n["node_attrs"]["div"] < min_div_per_cluster[cluster]:
                min_div_per_cluster[cluster] = n["node_attrs"]["div"]

    # modify clusters to each have divergence starting from 0
    for n in flat_tree:
        cluster = get_cluster(n)
        if "hidden" in n["node_attrs"]:
            n["node_attrs"]["div"] = 0
        elif cluster:
            n["node_attrs"]["div"] -= min_div_per_cluster[cluster]


if __name__ == '__main__':
    args = parse_args()

    with open(args.input, "rU") as fh:
        unified = json.load(fh)

    tree = unified["tree"]
    post_order_traversal_iterative(tree, mark_clusters_as_hidden)
    flat_tree = flatten_tree(tree)
    shift_hidden_div_and_time(flat_tree)

    unified["tree"] = tree
    with open(args.output, "w") as fh:
        json.dump(unified, fh, indent=2, sort_keys=True)
