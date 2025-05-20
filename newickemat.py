from collections import deque

def traverse_tree(tree):
    traverse_track = []
    if tree.size() == 0:
        return traverse_track

    stack = deque([(tree.root, -1)])
    

    while stack:
        node, ix = stack.pop()

        if ix >= 0:
            traverse_track.append((node, ix))
            continue
        
        children = tree.at(node).children
        num_children = len(children)
        stack.append((node, num_children))
        for i in reversed(range(num_children)):
            stack.append((children[i], -1))
            stack.append((node, i))  

    return traverse_track


    
def output_newick(tree, node_to_tip, os):
    for node, child_ix in traverse_tree(tree):
        node_ref = tree.at(node)

        if node == tree.root:
            parent_t = 0.0
        else:
            parent_t = tree.at(node_ref.parent).t

        if child_ix == 0 and node_ref.is_inner_node():
            os.write("(")

        if child_ix > 0 and child_ix < len(node_ref.children):
            os.write(",")

        if child_ix == len(node_ref.children):
            if node_ref.is_inner_node():
                os.write(")")
            elif node_ref.is_tip():
                os.write(node_to_tip.get(node, ""))

            label_entries = []

            if node_ref.mutations:
                curr_mutations =[]
                for m in node_ref.mutations:
                    curr_mutations.append(f"{m.from_base}{m.site + 1}{m.to_base},{m.t - parent_t:.6f}")
                label_entries.append("&mutations={" + ",".join(curr_mutations) + "}")

            if node_ref.missations:
                curr_missations = []
                for miss_map in node_ref.missations:
                    for site, from_base in miss_map.from_states.items():
                        curr_missations.append(f"{site}:{from_base}")
                label_entries.append("&missations={" + ",".join(curr_missations) + "}")

            if label_entries:
                os.write("[" + ",".join(label_entries) + "]")
            os.write(f":{node_ref.t - parent_t:.6f}")

    os.write(";")
