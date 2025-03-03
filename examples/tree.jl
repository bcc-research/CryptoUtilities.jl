using ..CryptoUtilities


leaves = ["1", "2", "3", "4", "5", "6", "7", "8"]
tree = build_merkle_tree(leaves)
depth = get_depth(tree)
root = get_root(tree)

queries = [0, 2, 5, 7]
proof = create_proof(tree, queries)

queried_leaves = [leaves[q + 1] for q in queries]
res = verify_proof(root, depth, queried_leaves, queries, proof)
@assert res == true