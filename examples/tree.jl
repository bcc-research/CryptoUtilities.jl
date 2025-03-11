using ..CryptoUtilities

leaves = [UInt16[1, 2, 3], UInt16[4, 5, 6], UInt16[7, 8, 9], UInt16[10, 11, 12], UInt16[13, 14, 15], UInt16[16, 17, 18], UInt16[19, 20, 21], UInt16[22, 23, 24]]
tree = build_merkle_tree(leaves)

# for layer in tree
#     println(layer)
# end

depth = get_depth(tree)
root = get_root(tree)

queries = [0, 2, 5, 7]
proof = create_proof(tree, queries)

queried_leaves = [leaves[q + 1] for q in queries]
res = verify_proof(root, depth, queried_leaves, queries, proof)
@assert res == true