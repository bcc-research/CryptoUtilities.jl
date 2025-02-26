using ..CryptoUtilities


data = ["block1", "block2", "block3", "block4"]
tree = build_merkle_tree(data)
for (i, layer) in enumerate(tree)
    println("Layer $i: ", layer)
end
merkle_root = compute_merkle_root(tree)
println("Merkle Root: ", merkle_root)