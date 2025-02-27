using ..CryptoUtilities


leaves = ["1", "2", "3", "4", "5", "6", "7", "8"]
tree = build_merkle_tree(leaves)
for (i, layer) in enumerate(tree)
    println("Layer $i: ", layer)
end
# merkle_root = get_root(tree)
# println("Merkle Root: ", merkle_root)

path = get_path(3, 3)
# for pp in path
#     println(pp)
# end

layers = batch_paths(3, [6, 7, 3])
# Layer 1: Dict(5 => 1, 4 => 3, 8 => 2)
# Layer 2: Dict(4 => 1, 3 => 2, 1 => 3)
# Layer 3: Dict(2 => 2, 1 => 1)
# for i in 1:3
#     println("Layer $i: ", layers[i])
# end

(shifts, proof_size) = flatten_layers(layers)

println(shifts)
println(proof_size)

batched_proof = batch_proofs(3, tree, [6, 7, 3])
println(batched_proof)

opened_leaves = [leaves[6], leaves[7], leaves[3]]
res = verify_batched_proof(3, get_root(tree), opened_leaves, [6, 7, 3], batched_proof)
println(res)