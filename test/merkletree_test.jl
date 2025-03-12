using CryptoUtilities
using Random

n = 20
K = 4
Q = 1000
N = 2^n

function sample(N::Int, Q::Int)
    if Q > N
        throw(ArgumentError("Q must be less than or equal to N"))
    end

    sampled_indices = Set{Int}()
    while length(sampled_indices) < Q
        index = rand(0:N-1)
        push!(sampled_indices, index)
    end

    sorted_indices = sort(collect(sampled_indices))
    return sorted_indices
end

function generate_random_leaves(N, K)
    return [rand(UInt16, K) for _ in 1:N]
end

leaves = generate_random_leaves(N, K)
tree = build_merkle_tree(leaves)

queries = sample(N, Q)

proof = create_proof(tree, queries)

queried_leaves = [leaves[q + 1] for q in queries]
depth = get_depth(tree)
root = get_root(tree)

res = verify_proof(root, depth, queried_leaves, queries, proof)
@assert res == true