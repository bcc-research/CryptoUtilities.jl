using StatsBase

n = 20
K = 4
Q = 1000
N = 2^n

function generate_random_leaves(N, K)
    return [rand(UInt16, K) for _ in 1:N]
end

@testset "Random leaves" begin
    leaves = generate_random_leaves(N, K)
    tree = build_merkle_tree(leaves)

    queries = sort(sample(1:N, Q, replace=false))

    proof = prove(tree, queries)

    queried_leaves = [leaves[q] for q in queries]
    depth = get_depth(tree)
    root = get_root(tree)

    @test verify(root, proof; depth, leaves=queried_leaves, leaf_indices=queries)
end