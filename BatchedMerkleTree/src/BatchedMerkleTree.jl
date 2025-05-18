module BatchedMerkleTree

using SHA
import Base: sizeof

export build_merkle_tree, get_root, get_depth
export CompleteMerkleTree, MerkleRoot, BatchedMerkleProof
export prove, verify

# XXX: make all proofs into reasonable objects with good ergonomics

# for now simply assume that isbitstype(T) = true
hash_leaf(leaf) = sha256(reinterpret(UInt8, leaf))
hash_siblings(left_hash::Vector{UInt8}, right_hash::Vector{UInt8}) = sha256(vcat(left_hash, right_hash))

# function hash_siblings_hex(left_hash::String, right_hash::String)
#     left_bytes = hex2bytes(left_hash)
#     right_bytes = hex2bytes(right_hash)
    
#     hash_bytes = sha256(vcat(left_bytes, right_bytes))
#     return bytes2hex(hash_bytes)
# end

function is_power_of_two(n::Int)
    return n > 0 && (n & (n - 1)) == 0
end

struct CompleteMerkleTree
    tree::Vector{Vector{Vector{UInt8}}}
end

function build_merkle_tree(leaves)
    @assert is_power_of_two(length(leaves))
    if isempty(leaves)
        return []
    end

    current_layer = [hash_leaf(leaf) for leaf in leaves]
    tree = [current_layer]

    while length(current_layer) > 1
        next_layer_size = length(current_layer) ÷ 2
        next_layer = [Vector{UInt8}(undef, 32) for _ in 1:next_layer_size]

        cnt = 1
        for i in 1:2:length(current_layer)
            left = current_layer[i]
            right = current_layer[i + 1]
            next_layer[cnt] = hash_siblings(left, right)
            cnt += 1
        end

        push!(tree, next_layer)
        current_layer = next_layer
    end

    return CompleteMerkleTree(tree)
end

struct MerkleRoot
    root::Union{Vector{UInt8}, Nothing}
end

sizeof(x::MerkleRoot) = sizeof(x.root)
get_root(c::CompleteMerkleTree) = MerkleRoot(isempty(c.tree) ? nothing : c.tree[end][1])
get_depth(c::CompleteMerkleTree) = length(c.tree) - 1

function ith_layer!(current_layer, queries_len, queries, proof)
    next_queries_len = 0
    i = 1

    while i <= queries_len
        query = queries[i]
        sibling = query ⊻ 1

        next_queries_len += 1
        queries[next_queries_len] = query >> 1

        if i == queries_len
            push!(proof, current_layer[sibling + 1])
            break
        end

        if query % 2 != 0
            push!(proof, current_layer[sibling + 1])
            i += 1
        else
            if queries[i + 1] != sibling
                push!(proof, current_layer[sibling + 1])
                i += 1
            else
                i += 2
            end
        end
    end

    return next_queries_len
end

struct BatchedMerkleProof
    proof::Vector{Vector{UInt8}}
end

sizeof(x::BatchedMerkleProof) = sum(length, x.proof)

function prove(tree::CompleteMerkleTree, queries)
    proof = Vector{Vector{UInt8}}()
    depth = length(tree.tree) - 1

    queries_buff = copy(queries)

    queries_buff .-= 1
    queries_cnt = length(queries)
    for i in 1:depth
        queries_cnt = ith_layer!(tree.tree[i], queries_cnt, queries_buff, proof)
    end

    return BatchedMerkleProof(proof)
end

function verify_ith_layer!(layer, queries, curr_cnt, proof, proof_cnt)
    next_cnt = 0
    i = 1

    while i <= curr_cnt
        query = queries[i]
        sibling = query ⊻ 1

        next_cnt += 1
        queries[next_cnt] = query >> 1

        if i == curr_cnt
            proof_cnt += 1
            pp = proof[proof_cnt]
            if query % 2 != 0
                layer[next_cnt] = hash_siblings(pp, layer[i])
            else
                layer[next_cnt] = hash_siblings(layer[i], pp)
            end
            break
        end

        if query % 2 != 0
            proof_cnt += 1
            pp = proof[proof_cnt]
            layer[next_cnt] = hash_siblings(pp, layer[i])
            i += 1
        else
            if queries[i + 1] != sibling
                proof_cnt += 1
                pp = proof[proof_cnt]
                layer[next_cnt] = hash_siblings(layer[i], pp)
                i += 1
            else
                layer[next_cnt] = hash_siblings(layer[i], layer[i + 1])
                i += 2
            end
        end
    end

    return (next_cnt, proof_cnt)
end

function verify(root::MerkleRoot, batched_proof::BatchedMerkleProof; depth, leaves, leaf_indices)
    if root.root === nothing
        return false
    end
    
    proof = copy(batched_proof.proof)
    layer = [hash_leaf(leaf) for leaf in leaves]
    queries = copy(leaf_indices)

    queries .-= 1

    curr_cnt = length(queries)
    proof_cnt = 0

    for _ in 1:depth
        (curr_cnt, proof_cnt) = verify_ith_layer!(layer, queries, curr_cnt, proof, proof_cnt)
    end

    return layer[1] == root.root
end

end # module BatchedMerkleTree
