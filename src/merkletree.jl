using SHA

# for now simply assume that isbitstype(T) = true
hash_leaf(leaf::Vector{T}) where T = sha256(reinterpret(UInt8, leaf))
hash_siblings(left_hash::Vector{UInt8}, right_hash::Vector{UInt8}) = sha256(vcat(left_hash, right_hash))

function hash_siblings_hex(left_hash::String, right_hash::String)
    left_bytes = hex2bytes(left_hash)
    right_bytes = hex2bytes(right_hash)
    
    hash_bytes = sha256(vcat(left_bytes, right_bytes))
    return bytes2hex(hash_bytes)
end

function is_power_of_two(n::Int)
    return n > 0 && (n & (n - 1)) == 0
end

function build_merkle_tree(leaves::Vector{Vector{T}}) where T
    @assert is_power_of_two(length(leaves))
    if isempty(leaves)
        return []
    end

    current_layer = [hash_leaf(leaf) for leaf in leaves]
    tree = [[bytes2hex(leaf) for leaf in current_layer]]

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

        push!(tree, [bytes2hex(node) for node in next_layer])
        current_layer = next_layer
    end

    return tree
end

function get_root(tree)
    return isempty(tree) ? "" : tree[end][1]
end

function get_depth(tree)
    return length(tree) - 1
end

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

function create_proof(tree_levels, queries)
    proof = String[]
    depth = length(tree_levels) - 1

    queries_buff = copy(queries)
    queries_cnt = length(queries)
    for i in 1:depth
        queries_cnt = ith_layer!(tree_levels[i], queries_cnt, queries_buff, proof)
    end

    return proof
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
                layer[next_cnt] = hash_siblings_hex(pp, layer[i])
            else
                layer[next_cnt] = hash_siblings_hex(layer[i], pp)
            end
            break
        end

        if query % 2 != 0
            proof_cnt += 1
            pp = proof[proof_cnt]
            layer[next_cnt] = hash_siblings_hex(pp, layer[i])
            i += 1
        else
            if queries[i + 1] != sibling
                proof_cnt += 1
                pp = proof[proof_cnt]
                layer[next_cnt] = hash_siblings_hex(layer[i], pp)
                i += 1
            else
                layer[next_cnt] = hash_siblings_hex(layer[i], layer[i + 1])
                i += 2
            end
        end
    end

    return (next_cnt, proof_cnt)
end

function verify_proof(root, depth, leaves, leaf_queries, batched_proof)
    proof = copy(batched_proof)
    layer = [bytes2hex(hash_leaf(leaf)) for leaf in leaves]
    queries = copy(leaf_queries)

    curr_cnt = length(queries)
    proof_cnt = 0

    for _ in 1:depth
        (curr_cnt, proof_cnt) = verify_ith_layer!(layer, queries, curr_cnt, proof, proof_cnt)
    end

    return layer[1] == root
end

export build_merkle_tree, get_root, get_depth
export create_proof, verify_proof