using SHA
function hash_pair(left::String, right::String)
    return bytes2hex(sha256(left * right))
end

function build_merkle_tree(data::Vector{String})
    if isempty(data)
        return []
    end

    current_layer = [bytes2hex(sha256(d)) for d in data]
    tree = [current_layer]

    while length(current_layer) > 1
        next_layer = []
        for i in 1:2:length(current_layer)
            left = current_layer[i]
            right = current_layer[i + 1]
            push!(next_layer, hash_pair(left, right))
        end
        push!(tree, next_layer)
        current_layer = next_layer
    end

    return tree
end

function get_root(tree::Vector{Vector{String}})
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
                layer[next_cnt] = hash_pair(pp, layer[i])
            else
                layer[next_cnt] = hash_pair(layer[i], pp)
            end
            break
        end

        if query % 2 != 0
            proof_cnt += 1
            pp = proof[proof_cnt]
            layer[next_cnt] = hash_pair(pp, layer[i])
            i += 1
        else
            if queries[i + 1] != sibling
                proof_cnt += 1
                pp = proof[proof_cnt]
                layer[next_cnt] = hash_pair(layer[i], pp)
                i += 1
            else
                layer[next_cnt] = hash_pair(layer[i], layer[i + 1])
                i += 2
            end
        end
    end

    return (next_cnt, proof_cnt)
end

function verify_proof(root, depth, leaves, leaf_queries, batched_proof)
    proof = copy(batched_proof)
    layer = [bytes2hex(sha256(l)) for l in leaves]
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