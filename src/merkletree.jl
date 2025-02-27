using SHA
function hash_children(left::String, right::String)
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
            push!(next_layer, hash_children(left, right))
        end
        push!(tree, next_layer)
        current_layer = next_layer
    end

    return tree
end

function get_root(tree::Vector{Vector{String}})
    return isempty(tree) ? "" : tree[end][1]
end

function get_path(depth::Int, query::Int)
    path = Vector{Int}(undef, depth)

    for i in 1:depth
        mod = query % 2
        path[i] = mod == 1 ? query + 1 : query - 1
        query = (query + mod) ÷ 2
    end

    return path
end

function batch_paths(depth::Int, queries::Vector{Int})
    layers = [Dict{Int, Int}() for _ in 1:depth]

    for q in queries
        path = get_path(depth, q)

        for i in 1:depth
            if !haskey(layers[i], path[i])
                layers[i][path[i]] = length(layers[i]) + 1
            end
        end
    end

    return layers
end

function flatten_layers(layers::Vector{Dict{Int, Int}})
    shifts = Vector{Int}(undef, length(layers))
    proof_size = 0

    for i in 1:length(layers)
        shifts[i] = proof_size
        proof_size += length(layers[i])
    end

    return (shifts, proof_size)
end

# PROVER
function batch_proofs(depth::Int, tree_levels::Vector{Vector{String}}, queries::Vector{Int})
    layers = batch_paths(depth, queries)
    (shifts, proof_size) = flatten_layers(layers)

    proof = Vector{String}(undef, proof_size)
    for i in 1:depth
        for (k, v) in layers[i]
            idx = v + shifts[i]
            proof[idx] = tree_levels[i][k]
        end
    end

    return proof
end

# VERIFIER 
function verify_merkle_proof(leaf::String, proof::Vector{String}, root::String, query::Int)
    current_hash = bytes2hex(sha256(leaf))
    for sibling_hash in proof
        mod = query % 2
        if mod == 1
            current_hash = hash_children(current_hash, sibling_hash)
        else
            current_hash = hash_children(sibling_hash, current_hash)
        end
        query = (mod + query) ÷ 2
    end
    return current_hash == root
end

function index_batched_proof(depth::Int, layers::Vector{Dict{Int, Int}}, shifts::Vector{Int}, path::Vector{Int}, batched_proof::Vector{String})
    reconstructed_proof = Vector{String}(undef, depth)

    for i in 1:depth
        layer_idx = layers[i][path[i]]
        abs_idx = layer_idx + shifts[i]
        reconstructed_proof[i] = batched_proof[abs_idx]
    end

    return reconstructed_proof
end

function verify_batched_proof(depth::Int, root::String, leaves::Vector{String}, queries::Vector{Int}, batched_proof::Vector{String})
    layers = batch_paths(depth, queries)
    (shifts, proof_size) = flatten_layers(layers)

    @assert length(batched_proof) == proof_size
    @assert length(leaves) == length(queries)

    # For each leaf
    for i in 1:length(leaves)
        # Get path for leaf
        path = get_path(depth, queries[i])
        # Verify opening proof from reconstructed proof
        reconstructed_proof = index_batched_proof(depth, layers, shifts, path, batched_proof)
        @assert verify_merkle_proof(leaves[i], reconstructed_proof, root, queries[i])
    end

    return true
end

export hash_children, build_merkle_tree, get_root, get_path, batch_paths, flatten_layers
export batch_proofs, verify_merkle_proof, index_batched_proof, verify_batched_proof