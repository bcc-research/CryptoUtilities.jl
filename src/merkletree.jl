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

function compute_merkle_root(tree::Vector{Vector{String}})
    return isempty(tree) ? "" : tree[end][1]
end

