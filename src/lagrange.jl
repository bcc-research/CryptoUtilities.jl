using BinaryFields

function evaluate_lagrange_basis(rs::Vector{T}) where T <: BinaryElem
    one_elem = T(one(T))
    current_layer = [one_elem + rs[1], rs[1]]
    len = 2
    for i in 2:length(rs)
        next_layer_size = 2 * len
        next_layer = Vector{T}(undef, next_layer_size)

        ri_p_one = one_elem + rs[i]
        for j in 1:len
            next_layer[2*j - 1] = current_layer[j] * ri_p_one
            next_layer[2*j]   = current_layer[j] * rs[i]
        end

        current_layer = next_layer
        len *= 2
    end

    return current_layer
end