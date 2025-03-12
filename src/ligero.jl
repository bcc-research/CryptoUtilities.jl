using BinaryFields

struct Message{T <: BinaryElem}
    data::Matrix{T}
end

function Message(data::Matrix{T}, column_length::Int) where T <: BinaryElem
    if size(data, 1) != column_length
        error("Invalid message: Expected each column to have length $column_length, but got $(size(data, 1)).")
    end
    Message{T}(data)
end
export Message

struct InterleavedCode{T <: BinaryElem}
    code::Matrix{T}
end
export InterleavedCode

function InterleavedCode(data::Matrix{T}, column_length::Int) where T <: BinaryElem
    if size(data, 1) != column_length
        error("Invalid code: Expected each column to have length $column_length, but got $(size(data, 1)).")
    end
    InterleavedCode{T}(data)
end

function encode(msg::Message{T}, rs::ReedSolomonEncoding{T}) where {T <: BinaryElem}
    @assert size(msg.data, 1) == message_length(rs) "Each message column must have length $(message_length(rs)), got $(size(msg.data, 1))"
    encoded_columns = [encode(rs, Vector(col)) for col in eachcol(msg.data)]
    encoded_matrix = hcat(encoded_columns...)
    return InterleavedCode(encoded_matrix)
end

function commit(code::InterleavedCode{T}) where {T <: BinaryElem}
    leaves = [map(BinaryFields.binary_val, Vector(row)) for row in eachrow(code.code)]
    return build_merkle_tree(leaves)
end

# TODO: convert from GF_2^16 to GF_2^128
function prove(msg::Message{T}, gr::Vector{T}) where {T <: BinaryElem}
    @assert size(msg.data, 2) == length(gr) "Number of columns must match length of gr"
    return msg.data * gr
end

export encode, commit, prove