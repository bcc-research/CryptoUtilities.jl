struct Message{T} where T <: BinaryElem
    data::Matrix{T}
end

function Message(data::Matrix{T}, column_length::Int) where T <: BinaryElem
    if size(data, 1) != column_length
        error("Invalid message: Expected each column to have length $column_length, but got $(size(data, 1)).")
    end
    new{T}(data)
end

export Message