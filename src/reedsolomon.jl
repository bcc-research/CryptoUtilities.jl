mutable struct ReedSolomonEncoding{T <: BinaryFieldElem}
    deg::UInt
    log_block_length::UInt
    parity_elem_size::UInt
    log_parity_elem_size::UInt
    twiddles_chunks::Matrix{T}
    set_twiddles::Bool

    function ReedSolomonEncoding{T}(deg, log_block_length) where T <: BinaryFieldElem
        parity_elem_size = 2^log_block_length - deg
        @assert is_pow_2(parity_elem_size)
        log_parity_elem_size = round(UInt, log2(parity_elem_size))

        return new{T}(
            deg,
            log_block_length,
            parity_elem_size,
            log_parity_elem_size,
            zeros(T, log_chunk_size, div(2^log_block_length, parity_elem_size)+1)
        )
    end
end

chunk_size(rs::ReedSolomonEncoding) = rs.chunk_size
log_chunk_size(rs::ReedSolomonEncoding) = rs.log_chunk_size
degree(rs::ReedSolomonEncoding) = rs.deg

function compute_twiddles!(rs::ReedSolomonEncoding{T}) where T
    for (i, twiddles_i) in enumerate(eachcol(rs.twiddles_chunks))
        beta = T((i-1)*chunk_size)
        compute_twiddles!(twiddles_i, beta, rs.log_chunk_size)
    end

    rs.set_twiddles = true
end

split_vector(rs::ReedSolomonEncoding{T}, v) where T <: BinaryFieldElem = Iterators.partition(v, rs.parity_elem_size)

function encode(rs::ReedSolomonEncoding{T}, message::Vector{T}) where T
    @assert length(message) == rs.deg
    @assert rs.set_twiddles
    message_copy = deepcopy(message)
    msg_chunks = split_vector(rs, message_copy)
    
    v_parity = zeros(T, rs.parity_elem_size)

    # Drops the first chunk, used later
    for (msg_chunk, twiddle_chunk) in Iterators.drop(zip(msg_chunks, eachcol(rs.twiddles_chunks)), 1)
        ifft_twiddles!(msg_chunk; twiddles=twiddle_chunk)
        v_parity .+= msg_chunk
    end

    fft!(v_parity, twiddles=rs.twiddles_chunks[1]);

    return [v_parity; message]
end

function encode!(vector::Vector{T}; inv_rate=2, twiddles=nothing) where T <: BinaryFieldElem
    @assert is_pow_2(length(vector))

    message_length = div(length(vector), inv_rate)
    if isnothing(twiddles)
        twiddles = compute_twiddles(zero(GF2_128Elem), round(UInt, log2(length(vector))))
    end
    ifft!((@view vector[1:message_length]); twiddles=(@view twiddles[1:message_length]))
    vector[message_length+1:end] .= zero(T)
    fft!(vector; twiddles)
end