mutable struct ReedSolomonEncoding{T <: BinaryElem}
    log_message_length::UInt
    log_block_length::UInt
    twiddles::Union{Nothing, Vector{T}}

    function ReedSolomonEncoding{T}(log_message_length, log_block_length; compute_twiddles=true) where T <: BinaryElem
        return new{T}(
            log_message_length,
            log_block_length,
            nothing
        )
    end
end

log_message_length(rs::ReedSolomonEncoding) = rs.log_message_length
message_length(rs::ReedSolomonEncoding) = 2^log_message_length(rs)

log_block_length(rs::ReedSolomonEncoding) = rs.log_block_length
block_length(rs::ReedSolomonEncoding) = 2^log_block_length(rs)


function compute_twiddles!(rs::ReedSolomonEncoding{T}; beta=T(0)) where T
    rs.twiddles = compute_twiddles(T, log_block_length(rs); beta)
end

# Converts twiddles from long vector to twiddles from shorter one
function short_from_long_tw(l_tw, n, k)
	s_tw = Vector{eltype(l_tw)}(undef, 2^k - 1)
	jump = 2^(n - k)
	s_tw[1] = l_tw[jump]

	idx = 2
	for i in 1:k-1
		jump *= 2
		take = 2^i

		@views @. s_tw[idx:idx + take - 1] = l_tw[jump:jump + take - 1]
		idx += take
	end

	return s_tw
end

function short_from_long_twiddles(rs::ReedSolomonEncoding{T}) where T 
    
    short_from_long_tw(rs.twiddles, log_block_length(rs), log_message_length(rs))
end

function encode(rs::ReedSolomonEncoding{T}, message::Vector{T}; verbose=false) where T
    message_coeffs = zeros(T, block_length(rs))

    message_coeffs_view = @view message_coeffs[1:message_length(rs)]
    message_coeffs_view .= message

    encode!(rs, message_coeffs; verbose)

    return message_coeffs
end

function encode!(rs::ReedSolomonEncoding{T}, v::Vector{T}; verbose=false, fill_zeros=true, short_twiddles=nothing) where T
    @assert length(v) == block_length(rs)
    @assert !isnothing(rs.twiddles)

    message_coeffs_view = @view v[1:message_length(rs)]
    if fill_zeros
        message_coeffs_view[message_length(rs):end] .= T(0)
    end

    s_tw = isnothing(short_twiddles) ? short_from_long_tw(rs.twiddles, log_block_length(rs), log_message_length(rs)) : short_twiddles

    ifft!(message_coeffs_view; twiddles=s_tw)
    fft!(v, twiddles=rs.twiddles; verbose);

    v
end

function reed_solomon(::Type{T}, message_length, block_length) where T
    @assert is_pow_2(block_length) && is_pow_2(message_length)
    @assert message_length < block_length

    log_message_length = round(Int, log2(message_length))
    log_block_length = round(Int, log2(block_length))

    rs = ReedSolomonEncoding{T}(log_message_length, log_block_length)
    compute_twiddles!(rs)

    return rs
end