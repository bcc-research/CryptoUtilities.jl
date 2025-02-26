is_pow_2(n) = 2^(round(Int, log2(n))) == n
unwrap_value(a::VecElement{T}) where T = a.value