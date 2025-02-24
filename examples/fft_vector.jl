using ..CryptoUtilities

n = 2^22

v = rand(GF2_128Elem, n)
v_clone = deepcopy(v)

fft!(v)
@assert v !== v_clone
ifft!(v)

@assert v == v_clone
