using CryptoUtilities
using BinaryFields

k = 5 
rs = [rand(BinaryElem16) for _ in 1:k]

evals = CryptoUtilities.evaluate_lagrange_basis(rs)