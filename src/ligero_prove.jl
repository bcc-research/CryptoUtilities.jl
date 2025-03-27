using CryptoUtilities
using BinaryFields

k = 10 
n = 12

rs = reed_solomon(BinaryElem16, 2^k, 2^n)

num_columns = 8
num_rows = 2^k
data = [rand(BinaryElem16) for _ in 1:(num_rows * num_columns)]
data = reshape(data, num_rows, num_columns)
msg = Message(data, 2^k)

code = encode(msg, rs)
tree = commit(code)

rs = [rand(BinaryElem16) for _ in 1:Int(log2(num_columns))]
gr = CryptoUtilities.evaluate_lagrange_basis(rs)
p = prove(msg, gr)