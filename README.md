# CryptoUtilities.jl

`CryptoUtilities.jl` is a Julia package providing utilities for succinct proofs,
including hardware-optimized implementations of binary fields, Reed-Solomon
encodings, and Merkle trees with multiopening proofs.

## Features

- **Binary Fields**: Efficient arithmetic operations over binary fields.
- **Reed-Solomon Encoding**: Encoding and decoding using Reed-Solomon codes w/ `O(n log n)` cost, implementing [this paper](https://ieee-focs.org/FOCS-2014-Papers/6517a316.pdf).
- **Merkle Trees**: Construction, proof generation, and verification for Merkle tree batched openings.

## Installation

To use `BinaryFields` or any of the other packages, simply add the package to your Julia environment:

```julia
using Pkg
Pkg.add("https://github.com/bcc-research/CryptoUtilities.jl.git", subdir="BinaryFields")
```
Replace `subdir="BinaryFields"` with `subdir="BinaryReedSolomon"` or
`subdir="MerkleTree"` to install the other packages.

(We are working on adding the packages to the Julia package registry, but for
now, you can install them directly from the GitHub repository.)

## Usage

### Notebooks

We provide examples and walkthroughs of the package (along with some Julia
basics) as [Pluto.jl](https://plutojl.org) notebooks in the
[CryptoUtilitiesNotebooks](https://github.com/bcc-research/CryptoUtilitiesNotebooks)
repository. We _highly_ recommend reading those, which should be all you need to
get started using this library. We also provide some (very basic) exmaples below, for
completeness.

### Binary Fields

```julia
using BinaryFields

a = BinaryElem16(0x1)
b = BinaryElem16(0x2)
c = a * b
println(c)
```

### Reed-Solomon Encoding

```julia
using BinaryReedSolomon

rs = reed_solomon(BinaryElem16, 256, 1024)
message = rand(BinaryElem16, 256)
encoded = encode(rs, message)
```

### Merkle Trees

```julia
using MerkleTree

leaves = [rand(UInt16, 4) for _ in 1:16]
tree = build_merkle_tree(leaves)
root = get_root(tree)
```

### Ligero Protocol

```julia
using Ligero

poly = rand(BinaryElem16, 2^10)
comm = commit(poly)
proof = prove(comm, rand(BinaryElem128, size(matrix(comm), 2)), [1, 2, 3])
```

## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

