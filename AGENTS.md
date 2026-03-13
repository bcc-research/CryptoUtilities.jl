# CryptoUtilities.jl Agent Guide

## Repo shape
- This is a small Julia monorepo with three standalone packages: `BinaryFields`, `BinaryReedSolomon`, and `BatchedMerkleTree`.
- There is no root `Project.toml`. Work inside the package you are changing.
- Each package carries its own `Project.toml`, `Manifest.toml`, `src/`, and `test/`.
- Keep changes package-local. Do not introduce repo-wide scaffolding unless the work genuinely spans packages.

## Style to preserve
- Keep the code minimal, direct, and close to the math. The repo favors plain Julia over framework-heavy abstractions.
- Use 4-space indentation and small, focused functions.
- Export lists live in the top-level module file. Keep exports explicit there.
- When extending Base, import the names explicitly near the module header, then define the methods in the relevant source file.
- Hot-path code should stay allocation-aware and type-stable. Existing code uses `@generated` functions, `@views`, broadcast fusion, and architecture-specific intrinsics only where they clearly pay for themselves.
- Comments are sparse and technical. Add them when they clarify an algorithm, invariant, or low-level trick, not for obvious line-by-line narration.
- Docstrings are selective, not exhaustive. Match the local tone: brief, factual, and focused on behavior.
- Avoid new dependencies unless they are clearly justified. This repo is intentionally light.
- If you touch low-level arithmetic in `BinaryFields`, preserve the current split between generic code and architecture-specific carryless multiply paths.
- `Manifest.toml` files are checked in per package. Do not churn them unless the dependency change is part of the task.

## Package notes
- `BinaryFields` is the low-level performance-sensitive package. It mixes generic polynomial/field code with explicit architecture branches and a `warmup.jl` include used to stabilize compilation behavior.
- `BinaryReedSolomon` builds on `BinaryFields` and prefers in-place transforms and precomputed twiddles.
- `BatchedMerkleTree` is intentionally simple and straightforward; do not overengineer it.

## Tests
- Tests are package-local and use `Test` directly.
- The usual pattern is a thin `test/runtests.jl` that includes one or more focused test files.
- Favor algebraic laws, round-trip properties, randomized loops with modest sample counts, and fixed regression vectors.
- Match existing test style: explicit `@testset`s, small helper functions only when they remove repetition, and no unnecessary test infrastructure.
- Standard commands:
  - `julia --project=BinaryFields -e 'using Pkg; Pkg.test()'`
  - `julia --project=BinaryReedSolomon -e 'using Pkg; Pkg.test()'`
  - `julia --project=BatchedMerkleTree -e 'using Pkg; Pkg.test()'`
- In this Codex environment, launching Julia may require elevated permissions because the local `julia` command is managed by `juliaup`.
- Current observed status on March 13, 2026:
  - `BinaryFields` tests pass.
  - `BinaryReedSolomon` tests pass.
  - `BatchedMerkleTree` tests currently fail immediately because `test/runtests.jl` imports `MerkleTree` instead of `BatchedMerkleTree`.
- The checked-in manifests were resolved with Julia 1.11.2; running tests under Julia 1.12.4 emits manifest-version warnings.

## Benchmarking guidance
- For Julia benchmarks, follow standard `BenchmarkTools` practice:
  - benchmark steady-state code, not compilation;
  - put the timed work in a function;
  - interpolate external values with `$`;
  - move input generation into `setup=` or precompute it outside the timed expression.
- If persistent benchmarks are added, keep them package-local. For this repo, that means `BinaryFields/benchmark/`, not a repo-root benchmark harness.
- Prefer the usual Julia package layout for benchmark suites:
  - `BinaryFields/benchmark/Project.toml` for benchmark-only dependencies;
  - `BinaryFields/benchmark/benchmarks.jl` defining a `BenchmarkGroup`/`SUITE`, which matches common `PkgBenchmark` usage.
- For arithmetic kernels, benchmark one operation at a time on concrete element types and fixed operand shapes. Avoid timing random number generation, conversions, or setup work unless that setup cost is the thing being studied.
- When benchmarking finite-field multiplication, keep the comparison honest:
  - warm up first;
  - benchmark per-type implementations separately;
  - include representative operand patterns, not just fully random inputs;
  - verify correctness against the current implementation before comparing timings.

## Editing expectations
- Prefer small diffs.
- Preserve the existing file layout unless there is a concrete reason to change it.
- When adding a new source file, wire it in through the package module file with an explicit `include`.
- When adding performance work, pair it with correctness tests before or alongside the optimization.
