# CUDA 13 upgrade — change summary

This documents the changes on `gpuOptimization` that get Vlasiator building and running
against **CUDA 13** on Roihu-GPU (GH200), via a new `nvc++`-based build architecture,
`roihu_gpu_nvhpc`. An earlier phase of this effort also tried out fsgrid's experimental
`std::execution::par`-based `experimental_for` feature (stdpar/RDC) on the same toolchain;
that approach was abandoned mid-effort due to a cascade of `nvc++`/RDC-specific issues
(nvlink duplicate-registration collisions, GCC15-hardening interactions, a runtime hang) and
fully reverted — nothing below relates to it. Everything here is what a "just get this
building and running on CUDA 13 / the `nvhpc` module" effort requires, using the same classic
`__global__`/`__device__` kernel style as the existing `roihu_gpu` (nvcc) build.

## Toolchain / build scaffolding (new arch `roihu_gpu_nvhpc`)

- **`MAKE/Makefile.roihu_gpu_nvhpc`** (new file) — build arch definition using `nvc++`/`mpic++`
  as `CMP`/`LNK` instead of classic `nvcc`, since `nvhpc/26.3` is the only module on Roihu that
  bundles CUDA 13 (13.1) alongside a GH200-compatible HPCX MPI. `-gpu=cc90` restricts codegen
  to Hopper/sm_90 only (without it, `nvc++` generates device code for every compute capability
  it knows about, ~12 architectures, which is pathologically slow to compile). `-cuda` is
  required in `LDFLAGS`: `-gpu=cc90` alone has no effect at link time unless paired with a mode
  flag (`nvc++` warns exactly this) — `-x cu` masks this at compile time, but the separate link
  step needs `-cuda` explicitly, or `nvc++` silently skips linking its own CUDA runtime and
  every `cuda*` symbol comes back undefined. No `--gcc-toolchain=` override is needed — this was
  tried, but confirmed unnecessary on the real system (an earlier concern that `nvc++` couldn't
  find a usable C++ standard library without one turned out to be specific to a sandboxed test
  environment, not this system).
- **`modules/roihu_gpu_nvhpc.sh`** (new file) — just `module load nvhpc/26.3`. Nothing else is
  needed: no separate `gcc` module (no toolchain override required), no separate `openmpi`
  module (`nvhpc` bundles its own HPCX MPI — loading another would mix two MPI/OpenMP runtime
  implementations in one binary), no `papi` module (PAPI isn't linked for this arch yet).
- **`build_fetched_libraries.sh`** — added a `-roihu_gpu_nvhpc` case for both build phases:
  points dependency-library builds at `mpic++` and at NVHPC's bundled CUDA 13.1 `nvtx3` include
  path (`.../nvhpc/Linux_aarch64/26.3/cuda/13.1/targets/sbsa-linux/include/nvtx3/`) instead of
  the classic `$CUDA_HOME/include/nvtx3`.

## CUDA 13 API compatibility

CUDA 13's CCCL/runtime headers dropped or changed a few APIs the code depended on. Both fixes
are version-gated on `CUDART_VERSION >= 13000` so the CUDA-12.9.1-based classic `roihu_gpu`
build path is unaffected.

- **`cub::Max()`/`cub::Min()` removed** — these were deprecated aliases for
  `cuda::maximum<>()`/`cuda::minimum<>()` in CUDA 12.x's CCCL and are gone outright in CUDA 13's.
  Replaced with the direct equivalents (present in both CUDA 12 and 13) in
  `arch/arch_device_cuda.h`'s reduction kernel (`max`/`min` reduction ops), plus added
  `#include <cuda/functional>`. Unconditional — the replacement works on both CUDA versions.
- **`cudaMemPrefetchAsync`/`cudaMemAdvise` signature change** — CUDA 13 dropped the
  `(ptr, count, int device, stream)` overloads in favor of a `cudaMemLocation`-based signature.
  `arch/arch_device_cuda.h` has a version-gated `gpuMemPrefetchAsync` wrapper that keeps call
  sites on the old int-device form and translates to `cudaMemLocation` internally when
  `CUDART_VERSION >= 13000`. The matching fix in `submodules/hashinator` (`split_gpuMemPrefetchAsync`/
  `split_gpuMemAdvise` in `include/splitvector/archMacros.h`) is already committed in that
  submodule's own history (`a188c5f cuda13 compatibility`), not part of this repo's working diff.

## General `nvc++` compiler compatibility

Needed to compile the codebase's existing GPU code under `nvc++` at all.

- **`hashinator`'s `__NVCC__` → `__NVCC__ || __NVCOMPILER` compiler-detection fixes**, and its
  own CUDA13 compatibility work, are committed directly in that submodule's own history — see
  `submodules/hashinator`'s log, not this repo's diff.
- **`arch/arch_device_cuda.h`**: added `#include <cstdio>`/`#include <cstdlib>` — the file has
  several pre-existing `printf`/`exit()` calls that apparently relied on these being pulled in
  transitively under `nvcc`; `nvc++`'s CUDA/CCCL headers don't provide them the same way.
- **Eigen NEON `Complex.h`** (`submodules/eigen/Eigen/src/Core/arch/NEON/Complex.h`, uncommitted
  working-tree fix): widened an existing Clang-specific workaround
  (`EIGEN_COMP_CLANG || EIGEN_COMP_CASTXML`, for a known Eigen issue where `vld1q_u64` expands to
  a statement expression invalid at file/global scope) to also cover
  `defined(__NVCOMPILER_LLVM__)` — `nvc++` is LLVM-based and hits the identical problem, but
  doesn't define `__clang__` so wasn't previously caught by the guard.
- **`VLASIATOR_IF_DEVICE`/`VLASIATOR_IF_DEVICE_ONLY`/`VLASIATOR_IF_HOST_ONLY`** (new macros in
  `arch/arch_device_api.h`, `NV_IF_TARGET`-based): `nvc++` does not reliably honor a raw
  `#if defined(__CUDA_ARCH__)` guard inside a `__host__ __device__` function — both branches can
  end up compiled into device code regardless of the guard, which is a hard compile error
  whenever the host branch touches something host-only (a global, a non-device-callable
  function). This is a genuine `nvc++`-vs-`nvcc` behavioral difference, not stdpar-specific — it
  reproduces under plain classic-style `-x cu -gpu=ccXX -cuda` compilation. Applied at every site
  that actually needs dual host/device dispatch: `spatial_cells/velocity_block_container.h`,
  `spatial_cells/velocity_mesh_gpu.h` (including `arch::buf<T>::operator[]`, which — before this
  fix — silently compiled and returned the *host* pointer in device code, a live correctness
  bug, not a compile failure), and `velocity_mesh_parameters.h`'s `getMeshWrapper()`. Falls
  through to the original `__CUDA_ARCH__`/`__HIP_DEVICE_COMPILE__` check under classic
  `nvcc`/`hipcc`, so the `roihu_gpu` and LUMI/HIP builds are unaffected.
- **Range-based `for` loop → index-based loop inside `#pragma omp parallel for`**
  (`vdf_compression/compression.cpp`, `vlasovsolver/gpu_trans_map_amr.cpp`): `nvc++`'s OpenMP
  frontend enforces strict canonical loop form and hard-rejects a range-based `for` under
  `#pragma omp parallel for` (`NVC++-S-0155: #pragma pfor does not match canonical form`); GCC's
  `libgomp` tolerates it, which is why the classic build compiles unchanged. Converted to an
  index-based loop at each site.
- **`meshWrapperDevInstance` internal linkage** (`velocity_mesh_parameters.h`): this
  `__device__ __constant__` variable is declared at namespace scope in a header with the stated
  intent that "each compilation unit will use its own" — but was missing the `static` needed to
  actually enforce that. Classic `nvcc`'s non-RDC device linker tolerates same-named `__device__`
  globals across translation units; `nvc++`'s device linker (`nvdd`) merges them as one external
  symbol and fails with `nvlink error: Multiple definition of ...` once enough TUs including this
  header get linked together. This is a pre-existing bug (present before any of this effort),
  just never triggered under the classic build. Fixed by adding `static`.
- **`grid.cpp`'s `transferInParts`**: guarded a `#pragma omp parallel for` with `if (nIncoming >
  0) { ... }` before entering it. `nvc++`'s `nvomp` OpenMP runtime showed inconsistent
  (crash/hang) behavior for a zero-trip-count canonical `for` loop in this codebase; skipping the
  region entirely when there's nothing to do sidesteps it regardless of the exact root cause.
