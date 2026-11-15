# Performance and Memory Analysis

The project remains a sequential C++17 baseline. Benchmark timings are collected with `std::chrono::steady_clock` and are exported to `benchmarks/`.

## Memory model

The main grid storage is contiguous:

- `u`, `uPrev`, `uWeight`: `3 * (NX + 1) * NY` floats
- `v`, `vPrev`, `vWeight`: `3 * NX * (NY + 1)` floats
- pressure buffers and divergence: `4 * NX * NY` floats
- cell types: `NX * NY` bytes
- particles: `4` floats per particle

The approximate grid memory is therefore:

```text
4 * (3*(NX+1)*NY + 3*NX*(NY+1) + 4*NX*NY) + NX*NY bytes
```

This excludes allocator overhead and temporary I/O buffers. The pressure solver uses its existing double buffers and does not allocate inside the Jacobi loop. Velocity extrapolation updates the existing arrays in place.

Approximate field plus particle storage for the standard benchmark sizes:

| Grid | Particles | Estimated storage |
|---|---:|---:|
| 64x64 | 3,225 | 216 KiB |
| 128x128 | 13,137 | 864 KiB |
| 256x256 | 53,504 | 3.38 MiB |

These are lower-bound estimates for the simulation arrays; the process working set will be larger because of the executable, C++ runtime, vector capacity, and operating-system allocations.

## Current benchmark baseline

The benchmark command is:

```powershell
.\build\fluid_sim.exe --benchmark
```

It writes one per-step CSV per grid and appends aggregate rows to `benchmarks/benchmark_summary.csv`. The most recent measured baseline should be read from that CSV rather than copied into documentation, because timings depend on CPU, compiler, build type, and operating-system load.

The pressure solve is expected to dominate runtime because each Jacobi iteration scans the fluid cells. Particle-to-grid and grid-to-particle stages scale with particle count. The current implementation is intentionally sequential and contains no OpenMP, SIMD, CUDA, or other parallel optimization.

## Recommended comparison procedure

1. Configure and build a Release binary.
2. Run the benchmark on an otherwise idle machine.
3. Repeat each resolution at least three times.
4. Compare total time, average step time, and stage percentages from the CSV files.
5. Record compiler version, CMake version, grid size, particle count, pressure iterations, tolerance, and operating system.

For memory profiling, use Windows Task Manager or a process profiler while running a fixed configuration. Compare peak working-set memory against the formula above and record the result with the benchmark CSV.
