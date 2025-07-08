# Particle Simulation Optimization

## Contributions

* **Write-Up**: All
* **Brandon**: Serial optimization design, implementations, and graphing; computation vs. synchronization
* **Kevin**: OpenMP implementation, further serial optimizations, graphing, and benchmarking parallel

---

## Introduction

This assignment focuses on optimizing a 2D particle simulation governed by a repulsive force equation. The project had two primary goals:

1. Optimize the naïve serial implementation, which originally ran in `O(n^2)` time.
2. Parallelize the optimized solution using OpenMP for shared-memory architectures.

---

## Serial Optimization Approach

The naïve version used a brute-force method that computed forces between all particle pairs, leading to `O(n^2)` complexity.

### Optimization Strategy

To improve performance, we:

* **Implemented spatial binning** to reduce computations to `O(n)`
* **Divided the simulation space into bins**, only considering particles within neighboring bins
* Introduced a `Bin` struct to group `particle_t*` by location
* Used **double buffering** with `current_frame` and `next_frame` to manage particle states across time steps

This binning approach ensures constant-time per-particle computation and allows easy parallelization.

### Performance

* **Naïve implementation** showed quadratic scaling with particle count
* **Optimized serial version** scaled linearly with the number of particles (`O(n)`), validating the effectiveness of binning

---

## Parallel Implementation with OpenMP

The optimized serial version was parallelized using OpenMP. Key strategies included:

* Adding `#pragma omp for` to loops over bins
* Splitting work for:

  * Force computation
  * Clearing next frame
  * Moving particles
  * Reassigning particles to bins

### Synchronization Strategy

* **No sync needed** when threads modify only local particle data
* **Locks required** when particles are reassigned across bins to avoid race conditions
* `#pragma omp single` used for frame-swapping at each time step

### Results

* **1-thread OpenMP** version performed worse due to overhead
* **64-thread version** performed well, except for small particle counts where overhead dominates
* **256-thread version** degraded due to context-switching overhead on AMD EPYC 7763 CPUs (64 cores, 128 threads)

---

## Scaling Behavior

### Strong Scaling

* Measures speedup for a fixed particle count as thread count increases
* Ideal slope: `-1` on a log-log plot
* Midrange thread counts (4–32) approached ideal
* Extremes (1–2 and 64–128 threads) showed overhead effects

### Weak Scaling

* Increases both particle count and thread count proportionally
* Ideal slope: `0` on log-log scale
* Midrange performance matched expectations; extremes diverged due to overhead and synchronization

---

## Computation vs Synchronization Time

* **Computation** scaled efficiently with threads (∼ `T/p`)
* **Synchronization** did not scale as well due to:

  * Lock contention during bin reassignment
  * Thread management overhead

### Insights

* There's an optimal balance between thread count and particle count
* Too many threads can **reduce performance** due to sync overhead
* Further improvements could include:

  * Fine-grained locking
  * Cache-aware data layout
  * Better bin iteration strategies

---

## Conclusion

We successfully:

* Reduced complexity from `O(n^2)` to `O(n)` using spatial binning
* Leveraged OpenMP for parallel speedup
* Identified trade-offs between computation and synchronization
* Analyzed strong and weak scaling behavior

This framework provides a scalable foundation for high-performance particle simulations.
