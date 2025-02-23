#include "common.h"
#include <omp.h>
#include <cmath>
#include <vector>
#include <cstdio>

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

// ================== Left Alone ==============================

// Define bin size
const double bin_size = 1.2 * cutoff;

// Bin struct to hold particles
struct Bin {
    std::vector<particle_t*> particles;
	omp_lock_t lock;
};

// 2D grid of bins for two frames
std::vector<std::vector<Bin>> bins_frame_1, bins_frame_2;
std::vector<std::vector<Bin>>* current_bins;
std::vector<std::vector<Bin>>* next_bins;
int bin_count_x, bin_count_y;

// Compute bin index for a given position
void get_bin_index(double x, double y, int& bx, int& by, double size) {
    // Apply the same reflection logic as move()
    if (x < 0) x = -x;
    if (x > size) x = 2 * size - x;
    if (y < 0) y = -y;
    if (y > size) y = 2 * size - y;

    // Compute bin indices
    bx = static_cast<int>(x / bin_size);
    by = static_cast<int>(y / bin_size);

    // Ensure bin indices stay within valid range
    bx = std::max(0, std::min(bx, bin_count_x - 1));
    by = std::max(0, std::min(by, bin_count_y - 1));
}


// Initialize bins based on the simulation domain size
void init_simulation(particle_t* parts, int num_parts, double size) {
    bin_count_x = static_cast<int>(size / bin_size);
    bin_count_y = static_cast<int>(size / bin_size);

    // Resize and initialize bins
    bins_frame_1.assign(bin_count_x, std::vector<Bin>(bin_count_y));
    bins_frame_2.assign(bin_count_x, std::vector<Bin>(bin_count_y));
    
    current_bins = &bins_frame_1;
    next_bins = &bins_frame_2;

    // Ensure all bins exist before assignment
    for (int i = 0; i < bin_count_x; ++i) {
        for (int j = 0; j < bin_count_y; ++j) {
            bins_frame_1[i][j] = Bin();
            bins_frame_2[i][j] = Bin();
        }
    }

	// Initialize locks
	for (int i = 0; i < bin_count_x; ++i) {
		for (int j = 0; j < bin_count_y; ++j) {
			omp_init_lock(&bins_frame_1[i][j].lock);
			omp_init_lock(&bins_frame_2[i][j].lock);
		}
	}

    // Assign particles to their bins in the initial frame
    for (int i = 0; i < num_parts; ++i) {
        int bx, by;
        get_bin_index(parts[i].x, parts[i].y, bx, by, size);
        (*current_bins)[bx][by].particles.push_back(&parts[i]);
    }
}

// Apply force with binning for all neighbors of particles
void apply_force_binning(particle_t& particle, Bin& bin) {
    for (auto neighbor : bin.particles) {
        if (neighbor != &particle) {
            apply_force(particle, *neighbor);
        }
    }
}

// Optimized force calculation considering only necessary bins
void compute_forces() {
	#pragma omp for collapse(2)
    for (int i = 0; i < bin_count_x; ++i) {
        for (int j = 0; j < bin_count_y; ++j) {
            Bin& bin = (*current_bins)[i][j];

            for (auto p : bin.particles) {
            // Reset acceleration before force accumulation
            p->ax = p->ay = 0;
                // Check self-bin and neighboring bins within cutoff
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        int ni = i + dx;
                        int nj = j + dy;
                        if (ni >= 0 && ni < bin_count_x && nj >= 0 && nj < bin_count_y) {
                            apply_force_binning(*p, (*current_bins)[ni][nj]);
                        }
                    }
                }
            }
        }
    }
}

// Simulate one step with binning using double buffering
void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Compute Forces using optimized binning
    compute_forces();

    // Reset next_bins
	#pragma omp for collapse(2)
    for (int i = 0; i < bin_count_x; ++i) {
        for (int j = 0; j < bin_count_y; ++j) {
            (*next_bins)[i][j].particles.clear();
        }
    }

	// Move particles and reassign to new bins
	#pragma omp for collapse(2)
	for (int bx = 0; bx < bin_count_x; ++bx) {
		for (int by = 0; by < bin_count_y; ++by) {
			// Process each particle in this bin.
			Bin& current_bin = (*current_bins)[bx][by];
			for (auto p : current_bin.particles) {
				move(*p, size);
				int new_bx, new_by;
				get_bin_index(p->x, p->y, new_bx, new_by, size);

				// Lock the next_bin and update
				omp_set_lock(&((*next_bins)[new_bx][new_by].lock));
				(*next_bins)[new_bx][new_by].particles.push_back(p);
				omp_unset_lock(&((*next_bins)[new_bx][new_by].lock));
			}
		}
	}

    // Swap frames
    #pragma omp single
    std::swap(current_bins, next_bins);
}
