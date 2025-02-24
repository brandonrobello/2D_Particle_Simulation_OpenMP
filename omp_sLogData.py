import numpy as np
import subprocess
import os

# Create output folder
out_folder = "omp_out"
os.makedirs(out_folder, exist_ok=True)

# Define the range for n_particles
n_min = 100
n_max = 100000
num_values = 10  # Number of values in the series

# Generate a semilog spaced series
n_particles_values = np.logspace(np.log10(n_min), np.log10(n_max), num=num_values, dtype=int)

# Define the list of thread counts to test
p_threads = [1, 2, 4, 8, 16, 32, 64, 96, 128, 160, 192, 224, 256]

for p_thread in p_threads:
    # Output file for storing simulation times and thread counts
    time_log_file = os.path.join(out_folder, f"simulation_times_{p_thread}.txt")

    # Run the commands directly using subprocess and capture output
    with open(time_log_file, "w") as log_file:
        # Write header for clarity
        log_file.write("n_particles threads simulation_time\n")
        
        for n in n_particles_values:
            # Prepare the command to run your executable
            command = ["./build_omp/openmp", "-n", str(n), "-s", "21"]

            # Set the environment variable OMP_NUM_THREADS for this run
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(p_thread)
            
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
            
            # Print stdout for debugging
            print(f"STDOUT for n={n} with {p_thread} threads:\n", result.stdout)
            
            # Initialize variables for extraction
            sim_time = None
            num_threads = None
            
            # Extract simulation time and thread count from stdout
            for line in result.stdout.split("\n"):
                if "Simulation Time" in line:
                    parts = line.split("=")
                    if len(parts) > 1:
                        sim_time = parts[1].split()[0]
                if "Number of threads" in line:
                    parts = line.split("=")
                    if len(parts) > 1:
                        num_threads = parts[1].split()[0]
                        
            # Fallback to the environment variable if not printed by the program
            if num_threads is None:
                num_threads = env.get("OMP_NUM_THREADS", "N/A")
            
            # Log the results (simulation time may be "N/A" if not found)
            if sim_time is not None:
                log_file.write(f"{n} {num_threads} {sim_time}\n")
            else:
                log_file.write(f"{n} {num_threads} N/A\n")

print("All jobs have been executed and logged.")
