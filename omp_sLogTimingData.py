import numpy as np
import subprocess
import os

# Create output folder
out_folder = "omp_out_timing"
os.makedirs(out_folder, exist_ok=True)

# Define the range for n_particles
n_min = 100
n_max = 100000
num_values = 10 

# Generate a semilog spaced series
n_particles_values = np.logspace(np.log10(n_min), np.log10(n_max), num=num_values, dtype=int)

# Define the list of thread counts to test
p_threads = [1, 2, 4, 8, 16, 32, 64, 96, 128, 160, 192, 224, 256]

# Output file for capturing the full timing data
timing_data_file = os.path.join(out_folder, "timing_data.txt")

with open(timing_data_file, "w") as f:
    # Write header
    f.write("n_particles threads simulation_time total_time computation_time synchronization_time communication_time\n")
    
    for p_thread in p_threads:
        for n in n_particles_values:
            # Prepare the command to run your executable
            command = ["./build_omp/openmp", "-n", str(n), "-s", "21"]

            # Set the environment variable OMP_NUM_THREADS for this run
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(p_thread)
            
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                    universal_newlines=True, env=env)
            
            # Debug print
            print(f"STDOUT for n={n} with {p_thread} threads:\n", result.stdout)
            
            # Initialize timing variables
            sim_time = None
            total_time = None
            comp_time = None
            sync_time = None
            comm_time = None
            
            # Process output lines to extract timing information
            for line in result.stdout.split("\n"):
                if "Simulation Time" in line:
                    # Example: "Simulation Time = 0.12345 seconds for 1000 particles."
                    parts = line.split("=")
                    if len(parts) > 1:
                        sim_time = parts[1].split()[0]
                if "time Data:" in line:
                    # Example: "time Data: total = 0.456, computation = 0.123, synchronization = 0.078, communcation = 0.255"
                    data_str = line.replace("time Data:", "").strip()
                    parts = data_str.split(",")
                    for part in parts:
                        key_value = part.strip().split("=")
                        if len(key_value) == 2:
                            key = key_value[0].strip()
                            value = key_value[1].strip()
                            if key == "total":
                                total_time = value
                            elif key == "computation":
                                comp_time = value
                            elif key == "synchronization":
                                sync_time = value
                            elif key == "communcation":  # note: using the same misspelling as printed
                                comm_time = value
            
            # Write the results to the timing data file
            if sim_time is not None:
                f.write(f"{n} {p_thread} {sim_time} {total_time} {comp_time} {sync_time} {comm_time}\n")
            else:
                f.write(f"{n} {p_thread} N/A N/A N/A N/A N/A\n")

print("All jobs have been executed and timing data logged.")
