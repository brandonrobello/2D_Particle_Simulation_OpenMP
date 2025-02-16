import numpy as np
import subprocess
import os

# Create output folder
## CHANGE OUTFOLDER TO DESIRED FOLDER
out_folder = "serial_out_2"
os.makedirs(out_folder, exist_ok=True)

# Define the range for n_particles
n_min = 100
n_max = 100000
num_values = 10  # Number of values in the series

# Generate a semilog spaced series
n_particles_values = np.logspace(np.log10(n_min), np.log10(n_max), num=num_values, dtype=int)

# Output file for storing simulation times
time_log_file = os.path.join(out_folder, "simulation_times.txt")

# Run the commands directly using subprocess and capture output
with open(time_log_file, "w") as log_file:
    for n in n_particles_values:
        command = ["./build_2/serial", "-n", str(n), "-o", f"{out_folder}/serial_{n}.out", "-s", "21"]
        result = subprocess.run(command, capture_output=True, text=True)
        
        # Print stdout for debugging
        print(f"STDOUT for n={n}:", result.stdout)
        
        # Extract simulation time from stdout
        for line in result.stdout.split("\n"):
            if "Simulation Time" in line:
                parts = line.split("=")
                if len(parts) > 1:
                    time_value = parts[1].split()[0]
                    log_file.write(f"{n} {time_value}\n")
                    break

print("All serial jobs have been executed and logged.")
