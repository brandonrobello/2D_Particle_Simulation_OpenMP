import numpy as np
import subprocess
import os

# Create output folder
out_folder = "results"
os.makedirs(out_folder, exist_ok=True)

# Define the range for n_particles
n_min = 1000
n_max = 1000000
num_values = 10  # Number of values in the series

# Generate a semilog spaced series
n_particles_values = np.logspace(np.log10(n_min), np.log10(n_max), num=num_values, dtype=int)

# Output file for storing simulation times
time_log_file = os.path.join(out_folder, "simulation_times.txt")

# Run the commands directly using subprocess and capture output
with open(time_log_file, "w") as log_file:
    log_file.write("Serial:\n")
    for n in n_particles_values:
        print(f"Serial {n}")

        command = ["./build/serial", "-n", str(n), "-s", "21"]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        # result = subprocess.run(command, capture_output=True, text=True)
        
        # Extract simulation time from stdout
        for line in result.stdout.split("\n"):
            if "Simulation Time" in line:
                parts = line.split("=")
                if len(parts) > 1:
                    time_value = parts[1].split()[0]
                    log_file.write(f"{n} {time_value}\n")
                    break
        
        print(f"Serial {n} {time_value}")

    thread_counts = [1, 4, 16, 64, 256]
    for threads in thread_counts:
        log_file.write(f"Parallel ({threads} threads):\n")
        for n in n_particles_values:
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(threads)

            print(f"Parallel {threads} {n}")

            command = ["./build/openmp", "-n", str(n), "-s", "21"]
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
            # result = subprocess.run(command, capture_output=True, text=True, env=env)
            
            # Extract simulation time from stdout
            for line in result.stdout.split("\n"):
                if "Simulation Time" in line:
                    parts = line.split("=")
                    if len(parts) > 1:
                        time_value = parts[1].split()[0]
                        log_file.write(f"{n} {time_value}\n")
                        break

            print(f"Parallel {threads} {n} {time_value}")


# Output file for storing simulation times
time_log_file = os.path.join(out_folder, "strong_weak_times.txt")

# Run the commands directly using subprocess and capture output
with open(time_log_file, "w") as log_file:
    for threads in [1, 2, 4, 8, 16, 32, 64, 128, 256]:
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(threads)

        # Strong
        for n in [1000, 10000, 100000, 1000000]:
            log_file.write(f"Strong ({threads} threads, {n} particles): ")

            print(f"Parallel {threads} {n}")
            
            command = ["./build/openmp", "-n", str(n), "-s", "21"]
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
            # result = subprocess.run(command, capture_output=True, text=True, env=env)
            
            # Extract simulation time from stdout
            for line in result.stdout.split("\n"):
                if "Simulation Time" in line:
                    parts = line.split("=")
                    if len(parts) > 1:
                        time_value = parts[1].split()[0]
                        log_file.write(f"{time_value}\n")
                        break
            
            print(f"Parallel {threads} {n} {time_value}")
        
        # Weak
        for parts_per_thread in [1000, 2000, 5000, 10000]:
            log_file.write(f"Weak ({threads} threads, {parts_per_thread} parts per thread): ")

            n = parts_per_thread * threads

            print(f"Parallel {threads} {parts_per_thread} {n}")

            command = ["./build/openmp", "-n", str(n), "-s", "21"]
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
            # result = subprocess.run(command, capture_output=True, text=True, env=env)
            
            # Extract simulation time from stdout
            for line in result.stdout.split("\n"):
                if "Simulation Time" in line:
                    parts = line.split("=")
                    if len(parts) > 1:
                        time_value = parts[1].split()[0]
                        log_file.write(f"{time_value}\n")
                        break
            
            print(f"Parallel {threads} {parts_per_thread} {n} {time_value}")

