import yt
import glob
import matplotlib.pyplot as plt

# 1. Get a sorted list of all your plotfiles
file_list = sorted(glob.glob("hdf5/Plt56_*.3d.hdf5"))

times = []
avg_lapses = []

print(f"Found {len(file_list)} files. Starting extraction...")

for f in file_list:
    # Load the dataset
    # use 'hint' to speed up loading if yt doesn't recognize the format immediately
    ds = yt.load(f)
    
    # Store the simulation time (usually ds.current_time)
    times.append(float(ds.current_time))
    
    # Calculate volume-weighted average of the lapse (alpha)
    ad = ds.all_data()
    mean_val = ad.quantities.weighted_average_quantity("alpha", "cell_volume")
    
    avg_lapses.append(float(mean_val))
    print(f"Processed {f}: Time = {ds.current_time:.2f}, Avg Lapse = {mean_val:.4f}")

# 3. Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(times, avg_lapses, marker='o', linestyle='-', color='b')
plt.xlabel('Simulation Time ($t$)')
plt.ylabel('Volume-Averaged Lapse $\langle \\alpha \\rangle$')
plt.title('Evolution of Average Lapse Over Time')
plt.grid(True)
plt.show()
