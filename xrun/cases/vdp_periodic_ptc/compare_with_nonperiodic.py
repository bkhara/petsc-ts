import numpy as np
import matplotlib.pyplot as plt

# Path to your CSV file
csv_file_periodic = './solution.txt'
csv_file_non_periodic = '/home/khara/projects/petsc-ts/xrun/results/vdp_bvp_ptc/solution.txt'

# Load the CSV data: assume 3 numeric columns, no header
data_periodic = np.loadtxt(csv_file_periodic, delimiter=',', skiprows=1)

# Extract columns
t = data_periodic[:, 0]
z = data_periodic[:, 1:3]

# Load the CSV data: assume 3 numeric columns, no header
data_non_periodic = np.loadtxt(csv_file_non_periodic, delimiter=',', skiprows=1)

# Extract columns
t_np = data_non_periodic[:, 0]
z_np = data_non_periodic[:, 1:3]

# Plot
# time series plot
fig, axs = plt.subplots(nrows=2,ncols=1,figsize=(16,8),subplot_kw={'aspect': 'auto'}, sharex=False, sharey=False, squeeze=True)
axs[0].plot(t,z[:,0],label=r'$y_1$_per')
axs[0].plot(t_np,z_np[:,0],label=r'$y_1$_nper')
axs[0].set_xlabel('time')
axs[0].legend()
axs[1].plot(t,z[:,1],label=r'$y_2$_per')
axs[1].plot(t_np,z_np[:,1],label=r'$y_2$_nper')
axs[1].set_xlabel('time')
axs[1].legend(loc='best')
plt.savefig('plot_time-series-comp.png')

# time series AND phase portrait
fig, axs = plt.subplots(nrows=2,ncols=2,figsize=(20,10),subplot_kw={'aspect': 'auto'}, sharex=False, sharey=False, squeeze=True)
axs[0,0].plot(t,z[:,0],label=r'$y_1$_per')
axs[0,0].plot(t_np,z_np[:,0],label=r'$y_1$_nper')
axs[0,0].set_xlabel('time')
axs[0,0].legend()
axs[0,1].plot(t,z[:,1],label=r'$y_2$_per')
axs[0,1].plot(t_np,z_np[:,1],label=r'$y_2$_nper')
axs[0,1].set_xlabel('time')
axs[0,1].legend(loc='best')
axs[1,0].plot(z[:,0], z[:,1], label='Ph. port.(per.)')
axs[1,0].plot(z_np[:,0], z_np[:,1], '--', label='Ph. port.(non-per.)')
axs[1,0].legend()
# axs[1,1].plot(freqs1, power1,'b*-',label='FFT: '+r'$y_1$')
# axs[1,1].plot(freqs2, power2,'r*--',label='FFT: '+r'$y_2$')
plt.savefig('plot_results-comp.png')