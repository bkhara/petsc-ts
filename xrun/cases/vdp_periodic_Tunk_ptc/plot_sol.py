import numpy as np
import matplotlib.pyplot as plt

# Path to your CSV file
csv_file = 'solution.txt'

# Load the CSV data: assume 3 numeric columns, no header
data = np.loadtxt(csv_file, delimiter=',', skiprows=1)

# Extract columns
t = data[:,0]
z = data[:,1:4]

# Plot
# time series plot
fig, axs = plt.subplots(nrows=2,ncols=1,figsize=(16,8),subplot_kw={'aspect': 'auto'}, sharex=False, sharey=False, squeeze=True)
axs[0].plot(t,z[:,0],'k-',label=r'$y_1$')
axs[0].set_xlabel('time')
axs[0].legend()
axs[1].plot(t,z[:,1],'m--',label=r'$y_2$')
axs[1].set_xlabel('time')
axs[1].legend(loc='best')
plt.savefig('plot_time-series.png')

# time series AND phase portrait
fig, axs = plt.subplots(nrows=2,ncols=2,figsize=(20,10),subplot_kw={'aspect': 'auto'}, sharex=False, sharey=False, squeeze=True)
axs[0,0].plot(t,z[:,0],'k-',label=r'$y_1$')
axs[0,0].set_xlabel('time')
axs[0,0].legend()
axs[0,1].plot(t,z[:,1],'m--',label=r'$y_2$')
axs[0,1].set_xlabel('time')
axs[0,1].legend(loc='best')
axs[1,0].plot(z[:,0], z[:,1], 'g-', label='Phase portrait')
axs[1,0].legend()
axs[1,1].plot(t,z[:,2],'b--',label=r'$y_3$')
axs[1,1].set_xlabel('time')
axs[1,1].legend(loc='best')
# axs[1,1].plot(freqs1, power1,'b*-',label='FFT: '+r'$y_1$')
# axs[1,1].plot(freqs2, power2,'r*--',label='FFT: '+r'$y_2$')
plt.savefig('plot_results.png')