import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.optimize import curve_fit
import math

#Initialisation
result_boxsize_1core = []
result_time_1core = []


#Normalisation using boxsize = 100.0
N_tot_avg = 1475
c = N_tot_avg/(100.*100.)


#Data from using 1 core
boxsize_1core = [100.0, 110.00000000000001, 120.00000000000001, 130.0, 140.0, 150.0, 160.0, 170.00000000000003, 180.0, 190.0]

time_1core = [176.43241691589355, 206.67291402816772, 242.3817000389099, 282.6246249675751, 411.8980619907379, 435.7475280761719, 502.293270111084, 582.0376679897308, 659.9761738777161, 817.1540269851685]


#Calculating number of particles
for i in boxsize_1core:
        result_boxsize_1core.append(c*i*i)

#Calculating time for 1 orbit
for i in time_1core:
        result_time_1core.append(i/10.0)


#Data from using 4 cores
boxsize_4cores = [1475.0, 1784.7500000000005, 2124.0000000000005, 2492.75, 2891.0000000000005, 3318.75, 3776.0, 4262.750000000001, 4779.0, 5324.75]

time_4cores = [8.629414200782776, 10.309395098686219, 12.731140494346619, 15.663812708854675, 17.607955598831175, 20.78720109462738, 25.354449701309203, 28.570127391815184, 30.524536085128783, 35.46618208885193]


#writing to a file
a1 = np.array(result_boxsize_1core)
b1 = np.array(result_time_1core)
a4 = np.array(boxsize_4cores)
b4 = np.array(time_4cores)

output = open('plot_data_1core_4cores.data', 'w')
for i in np.array(xrange(len(a1))):
	output.write("%f\t%f\t%f\t%f\n" %(a1[i], b1[i], a4[i], b4[i]))

output.close()

#Plot of N vs t (seconds/orbit) for 1 core & 4 cores 
plt.plot(result_boxsize_1core, result_time_1core, 'r^')
plt.plot(boxsize_4cores, time_4cores, 'g^')
plt.xlabel("Number of particles")
plt.ylabel("Time (seconds/orbit)")
plt.title("Plot of computation time for 1 orbit vs Number of Particles \n red: 1 core \n green: 4 cores")
plt.show()

