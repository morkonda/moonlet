import matplotlib.pyplot as plt

with open('position_moonlet.txt','r') as particle_data_file:
	data = particle_data_file.readlines()

with open('position_x.txt','r') as file_to_find_number_of_particles:
	number_data = file_to_find_number_of_particles.readlines()

n				= len(number_data)
x_data 				= []
y_data 				= []
z_data 				= []
time_data 			= [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.]

for line in data:
	position_velocity = line.split()
	x_data.append(float(position_velocity[0]))
	y_data.append(float(position_velocity[1]))
	z_data.append(float(position_velocity[2]))

plt.plot(time_data, x_data, 'r', label="x position")
plt.plot(time_data, y_data, 'ro', label="y position")
plt.plot(time_data, z_data, 'r^', label="z position")
plt.title("Plot of x, y and z positions of the moonlet versus time "+str(n)+" particles")
plt.xlabel("Time (orbits)")
plt.ylabel("Position (m)")
plt.legend(loc='lower left')
plt.show()
