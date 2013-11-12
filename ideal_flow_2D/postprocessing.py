#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys

pressure = sys.argv[1]

def load(fname):
	''' load file using std open '''
        f = open(fname, 'r')

	data = []
        for line in f.readlines():
		data.append(line.replace('\n','').replace(',', ' ').strip(' ').split(' '))

        f.close()

	return np.array(np.concatenate(data), dtype=np.float64)

def main():
	#output files will be loaded into 1-D arrays left to right, top to bottom
	#grid_shape_data = np.loadtxt(boundary_conditions_filename, delimiter = ' ,', unpack = True)

	pressure_values = load(pressure);	

	print pressure_values
	print " "

	print pressure_values[0:pressure_values.size:3].size 
        print pressure_values[1:pressure_values.size:3].size 
        print pressure_values[2:pressure_values.size:3].size

	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.scatter(pressure_values[0:66:3], pressure_values[1:66:3], pressure_values[2:66:3], label='parametric curve')
	ax.legend()
	plt.show()

if __name__ == "__main__":
	main()
