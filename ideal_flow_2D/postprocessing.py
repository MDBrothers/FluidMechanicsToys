#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys

node_designation_filename = sys.argv[1]
stream_function_output_filename = sys.argv[2]
boundary_conditions_filename = sys.argv[3]

def load(fname):
	''' load file using std open '''
        f = open(fname, 'r')

	data = []
        for line in f.readlines():
		data.append(line.replace('\n','').replace(',', ' ').strip(' ').split(' '))

        f.close()

	return np.concatenate(data)

def main():
	#output files will be loaded into 1-D arrays left to right, top to bottom
	#grid_shape_data = np.loadtxt(boundary_conditions_filename, delimiter = ' ,', unpack = True)

	grid_shape_data = load(boundary_conditions_filename)

	print grid_shape_data
	stream_function_values = np.array(load(stream_function_output_filename));
	evaluation_nodes = np.array(load(node_designation_filename));	

	uniform_velocity = grid_shape_data[0]
	grid_spacing = grid_shape_data[1]
	numrows = grid_shape_data[2]
	numcols = grid_shape_data[3]

	u_velocity_values = np.zeros(stream_function_values.size())
	v_velocity_values = np.zeros(stream_function_values.size())
	pressure_values = np.zeros(stream_function_values.size())
	
	for entry in range(boundary_nodes.size()):
		if boundary_nodes[entry] == 1.0:
			if stream_function_values[entry] != stream_function_values[entry - 1]:
				v_velocity_values[entry] = 3214234.987

if __name__ == "__main__":
	main()
