//Author: Michael Brothers
//Date: Nov 5 2013
//Class: ME5613: Advanced Fluid Mechanics
//Title: Problem 4: 2D ideal flow numeric simulation

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void load_boundary_conditions(double* boundary_grid,
                              double* grid_mask,
                              char* bc_input_filename,
                              char* interp_input_filename,
                              int* number_of_rows,
                              int* number_of_columns);


void calculate_linear_interpolation(const int length_of_interpolation,
                                    double** interpolation,
                                    const double min_value,
                                    const double max_value);

void compute_stream_values_in_rectangular_region(const int first_column,
                                                 const int last_column,
                                                 const int first_row,
                                                 const int last_row,
                                                 const int number_of_rows,
                                                 const int number_of_columns,
                                                 const double* problem_grid_old,
                                                 double* problem_grid_new,
                                                 const double* boundary_grid,
                                                 const double* grid_mask);

double compute_root_mean_square_error(const double* problem_grid_old,
                                    const double* problem_grid_new,
                                    const double number_of_interior_points,
                                    const int number_of_rows,
                                    const int number_of_columns);




