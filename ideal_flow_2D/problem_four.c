//Author: Michael Brothers
//Date: Nov 5 2013
//Class: ME5613: Advanced Fluid Mechanics
//Title: Problem 4: 2D ideal flow numeric simulation

#include "./problem_four.h"
const double tolerance = 1.E-6;
const int maxiters = 100;
const double rho = 1.2041;
const double inlet_pressure = 101.25*1000.0;

int main(int argc, char* argv[])
{
    double* grid_mask;
    double* boundary_grid;
    double* old_grid;
    double* new_grid;
    double** interpolation;
    double grid_spacing = 0.0;
    double inlet_velocity = 0.0;
    double pos_d_y, pos_d_x, neg_d_y, neg_d_x;
    double pressure, designation, u_vel, v_vel;

    int number_of_rows = 0;
    int number_of_columns = 0;
    int number_of_interior_points = 0;
    int column, row;

    int catcher;
    

    //printf("%d%s", argc, "\n");
    //printf("%s%s", argv[1], "\n");
    //printf("%s%s", argv[2], "\n");

    char *buf = malloc(256*sizeof(char));
    FILE *infile = fopen(argv[1], "r");
    char* token;

    if (infile == NULL) 
    {
         printf("%s", "failed to open file\n");
    }
    
    //Sucessfull opened a file with name matching argv[1] for loading boundary conditions
    else 
    {
	catcher = fscanf(infile, "%lf", &inlet_velocity);
	catcher = fscanf(infile, "%lf", &grid_spacing);
        catcher = fscanf(infile, "%d %d", &number_of_rows, &number_of_columns);

        //printf("%d %s", number_of_rows, "rows\n");
        //printf("%d %s", number_of_columns, "columns\n");

        //dynamically resize arrays to fit problem specified in input file
        boundary_grid = malloc(number_of_rows*number_of_columns*sizeof(double));
        grid_mask = malloc(number_of_rows*number_of_columns*sizeof(double));
        old_grid = malloc(number_of_rows*number_of_columns*sizeof(double));
        new_grid = malloc(number_of_rows*number_of_columns*sizeof(double));

        //We address the array in (row, column) order
        //When results are plotted row is the y axis
        //and column is the x axis, this makes it easier 
        //for me to understand.
        int row = number_of_rows -1;
        column = 0;
       
        //obtain a line from the boundary condition file
        while(fscanf(infile, "%s", buf) != EOF)
        {
            //we need to tokenize the line so that we can parse out the numerical information
            //we define a token as a floating point number, delimited by a comma
            token = strtok(buf, ",");
            column = 0;

            //while not out of tokens, do
            while(token != NULL)
            {
                //process the token as a floating point number
                sscanf(token, "%lf", &boundary_grid[number_of_columns*row + column]);

                old_grid[number_of_columns*row + column] = 1.0;
                new_grid[number_of_columns*row + column] = 0.0;

                //set the grid mask to the appropriate value
                //grid mask de-activates the modification of boundary nodes
                if(boundary_grid[number_of_columns*row + column] < 0.0)
                {
                    grid_mask[number_of_columns*row + column] = 1.0;
                    boundary_grid[number_of_columns*row + column] = 0.0;
                    number_of_interior_points++;
                }
                else
                    grid_mask[number_of_columns*row + column] = 0.0;

                //printf("Row %d, Col %d\n", row, column);

                //get next token in the line we just read from the file
                token = strtok(NULL, ",");

                column++;
            }
            row--;
        }
    }

    //use stream infile to look at another file, the interpolation directions
    fclose(infile);
    infile = fopen(argv[2], "r");

    if (infile == NULL) 
    {
         fprintf(stderr, "%s", buf);
         exit(1);
    }

    //We sucessfully opened a file that matches the name stored in argv[2], the interpolation
    //directions file
    else 
    {
        //these variables will store information about
        //the first and last node in a linear interpolated set
        int first_row = 0;
        int first_column = 0;
        int last_row = 0;
        int last_column = 0;
        int num_interp_rows = 0;
        int num_interp_columns = 0;
        double min_value = 0.0;
        double max_value = 0.0;

        int entry = 0;

        //resize an array of pointers to addresses in boundary_grid
        //modifying these is just like modifying the original array
        if(number_of_columns > number_of_rows)
            interpolation = malloc(number_of_columns*sizeof(double*));
        else
            interpolation = malloc(number_of_rows*sizeof(double*));

        //read the interpolation directions file line by line and produce the interpolations
        while(fscanf(infile, "%d %d %d %d", &first_row, &first_column, &last_row, &last_column) != EOF){
            //These lines really must be vertical or horizontal
            min_value = boundary_grid[number_of_columns*first_row + first_column];
            max_value = boundary_grid[number_of_columns*last_row + last_column];
            num_interp_rows = last_row - first_row;
            num_interp_columns = last_column - first_column;

            //Since the interpolation node sets are either vertical or horizontal, were increment in
            //only one of those two directions
            if(num_interp_columns > 0)
            {
                for(entry = 0; entry < num_interp_columns; entry ++)
                {
                    interpolation[entry] = &boundary_grid[number_of_columns*first_row + first_column + entry];
                }
                    calculate_linear_interpolation(num_interp_columns, interpolation, min_value, max_value);
            }

            if(num_interp_rows > 0)
            {
                for(entry = 0; entry < num_interp_rows; entry ++)
                {
                    interpolation[entry] = &boundary_grid[number_of_columns*(first_row + entry) + first_column];
                }
                    calculate_linear_interpolation(num_interp_rows, interpolation, min_value, max_value);
            }
        }
    }

    fclose(infile);
    
    //perform the simulation according to a given convergence criterion specified at the top of this
    //file
    double error = 999.0;
    int numiters = 0;
    int entry = 0;
    int maxentry = number_of_rows*number_of_columns;

    while((error > tolerance) && (numiters < maxiters))
    {
        compute_stream_values_in_rectangular_region(0, //const int first_column,
                                                 number_of_columns, //const int last_column,
                                                 0, //const int first_row,
                                                 number_of_rows, //const int last_row,
                                                 number_of_rows, //const int number_of_rows,
                                                 number_of_columns, //const int number_of_columns,
                                                 old_grid, //const double* problem_grid_old,
                                                 new_grid, //double* problem_grid_new,
                                                 boundary_grid, //const double* boundary_grid,
                                                 grid_mask); //const double* grid_mask)

        error = compute_root_mean_square_error(old_grid, //const double* problem_grid_old,
                                       new_grid, //const double* problem_grid_new,
                                       number_of_interior_points, //const double number_of_interior_points,
                                       number_of_rows, //const int number_of_rows,
                                       number_of_columns);//const int number_of_columns)

        for(entry = 0; entry < maxentry; entry ++)
        {
            old_grid[entry] = new_grid[entry];
        }
                               
        numiters++;
   }

    printf("Iteration %d complete!\nThe error is %lf\n", numiters, error);
    printf("See \"stream_output.txt\" for stream function values\n");
    infile = fopen("./stream_output.txt", "w");

    //display our stream function values at X = column, Y=row in a (row, column) order
    //print the same information to a file called "output.txt".
    //These positions on each of these text representations will match that in the book
    row = 0;
    column = 0;
    for(row = number_of_rows - 1; row >= 0; row--)
    {
        for(column =0; column < number_of_columns; column++)
        {
            //printf("%.1f (%d, %d), ", new_grid[number_of_columns*row + column], row, column);
            fprintf(infile, "%lf ", new_grid[number_of_columns*row + column]);
        }
        fprintf(infile, "\n");
        //printf("\n");
    }

    fclose(infile);

    printf("See \"pressure_output.txt\" for pressure_distribution values\ncorresponding to the top and bottom walls\n");
    FILE * ofile = fopen("./pressure_output.txt", "w");

    infile = fopen(argv[3], "r");

    if ((infile == NULL) || (ofile == NULL)) 
    {
         printf("%s", "failed to open file\n");
    }

    entry = 0;
       //obtain a line from the boundary condition file
    while(fscanf(infile, "%s", buf) != EOF)
    {
            //we need to tokenize the line so that we can parse out the numerical information
            //we define a token as a floating point number, delimited by a comma
            token = strtok(buf, ",");
            column = 0;

            //while not out of tokens, do
            while(token != NULL)
            {
                //process the token as a floating point number
                sscanf(token, "%lf", &designation);

		if(designation == 1.0)
		{
			//Am I at the top, bottom edges, vertical wall, 'interior' horizontal boundary?
			if((entry/number_of_columns) < 1)
				u_vel = (new_grid[entry] - new_grid[entry + number_of_columns])/grid_spacing;
			else if (entry > (number_of_rows*number_of_columns - number_of_columns))
                                u_vel = (new_grid[entry-number_of_columns] - new_grid[entry])/grid_spacing;
			else if (entry%number_of_columns == 0)
			{
				u_vel = 0.0;
				v_vel = (new_grid[entry+1] - new_grid[entry])/(-grid_spacing);
			}
			else if (((entry+1)%number_of_columns) == 0)
			{
				u_vel = 0.0;
				v_vel = (new_grid[entry-1] - new_grid[entry])/grid_spacing;
			}
			else
			{
				pos_d_y = new_grid[entry-number_of_columns] - new_grid[entry];
				neg_d_y = new_grid[entry+number_of_columns] - new_grid[entry];
				pos_d_x = new_grid[entry+1] - new_grid[entry];
				neg_d_x = new_grid[entry-1] - new_grid[entry];

				if (pos_d_y == 0.0)
					u_vel = -neg_d_y/grid_spacing;
				else if (neg_d_y == 0.0)
					u_vel = pos_d_y/grid_spacing;
				else
					u_vel = (pos_d_y - neg_d_y)/(2.0*grid_spacing);

				if (pos_d_x == 0.0)
					v_vel = neg_d_x/grid_spacing;	
				else if (neg_d_x == 0.0)
					v_vel = -pos_d_x/grid_spacing;
				else
					v_vel = (neg_d_x - pos_d_x)/(2.0*grid_spacing);
			}

			//compute pressure using Bernouli equation
		 	double vel_mag_squared = pow(u_vel, 2.0) + pow(v_vel, 2.0);
			pressure = -.5*rho*(vel_mag_squared) + inlet_pressure + .5*rho*(inlet_velocity);
			int my_column = entry%number_of_columns;
			int my_row = (entry - my_column)/number_of_columns;
			int y_coordinate = my_row*-1 + number_of_rows;  
			fprintf(ofile, "%d %d %lf\n", my_column, y_coordinate, pressure);
		
		}
		
                entry++;
                token = strtok(NULL, ",");
            }
    }

    fclose(infile);
    fclose(ofile);

    free(boundary_grid);
    boundary_grid = NULL;

    free(buf);
    buf = NULL;

    free(old_grid);
    old_grid = NULL;

    free(new_grid);
    new_grid = NULL;

    //free(*interpolation);
    free(interpolation);
    interpolation = NULL;

    return 0;
}

void calculate_linear_interpolation(int length_of_interpolation,
                                    double** interpolation,
                                    const double min_value,
                                    const double max_value)
{
    double interval = (max_value - min_value)/((double)length_of_interpolation);
    int coordinate = 0;

    for(coordinate = 0; coordinate < length_of_interpolation; coordinate++)
    {
        *interpolation[coordinate] = min_value + interval*((double)coordinate);
    }
}

void compute_stream_values_in_rectangular_region(const int first_column,
                                                 const int last_column,
                                                 const int first_row,
                                                 const int last_row,
                                                 const int number_of_rows,
                                                 const int number_of_columns,
                                                 const double* problem_grid_old,
                                                 double* problem_grid_new,
                                                 const double* boundary_grid,
                                                 const double* grid_mask)
{
    int row = 0;
    int column = 0;

    for(row = first_row; row < last_row; row++)
    {
        for(column = first_column; column < last_column; column++)
        {
            problem_grid_new[number_of_columns*row + column] = boundary_grid[number_of_columns*row + column] +
                                                                 grid_mask[number_of_columns*row + column]*
                                                                 .25*(problem_grid_old[number_of_columns*row + column - 1] +
                                                                      problem_grid_old[number_of_columns*row + column + 1] +
                                                                      problem_grid_old[number_of_columns*(row -1) + column] +
                                                                      problem_grid_old[number_of_columns*(row + 1) + column]);
        }
    }


}

double compute_root_mean_square_error(const double* problem_grid_old,
                                    const double* problem_grid_new,
                                    const double number_of_interior_points,
                                    const int number_of_rows,
                                    const int number_of_columns)
{
    double running_total = 0.0;
    int row = 0;
    int column = 0;

    for(row = 0; row < number_of_rows; row++)
    {
        for(column = 0; column < number_of_columns; column++)
        {
            running_total += pow(problem_grid_old[number_of_columns*row + column] - 
                                 problem_grid_new[number_of_columns*row + column] , 2.0);
        }
    }

    running_total /= number_of_interior_points;

    return running_total;
}




