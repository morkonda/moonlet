/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem uses shearing sheet boundary
 * conditions. Particle properties resemble those found in 
 * Saturn's rings. 
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <unistd.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "input.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "display.h"

#define __USE_C99_MATH

extern double OMEGA;
int number_of_particles = 0;
double a;
double M_saturn;
double a_0;
double total_mass;
double total_mass_1;
double total_mass_2;
double total_mass_3;
double total_mass_4;
double total_mass_5;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v);
bool check(double x, double y);
bool check_to_add();
extern double opening_angle2;

//bool check_to_add()
//{
  //remove("position_moonlet.txt");
//FILE * fp_temp = fopen("position.txt", "r");
//}

void problem_init(int argc, char* argv[])
{
	remove("position_moonlet.txt");
	printf("test");
	//FILE *fp_moonlet;
	//fp_moonlet = fopen("position_moonlet.txt", "a");
	//fclose(fp_moonlet);
        for(int i=0; i<argc;i++)
	{
		printf("%d\t%s\n",i,argv[i]);
	}
	// Setup constants
#ifdef GRAVITY_TREE
	opening_angle2	= .5;
#endif // GRAVITY_TREE
	OMEGA 				= 0.00013143527;	// 1/s
	tmax                            = 10.* 2.*M_PI/OMEGA;  
	G 				= 6.67428e-11;		// N / (1e-5 kg)^2 m^2
	M_saturn 			= 568.36e24;		// kg
	a 				= pow(((G*M_saturn)/(OMEGA*OMEGA)),(1/3));		// m
	a_0 				= 134912000. ;		// m 
	softening 			= 0.1;			// m
	dt 				= 1e-3*2.*M_PI/OMEGA;	// s
#ifdef OPENGL
	display_rotate_z		= 20;			// Rotate the box by 20 around the z axis, then 
	display_rotate_x		= 60;			// rotate the box by 60 degrees around the x axis	
#ifdef LIBPNG
	system("mkdir png");
#endif // LIBPNG
#endif // OPENGL
	root_nx = 2; root_ny = 2; root_nz = 1;
	nghostx = 2; nghosty = 2; nghostz = 0; 			// Use two ghost rings

//	double surfacedensity		= 400; 			// kg/m^2

	double particle_density		= 400;			// kg/m^3
	double particle_radius_min 	= 1;			// m
	double particle_radius_max 	= 4;			// m
	double particle_radius_slope 	= -3;	
	//double increment = input_get_double(argc,argv,"increment",1.1);
///	printf("%f\n",input_get_double(argc,argv,"a",123));
	double increment 		= 1.9;
	boxsize 			= 100.*increment;	// m
///	if (argc>1){						//Try to read boxsize from command line
	//	boxsize = atof(argv[1]);
	//}
///

init_box();
	
	// Initial conditions

	//printf("Toomre wavelength: %f\n",2.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*G);

	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear

	//double total_mass = surfacedensity*boxsize_x*boxsize_y;
	bool result 				= false;
	double mass 				= 0.;
	double temp_x 				= 0.;
	double temp_sigma 			= 0.;
	double number_density 			= 0.;
	double surface_area 			= 0.;
	double x_left_1		 		= -32.5*increment;
	double x_left_2				= 0.;
        //double x_left_2		 	= -28.0*increment;
	double x_right_1 	 		= 27.5*increment;
	//double x_right_2	 		= 46.0*increment;
	double x_right_2 			= 0.;
	double slope_left 			= -1.e100;
	double slope_right 			= 1.e100;
	double surface_density_upper 		= 338.;
	double surface_density_lower 		= 3.;

	x_left_2 = (1/slope_left)*(surface_density_lower-surface_density_upper+(slope_left*x_left_1));
	x_right_2 = (1/slope_right)*(surface_density_upper-surface_density_lower+(slope_right*x_right_1));

	FILE *fp;
	fp = fopen("position_x.txt", "w");
	surface_area = (boxsize_x/2.)*(boxsize_y/2);

	fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n", x_left_1, x_left_2, x_right_1, x_right_2, surface_density_lower, surface_density_upper, slope_left, slope_right);

////
bool check(double x, double y)
{
        if(x < x_left_1)
         {
                if(y <= surface_density_upper)
                 {
                        return true;
                 }
         }

        if(x >= x_left_1 && x < x_left_2)
         {
                if(y <= ((slope_left*x)+(surface_density_upper)-(slope_left*x_left_1)))
                 {
                        return true;
                 }
         }

        if(x >= x_left_2 && x < x_right_1)
         {
                if(y <= surface_density_lower)
                 {
                        return true;
                 }
         }

        if(x >= x_right_1 && x < x_right_2)
         {
                if(y <= ((slope_right*x)+(surface_density_lower)-(slope_right*x_right_1)))
                 {
                        return true;
                 }
         }

        if(x >= x_right_2)
         {
                if(y <= surface_density_upper)
                 {
                        return true;
                 }
         }

        return false;
}


	double total_mass_1 = (surface_density_upper)*((boxsize_x/2)+x_left_1)*boxsize_y;
	double total_mass_2 = ((x_left_2-x_left_1)*surface_density_lower)+(0.5*(x_left_2-x_left_1)*(surface_density_upper-surface_density_lower));
	double total_mass_3 = (surface_density_lower)*(x_right_1-x_left_2)*boxsize_y;
	double total_mass_4 = ((x_right_2-x_right_1)*surface_density_lower)+(0.5*(x_right_2-x_right_1)*(surface_density_upper-surface_density_lower));
	double total_mass_5 = (surface_density_upper)*((boxsize_x/2.)-x_right_2)*boxsize_y;

	total_mass = total_mass_1+total_mass_2+total_mass_3+total_mass_4+total_mass_5;


     while(mass<total_mass)
{
	temp_x = tools_uniform(-boxsize_x/2.,boxsize_x/2.);
	temp_sigma = tools_uniform(0, (surface_density_upper+50.));
	result = check(temp_x, temp_sigma);
	if(result)
	 {
        	struct particle pt;
        	pt.x           		= temp_x;
        	pt.y           		= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
        	pt.z           		= tools_normal(1.);                     
        	pt.vx          		= 0;
        	pt.vy          		= -1.5*pt.x*OMEGA;
        	pt.vz          		= 0;
        	pt.ax          		= 0;
        	pt.ay          		= 0;
        	pt.az          		= 0;
        	double radius  		= tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
        	//#ifndef COLLISIONS_NONE
		pt.r           		= radius;                          
		//#endif
        	double particle_mass    = particle_density*4./3.*M_PI*radius*radius*radius;
        	pt.m           		= particle_mass;    
        	particles_add(pt);
        	mass += particle_mass;
		number_of_particles++;
		number_density 		= total_mass/(particle_mass*surface_area);
		fprintf(fp, "%f \t %f \t %f \n", pt.x, temp_sigma, number_density);
		
	 }
}
	struct particle pt;
	pt.x				= 0;
	//pt.x 				= (1./3.)*(x_left_1);
	pt.y				= tools_normal(1.);
	pt.z				= tools_normal(1.);
	pt.vx				= 0;
 	pt.vy				= -1.5*pt.x*OMEGA;
	pt.vz				= 0;
	pt.ax				= 0;
	pt.ay				= 0;
	pt.az				= 0;
	double radius 			= 100.;
	pt.r 				= radius;
	double particle_mass 		= particle_density*4./3.*M_PI*radius*radius*radius;
	pt.m				= particle_mass;
	particles_add(pt);
	number_of_particles++;

	fclose(fp);

}                
                 
double coefficient_of_restitution_bridges(double v)
{
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void problem_inloop()
{
}

void problem_output()
{
	if (output_check(1e-3*2.*M_PI/OMEGA))
	{
		output_timing();
		//output_append_velocity_dispersion("veldisp.txt");
	}
	if (output_check(2.*M_PI/OMEGA))
	{
		char fp_position_data[1024];
		sprintf(fp_position_data,"position_%2.1f.txt",t/(2.*M_PI/OMEGA));
		output_ascii(fp_position_data);

	}

	if (output_check(2.*M_PI/OMEGA))
	{
		FILE* fp_moonlet_data = fopen("position_moonlet.txt","a+");
		for(int i=0;i<number_of_particles;i++)
		 {
			if (particles[i].r == 100.)
		         { 
		     		struct particle ml = particles[i];
				fprintf(fp_moonlet_data,"%e\t%e\t%e\t%e\t%e\t%e\n", ml.x, ml.y, ml.z, ml.vx, ml.vy, ml.vz);
				fclose(fp_moonlet_data);
			 }
		 }
	}
}

void problem_finish()
{
}
