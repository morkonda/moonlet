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
double r_h;
double toomre_wavelength;
double total_mass;
double total_mass_1;
double total_mass_2;
double total_mass_3;
double total_mass_4;
double total_mass_5;
double moonlet_radius;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v);
bool check(double x, double y);
double toomre_wavelength_function(double x);
double sigma_function(double x);
double deltaJ_function(double x);
extern double opening_angle2;


void problem_init(int argc, char* argv[])
{
	remove("position_moonlet.txt");
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
	a 				= pow(((G*M_saturn)/(OMEGA*OMEGA)),(1./3.));		//m
	a_0 				= 134912000. ;		// m 
	moonlet_radius                  = 100.;			// m
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

	double surfacedensity		= 400; 			// kg/m^2

	double particle_density		= 400;			// kg/m^3
	//--> double particle_radius_min 	= 1;
	double particle_radius_min	= 1;			// m
	//--> double particle_radius_max 	= 4;
	double particle_radius_max	= 4;			// m
	double particle_radius_slope 	= -3;	
	//double increment = input_get_double(argc,argv,"increment",1.1);
///	printf("%f\n",input_get_double(argc,argv,"a",123));
	double increment 		= 5.;
	double moonlet_mass             = particle_density*4./3.*M_PI*moonlet_radius*moonlet_radius*moonlet_radius;

	boxsize 			= 100.*increment;	// m
///	if (argc>1){						//Try to read boxsize from command line
	//	boxsize = atof(argv[1]);
	//}
///

init_box();
	
	// Initial conditions

	//printf("Toomre wavelength: %f m \n",2.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*G);

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
	//--> double x_right_1 	 		= 27.5*increment;
	double x_right_1			= 32.5*increment;
	double x_right_2 			= 0.;
	double slope_left 			= -1.e100;
	double slope_right 			= 1.e100;
	//--> double surface_density_upper 	= 338.;
	double surface_density_upper		= 338./2.;
	double surface_density_lower 		= 3.;
	double alpha				= 2.46;
	double K_0				= 0.69676999;
	double K_1				= 0.45731925;
	double n_0				= OMEGA;		// 1/s
	double J_m				= a*a*n_0;
	r_h                                     = a*(pow((moonlet_mass/M_saturn),(1./3.)));
	x_left_2 				= (1/slope_left)*(surface_density_lower-surface_density_upper+(slope_left*x_left_1));
	x_right_2 				= (1/slope_right)*(surface_density_upper-surface_density_lower+(slope_right*x_right_1));

	FILE *fp;
	fp = fopen("position_x.txt", "w");
	surface_area = (boxsize_x/2.)*(boxsize_y/2);

	fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n", x_left_1, x_left_2, x_right_1, x_right_2, surface_density_lower, surface_density_upper, slope_left, slope_right);

	printf("Hill radius of the moonlet = %e m \n",r_h);

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


double result_left_1(double x)
{
	double lower_limit = -boxsize_x/2.;
	double upper_limit = -2.5*r_h;
	double constant_neg_mdr = -3.*n_0*r_h*r_h*((64.*pow(G*moonlet_mass,2.)*a)/(243.*pow(OMEGA,3.)))*(pow(((2.*K_0)+(K_1)),2.));

	double integral_left_1_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_upper_value = ((-1./(3.*pow(-integral_upper_limit,3)))-(alpha/(2.*a*pow(-integral_upper_limit,2))));
		double integral_lower_value = ((-1./(3.*pow(-integral_lower_limit,3)))-(alpha/(2.*a*pow(-integral_lower_limit,2))));
		return constant_neg_mdr*integral_sigma*(integral_upper_value-integral_lower_value);
	 }
	
		
	double integral_left_1_for_slope_sigma(double integral_lower_limit, doublt integral_upper_limit)
	 {
		double integral_upper_value = (-((slope_left*r_h)/(2.*pow(-integral_upper_limit,2)))-((surface_density_upper-(slope_left*x_left_1))/(3.*pow(-integral_upper_limit,3)))-((alpha*(slope_left*r_h))/(a*-integral_upper_limit))-((alpha*(surface_density_upper-(slope_left*x_left_1)))/(2.*a*pow(-integral_upper_limit,2))));
		double integral_lower_value = (-((slope_left*r_h)/(2.*pow(-integral_upper_limit,2)))-((surface_density_upper-(slope_left*x_left_1))/(3.*pow(-integral_lower_limit,3)))-((alpha*(slope_left*r_h))/(a*-integral_lower_limit))-((alpha*(surface_density_upper-(slope_left*x_left_1)))/(2.*a*pow(-integral_lower_limit,2))));
		return constant_neg_mdr*(integral_upper_value-integral_lower_value);
	 }


//x_left_1 is to the right of (-2.5*r_h)
	if(upper_limit <= x_left_1)
	 {
		return integral_left_1_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 }
	
//x_left_1 is to the left of (-2.5*r_h) and x_left_2 is to the right of (-2.5*r_h)
	else if(upper_limit > x_left_1 && upper_limit < x_left_2) 
	 {
		double temp_result_1 = integral_left_1_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);		
		double temp_result_2 = integral_left_1_for_slope_sigma(x_left_1, upper_limit);	
		return temp_result_1 + temp_result_2;
	 }
	
//x_left_1 and x_left_2 are to the left of (-2.5*r_h)
	else if(upper_limit > x_left_2)
	 {
		double temp_result_1 = integral_left_1_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_1_for_slope_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_1_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }
}

double result_left_2(double x)
{
	double lower_limit = -2.5*r_h;
	double upper_limit = -1.8*r_h;
	double constant_neg_cdr = -3.*n_0*r_h*r_h*(r_h/(2.*a))*J_m;
	
	double integral_left_2_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_upper_value = (1./3.)*pow(-integral_upper_limit,3);
		double integral_lower_value = (1./3.)*pow(-integral_lower_limit,3);
		return constant_neg_cdr*integral_sigma*(integral_upper_value - integral_lower_value);
	 }

	double integral_left_2_for_slope_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double integral_upper_value = (((slope_left*r_h*pow(-integral_upper_limit,4))/4.)+(((surface_density_upper-(slope_left*x_left_1))*pow(-integral_upper_limit,3))/3.));
		double integral_lower_value = (((slope_left*r_h*pow(-integral_lower_limit,4))/4.)+(((surface_density_upper-(slope_left*x_left_1))*pow(-integral_lower_limit,3))/3.));
		return constant_neg_cdr*(integral_upper_value - integral_lower_value);
	 }

//x_left_1 is to the right of (-1.8*r_h)
	if(upper_limit <= x_left_1)
	 {
		return integral_left_2_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 }

//x_left_1 is in between (-2.5*r_h) and (-1.8*r_h) and x_left_2 is to the right of (1.8*r_h)
	else if(lower_limit < x_left_1 && upper_limit < x_left_2)
	 {
		double temp_result_1 = integral_left_2_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_2_for_slope_sigma(x_left_1, upper_limit);
		return temp_result_1 + temp_result_2;
	 }

//x_left_1 and x_left_2 are in between (-2.5*r_h) and (-1.8*r_h)
	else if(lower_limit < x_left_1 && upper_limit < x_left_2)
	 {
		double temp_result_1 = integral_left_2_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_2_for_slope_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_2_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2 + temp_result_3;
	}
	
//x_left_1 is to the left of (-2.5*r_h) and x_left_2 is in between (-2.5*r_h) and (-1.8*r_h)
	else if(lower_limit > x_left_1 && upper_limit > x_left_2)
	 {
		double temp_result_1 = integral_left_2_for_slope_sigma(lower_limit, x_left_2);
		double temp_result_2 = integral_left_2_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);		
		return temp_result_1 + temp_result_2;
	}

//both x_left_1 and x_left_2 are to the left of (-2.5*r_h)
	else if(lower_limit > x_left_2)
	 {
		return integral_left_2_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
 	 }	
}


double result_left_3(double x)
{
	double lower_limit = -1.8*r_h;
	double upper_limit = 0;
	double constant_neg_hr = -3.*n_0*r_h*r_h*(r_h/a.)*J_m;

	double integral_left_3_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double integral_upper_value = (1./3.)*pow(integral_upper_limit,3);
		double integral_lower_value = (1./3.)*pow(integral_lower_limit,3);
		return constant_neg_hr*integral_sigma*(integral_upper_value - integral_lower_value);
	 }
	
	double integral_left_3_for_slope_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double integral_upper_value = (((slope_left*r_h*pow(-integral_upper_limit,4))/4.)+(((surface_density_upper-(slope_left*x_left_1))*pow(-integral_upper_limit,3))/3.));
		double integral_lower_value = (((slope_left*r_h*pow(-integral_lower_limit,4))/4.)+(((surface_density_upper-(slope_left*x_left_1))*pow(-integral_lower_limit,3))/3.));
		return constant_neg_hr*(integral_upper_value - integral_lower_value);
	 }

//both x_left_1 and x_left_2 are to the right of (-1.8*r_h)	
	if(x_left_1 > lower_limit)
	 {
		double temp_result_1 = integral_left_3_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);		
		double temp_result_2 = integral_left_3_for_slope_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_3_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//x_left_1 is to the left of (-1.8*r_h) and x_left_2 is to the right of (-1.8*r_h)
	else if(x_left_2 > lower_limit)
	 {
		double temp_result_1 = integral_left_3_for_slope_sigma(lower_limit, x_left_2);		
		double temp_result_2 = integral_left_3_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2;
	 }

//both x_left_1 and x_left_2 are to the left of (-1.8*r_h)	
	else
	 {
		return integral_left_3_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
	 }
}


double torque_function(double x)
{
		

}


double sigma_function(double x)
{
	if(x < x_left_1)
	 {
		return surface_density_upper;
	 }
	
	if(x >= x_left_1 && x < x_left_2)
	 {
		return ((slope_left*x)+(surface_density_upper)-(slope_left*x_left_1));
	 }

	if(x >= x_left_2 && x < x_right_1)
	 {
		return surface_density_lower;
	 }

	if(x >= x_right_1 && x < x_right_2)
	 {
		return ((slope_right*x)+(surface_density_lower)-(slope_right*x_right_1));
	 }
	
	if(x >= x_right_2)
	 {
		return surface_density_upper;
	 }
}


double deltaJ_function(double x)
{
	double temp_b = x/r_h;

	//negative Moderate Deflection Region	
	if(temp_b < -2.5)
	 {
		double temp_constant = (64.*(pow(G*moonlet_mass,2.))*a)/(243.*pow(OMEGA,3.)*pow(fabs(temp_b),5.));
		double temp_value    = (pow(((2.*K_0)+K_1),2.))*(1.+(alpha*(fabs(temp_b)/a)));
		return -temp_constant*temp_value;
	 }
	
	//negative Chaotic Deflection Region
	if(temp_b >= -2.5 && temp_b < -1.8)
	 {
		return -(fabs(temp_b)*r_h*J_m)/(2.*a);
	 }

	//negative Horseshoe Region
	if(temp_b >= -1.8 && temp_b < 0)
	 {
		return -(fabs(temp_b)*r_h*J_m)/a;
	 }

	//positive Horseshoe Region
	if(temp_b >= 0 && temp_b < 1.8)
	 {
		return (temp_b*r_h*J_m)/a;
	 }
	
	//positive Chaotic Deflection Region
	if(temp_b >= 1.8 && temp_b < 2.5)
	 {
		return (temp_b*r_h*J_m)/(2.*a);
	 }

	//positive Moderate Deflection Region
	if(temp_b >= 2.5)
	 {
		double temp_constant = (64.*(pow((G*moonlet_mass),2.))*a)/(243.*pow(OMEGA,3.)*pow(temp_b,5.));
		double temp_value    = (pow(((2.*K_0)+K_1),2.))*(1.+(alpha*(temp_b/a)));
		return temp_constant*temp_value;
	 }
}


double toomre_wavelength_function(double x)
{
        return (2.*M_PI*M_PI*(sigma_function(x))/OMEGA/OMEGA*G);
}

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
		toomre_wavelength 	= toomre_wavelength_function(pt.x);
		fprintf(fp, "%f \t %f \t %f \t %f \n", pt.x, temp_sigma, number_density, toomre_wavelength);
		
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
	pt.r 				= moonlet_radius;
	pt.m				= moonlet_mass;
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
	if (output_check((2.*M_PI/OMEGA)/10.))
	{
		char fp_position_data[1024];
		sprintf(fp_position_data,"Data/position_%2.1f.txt",t/(2.*M_PI/OMEGA));
		output_ascii(fp_position_data);

	}

	if (output_check(2.*M_PI/OMEGA))
	//if(0)
	{
		FILE* fp_moonlet_data = fopen("position_moonlet.txt","a+");
		for(int i=0;i<number_of_particles;i++)
		 {
			if (particles[i].r == moonlet_radius)
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
