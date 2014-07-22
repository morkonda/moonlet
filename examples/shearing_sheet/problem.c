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
FILE *fp;
double a;
double M_saturn;
double a_0;
double J_m;
double r_h;
double toomre_wavelength;
double total_mass;
double total_mass_1;
double total_mass_2;
double total_mass_3;
double total_mass_4;
double total_mass_5;
double moonlet_x;
double moonlet_radius;
double moonlet_mass;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v);
bool check(double x, double y);
double result_left_1();
double result_left_2();
double result_left_3();
double result_right_1();
double result_right_2();
double result_right_3();
double torque_function();
double migration_rate_function();
double sigma_function(double x);
double toomre_wavelength_function(double x);
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
	OMEGA 					= 0.00013143527;				// 1/s
	tmax                           		= 10.*2.*M_PI/OMEGA;  
	G 					= 6.67428e-11;					// N / (1e-5 kg)^2 m^2
	//M_saturn 				= 568.36e24;					// kg
	M_saturn				= 5.e26;					// kg
	//a_0 					= pow(((G*M_saturn)/(OMEGA*OMEGA)),(1./3.));	// m
	//a_0 					= 134912000. ;					// m
	softening 				= 0.1;						// m
	dt 					= 1e-3*2.*M_PI/OMEGA;				// s
#ifdef OPENGL
	display_rotate_z			= 20;						// Rotate the box by 20 around the z axis, then 
	display_rotate_x			= 60;						// rotate the box by 60 degrees around the x axis	
#ifdef LIBPNG
	system("mkdir png");
#endif // LIBPNG
#endif // OPENGL
	root_nx = 2; root_ny = 2; root_nz = 1;
	nghostx = 2; nghosty = 2; nghostz = 0; 							// Use two ghost rings

	//double surfacedensity			= 400; 						// kg/m^2

	double particle_density			= 400;						// kg/m^3
	double moonlet_density			= 400;						// kg/m^3
	//--> double particle_radius_min 	= 1;
	double particle_radius_min		= 1;						// m
	//--> double particle_radius_max 	= 4;
	double particle_radius_max		= 4;						// m
	double particle_radius_slope 		= -3;	
	//double increment = input_get_double(argc,argv,"increment",1.1);
///	printf("%f\n",input_get_double(argc,argv,"a",123));
	double increment 			= 100.;

	boxsize 				= 100.*increment;				// m
///	if (argc>1){										//Try to read boxsize from command line
	//	boxsize = atof(argv[1]);
	//}
///

init_box();
	
	// Initial conditions

	//printf("Toomre wavelength: %f m \n",2.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*G);

	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  				// small fraction of the shear

	//double total_mass = surfacedensity*boxsize_x*boxsize_y;
	bool result 				= false;
	double mass 				= 0.;
	double temp_x 				= 0.;
	double temp_sigma 			= 0.;
	double number_density 			= 0.;
	double surface_area 			= 0.;
	double x_left_1			 	= -32.4*increment;
	//double x_right_1 	 		= 27.5*increment;
	//double x_left_1			= -35.5*increment;
	//double x_left_1			= -1500.;
	//double x_right_1			= 27.5*increment;
	//double slope_left			= -0.076;
	//double slope_right			= 0.0205;
	double slope_left 			= -0.048;		//real slope		// kg m^-2 m^-1
	double slope_right	 		= 0.0144;		//real slope		// kg m^-2 m^-1
	//slope_right				= 0.048;
	//slope_left				= -10.48;
	//slope_right				= 10.48;
	//double slope_left			= -1.e7;
	//double slope_right			= 1.e7;
	//double slope_left			= -0.076;
	//double slope_right			= 0.0205;
	//--> double surface_density_upper 	= 338.;
	double surface_density_upper		= 338.;
	//double surface_density_lower		= 300.;
	double surface_density_lower	 	= 338.;
	double alpha				= 2.46;
	double K_0				= 0.69676999;
	double K_1				= 0.45731925;
	double n_0				= OMEGA;					// 1/s
	//n_0					= 10.7648/(24*60*60);
	double x_left_2	 			= (1./slope_left)*(surface_density_lower-surface_density_upper+(slope_left*x_left_1));
	double x_right_1			= -x_left_2;
	double x_right_2 			= (1./slope_right)*(surface_density_upper-surface_density_lower+(slope_right*x_right_1));
	//double x_left_2 			= (((surface_density_lower)/(slope_left)-((surface_density_upper)/(slope_left))+x_left_1);
	////double x_left_2 			= ((surface_density_lower/slope_left)-(surface_density_upper/slope_left)+(x_left_1));
	//double x_right_2			= (((surface_density_upper-surface_density_lower)/(slope_right))+x_right_1);
	////double x_right_2 			= ((surface_density_upper/slope_right)-(surface_density_lower/slope_right)+(x_right_1));	

	//x_right_2 				= -x_left_1;
	//x_right_1 				= -x_left_2;	

	printf("1/s_l = %f\n", (1./slope_left)*(surface_density_lower-surface_density_upper));	

	moonlet_radius				= 1000.;					// m
	moonlet_mass 			        = moonlet_density*(4./3.)*M_PI*moonlet_radius*moonlet_radius*moonlet_radius;
	//moonlet_mass				= 2.e12;					// kg
	a_0					= pow(((G*(M_saturn+moonlet_mass))/(OMEGA*OMEGA)),(1./3.));
	//moonlet_x 				= 0.;						// m
	//moonlet_x				= (1./10.)*x_right_1;				// m 
	moonlet_x				= -300.;					// m
	a				        = a_0 + moonlet_x;                              // m		//Actual semi-major axis
	J_m 					= a_0*a_0*n_0;
	r_h                                     = a_0*(pow((moonlet_mass/(3.*M_saturn)),(1./3.)));
	fp 					= fopen("position_x.txt", "w");
	surface_area 				= (boxsize_x/2.)*(boxsize_y/2);

	fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n", x_left_1, x_left_2, x_right_1, x_right_2, surface_density_lower, surface_density_upper, slope_left, slope_right);

	printf("x_left_1, x_left_2 = %16.16f\t%16.16f\n",x_left_1,x_left_2);
	printf("x_right_1, x_right_2 = %16.16f\t%16.16f\n",x_right_1,x_right_2);

	//printf("moonlet_mass = %f\n",moonlet_mass);
	printf("a_0 = %f\n", a_0);
	//printf("a_moonlet = %f\n", a);
	//printf("J_m = %f\n", J_m);
	printf("Hill radius of the moonlet = %e m \n",r_h);
	printf("Moonlet_x = %f m \n", moonlet_x);

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

	total_mass_1 = (surface_density_upper)*((boxsize_x/2.)+x_left_1)*boxsize_y;
	total_mass_2 = ((x_left_2-x_left_1)*surface_density_lower)+(0.5*(x_left_2-x_left_1)*(surface_density_upper-surface_density_lower));
	total_mass_3 = (surface_density_lower)*(x_right_1-x_left_2)*boxsize_y;
	total_mass_4 = ((x_right_2-x_right_1)*surface_density_lower)+(0.5*(x_right_2-x_right_1)*(surface_density_upper-surface_density_lower));
	total_mass_5 = (surface_density_upper)*((boxsize_x/2.)-x_right_2)*boxsize_y;
	total_mass = total_mass_1+total_mass_2+total_mass_3+total_mass_4+total_mass_5;


double sigma_left(double x)
{
	return (slope_left*x)-(slope_left*x_left_1)+surface_density_upper;
}

double sigma_right(double x)
{
	return (slope_right*x)-(slope_right*x_right_1)+surface_density_lower;
}

	//printf("y_left_upper, y_right_upper = %16.16e\t%16.16e\n", sigma_left(x_left_1),sigma_right(x_right_2));
	//printf("y_left_lower, y_right_lower = %16.16e\t%16.16e\n", sigma_left(x_left_2),sigma_right(x_right_1));

//function that returns the value of b for a given x value
double check_limit_sign(double check_limit)
{
	if(check_limit <= 0.)
	 {
		return -(fabs(check_limit));
	 }
	else if(check_limit > 0.)
	 {
		return (fabs(check_limit));
	 }
}


double result_left_1()
{
	double lower_limit = -boxsize_x/2.;
	double upper_limit = -2.5*r_h;
        double constant_neg_mdr = 3.*n_0*r_h*r_h*((64.*pow(G*moonlet_mass,2.)*a)/(243.*pow(OMEGA,3.)))*(pow(((2.*K_0)+(K_1)),2.));

	double integral_left_1_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {

		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);

		double integral_lower_value = ((1./(3.*pow(lower_limit_local,3.)))-(alpha/(2.*a*pow(lower_limit_local,2.))));	
		double integral_upper_value = ((1./(3.*pow(upper_limit_local,3.)))-(alpha/(2.*a*pow(upper_limit_local,2.))));

		if(integral_lower_limit == lower_limit)
		 {
			integral_lower_value = 0.;
		 }

		double temp_result = (1./(r_h*r_h))*constant_neg_mdr*integral_sigma*(integral_upper_value-integral_lower_value);

		//printf("int_upper_value = %f\n", integral_upper_value);
		//printf("int_lower_value = %f\n", integral_lower_value);

		//printf("CHECK LEFT CONST = %16.16e\n", temp_result);
		return temp_result;
	 }
	
		
	double integral_left_1_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		//printf("-------------------------->Check\n");
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);	

		double integral_lower_value = (+((slope_left)/(2.*pow(lower_limit_local,2.)))+((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))/(3.*pow(lower_limit_local,3.)))-((alpha*slope_left)/(a*lower_limit_local))-((alpha*(-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x)))/(2.*a*pow(lower_limit_local,2.))));
		double integral_upper_value = (+((slope_left)/(2.*pow(upper_limit_local,2.)))+((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))/(3.*pow(upper_limit_local,3.)))-((alpha*slope_left)/(a*upper_limit_local))-((alpha*(-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x)))/(2.*a*pow(upper_limit_local,2.))));

		if(integral_lower_limit == lower_limit)
		 {
			integral_lower_value = 0.;
		 }		

		double temp_result = (1./(r_h*r_h))*constant_neg_mdr*(integral_upper_value-integral_lower_value);

		return temp_result;
	 }


//x_left_1 is to the right of (-2.5*r_h)
	if(upper_limit <= x_left_1)
	 {
		//printf("..............1\n");
		return integral_left_1_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 }
	
//x_left_1 is to the left of (-2.5*r_h) and x_left_2 is to the right of (-2.5*r_h)
	else if(upper_limit > x_left_1 && upper_limit < x_left_2) 
	 {
		//printf("..............2\n");
		double temp_result_1 = integral_left_1_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);		
		double temp_result_2 = integral_left_1_for_function_sigma(x_left_1, upper_limit);	
		return temp_result_1 + temp_result_2;
	 }
	
//x_left_1 and x_left_2 are to the left of (-2.5*r_h)
	else if(upper_limit > x_left_2)
	 {
		//printf("..............3\n");
		double temp_result_1 = integral_left_1_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_1_for_function_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_1_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);

		//printf("Temp_result_1_left = %f\n", temp_result_1);
		//printf("Temp_result_2_left = %f\n", temp_result_2);
		//printf("Temp_result_3_left = %f\n", temp_result_3);

		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_left_1");
		return 0.;
	 }
}

double result_left_2()
{
	double lower_limit = -2.5*r_h;
	double upper_limit = -1.8*r_h;
	double constant_neg_cdr = 3.*n_0*r_h*r_h*(r_h/(2.*a_0))*J_m;

	double integral_left_2_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);		
	
		//printf("Left Constant: lower_limit, upper_limit = %16.16f\t%16.16f\n",lower_limit_local, upper_limit_local);
		//printf("Constant sigma: left_lower_limit, left_upper_limit = %f\t%f\n", lower_limit_local, upper_limit_local);

		double integral_upper_value = -(1./3.)*pow(upper_limit_local,3.);
		double integral_lower_value = -(1./3.)*pow(lower_limit_local,3.);
		double temp_result = (1./(r_h*r_h*r_h))*constant_neg_cdr*integral_sigma*(integral_upper_value - integral_lower_value);
	
		//printf("result_left_2_const = %f\n", temp_result);		

		return temp_result;
	 }

	double integral_left_2_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);		

		//printf("Left Function: lower_limit, upper_limit = %16.16f\t%16.16f\n", lower_limit_local, upper_limit_local);
		//printf("Function sigma: left_lower_limit, left_upper_limit = %f\t%f\n", lower_limit_local, upper_limit_local);
		
		double integral_upper_value = -(((slope_left*pow(upper_limit_local,4.))/4.)+(((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))*pow(upper_limit_local,3.))/3.));
		double integral_lower_value = -(((slope_left*pow(lower_limit_local,4.))/4.)+(((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))*pow(lower_limit_local,3.))/3.));
		double temp_result = (1./(r_h*r_h*r_h))*constant_neg_cdr*(integral_upper_value - integral_lower_value);
		return temp_result;
	 }

//x_left_1 is to the right of (-1.8*r_h)
	if(upper_limit <= x_left_1)
	 {
		//printf("CHECK_LEFT_1\n");
		return integral_left_2_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 }

//x_left_1 is in between (-2.5*r_h) & (-1.8*r_h) and x_left_2 is to the right of (1.8*r_h)
	else if(lower_limit < x_left_1 && upper_limit > x_left_1 && upper_limit < x_left_2)
	 {
		//printf("CHECK_LEFT_2\n");
		double temp_result_1 = integral_left_2_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_2_for_function_sigma(x_left_1, upper_limit);
		return temp_result_1 + temp_result_2;
	 }

//x_left_1 and x_left_2 are in between (-2.5*r_h) and (-1.8*r_h)
	else if(lower_limit < x_left_1 && upper_limit > x_left_2)
	 {
		//printf("CHECK_LEFT_3\n");
		double temp_result_1 = integral_left_2_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);
		double temp_result_2 = integral_left_2_for_function_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_2_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 } 
	
//x_left_1 is to the left of (-2.5*r_h) and x_left_2 is in between (-2.5*r_h) & (-1.8*r_h)
	else if(lower_limit > x_left_1 && lower_limit < x_left_2 && upper_limit > x_left_2)
	 {
		//printf("CHECK_LEFT_4\n");
		double temp_result_1 = integral_left_2_for_function_sigma(lower_limit, x_left_2);
		double temp_result_2 = integral_left_2_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);

		//printf("Integral_function_left = %f\n", temp_result_1);
		//printf("Integral_constant_left = %f\n", temp_result_2);		
		
		return temp_result_1 + temp_result_2;
	 }

//both x_left_1 and x_left_2 are to the left of (-2.5*r_h)
	else if(lower_limit > x_left_2)
	 {
		//printf("CHECK_LEFT_5\n");
		return integral_left_2_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
 	 }

//x_left_1 is to the left of (-2.5*r_h) and x_left_2 is to the right of (-1.8*r_h)
	else if(lower_limit > x_left_1 && upper_limit < x_left_2)
	 {
		//printf("CHECK_LEFT_6\n");
		return integral_left_2_for_function_sigma(lower_limit, upper_limit);
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_left_2");
		return 0.;
	 }	
}


double result_left_3()
{
	double lower_limit = -1.8*r_h;
	double upper_limit = 0.;
	double constant_neg_hr = 3.*n_0*r_h*r_h*(r_h/a_0)*J_m;

	double integral_left_3_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);

		double integral_upper_value = -(1./3.)*pow(upper_limit_local,3.);
		double integral_lower_value = -(1./3.)*pow(lower_limit_local,3.);

		//printf("Left_lower_limit = %f\n", lower_limit_local);
		//printf("Left_upper_limit = %f\n", upper_limit_local);

		double temp_result = (1./(r_h*r_h*r_h))*constant_neg_hr*integral_sigma*(integral_upper_value - integral_lower_value);	
		return temp_result; 
	 }
	
	double integral_left_3_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);

		double integral_upper_value = -(((slope_left*pow(upper_limit_local,4.))/4.)+(((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))*pow(upper_limit_local,3.))/3.));
		double integral_lower_value = -(((slope_left*pow(lower_limit_local,4.))/4.)+(((-(slope_left*x_left_1)+surface_density_upper+(slope_left*moonlet_x))*pow(lower_limit_local,3.))/3.));
		double temp_result = (1./(r_h*r_h*r_h))*constant_neg_hr*(integral_upper_value - integral_lower_value);
		
		//printf("result_left_function = %f\n", temp_result);
	
		return temp_result; 
	 } 

//both x_left_1 and x_left_2 are to the right of (-1.8*r_h)	
	if(x_left_1 > lower_limit)
	 {
		printf("---------------Check_3\n");
		double temp_result_1 = integral_left_3_for_constant_sigma(lower_limit, x_left_1, surface_density_upper);		
		double temp_result_2 = integral_left_3_for_function_sigma(x_left_1, x_left_2);
		double temp_result_3 = integral_left_3_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//x_left_1 is to the left of (-1.8*r_h) and x_left_2 is to the right of (-1.8*r_h)
	else if(x_left_1 < lower_limit && x_left_2 > lower_limit)
	 {
		printf("---------------Check_2\n");
		double temp_result_1 = integral_left_3_for_function_sigma(lower_limit, x_left_2);		
		double temp_result_2 = integral_left_3_for_constant_sigma(x_left_2, upper_limit, surface_density_lower);
		return temp_result_1 + temp_result_2;
	 }

//both x_left_1 and x_left_2 are to the left of (-1.8*r_h)	
	else if(x_left_2 < lower_limit)
	 {
		//printf("---------------Check_1\n");

		double temp_result = integral_left_3_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
		return temp_result;
		
		//double temp_1 = -(1./(r_h*r_h*r_h))*(constant_neg_hr*surface_density_lower)*(1./3.)*pow(-300.-(-2625.993411),3.);
		//double temp_2 = -(1./(r_h*r_h*r_h))*(constant_neg_hr*surface_density_lower)*(1./3.)*pow(0.-(-300.),3.);
		//return temp_1-temp_2;

		//printf("Left_result = %f\n", temp_result);			
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_left_3");
		return 0.;
	 }
}


//Clockwise: -ve torque
double result_right_1()
{
	double lower_limit = 0.;
	double upper_limit = 1.8*r_h;
	double constant_pos_hr = 3.*n_0*r_h*r_h*(r_h/a_0)*J_m;

	double integral_right_1_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);

		double integral_upper_value = (1./3.)*pow(upper_limit_local,3.);
		double integral_lower_value = (1./3.)*pow(lower_limit_local,3.);

		//printf("Right_lower_limit = %f\n", lower_limit_local);
		//printf("Right_upper_limit = %f\n", upper_limit_local);

		double temp_result = (1./(r_h*r_h*r_h))*constant_pos_hr*integral_sigma*(integral_upper_value - integral_lower_value);
		return temp_result; 
	 }

	double integral_right_1_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);

		double integral_upper_value = (((slope_right*pow(upper_limit_local,4.))/4.)+(((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))*pow(upper_limit_local,3.))/3.));
		double integral_lower_value = (((slope_right*pow(lower_limit_local,4.))/4.)+(((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))*pow(lower_limit_local,3.))/3.));
		double temp_result = (1./(r_h*r_h*r_h))*constant_pos_hr*(integral_upper_value - integral_lower_value);
	
		//printf("result_right_function = %f\n", temp_result);

		return temp_result; 
	}

//x_right_1 is to the right of (1.8*r_h)
	if(x_right_1 > upper_limit)
	 {
		//printf("Check_1------------\n");

		double temp_result = integral_right_1_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
		
		//double temp_1 = (1./(r_h*r_h*r_h))*(constant_pos_hr*surface_density_lower)*(1./3.)*pow(300.-0.,3.);
		//double temp_2 = integral_right_1_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
		//return temp_1 + temp_2;

		//printf("Right_result = %f\n", temp_result);
		return temp_result;
	 }

//x_right_1 is to the left of (1.8*r_h) and x_right_2 is to the right of (1.8*r_h)
	else if(x_right_1 < upper_limit && x_right_2 > upper_limit)
	 {
		printf("Check_2------------\n");
		double temp_result_1 = integral_right_1_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_1_for_function_sigma(x_right_1, upper_limit);
		return temp_result_1 + temp_result_2;
	 }

//both x_right_1 and x_rigth_2 are to the left of (1.8*r_h)
	else if(x_right_2 < upper_limit)
	 {
		printf("Check_3------------\n");
		double temp_result_1 = integral_right_1_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_1_for_function_sigma(x_right_1, x_right_2);
		double temp_result_3 = integral_right_1_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_right_1");
		return 0.;
	 }	
}


double result_right_2()
{
	double lower_limit = 1.8*r_h;
	double upper_limit = 2.5*r_h;
	double constant_pos_cdr = 3.*n_0*r_h*r_h*(r_h/(2.*a_0))*J_m;

	double integral_right_2_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);

		//printf("Right Constant: lower_limit, upper_limit = %16.16f\t%16.16f\n", lower_limit_local, upper_limit_local);
		//printf("Constant sigma: right_lower_limit, right_upper_limit = %f\t%f\n", lower_limit_local, upper_limit_local);

		double integral_upper_value = (1./3.)*pow(upper_limit_local,3.);
		double integral_lower_value = (1./3.)*pow(lower_limit_local,3.);
		double temp_result = (1./(r_h*r_h*r_h))*constant_pos_cdr*integral_sigma*(integral_upper_value - integral_lower_value);
		
		//printf("result_right_2_const = %f\n", temp_result);
	
		return temp_result; 
	 }
	
	double integral_right_2_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);
	
		//printf("Right function: lower_limit, upper_limit = %16.16f\t%16.16f\n", lower_limit_local, upper_limit_local);
		//printf("Function sigma: right_lower_limit, right_upper_limit = %f\t%f\n", lower_limit_local, upper_limit_local);

		double integral_upper_value = (((slope_right*pow(upper_limit_local,4.))/4.)+(((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))*pow(upper_limit_local,3.))/3.));
		double integral_lower_value = (((slope_right*pow(lower_limit_local,4.))/4.)+(((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))*pow(lower_limit_local,3.))/3.));
		double temp_result = (1./(r_h*r_h*r_h))*constant_pos_cdr*(integral_upper_value - integral_lower_value);
		return temp_result;
	 }

//both x_right_1 and x_right_2 are to the right of (2.5*r_h)
	if(x_right_1 > upper_limit)
	 {
		printf("CHECK_RIGHT_1\n");
		return integral_right_2_for_constant_sigma(lower_limit, upper_limit, surface_density_lower);
	 }

//x_right_1 is in between (1.8*r_h) & (2.5*r_h) and x_right_2 is to the right of (2.5*r_h)
	else if(x_right_1 > lower_limit && x_right_1 < upper_limit && x_right_2 > upper_limit)
	 {
		//printf("CHECK_RIGHT_2\n");
		double temp_result_1 = integral_right_2_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_2_for_function_sigma(x_right_1, upper_limit);

		//printf("Integral_constant_right = %f\n", temp_result_1);
		//printf("Integral_function_right = %f\n", temp_result_2);

		return temp_result_1 + temp_result_2;
	 }

//both x_right_1 and x_right_2 are in between (1.8*r_h) & (2.5*r_h)
	else if(x_right_1 > lower_limit && x_right_2 < upper_limit)
	 {
		printf("CHECK_RIGHT_3\n");
		double temp_result_1 = integral_right_2_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_2_for_function_sigma(x_right_1, x_right_2);
		double temp_result_3 = integral_right_2_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//x_right_1 is to the left of (1.8*r_h) and x_right_2 is in between (1.8*r_h) & (2.5*r_h)
	else if(x_right_1 < lower_limit && x_right_2 > lower_limit && x_right_2 < upper_limit)
	 {
		printf("CHECK_RIGHT_4\n");
		double temp_result_1 = integral_right_2_for_function_sigma(lower_limit, x_right_2);
		double temp_result_2 = integral_right_2_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2;
	 }
//both x_right_1 and x_right_2 are to the left of (1.8*r_h)
	else if(x_right_2 < lower_limit)
	 {
		printf("CHECK_RIGHT_5\n");
		return integral_right_2_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 }

//x_right_1 is to the left of (1.8*r_h) and x_right_2 is to the right of (2.5*r_h)
	else if(x_right_1 < lower_limit && x_right_2 > upper_limit)
	 {
		printf("CHECK_RIGHT_6\n");
		return integral_right_2_for_function_sigma(lower_limit, upper_limit);
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_right_2");
		return 0.;
	 }
}


double result_right_3()
{
	double lower_limit = 2.5*r_h;
	double upper_limit = boxsize_x/2.;
	double constant_pos_mdr = 3.*n_0*r_h*r_h*((64.*pow(G*moonlet_mass,2.)*a)/(243.*pow(OMEGA,3.)))*(pow(((2.*K_0)+(K_1)),2.));

	double integral_right_3_for_constant_sigma(double integral_lower_limit, double integral_upper_limit, double integral_sigma)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);

		double integral_lower_value = (-(1./(3.*pow(lower_limit_local,3.)))-(alpha/(2.*a*pow(lower_limit_local,2.))));
		double integral_upper_value = (-(1./(3.*pow(upper_limit_local,3.)))-(alpha/(2.*a*pow(upper_limit_local,2.))));

		if(integral_upper_limit == upper_limit)
		 {
			integral_upper_value = 0.;
		 }

		double temp_result = (1./(r_h*r_h))*constant_pos_mdr*integral_sigma*(integral_upper_value-integral_lower_value);
	 	return temp_result;
	 }

	double integral_right_3_for_function_sigma(double integral_lower_limit, double integral_upper_limit)
	 {
		double lower_limit_local = check_limit_sign(integral_lower_limit);
		double upper_limit_local = check_limit_sign(integral_upper_limit);	

		double integral_lower_value = (-((slope_right)/(2.*pow(lower_limit_local,2.)))-((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))/(3.*pow(lower_limit_local,3.)))-((alpha*slope_right)/(a*lower_limit_local))-((alpha*(-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x)))/(2.*a*pow(lower_limit_local,2.))));
		double integral_upper_value = (-((slope_right)/(2.*pow(upper_limit_local,2.)))-((-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x))/(3.*pow(upper_limit_local,3.)))-((alpha*slope_right)/(a*upper_limit_local))-((alpha*(-(slope_right*x_right_1)+surface_density_lower+(slope_right*moonlet_x)))/(2.*a*pow(upper_limit_local,2.))));

		if(integral_upper_limit == upper_limit)
		 {
			integral_upper_value = 0.;
		 }

		double temp_result = (1./(r_h*r_h))*constant_pos_mdr*(integral_upper_value - integral_lower_value);
		return temp_result; 
	}

//both x_right_1 and x_right_2 are to the left of (2.5*r_h)
	if(x_right_2 <= lower_limit)
	 {
		//printf("1.............\n");
		double temp_1 = integral_right_3_for_constant_sigma(lower_limit, upper_limit, surface_density_upper);
	 	return temp_1;
	 }
	
//x_right_1 is to the left of (2.5*r_h) and x_right_2 is to the right of (2.5*r_h)
	else if(x_right_1 < lower_limit && x_right_2 > lower_limit)
	 {
		//printf("2.............\n");
		double temp_result_1 = integral_right_3_for_function_sigma(lower_limit, x_right_2);
		double temp_result_2 = integral_right_3_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);
		return temp_result_1 + temp_result_2;
	 }

//both x_right_1 and x_right_2 are to the right of (2.5*r_h)
	else if(x_right_1 > lower_limit)
	 {
		//printf("3.............\n");
		double temp_result_1 = integral_right_3_for_constant_sigma(lower_limit, x_right_1, surface_density_lower);
		double temp_result_2 = integral_right_3_for_function_sigma(x_right_1, x_right_2);
		double temp_result_3 = integral_right_3_for_constant_sigma(x_right_2, upper_limit, surface_density_upper);

		//printf("Temp_result_1_right = %f\n", temp_result_1);
		//printf("Temp_result_2_right = %f\n", temp_result_2);
		//printf("Temp_result_3_right = %f\n", temp_result_3);

		return temp_result_1 + temp_result_2 + temp_result_3;
	 }

//none of the above is true
	else
	 {
		printf("Missing in result_right_3");
		return 0.;
	 }
}


double torque_function()
{
	double temp_result_left_1 = fabs(result_left_1());
	double temp_result_left_2 = -fabs(result_left_2());
	double temp_result_left_3 = -fabs(result_left_3());
	double result_torque_left = temp_result_left_1 + temp_result_left_2 + temp_result_left_3;
	double temp_result_right_1 = fabs(result_right_1());
	double temp_result_right_2 = fabs(result_right_2());
	double temp_result_right_3 = -fabs(result_right_3());
	double result_torque_right = temp_result_right_1 + temp_result_right_2 + temp_result_right_3;
	double result_torque = result_torque_left + result_torque_right;	

	printf("Torque_1 = %16.16e Nm \n", temp_result_left_1);
	printf("Torque_2 = %16.16e Nm \n", temp_result_left_2);
	printf("Torque_3 = %16.16e Nm \n", temp_result_left_3);
	printf("Torque_4 = %16.16e Nm \n", temp_result_right_1);
	printf("Torque_5 = %16.16e Nm \n", temp_result_right_2);
	printf("Torque_6 = %16.16e Nm \n", temp_result_right_3);
	printf("Torque (-ve x) = %16.16e Nm\n", result_torque_left);
	printf("Torque (+ve x) = %16.16e Nm\n", result_torque_right);
	return result_torque;		
}


double migration_rate_function();
{
	double result_torque = torque_function();
	double result_migration_rate = (2.*result_torque)/(moonlet_mass*n_0*a_0);		// m/s
	result_migration_rate = (result_migration_rate/1000.)*(365.*24.*60.*60.);
	printf("Calculated Torque = %e Nm \n", result_torque);
	printf("Calculated Migration rate = %16.16f km/yr\n", result_migration_rate);
	//return result_migration_rate;
}

exit(0);

while(mass<total_mass)
{
        temp_x = tools_uniform(-boxsize_x/2.,boxsize_x/2.);
        temp_sigma = tools_uniform(0, (surface_density_upper+50.));
        result = check(temp_x, temp_sigma);
        if(result)
         {
                struct particle pt;
                pt.x                            = temp_x;
                pt.y                            = tools_uniform(-boxsize_y/2.,boxsize_y/2.);
                pt.z                            = tools_normal(1.);
                pt.vx                           = 0;
                pt.vy                           = -1.5*pt.x*OMEGA;
                pt.vz                           = 0;
                pt.ax                           = 0;
                pt.ay                           = 0;
                pt.az                           = 0;
                double radius                   = tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
		//#ifndef COLLISIONS_NONE
                pt.r                            = radius;
		//#endif
                double particle_mass            = particle_density*4./3.*M_PI*radius*radius*radius;
                pt.m                            = particle_mass;
                particles_add(pt);
                mass                           += particle_mass;
                number_of_particles++;
                number_density                  = total_mass/(particle_mass*surface_area);
                fprintf(fp, "%f \t %f \t %f \n", pt.x, temp_sigma, number_density);
         }
}

        struct particle pt;
        pt.x                                    = moonlet_x;
        pt.y                                    = tools_normal(1.);        
        pt.z                                    = tools_normal(1.);
        pt.vx                                   = 0;
        pt.vy                                   = -1.5*pt.x*OMEGA;
        pt.vz                                   = 0;
        pt.ax                                   = 0;
        pt.ay                                   = 0;
        pt.az                                   = 0;
        pt.r                                    = moonlet_radius;
        pt.m                                    = moonlet_mass;
        particles_add(pt);
        number_of_particles++;

        fclose(fp);


double sigma_function(double x)
{
	if(x < x_left_1)
	 {
		return surface_density_upper;
	 }
	
	if(x >= x_left_1 && x < x_left_2)
	 {
		return ((slope_left*x)-(slope_left*x_left_1)+surface_density_upper);
	 }

	if(x >= x_left_2 && x < x_right_1)
	 {
		return surface_density_lower;
	 }

	if(x >= x_right_1 && x < x_right_2)
	 {
		return ((slope_right*x)-(slope_right*x_right_1)+surface_density_lower);
	 }
	
	if(x >= x_right_2)
	 {
		return surface_density_upper;
	 }
}


double toomre_wavelength_function(double x)
{
        return (2.*M_PI*M_PI*(sigma_function(x))/OMEGA/OMEGA*G);
}

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
