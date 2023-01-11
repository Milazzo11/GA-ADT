// File: test.c
//
// Description: creates GA instance and solves basic test problem
//
// @author Maximus Milazzo
//
///////////////////////////////////////////////////////////



#include <stdio.h>
#include <float.h>
#include <math.h>
#include "GA.h"



/// Fitness function to evaluate variables in the equation: 3x + 5y = 1.
/// The closer this equation is to being true given values of "x" and "y"
/// derived from the "entity" parameter, the larger the returned "fitness" value.
///
/// @param entity - the individual in the genetic algorithm's population whose
///    fitness value is being determined
///
/// @return the entity's fitness value

double eq_fitness(double * entity) {

  double res = 3 * entity[0] + 5 * entity[1] - 1;
  // finds result of expression: "3x + 5y - 1"

  if (res == 0)
    return DBL_MAX;
    // returns the maximum possible value if equation is solved
  
  return fabs(1 / res);
  // the reciprocal of the expression result is taken so that values closer to
  // solving the equation result in a higher fitness value being returned
  
}



/// Main entry point for GA test.
///
/// @return 0 on success or failure code on failure

int main(void) {
  
  GA eq_solver =  GA_default(eq_fitness);
  double * res = GA_run(eq_solver, 100);
  // creates and runs GA instance

  printf("x = %f\n", res[0]);
  printf("y = %f\n\n", res[1]);
  printf("3x + 5y = %f\n\n", res[0] * 3 + res[1] * 5);
  // displays results

  GA_destroy(eq_solver);
  // deallocates GA memory
  
  return 0;
  
}
