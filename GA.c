// File: GA.c
//
// Description: GA (genetic algorithm) ADT implementation
//
// @author Maximus Milazzo
//
///////////////////////////////////////////////////////////



#include "GA.h"



/// Defines the genetic algorithm struct
struct GA_s {

  double ** parents;
  double ** children;
  double * fitness;
  // population data

  double (*fitness_calc)(double * entity);
  // fitness calculation function

  void (*crossover)(GA optimizer, double * parent_1, double * parent_2, 
    int child_index);
    // function pointer to be set to crossover function corresponding with
    // selected crossover mode

  double gene_combination_chance;
  double mutation_chance;
  double mutation_factor;
  double gene_swap_chance;
  // double-type schema

  int num_parents;
  int num_children;
  int num_genes;
  // integer-type schema
  
};



/// Swaps a parent with the next indexed parent in the population list, as well
/// as associated fitness values.
///
/// @param optimizer - the GA instance
/// @param index - index of parent to swap

static void swap_parents(GA optimizer, int index) {

  double temp_fit = optimizer->fitness[index];
  optimizer->fitness[index] = optimizer->fitness[index + 1];
  optimizer->fitness[index + 1] = temp_fit;
  // swaps fitness values

  for (int i = 0; i < optimizer->num_genes; ++i) {
    
    double temp_gene = optimizer->parents[index][i];
    optimizer->parents[index][i] = optimizer->parents[index + 1][i];
    optimizer->parents[index + 1][i] = temp_gene;
    // swaps all genes between parents

  }
  
}



/// Shifts the parent list so that a newly added parent to position 0 will be
/// in the correcltly ordered place based on fitness values.
///
/// @param optimizer - the GA instance

static void shift_parents(GA optimizer) {

  for (int i = 0; i < optimizer->num_parents - 1; ++i) {

    if (optimizer->fitness[i] > optimizer->fitness[i + 1]) {
    
      swap_parents(optimizer, i);
      // swaps parents so that they are fitness ordered --
      // larger fitness values are said to be more "fit"

    } else {

      return;
      // once no more swapping is needed, the function returns
    
    }

  }
  
}



/// Performs a uniform type genetic crossover for breeding.
///
/// @param optimizer - the GA instance
/// @param parent_1 - the first breeding parent
/// @param parent_2 - the second breeding parent
/// @param child_index - the resuling child index

static void crossover_uniform(GA optimizer, double * parent_1,
  double * parent_2, int child_index) {

  for (int j = 0; j < optimizer->num_genes; ++j ) {
  
    if (optimizer->gene_combination_chance > ((double) rand()) / 
      ((double) RAND_MAX)) {
      // gene combination case

      optimizer->children[child_index][j] = (parent_1[j] + parent_2[j]) / 2;
        
    } else if (0.5 > ((double) rand()) / ((double) RAND_MAX)) {
      // parent 1 dominance case

      optimizer->children[child_index][j] = parent_1[j];
          
    } else {
      // parent 2 dominance case

      optimizer->children[child_index][j] = parent_2[j];
  
    }
      
  }
  
}



/// Performs a single point type genetic crossover for breeding.
///
/// @param optimizer - the GA instance
/// @param parent_1 - the first breeding parent
/// @param parent_2 - the second breeding parent
/// @param child_index - the resuling child index

static void crossover_single_point(GA optimizer, double * parent_1, 
  double * parent_2, int child_index) {

  int cross_point_index = rand() % optimizer->num_genes;
  // random cross point

  for (int i = 0; i < cross_point_index; ++i)
    optimizer->children[child_index][i] = parent_1[i];
    // sets values up until cross point to parent 1 genes

  for (int j = cross_point_index; j < optimizer->num_genes; ++j)
    optimizer->children[child_index][j] = parent_2[j];
    // sets values at and after cross point to parent 2 genes
  
}



/// Performs a double point type genetic crossover for breeding.
///
/// @param optimizer - the GA instance
/// @param parent_1 - the first breeding parent
/// @param parent_2 - the second breeding parent
/// @param child_index - the resuling child index

static void crossover_double_point(GA optimizer, double * parent_1, 
  double * parent_2, int child_index) {

  int cross_point_upper_index;
  int cross_point_lower_index;
  // variables to hold cross point indices
  
  int cross_point_1_index = rand() % optimizer->num_genes;
  int cross_point_2_index = rand() % optimizer->num_genes;
  // generates two ransom cross points

  if (cross_point_1_index > cross_point_2_index) {

    cross_point_upper_index = cross_point_1_index;
    cross_point_lower_index = cross_point_2_index;
    
  } else {
    
    cross_point_upper_index = cross_point_2_index;
    cross_point_lower_index = cross_point_1_index;
    
  }
  // sets upper and lower index values

  for (int i = 0; i < cross_point_lower_index; ++i)
    optimizer->children[child_index][i] = parent_1[i];
    // sets values up until lower cross point to parent 1 genes

  for (int j = cross_point_lower_index; j < cross_point_upper_index; ++j)
    optimizer->children[child_index][j] = parent_2[j];
    // sets values at and after lower cross point to parent 2 genes up until
    // upper cross point
    

  for (int k = cross_point_upper_index; k < optimizer->num_genes; ++k)
    optimizer->children[child_index][k] = parent_1[k];
    // sets values at and after upper cross point to parent 1 genes

}



/// Sets GA crossover function to reflect selected crossover mode.
///
/// @param optimizer - the GA instance
/// @param crossover_mode - the selected crossover mode

static void set_crossover_mode(GA optimizer, int crossover_mode) {

  switch (crossover_mode) {
    case 0:
      optimizer->crossover = crossover_uniform;
      break;
      // selects uniform crossover

    case 1:
      optimizer->crossover = crossover_single_point;
      break;
      // selects single point crossover

    default:
      optimizer->crossover = crossover_double_point;
      break;
      // selects double point crossover

  }

}



/// Initializes a new genetic algorithm (GA).

GA GA_init(int num_parents, int num_children, int num_genes, int crossover_mode,
  double gene_combination_chance, double mutation_chance, double mutation_factor, 
  double gene_swap_chance, double ** start_parents, int parent_val_range, 
  double (*fitness_calc)(double * entity)) {

  srand(time(NULL));
  // seeds random number generator

  GA optimizer = (GA) malloc(sizeof(struct GA_s));
  // allocates memory for GA instance

  optimizer->num_parents = num_parents;
  optimizer->num_children = num_children;
  optimizer->num_genes = num_genes;
  // sets integer-typed schema values
  
  optimizer->gene_combination_chance = gene_combination_chance;
  optimizer->mutation_chance = mutation_chance;
  optimizer->mutation_factor = mutation_factor;
  optimizer->gene_swap_chance = gene_swap_chance;
  // sets double-typed sceham values

  optimizer->parents = (double **) malloc(sizeof(double *) * num_parents);
  optimizer->children = (double **) malloc(sizeof(double *) * num_children);
  optimizer->fitness = (double *) calloc(sizeof(double), num_parents);
  // allocates memory for population data

  for (int i = 0; i < num_parents; ++i)
    optimizer->parents[i] = (double *) malloc(sizeof(double) * num_genes);
    // allocates memory for parents
  
  for (int j = 0; j < num_parents; ++j) {

    if (start_parents == NULL) {

      for (int k = 0; k < num_genes; ++k)
        optimizer->parents[0][k] = ((double) (rand() % parent_val_range))
          * ((((double) rand()) - ((double) (RAND_MAX / 2))) / 
          ((double) (RAND_MAX / 2)));
          // parent genes set to random values in range (-RANGE, +RANGE)

    } else {
      
      for (int k = 0; k < num_genes; ++k)
        optimizer->parents[0][k] = start_parents[j][k];
        // start parent values copied to GA instance parent values

    }
    
    optimizer->fitness[0] = fitness_calc(optimizer->parents[0]);//
    shift_parents(optimizer);
    // parents are properly ordered in terms of initial fitness
    
  }

  for (int k = 0; k < num_children; ++k)
    optimizer->children[k] = (double *) malloc(sizeof(double) * num_genes);
    // allocates memory for children
  
  optimizer->fitness_calc = fitness_calc;
  // sets fitness calculation function

  set_crossover_mode(optimizer, crossover_mode);
  // sets crossover mode

  return optimizer;
  
}



/// Initializes a new genetic algorithm (GA) with default parameters.

GA GA_default(double (*fitness_calc)(double * entity)) {

  return GA_init(DEFAULT_GA_PARAMS, fitness_calc);

}



/// Breeds parents to create new children.
///
/// @param optimizer - the GA instance

static void breed(GA optimizer) {

  double * parent_1;
  double * parent_2;
  // holds selected breeding parents
  
  for (int i = 0; i < optimizer->num_children; ++i) {

    parent_1 = optimizer->parents[rand() % optimizer->num_parents];
    parent_2 = optimizer->parents[rand() % optimizer->num_parents];
    // picks random parents

    optimizer->crossover(optimizer, parent_1, parent_2, i);
    // performs crossover based on selected mode
    
  }
  
}



/// Performs a gene swap between two children.
///
/// @param optimizer - the GA instance
/// @param child_1_index - the index of the first child
/// @param child_2_index - the index of the second child
/// @param gene_index - the index of the allele to swap between children

static void swap_genes(GA optimizer, int child_1_index, int child_2_index,
  int gene_index) {

  double temp_gene = optimizer->children[child_1_index][gene_index];
  
  optimizer->children[child_1_index][gene_index] = 
    optimizer->children[child_2_index][gene_index];

  optimizer->children[child_2_index][gene_index] = temp_gene;

}



/// Carries out gene mutations on children.
///
/// @param optimizer - the GA instance

static void mutate(GA optimizer) {

  for (int i = 0; i < optimizer->num_children; ++i) {

    for (int j = 0; j < optimizer->num_genes; ++j) {

      if (optimizer->mutation_chance > ((double) rand()) / 
        ((double) RAND_MAX)) {

        optimizer->children[i][j] += optimizer->children[i][j] * 
          optimizer->mutation_factor * ((((double) rand()) - 
          ((double) (RAND_MAX / 2))) / ((double) (RAND_MAX / 2)));
          // shifts gene value up or down by mutation factor
          // multiplied by a random double value in range -1 to 1
      
      }

      if (optimizer->gene_swap_chance > ((double) rand()) / 
        ((double) RAND_MAX)) {

        swap_genes(optimizer, i, rand() % optimizer->num_children, j);
        // swaps genes between current child and another randomly selected child
        
      }
    
    }

  }
  
}



/// Selects surviving members of a generation that will become parents in the
/// next generation.
///
/// @param optimizer - the genetic algorithm optimizer

static void select_parents(GA optimizer) {
  
  double cur_fit;

  for (int i = 0; i < optimizer->num_children; ++i) {
    
    cur_fit = optimizer->fitness_calc(optimizer->children[i]);

    if (cur_fit > optimizer->fitness[0]) {

      optimizer->fitness[0] = cur_fit;

      for(int j = 0; j < optimizer->num_genes; ++j)
        optimizer->parents[0][j] = optimizer->children[i][j];
        // sets new parent to the first index (lowest "fitness" spot)

      shift_parents(optimizer);
      
    }

  }
  
}



/// Runs the genetic algorithm for a single generation.
///
/// @param optimizer - the GA instance

static void run_generation(GA optimizer) {
  
  breed(optimizer);
  mutate(optimizer);
  select_parents(optimizer);

}



/// Runs the genetic algorithm.

double * GA_run(GA optimizer, int generations) {

  for (int i = 0; i < generations; ++i)
    run_generation(optimizer);

  return optimizer->parents[optimizer->num_parents - 1];
  // returns the most "fit" parent from the last generation

}



/// Deallocates memory associated with GA instance.

void GA_destroy(GA optimizer) {

  for (int i = 0; i < optimizer->num_parents; ++i)
    free(optimizer->parents[i]);
    // deallocates individual parents

  for (int j = 0; j < optimizer->num_children; ++j)
    free(optimizer->children[j]);
    // deallocates individual children

  free(optimizer->parents);
  free(optimizer->children);
  free(optimizer->fitness);
  // deallocates population data
  
  free(optimizer);
  // deallocates GA instance
  
}
