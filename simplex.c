// Simplex solver in C. Uses 2 phase simplex. Start reading at [solve](#section-39).

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <curses.h>

#define DATATYPE float

#include "parse.h"

// Verbose mode, will be set with a cmd line flag
int verbosemode = true;
int interactivemode = false;
#define VP(x, ...) { if (verbosemode) printf(x, __VA_ARGS__); }

// Type of solution
typedef enum {
  UNBOUNDED,
  BOUNDED,
  INFEASIBLE
} solution_type_t;

// An optimum solution which contains the point that is optimal
typedef struct {
  solution_type_t type;
  DATATYPE * x;
} optimum_solution_t;

// An unbounded solution returns a cone
typedef struct {
  solution_type_t type;
  DATATYPE * x;
  DATATYPE d;
} unbounded_solution_t;

// A solution, which is either optimum, unbounded, or infeasible
typedef union {
  solution_type_t type;
  optimum_solution_t optimum;
  unbounded_solution_t unbounded;
} solution_t;

// Table representing simplex state
typedef struct {
  DATATYPE ** values;
  int * basic;
  int rows;
  int cols;
  int artificial;
  DATATYPE * c;
} tableau_t;

// ## check_unbounded
// Checks to see if a table is currently representing an unbounded solution
bool
check_unbounded (tableau_t * tableau) {
  
}

// ## check_infeasible
// Checks to see if a table is currently showing that the problem is infeasible
bool
check_infeasible (tableau_t * tableau, int entering_var) {
  int i;
  for (i = 1; i < tableau->rows - 1; i++) {
    if (tableau->values[i][entering_var] >= 0) return false;
  }
  return true;
}

// ## check_optimum
// Returns the index of the most negative cost or -1 if the problem is optimal
int
check_optimum (tableau_t * tableau) {
  DATATYPE min = tableau->values[0][0];
  int index = 0;
  int i;
  for (i = 0; i < tableau->cols - 1; i++) {
    if (tableau->values[0][i] < min) {
      min = tableau->values[0][i];
      index = i;
    }
  }
  if (min >= 0) return -1;
  else return index;
}

// ## mrt
// Runs the minimum ratio test and returns the index of the leaving var
int
mrt (tableau_t * tableau, int entering_var) {
  int i;
  DATATYPE min = INFINITY;
  int index = 1;
  for (i = 1; i < tableau->rows; i++) {
    DATATYPE ratio = tableau->values[i][tableau->cols - 1] / tableau->values[i][entering_var];
    if (ratio < 0) continue;
    if (ratio < min) {
      min = ratio;
      index = i;
    }
  }
  return index - 1;
}

// ## normalize_row
// Normalizes a row of the table
void
normalize_row (tableau_t * tableau, int row, int col) {
  DATATYPE factor = tableau->values[row][col];
  int i;
  for (i = 0; i < tableau->cols; i++) {
    tableau->values[row][i] = tableau->values[row][i] / factor;
  }
}

// ## gaussian_col
// Makes all of the other cols 0
void
gaussian_col (tableau_t * tableau, int row, int col) {
  int i, j;
  // For each row, if it's not the row of the leaving var, make it a 0 by
  // subtracting a multiple of the row of the leaving var.
  for (i = 0; i < tableau->rows; i++) {
    if (i != row) {
      DATATYPE factor = tableau->values[i][col];
      for (j = 0; j < tableau->cols; j++) {
        tableau->values[i][j] = tableau->values[i][j] - (factor * tableau->values[row][j]);
      }
    }
  }
}

// ## print_tableau
// Prints a tableau
void
print_tableau (tableau_t * tableau) {
  if (!verbosemode) return;
  int i, j;
  for (i = 0; i < tableau->rows; i++) {
    for (j = -1; j < tableau->cols; j++) {
      if (i == 0 && j == -1) printf("\t");
      if (j == -1 && i > 0) printf("x%d\t", tableau->basic[i-1] + 1);
      else printf("");
      if (j != -1) printf("%3g\t", tableau->values[i][j]);
    }
    printf("\n");
  }
}

// ## update_basic_vars
// Replaces a basic var with a different basic var
void
update_basic_vars(tableau_t * tableau, int leaving, int entering) {
  int i;
  for (i = 0; i < tableau->rows - 1; i++) {
    if (tableau->basic[i] == tableau->basic[leaving]) {
      tableau->basic[i] = entering;
    }
  }
}

// ## simplex_iteration
// Single simplex iteration. If soln is optimum, will return same table
tableau_t *
simplex_iteration (tableau_t * tableau) {
  int entering_var = check_optimum(tableau);
  int leaving_idx = mrt(tableau, entering_var);

  VP("leaving idx: %d var: %d, entering var: %d\n", leaving_idx, tableau->basic[leaving_idx] + 1, entering_var + 1);
  if (verbosemode) {
    printf("Current basis: ");
    int i;
    for (i = 0; i < tableau->rows - 1; i++) printf(" x%d", tableau->basic[i] + 1);
    printf("\n");
  }

  // make the leaving var, entering var element a 1
  normalize_row(tableau, leaving_idx + 1, entering_var);

  // make the columns of the entering var 0
  gaussian_col(tableau, leaving_idx + 1, entering_var);

  update_basic_vars(tableau, leaving_idx, entering_var);

  return tableau;
}

// ## simplex
// Runs simplex on the given table and returns the table itself. Must have an
// initial basis setup
tableau_t *
simplex (tableau_t * tableau) {
  while (true) {
    int entering_var;
    if ((entering_var = check_optimum(tableau)) == -1) break;
    if (check_infeasible(tableau, entering_var)) break;

    tableau = simplex_iteration(tableau);

    if (verbosemode) {
      print_tableau(tableau);
      char buf[2];
      if (interactivemode) fgets(buf, 2, stdin);
    }
  }

  return tableau;
}

// ## is_ident
// Checks to see if the index row of the table looks like the identity. If so
// returns what row of the identity. if not returns -1
int
is_ident (tableau_t * tableau, int index) {
  int found_index = -1, i;
  for (i = 1; i < tableau->rows; i++) {
    DATATYPE val = tableau->values[i][index];
    if (val == 1) {
      if (found_index != -1) return -1;
      else found_index = i - 1;
    }
    else if (val == 0) continue;
    else return -1;
  }
  return found_index;
}

// ## find_basis
// Returns the necessary # of artificial vars
int
find_basis (tableau_t * tableau) {
  // tableau->basic should be of the right size, but all -1
  // Loop over all of the columns and look for something that looks like the
  // identity.

  int i, num_bas = 0;
  for (i = 0; i < tableau->cols - 1; i++) {

    // If a col looks like the identity, num_bas++ and add it to the correct
    // spot in the basic vars array.

    int index = is_ident(tableau, i);
    if (index != -1) {
      num_bas++;
      tableau->basic[index] = i;
    }
  }

  // VP("num vars found for basis: %d\n", num_bas);

  // If we end up with num_bas == rows - 1, we've found a basis, return TRUE.
  if (num_bas == tableau->rows - 1) return 0;

  // Otherwise we haven't found a basis. Return the # of artifical vars that
  // will now be necessary in order to find a basis for phase-1 of two phase
  // simplex.
  else return tableau->rows - 1 - num_bas;
}

// ## artificialize
// Checks for a basis. Returns TWOPHASE if two phases must be run and ONEPHASE
// if only one phase must be run (ie, there's a basis already).
enum {TWOPHASE, ONEPHASE}
artificialize (tableau_t * tableau) {
  int num_artificial = find_basis(tableau);
  VP("Need %d artificial variables.\n", num_artificial);

  // If we need no artifical vars, it's a single phase problem.
  if (num_artificial == 0) return ONEPHASE;
  else {
    tableau->artificial = num_artificial;
    // realloc tableau->values to be of the right size. It'll need
    // num_artifical more cols.

    int i;
    for (i = 0; i < tableau->rows; i++) {
      tableau->values[i] = realloc(
        tableau->values[i], 
        sizeof(DATATYPE) * (tableau->cols + num_artificial)
      );

      // Copy right hand side over
      tableau->values[i][num_artificial + tableau->cols - 1] = tableau->values[i][tableau->cols - 1];

      // initialize the columns corresponding to the artificial vars to 0
      int j;
      for (j = tableau->cols - 1; j < num_artificial + tableau->cols - 1; j++) {
        tableau->values[i][j] = 0;
      }
    }

    // Reset objective row.
    for (i = 0; i < tableau->cols; i++) tableau->values[0][i] = 0;

    int curr_artificial = tableau->cols - 1;
    // Loop over basic vars to find what cols of the identity we need so
    // we can figure out where the artifical vars need to be added.
    for (i = 0; i < tableau->cols - 1; i++) { // this is wrong
      if (tableau->basic[i] == -1) {
        tableau->basic[i] = curr_artificial;
        // Means we need a col of the ident with the 1 in position i
        tableau->values[i + 1][curr_artificial] = 1;
        curr_artificial += 1;

        // Add the contents of this row to the objective row since it 
        // corresponds to a row of the matrix that corresponds to an artificial
        // variable.
        int j;
        for (j = 0; j < tableau->cols + num_artificial; j++) 
          tableau->values[0][j] += tableau->values[i + 1][j];
      }
    }

    for (i = 0; i < tableau->cols + num_artificial; i++) 
      tableau->values[0][i] = -tableau->values[0][i];

    for (i = tableau->cols - 1; i < tableau->cols + num_artificial - 1; i++)
      tableau->values[0][i] = 0;

    // Fix table state
    tableau->cols += num_artificial;

    return TWOPHASE;
  }
}

// ## deartificialize
// Undoes the artificialization, reconstructing the objective row
void
deartificialize (tableau_t * tableau) {
  tableau->cols -= tableau->artificial;
  int i, j;
  DATATYPE product;

  // Move RHS over
  for (i = 0; i < tableau->rows; i++)
    tableau->values[i][tableau->cols - 1] = tableau->values[i][tableau->cols + tableau->artificial - 1];

  // Reconstruct objective row
  for (i = 0; i < tableau->cols - 1; i++) {
    product = 0;
    for (j = 1; j < tableau->rows - 1; j++)
      product += (tableau->c[tableau->basic[j - 1]] * tableau->values[j][i]);
    product -= tableau->c[i];
    tableau->values[0][i] = product;
  }

  product = 0;
  for (j = 1; j < tableau->rows - 1; j++)
    product += (tableau->c[tableau->basic[j - 1]] * tableau->values[j][tableau->cols - 1]);
  tableau->values[0][tableau->cols - 1] = product;
}

bool
verify_solution (tableau_t * tableau, DATATYPE ** A, DATATYPE * b, int num_vars, int num_constraints) {
  int i;
  DATATYPE * x = calloc(num_vars, sizeof(DATATYPE));
  for (i = 0; i < tableau->rows - 1; i++) {
    x[tableau->basic[i]] = tableau->values[i + 1][tableau->cols - 1];
  }

  DATATYPE accum;
  int j;
  for (i = 0; i < num_constraints; i++) {
    accum = 0;
    for (j = 0; j < num_vars; j++) {
      accum += A[i][j] * x[j];
    }
    if (accum != b[i]) printf("Solution FAILS constraint %d! got %g, expected %g!\n", i, accum, b[i]);
  }

  accum = 0;
  for (i = 0; i < num_vars; i++) {
    accum += tableau->c[i] * x[i];
  }
  if (accum != tableau->values[0][tableau->cols - 1]) 
    printf("INCORRECT objective value in table. got %g, expected %g", tableau->values[0][tableau->cols - 1], accum);
}

// ## solve
// Actually solves an LP problem
solution_t *
solve (DATATYPE ** A, DATATYPE * b, DATATYPE * c, int num_vars, int num_constraints) {
  tableau_t * tableau = malloc(sizeof(tableau_t));
  tableau->values = calloc(num_constraints + 1, sizeof(DATATYPE *));
  tableau->basic = calloc(num_constraints, sizeof(DATATYPE));
  tableau->rows = num_constraints + 1;
  tableau->cols = num_vars + 1;
  tableau->artificial = 0;
  tableau->c = c;

  int i, x;

  for (i = 0; i < num_constraints + 1; i++) {
    tableau->values[i] = calloc(num_vars + 1, sizeof(DATATYPE));
  }

  for (i = 0; i < num_constraints; i++)
    tableau->basic[i] = -1;

  // setup simplex table
  for (i = 0; i < num_vars; i++) {
    tableau->values[0][i] = -c[i];
  }

  int j;
  for (i = 0; i < num_constraints; i++)
    for (j = 0; j < num_vars; j++)
      tableau->values[i + 1][j] = A[i][j];

  for (i = 0; i < num_constraints; i++)
    tableau->values[i + 1][num_vars] = b[i];

  VP("initial tableau:\n", "");
  print_tableau(tableau);

  // ### The important part.
  if (artificialize(tableau) == TWOPHASE) {
    VP("tableau after artificialized:\n", "");
    print_tableau(tableau);
    simplex(tableau);
    VP("tableau after phase 1:\n", "");
    print_tableau(tableau);
    deartificialize(tableau);
    VP("tableau after deartificialization:\n", "");
    print_tableau(tableau);
  }

  simplex(tableau);

  int ent;
  if ((ent = check_optimum(tableau)) == -1) {
    for (i = 0; i < tableau->rows - 1; i++) { 
      printf("x%d\t%g\n", tableau->basic[i] + 1, tableau->values[i + 1][tableau->cols - 1]);
    }
    printf("OBJECTIVE: %g\n", tableau->values[0][tableau->cols - 1]);
    verify_solution(tableau, A, b, num_vars, num_constraints);
  } else if (check_infeasible(tableau, ent)) {
    printf("INFEASIBLE!");
  }
}

int
main(int argc, const char *argv[]) {

/*
  float a1[] = {1, 1, 1, 1, 0, 0};
  float a2[] = {1, 3, 0, 0, 1, 0};
  float a3[] = {0, 0, 1, 0, 0, 1};
  float * A[] = {a1, a2, a3};

  float c[] = {1, 2, 3, 0, 0, 0};
  float b[] = {5, 6, 3};
  solve(A, b, c, 6, 3);
*/
/*
  float a1[] = {1, 1};
  float a2[] = {2, 1};
  float * A[] = {a1, a2};
  float c[] = {1, 1};
  float b[] = {6, 8};
  solve(A, b, c, 2, 2);
*/

/*
  float a1[] = {1, 1, -1, 0, 0};
  float a2[] = {0, 1, 0, -1, 0};
  float a3[] = {1, 0, 0, 0, 1};
  float * A[] = {a1, a2, a3};
  float c[] = {1, -1, 0, 0, 0};
  float b[] = {3, 2, 4};

  solve(A, b, c, 5, 3);
  */


  float ** A, * b, * c;
  char * A_str, * b_str, * c_str;
  int num_vars, num_constraints;

  get_file("A.txt", &A_str);
  get_file("b.txt", &b_str);
  get_file("c.txt", &c_str);
  parse_matrix(&num_constraints, &num_vars, &A, A_str);

  parse_rows(&b, b_str, "\n");
  parse_rows(&c, c_str, "\n");

  printf("vars: %d, constraints: %d\n", num_vars, num_constraints);
  solve(A, b, c, num_vars, num_constraints);

  return 0;
}

