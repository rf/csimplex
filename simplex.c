#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <curses.h>

#define DATATYPE float

#define DEBUGMODE true

#define DP(x, ...) { if (DEBUGMODE) printf(x, __VA_ARGS__); }

typedef enum {
  UNBOUNDED,
  BOUNDED,
  INFEASIBLE
} solution_type_t;

typedef struct {
  solution_type_t type;
  DATATYPE * x;
} optimum_solution_t;

typedef struct {
  solution_type_t type;
  DATATYPE x;
  DATATYPE d;
} unbounded_solution_t;

typedef union {
  solution_type_t type;
  optimum_solution_t optimum;
  unbounded_solution_t unbounded;
} solution_t;

typedef struct {
  DATATYPE ** values;
  int * basic;
  int rows;
  int cols;
} tableau_t;

bool
check_unbounded (tableau_t * tableau) {
  
}

bool
check_infeasible (tableau_t * tableau) {

}

// Returns the index of the most negative cost or -1 if the problem is optimal
int
check_optimum (tableau_t * tableau) {
  DATATYPE min = tableau->values[0][0];
  int index = 0;
  int i;
  for (i = 0; i < tableau->rows; i++) {
    if (tableau->values[0][i] < min) {
      min = tableau->values[0][i];
      index = i;
    }
  }
  if (min >= 0) return -1;
  else return index;
}

int
mrt (tableau_t * tableau, int entering_var) {
  int i;
  DATATYPE min = tableau->values[1][tableau->cols - 1] / tableau->values[1][entering_var];
  int index = 1;
  for (i = 1; i < tableau->rows; i++) {
    DATATYPE ratio = tableau->values[i][tableau->cols - 1] / tableau->values[i][entering_var];
    DP("got ratio %f for row %d\n", ratio, i);
    if (ratio < 0) continue;
    if (ratio < min) {
      DP("new idx is %d\n", i);
      min = ratio;
      index = i;
    }
  }
  return index - 1;
}

void
normalize_row (tableau_t * tableau, int row, int col) {
  DATATYPE factor = tableau->values[row][col];
  int i;
  for (i = 0; i < tableau->cols; i++) {
    tableau->values[row][i] = tableau->values[row][i] / factor;
  }
}

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

void
print_tableau (tableau_t * tableau) {
  int i, j;
  for (i = 0; i < tableau->rows; i++) {
    for (j = -1; j < tableau->cols; j++) {
      if (j == -1 && i > 0) printf("x%d\t", tableau->basic[i-1] + 1);
      else printf("\t");
      if (j != -1) printf("%f ", tableau->values[i][j]);
    }
    printf("\n");
  }
}

void
update_basic_vars(tableau_t * tableau, int leaving, int entering) {
  int i;
  for (i = 0; i < tableau->rows - 1; i++) {
    if (tableau->basic[i] == tableau->basic[leaving]) {
      tableau->basic[i] = entering;
    }
  }
}

// Single simplex iteration. If soln is optimum, will return same table
tableau_t *
simplex_iteration (tableau_t * tableau) {
  int entering_var = check_optimum(tableau);
  int leaving_idx = mrt(tableau, entering_var);

  DP("leaving idx: %d var: %d, entering var: %d\n", leaving_idx, tableau->basic[leaving_idx] + 1, entering_var + 1);

  // make the leaving var, entering var element a 1
  normalize_row(tableau, leaving_idx + 1, entering_var);

  // make the columns of the entering var 0
  gaussian_col(tableau, leaving_idx + 1, entering_var);

  update_basic_vars(tableau, leaving_idx, entering_var);

  return tableau;
}

// Runs simplex on the given table and returns the table itself. Must have an
// initial basis setup
tableau_t *
simplex (tableau_t * tableau) {
  while (true) {
    if (check_optimum(tableau) == -1) {
      break;
    }

    tableau = simplex_iteration(tableau);

    #ifdef DEBUGMODE
    print_tableau(tableau);
    char buf[2];
    fgets(buf, 2, stdin);
    #endif
  }

  return tableau;
}

// Checks to see if the index row of the table looks like the identity. If so
// returns what row of the identity. if not returns -1
int
is_ident (tableau_t * tableau, int index) {
  int found_index = -1;
  for (int i = 1; i < tableau->rows; i++) {
    DATATYPE val = tableau->values[i][index];
    if (val == 1) {
      if (found_index != -1) return -1;
      else found_index = i - 1;
    }
    else if (val == 0) continue;
    else return -1;
  }
}

// Returns the necessary # of artificial vars
int
find_basis (tableau_t * tableau) {
  // tableau->basic should be of the right size, but all -1
  // Loop over all of the columns and look for something that looks like the
  // identity.

  int i, num_bas;
  for (i = 0; i < tableau->cols - 1; i++) {

    // If a col looks like the identity, num_bas++ and add it to the correct
    // spot in the basic vars array.

    int index = is_ident(tableau, i);
    if (index != -1) {
      num_bas++;
      tableau->basic[index] = i;
    }
  }

  // If we end up with num_bas == rows - 1, we've found a basis, return TRUE.
  if (num_bas == tableau->rows - 1) return 0;

  // Otherwise we haven't found a basis. Return the # of artifical vars that
  // will now be necessary in order to find a basis for phase-1 of two phase
  // simplex.
  else return tableau->rows - 1 - num_bas;
}

// Checks for a basis. Returns TWOPHASE if two phases must be run and ONEPHASE
// if only one phase must be run (ie, there's a basis already).
enum {TWOPHASE, ONEPHASE}
artificialize (tableau_t * tableau) {
  int num_artificial = find_basis(tableau);

  // If we need no artifical vars, it's a single phase problem.
  if (num_artificial == 0) return ONEPHASE;
  else {
    // realloc tableau->values to be of the right size. It'll need
    // num_artifical more cols.

    // Loop over basic vars to find what cols of the identity we need so
    // we can figure out where the artifical vars need to be added.
  }
}

solution_t *
solve (DATATYPE ** A, DATATYPE * b, DATATYPE * c, int num_vars, int num_constraints) {
  tableau_t * tableau = malloc(sizeof(tableau_t));
  tableau->values = calloc(num_constraints + 1, sizeof(DATATYPE *));
  tableau->basic = calloc(num_constraints, sizeof(DATATYPE));
  tableau->rows = num_constraints + 1;
  tableau->cols = num_vars + 1;

  int i, x;

  for (i = 0; i < num_constraints + 1; i++) {
    tableau->values[i] = calloc(num_vars + 1, sizeof(DATATYPE));
  }

  // assume we have a basis for now
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

  #ifdef DEBUGMODE
  printf("initial tableau:\n");
  print_tableau(tableau);
  #endif

  simplex(tableau);

  printf("optimum\n");
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

  float a1[] = {1, 1};
  float a2[] = {2, 1};
  float * A[] = {a1, a2};
  float c[] = {1, 1};
  float b[] = {6, 8};

  solve(A, b, c, 2, 2);

  return 0;
}
