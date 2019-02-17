#ifndef CONFIG_H
#define CONFIG_H

/* dimension of matrices */
#define N 48U
/* number of NxN matrices at available at one time */
#define NxN_MAX 5
/* number of Nx1 matrices at available at one time */
#define Nx1_MAX 10

/* how many inversion test to run */
#define INVERSION_TEST_PRECISION 1e-6
#define NUM_INVERSION_TEST_RUNS 8

/* bounds on entries of p matrix in quadratic forms */
#define P_RAND_ENTRY_MIN -1e3
#define P_RAND_ENTRY_MAX 1e3
/* bounds on entries of q matrix in quadratic forms */
#define Q_RAND_ENTRY_MIN -1e3
#define Q_RAND_ENTRY_MAX 1e3
/* bounds on entries of initial state */
#define X_RAND_ENTRY_MIN -1e3
#define X_RAND_ENTRY_MAX 1e3

/* number of optimization test runs */
#define NUM_OPT_TEST_RUNS 16

/* constraints for admm */
#define ADMM_BOX_CONSTRAINT_MAX 1e12
#define ADMM_BOX_CONSTRAINT_MIN -1e12

/* python command on ur system */
#define PYTHON_COMMAND "python"

/* specify iterations each opt algo will take in worst case */
#define GRAD_ITERATIONS 1e4
#define HESS_ITERATIONS 1e1
#define ADMM_ITERATIONS 1e4

/* define or do not define to run test or not */
#define INV_TEST
#define PROD_TEST
#define GRAD_TEST
#define HESS_TEST
#define ADMM_TEST
#define REF_TEST

#endif
