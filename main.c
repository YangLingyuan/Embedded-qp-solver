#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "test.h"

int main(void)
{
	/* important */
	kmalloc_init();
	srand(time(0));

	/* some test for matrix lib */
	inversion_test();
	scalar_prod_test();

	/* prepare quadratic cost function */
	struct _matrix * p = matrix_alloc(NxN);
	struct _matrix * q = matrix_alloc(Nx1);
	double r = 0;
	struct _quadratic_form * qf = quadratic_form_alloc(p, q, r);
	if (!qf) {
		fprintf(stderr, "invalid quadratic "
				"from quadratic_form_alloc in test\n");
		return EXIT_FAILURE;
	}

	/* to hold initial state */
	struct _matrix * x0 = matrix_alloc(Nx1);

	/* test optimizers */
	for (unsigned i = 0; NUM_OPT_TEST_RUNS > i; i++) {
		matirx_random_pos_def(p, P_RAND_ENTRY_MIN, P_RAND_ENTRY_MAX);
		matrix_random(q, Q_RAND_ENTRY_MIN, Q_RAND_ENTRY_MAX);
		matrix_random(x0, X_RAND_ENTRY_MIN, X_RAND_ENTRY_MAX);

#ifdef GRAD_TEST
		test_optimizer(qf, GRAD_ITERATIONS, x0,
			       "gradient_descent_with_line_search",
				gradient_descent_with_line_search);
#endif
#ifdef HESS_TEST
		test_optimizer(qf, HESS_ITERATIONS, x0,
			       "newton_method_with_line_search",
				newton_method_with_line_search);
#endif
#ifdef ADMM_TEST
		test_optimizer(qf, ADMM_ITERATIONS, x0, "admm", admm);
#endif
#ifdef REF_TEST
		test_reference(qf);
#endif

		printf("\v");
	}

	/* clean up */
	matrix_free(x0);
	matrix_free(p);
	matrix_free(q);
	quadratic_form_free(qf);

	return EXIT_SUCCESS;
}
