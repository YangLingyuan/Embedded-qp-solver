#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#include "test.h"

#define TMP_FILENAME "tmp_test_file"
#define SHELL_REF_TEST_CMD PYTHON_COMMAND " test/qp_ref.py " TMP_FILENAME
#define ABS(x) (x < 0 ? -x : x)

static unsigned almost_identity(struct _matrix * m, double precision)
{
	double x = 0;
	unsigned nrows = MATRIX_GET_ROW(m);
	for (unsigned i = 0; nrows > i; i++)
		for (unsigned j = 0; nrows > j; j++) {
			x = matrix_get_entry(m, ME(i, j));
			if (isnan(x)) {
				fprintf(stderr,
					"entry (%u, %u) is nan\n", i, j);
				return 0;
			} else if (i == j && precision < ABS(1 - x)) {
				fprintf(stderr, "entry (%u, %u) is %.12lf "
						"but should be 1\n", i, j, x);
				return 0;
			} else if (i != j && precision < ABS(x)) {
				fprintf(stderr, "entry (%u, %u) is %.12lf "
						"but should be 0\n", i, j, x);
				return 0;
			}
		}
	return 1U;
}

void inversion_test(void)
{
	unsigned i = 0;
	for (; NUM_INVERSION_TEST_RUNS > i; i++) {
		struct _matrix * a = matrix_alloc(NxN);
		struct _matrix * b = matrix_alloc(NxN);
		struct _matrix * c = matrix_alloc(NxN);

		matirx_random_pos_def(a, P_RAND_ENTRY_MIN, P_RAND_ENTRY_MAX);
		matrix_copy(b, a);
		matrix_invert(a);
		matrix_mult(c, b, a);
		if (!almost_identity(c, INVERSION_TEST_PRECISION))
			fprintf(stderr, "error\n");
		matrix_free(a);
		matrix_free(b);
		matrix_free(c);
	}

	printf("%u inversion tests completed\n", i);
}

void scalar_prod_test(void)
{
	struct _matrix * a = matrix_alloc(NxN);
	for (unsigned i = 0; N_DIM > i; i++)
		matrix_set_entry(a, ME(i, 0), i);

	struct _matrix * b = matrix_alloc(Nx1);
	matrix_set_entry(b, ME(0, 0), 1);

	struct _matrix * c = matrix_alloc(Nx1);
	matrix_mult(c, a, b);

	struct _matrix * d = matrix_alloc(Nx1);
	for (unsigned i = 0; N_DIM > i; i++)
		matrix_set_entry(d, ME(i, 0), i);
	matrix_trans(c);
	matrix_trans(d);

	double should_be = (N_DIM*(N_DIM+1)*(2*N_DIM+1)/6) - N_DIM*N_DIM;
	if (should_be != matrix_scalar_prod(d, c))
		fprintf(stderr, "error in mult test\n");

	matrix_free(a);
	matrix_free(b);
	matrix_free(c);
	matrix_free(d);

	printf("scalar_prod test completed\n");
}

void
test_optimizer(struct _quadratic_form * qf, unsigned iterations,
	       struct _matrix * x0, const char * opt_name,
	       struct _matrix * (* opt)(struct _matrix *,
		                        unsigned,
		                        struct _quadratic_form *))
{
	time_t start = time(0);
	struct _matrix * x = opt(x0, iterations, qf);
	time_t end = time(0);

	printf("%s gave:\n", opt_name);
	printf("f(x_min) = %e (took around %li s)\n",
	       quadratic_form_eval(qf, x), end - start);
	fflush(stdout);

	matrix_free(x);
}

void test_reference(struct _quadratic_form * qf)
{
	FILE * fp = fopen(TMP_FILENAME, "w");

	double n = MATRIX_GET_ROW(qf->p);
	fwrite((const void *)&n,
	       sizeof(n), 1, fp);

	struct _matrix * ms[2] = {qf->p, qf->q,};
	for (unsigned i = 0; 2 > i; i++) {
		struct _matrix * m = ms[i];
		double * elements = m->elements;
		unsigned rows = MATRIX_GET_ROW(m);
		unsigned cols = MATRIX_GET_COL(m);
		fwrite((const void *)m->elements,
		       sizeof(*elements), rows * cols, fp);
	}

	fclose(fp);

	system(SHELL_REF_TEST_CMD);
	unlink(TMP_FILENAME);
}
