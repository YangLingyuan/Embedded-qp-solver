#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "test.h"

/* Timer used for bringing the system back to EM0. */
RTCDRV_TimerID_t xTimerForWakeUp;

int main(void)
{
	/* Chip revision alignment and errata fixes */
	CHIP_Init();

	/* If first word of user data page is non-zero, */
	/* enable Energy Profiler trace */
	BSP_TraceProfilerSetup();
	/* Initialize RTC timer. */
	RTCDRV_Init();
	RTCDRV_AllocateTimer(&xTimerForWakeUp);

	/* important */
	kmalloc_init();
	srand(time(0));

	/* some test for matrix lib */
#ifdef INV_TEST
	inversion_test();
#endif
#ifdef PROD_TEST
	scalar_prod_test();
#endif

	/* prepare quadratic cost function */
	struct _matrix *p = matrix_alloc(NxN);
	struct _matrix *q = matrix_alloc(Nx1);
	double r = 0;
	struct _quadratic_form *qf = quadratic_form_alloc(p, q, r);
	if (!qf)
	{

		return EXIT_FAILURE;
	}

	/* to hold initial state */
	struct _matrix *x0 = matrix_alloc(Nx1);

	/* 2 seconds of sleep mode before starting */
	RTCDRV_StartTimer(xTimerForWakeUp, rtcdrvTimerTypeOneshot, 2000, 0, 0);
	EMU_EnterEM3(true);

	/* task period 10 seconds, generate interrupt */
	/* periodically to wake the device up from EM3 */
	RTCDRV_StartTimer(xTimerForWakeUp, rtcdrvTimerTypePeriodic, 10000, 0, 0);
	for (unsigned i = 0; NUM_OPT_TEST_RUNS > i; i++)
	{
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

		/* enter sleep mode */
		EMU_EnterEM3(true);
	}

	/* clean up */
	matrix_free(x0);
	matrix_free(p);
	matrix_free(q);
	quadratic_form_free(qf);

	/* do nothing */
	while (1)
		EMU_EnterEM3(true);

	return 0;
}
