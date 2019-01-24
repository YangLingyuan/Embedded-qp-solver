#ifndef QP_SOLVERS
#define QP_SOLVERS

struct _matrix *
gradient_descent_with_line_search(struct _matrix * x0,
		                  unsigned iterations,
		                  struct _quadratic_form * qf);

struct _matrix *
newton_method_with_line_search(struct _matrix * x0,
		               unsigned iterations,
		               struct _quadratic_form * qf);

struct _matrix *
admm(struct _matrix * x0, unsigned iterations,
		          struct _quadratic_form * qf);

#endif
