#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include "matrix_type.h"
#include "kmalloc.h"

/* unary operations on matrices */
void matrix_trans(struct _matrix * m);
void matrix_neg(struct _matrix * m);
void matrix_invert(struct _matrix * m);
double matrix_norm(struct _matrix * m);
/* binary operations */
void matrix_add(struct _matrix * sum,
		struct _matrix * a, struct _matrix * b);
void matrix_sub(struct _matrix * diff,
                struct _matrix * a, struct _matrix * b);
void matrix_scalar_mult(struct _matrix * m, double s);
void matrix_mult(struct _matrix * prod,
		 struct _matrix * a, struct _matrix * b);
double matrix_scalar_prod(struct _matrix * a, struct _matrix * b);
/* misc operations */
void matrix_print(struct _matrix * m);
void matrix_zero_up(struct _matrix * m);
void matrix_copy(struct _matrix * a, struct _matrix * b);
double matrix_get_entry(struct _matrix * m, struct _matrix_entry);
void matrix_set_entry(struct _matrix * m, struct _matrix_entry, double val);
void matrix_random(struct _matrix * m, double min, double max);
void matirx_random_pos_def(struct _matrix * a, double min, double max);
struct _matrix * matrix_alloc(enum kmalloc_type type);
void matrix_free(struct _matrix * m);

#endif
