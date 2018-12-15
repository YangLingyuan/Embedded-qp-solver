#ifndef MATRIX_H
#define MATRIX_H

struct _matrix {
	unsigned nrows;
	unsigned ncols;
	double * elements;
	unsigned * permutation;
};

/* unary operations on matrices */
void matrix_trans(struct _matrix m);
void matrix_neg(struct _matrix m);
void matrix_invert(struct _matrix m,
		   struct _matrix tmp, double * a, double * b);
void matrix_lup_decompose(struct _matrix m);
/* binary operations */
void matrix_add(struct _matrix sum,
		struct _matrix a, struct _matrix b);
void matrix_sub(struct _matrix diff,
                struct _matrix a, struct _matrix b);
void matrix_mult(struct _matrix prod,
		 struct _matrix a, struct _matrix b);
/* misc operations */
void matrix_print(struct _matrix m);
void matrix_zero_up(struct _matrix m);
void matrix_copy(struct _matrix a, struct _matrix b);
void matrix_set_up(struct _matrix * m, unsigned nrows,
		unsigned ncols, double * elements, unsigned * permutation);
double matrix_get_entry(struct _matrix m, unsigned row, unsigned col);
void matrix_set_entry(struct _matrix m,
                      unsigned row, unsigned col, double val);


#endif /* MATRIX_H */
