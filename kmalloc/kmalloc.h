#ifndef KMALLOC_H
#define KMALLOC_H

#include "bits.h"

#define KM_ZERO_BIT 31
#define KM_ZERO BIT32(KM_ZERO_BIT)

enum kmalloc_type {
	NxN,
	Nx1,
	QUADRATIC_FORM,
	KMALLOC_TYPE_END,
};

void * kmalloc(enum kmalloc_type type, unsigned flags);
void kfree(void * me, enum kmalloc_type type);
void kmalloc_init(void);

#endif
