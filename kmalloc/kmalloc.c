#include <string.h>

#include "kmalloc.h"
#include "klist.h"

#include "matrix_type.h"
#include "qp.h"

/* to be appended to any type */
struct _km_header {
	struct _klist link;
};

/* supported types */
struct _km_NxN {
	struct _km_header h;
	struct _matrix val;
	double elements[N * N];
};
struct _km_Nx1 {
	struct _km_header h;
	struct _matrix val;
	double elements[N * 1];
};
struct _km_QUADRATIC_FORM {
	struct _km_header h;
	struct _quadratic_form val;
};

unsigned km_size[KMALLOC_TYPE_END];

#define KMALLOC_VAL(link) \
        ((void *)(((char *)link) + (sizeof(struct _km_header))))
#define KMALLOC_HEADER(val) \
        ((struct _km_header *)(((char *)val) - (sizeof(struct _km_header))))

static struct _km_NxN NxN_bucket[NxN_MAX];
static struct _km_Nx1 Nx1_bucket[Nx1_MAX];
static struct _km_QUADRATIC_FORM QUADRATIC_FORM_bucket[QUADRATIC_FORM_MAX];

static struct _klist * km_used[KMALLOC_TYPE_END];
static struct _klist * km_unused[KMALLOC_TYPE_END];

static void km_insert(struct _klist * me, struct _klist ** list);
static void km_remove(struct _klist * me, struct _klist ** list);
static struct _klist * km_pop(struct _klist ** list);

void kmalloc_init(void)
{
	/* ugly but meh.. */
	km_size[NxN] = sizeof(struct _km_NxN);
	for (unsigned i = 0; NxN_MAX > i; i++) {
		struct _km_NxN * bucket_element = &NxN_bucket[i];
		km_insert(&bucket_element->h.link, &km_unused[NxN]);
		bucket_element->val.elements = bucket_element->elements;
	}
	km_size[Nx1] = sizeof(struct _km_Nx1);
	for (unsigned i = 0; Nx1_MAX > i; i++) {
		struct _km_Nx1 * bucket_element = &Nx1_bucket[i];
		km_insert(&bucket_element->h.link, &km_unused[Nx1]);
		bucket_element->val.elements = bucket_element->elements;
	}
	km_size[QUADRATIC_FORM] = sizeof(struct _km_QUADRATIC_FORM);
	for (unsigned i = 0; QUADRATIC_FORM_MAX > i; i++)
		km_insert(&QUADRATIC_FORM_bucket[i].h.link,
				        &km_unused[QUADRATIC_FORM]);
}

static void km_insert(struct _klist * me, struct _klist ** list)
{
	if (!me || !list)
		return;

	me->np[KLIST_RIGHT] = *list;
	me->np[KLIST_LEFT] = 0;
	*list = me;
}

static void km_remove(struct _klist * me, struct _klist ** list)
{
	if (!me || !list)
		return;

	if (me == *list)
		*list = me->np[KLIST_RIGHT];
	klist_delete(me);
}

static struct _klist * km_pop(struct _klist ** list)
{
	if (!list)
		return 0;

	struct _klist * tmp = *list;
	if (tmp)
		*list = tmp->np[KLIST_RIGHT];
	klist_delete(tmp);

	return tmp;
}

void * kmalloc(enum kmalloc_type type, unsigned flags)
{
	if (KMALLOC_TYPE_END <= type)
		/* unkown type */
		return 0;

	struct _klist * tmp = km_pop(&km_unused[type]);
	if (!tmp)
		/* no elements of this type left */
		return tmp;

	if (KM_ZERO & flags)
		memset((void *)tmp, 0, km_size[type]);

	km_insert(tmp, &km_used[type]);

	return KMALLOC_VAL(tmp);
}

void kfree(void * me, enum kmalloc_type type)
{
	if (!me || KMALLOC_TYPE_END <= type)
		/* nothing to do */
		return;

	struct _klist * tmp = &KMALLOC_HEADER(me)->link;
	km_remove(tmp, &km_used[type]);
	km_insert(tmp, &km_unused[type]);
}
