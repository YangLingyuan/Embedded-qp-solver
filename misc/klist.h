#ifndef KLIST_H
#define KLIST_H

#define KLIST_LEFT 0
#define KLIST_RIGHT 1

#define KLIST_VAL(type, link) \
        ((type *)(((char *)link) - (sizeof(type) - sizeof(struct _klist))))

/* has to be embedded into some other structure, only provides the */
/* basic list mechanism nothing more (not even head tail) */
struct _klist {
	struct _klist * np[2];
};

/* inserts a klist node at the left or right of the klist node at */
/* x `elem` [KLIST_LEFT, KLIST_RIGHT] */
unsigned klist_insert(struct _klist * me, struct _klist * at, unsigned x);

/* deletes an arbitrary klist node form an existing list */
struct _klist * klist_delete(struct _klist * me);

#endif
