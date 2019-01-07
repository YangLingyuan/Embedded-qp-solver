#include "klist.h"

unsigned klist_insert(struct _klist * me, struct _klist * at, unsigned x)
{
	if (!at || !me || (KLIST_LEFT != x && KLIST_RIGHT != x))
		return 0;

	struct _klist * tmp = at->np[x];
	me->np[!x] = at;
	at->np[x] = me;
	me->np[x] = tmp;
	if (tmp)
		tmp->np[!x] = me;

	return 1;
}

struct _klist * klist_delete(struct _klist * me)
{
	if (!me)
		return me;

	struct _klist * tmp[2] = {me->np[KLIST_LEFT], me->np[KLIST_RIGHT]};
	for (unsigned x = KLIST_LEFT; KLIST_RIGHT >= x; x++)
		if (tmp[x])
			tmp[x]->np[!x] = tmp[!x];

	me->np[KLIST_LEFT] = me->np[KLIST_RIGHT] = 0;
	return me;
}
