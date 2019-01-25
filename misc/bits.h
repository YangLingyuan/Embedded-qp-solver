#ifndef BITS_H
#define BITS_H 1

#define BIT32(n) (1U   << (n))
#define BIT64(n) (1ULL << (n))

#define MASK32(n, m) (((~0U)   << (n)) & (~0U   >> (32 - 1 - (m))))
#define MASK64(n, m) (((~0ULL) << (n)) & (~0ULL >> (64 - 1 - (m))))

#define SETB32(n, word) ((word) |= BIT32(n))
#define SETB64(n, word) ((word) |= BIT64(n))
#define SETM32(n, m, word, what)                             \
	do {                                                 \
		CLRM32(n, m, word);                          \
		((word) |= (MASK32(n, m) & ((what) << (n))));\
	} while (0)
#define SETM64(n, m, word, what)                             \
	do {                                                 \
		CLRM64(n, m, word);                          \
		((word) |= (MASK64(n, m) & ((what) << (n))));\
	} while (0)
#define CLRB32(n, word) ((word) &= ~BIT32(n))
#define CLRB64(n, word) ((word) &= ~BIT64(n))
#define CLRM32(n, m, word) ((word) &= (~MASK32(n, m)))
#define CLRM64(n, m, word) ((word) &= (~MASK64(n, m)))

#define GETB32(n, word) ((word) & BIT32(n))
#define GETB64(n, word) ((word) & BIT64(n))
#define GETM32(n, m, word) ((word) & MASK32(n, m))
#define GETM64(n, m, word) ((word) & MASK64(n, m))

/* for 32 bit */
#define ALIGN8(x) (7U & (x) ? ((x) + 8U) & ~7U : (x))

#endif
