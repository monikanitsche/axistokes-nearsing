#define THIRD
c#define FOURTH

#define NOROUNDOFFFIX
c#define ROUNDOFFFIX

#define DBLE
c#define QUAD

#ifdef DBLE
#define MODE real*8
#endif

#ifdef QUAD
#define MODE real*16
#endif
