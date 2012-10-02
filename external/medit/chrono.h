#ifdef __cplusplus
extern "C" {
#endif

#include <time.h>

#ifndef  ON
#define  RESET  0
#define  ON     1
#define  OFF    2
#endif

#define  TIMEMAX   16
#define  MAXCLK    ( 1073741823. / (double)CLOCKS_PER_SEC )


typedef struct mytime {
  double    ctim,dtim;
  time_t    ptim;
  short     call;
} mytime;


/* prototypes */
void   chrono(int cmode,mytime *ptt);
double gttime(mytime t);
void   tminit(mytime *t,int maxtim);


#ifdef __cplusplus
}
#endif
