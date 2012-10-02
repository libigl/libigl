/* 
 *  simulation of a chronograph
 *  in : tim
 *  out: tim.dtim = elapsed time in micro-secs
 *       tim.ptim = elapsed time in secs
 *       tim.call = number of calls
 *
 *  Written by Pascal J. Frey
 *  email: Pascal.Frey@inria.fr, 1999
*/

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "chrono.h"


/* return elapsed time in secs. */
static double diftim(time_t t2,time_t t1) {
  struct  tm  *ptm;
  double  tim;
  int     hh1,mm1,ss1,hh2,mm2,ss2;

  ptm = localtime(&t1);
  hh1 = ptm->tm_hour;
  mm1 = ptm->tm_min;
  ss1 = ptm->tm_sec;

  ptm = localtime(&t2);
  hh2 = ptm->tm_hour;
  mm2 = ptm->tm_min;
  ss2 = ptm->tm_sec;
  if ( hh2 < hh1 )  hh2 += 24;
  
  tim  = 3600.0*(hh2-hh1);
  tim += 60.0*(mm2-mm1);
  tim += ss2-ss1;

  return(tim);
}


/* get system and user times in micro-seconds */
void  chrono(int cmode,mytime *ptt) {
  time_t tt;

  if ( cmode == RESET ) {
    ptt->dtim  = clock();
    ptt->ctim  = 0.0f;
    ptt->ptim  = 0;
    ptt->call  = 0;
  }
  else {
    ptt->dtim = difftime(clock(),ptt->dtim);  /* in secs */
    if ( cmode == ON ) {
      ptt->ptim = time(NULL);
      ptt->call++;
    }
    else if ( cmode == OFF ) {
      tt = time(NULL);
      ptt->ctim += diftim(tt,ptt->ptim);
      ptt->ptim  = 0;
    }
  }
}


/* return time (converted in secs */
double gttime(mytime t) {

  if ( t.ctim < MAXCLK )
    return(t.dtim / (double)CLOCKS_PER_SEC);
  else
    return(t.ctim);
}


/* initialize time table */
void  tminit(mytime *t,int maxtim) {
  int     k;

  for (k=0; k<maxtim; k++) {
    t[k].dtim = clock();
    t[k].ptim = 0;
    t[k].ctim = 0.0;
    t[k].call = 0;
  }
}


#ifdef __cplusplus
}
#endif

