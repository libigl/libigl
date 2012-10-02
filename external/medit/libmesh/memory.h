#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <assert.h>

/* prototype (re)definitions */
void  *M_malloc(size_t size,char *call);
void  *M_calloc(size_t nelem,size_t elsize,char *call);
void  *M_realloc(void *ptr, size_t size,char *call);
void   M_free(void *ptr);

/* ptototypes : tools */
int    M_memLeak();
void   M_memDump();
size_t M_memSize();


#ifdef __cplusplus
}
#endif
