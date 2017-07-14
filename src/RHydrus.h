#include <f2c.h>
#define MAIN__  hydrus_MAIN__

extern int MAIN__(void);

extern int xargc;
extern char **xargv;


#include <Rdefines.h>

int R_direct_MAIN__(const char *options_in, const char *selector_in,  const char *atmosph_in,
		    const char *profile_dat, const char *outdir, int printOut, SEXP r_profile);


int (*RCallbackRoutine)(integer , integer, logical, logical, real *, real *, real *, real *, integer *tmax);
