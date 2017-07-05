
#include <Rdefines.h>


#include "RHydrus.h"



SEXP
R_hydrus_standalone(SEXP r_args)
{

    xargc = Rf_length(r_args);
    xargv = (char **) R_alloc(sizeof(char *), xargc);
    for(int i = 0; i < xargc; i++)
	xargv[i] = CHAR(STRING_ELT(r_args, i));
    MAIN__();
    return(R_NilValue);
}
