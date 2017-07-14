
#include <Rdefines.h>


#include "RHydrus.h"

static SEXP R_CallbackFun;

int R_IterationCallback(integer iteration, integer numnp, logical lwat, logical lchem, real *vold, real *vnew, real *thvold, real *thvnew);

SEXP
R_hydrus_standalone(SEXP r_args)
{

    xargc = Rf_length(r_args);
    xargv = (char **) R_alloc(sizeof(char *), xargc);
    for(int i = 0; i < xargc; i++)
	xargv[i] = strdup(CHAR(STRING_ELT(r_args, i)));
    MAIN__();
    return(R_NilValue);
}


SEXP
R_hydrus_inputs(SEXP r_options, SEXP r_select, SEXP r_atmosph, SEXP r_profile, SEXP r_outdir, SEXP r_print, SEXP callback)
{
    char *args[] = {"RdummyName"};
    xargc = 1;
    xargv = args;

    if(Rf_length(callback))  {
	if(TYPEOF(callback) == EXTPTRSXP) {
	    RCallbackRoutine = R_ExternalPtrAddr(callback);
	} else {
	    R_CallbackFun = callback;
	    RCallbackRoutine = R_IterationCallback;
	}
    }
    
    R_MAIN__(CHAR(STRING_ELT(r_options, 0)), CHAR(STRING_ELT(r_select, 0)), 
             CHAR(STRING_ELT(r_atmosph, 0)), CHAR(STRING_ELT(r_profile, 0)), 
	     CHAR(STRING_ELT(r_outdir, 0)), INTEGER(r_print)[0]);

    R_CallbackFun = NULL;
    RCallbackRoutine = NULL;

    return(R_NilValue);
}

SEXP
R_hydrus_direct(SEXP r_options, SEXP r_select, SEXP r_atmosph, SEXP r_profile, SEXP r_outdir, SEXP r_print, SEXP r_profile_data)
{

    char *args[] = {"RdummyName"};
    xargc = 1;
    xargv = args;

    R_direct_MAIN__(CHAR(STRING_ELT(r_options, 0)), CHAR(STRING_ELT(r_select, 0)), 
		    CHAR(STRING_ELT(r_atmosph, 0)), CHAR(STRING_ELT(r_profile, 0)), 
                    CHAR(STRING_ELT(r_outdir, 0)), INTEGER(r_print)[0], r_profile_data);

    return(R_NilValue);
}



int
R_TestCallback(integer iteration, integer numnp, logical lwat, logical lchem, real *vold, real *vnew, real *thvold, real *thvnew)
{
    static int ctr = 0;
    if(iteration == 1)
	ctr = 0;

    ctr++;
    return( ctr == 10 ? -1 : 0);
}

int
R_IterationCallback(integer iteration, integer numnp, logical lwat, logical lchem, real *vold, real *vnew, real *thvold, real *thvnew)
{
    if(R_CallbackFun) {
	SEXP e, p;
	double *dold, *dnew;
	PROTECT(e = allocVector(LANGSXP, 6));
	p = e;
	SETCAR(p, R_CallbackFun); p= CDR(p);
	SETCAR(p, ScalarInteger(iteration)); p= CDR(p);
	SETCAR(p, ScalarInteger(numnp)); p= CDR(p);
	SETCAR(p, ScalarLogical(lwat)); p= CDR(p);
	SETCAR(p, dold = NEW_NUMERIC(iteration)); p= CDR(p);
	SETCAR(p, dnew = NEW_NUMERIC(iteration)); 
	for(int i = 0; i < iteration; i++) {
	    dold[i] = vold[i];
	    dnew[i] = vnew[i];
	}

	SEXP r_ans = Rf_eval(e, R_GlobalEnv);
	int ans = 0;
	if(TYPEOF(r_ans) == INTSXP)
	    ans = INTEGER(r_ans)[0];

	UNPROTECT(1);
	return(ans);
    } else
	return(0);
}


