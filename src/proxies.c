int xargc = 0;
char **xargv = 0;


#include "f2c.h"

// get the date year, month, day of month
int 
getdat_(integer *i, integer *month, integer *day)
{
    return(0);
}

// get the time hour, min, sec, millisecond(?)
int gettim_(shortint *a, shortint  *b, shortint *c, shortint *d)
{
    return(0);
}


integer len_trim__(char *x, ftnlen n)
{
    char *p = x + n-1;
    while(*p  == ' ')
	p--;
    return(p-x + 1);
//    return(strlen(x));
}

int nargs_()
{
    return(xargc);
}


#include <Rdefines.h>

/* Override any other s_stop() routine and just raise an error with R. */
int s_stop(char *msg, ftnlen len)
{
    if(!msg  || !msg[0])
	return(0);

    PROBLEM  "[s_stop] stopping Hydrus: %s", msg
      ERROR;
    return(0);
}

void sig_die(char *msg, int n)
{
    PROBLEM  "[sig_die] stopping Hydrus: %s", msg
      ERROR;
}



SEXP
R_len_trim(SEXP r_str)
{
    const char *str = CHAR(STRING_ELT(r_str, 0));
    int ans = len_trim__(str, strlen(str));
    return(ScalarInteger(ans));
}
