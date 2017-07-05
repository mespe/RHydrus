int xargc = 0;
char **xargv = 0;


#include "f2c.h"

int Size = sizeof(ftnlen);
const ftnlen val = 260;
int vali = val;
int 
getdat_(integer *i, integer *month, integer *day)
{
    return(0);
}

int gettim_(shortint *a, shortint  *b, shortint *c, shortint *d)
{
    return(0);
}


integer len_trim__(char *x, ftnlen n)
{
    return(strlen(x));
}

int nargs_()
{
    return(xargc);
}


#include <Rdefines.h>

/* Override any other s_stop() routine and just raise an error with R. */
int s_stop(char *msg, ftnlen len)
{
    PROBLEM  "[s_stop] stopping Hydrus: %s", msg
      ERROR;
    return(0);
}

void sig_die(char *msg, int n)
{
    PROBLEM  "[sig_die] stopping Hydrus: %s", msg
      ERROR;
}
