/* INPUT.f -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h> /* For exit() */
#include <f2c.h>

/* Table of constant values */

static integer c__30 = 30;
static integer c__1 = 1;
static integer c__8 = 8;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__9 = 9;
static integer c__32 = 32;
static integer c__0 = 0;
static real c_b678 = 1e-8f;
static real c_b682 = 1.f;
static integer c__2 = 2;
static real c_b699 = -1.f;
static doublereal c_b710 = 10.;
static integer c__5 = 5;
static integer c__31 = 31;
static integer c__33 = 33;

/* Source file INPUT.FOR |||||||||||||||||||||||||||||||||||||||||||||||| */
/* Subroutine */ int basinf_(real *cosalf, integer *maxit, real *tolth, real *
	tolh, logical *topinf, logical *botinf, logical *shorto, logical *
	lwat, logical *lchem, logical *sinkf, logical *wlayer, logical *qgwlf,
	 logical *freed, logical *seepf, logical *atmbc, integer *kodtop, 
	integer *kodbot, real *rtop, real *rroot, real *rbot, real *hcrits, 
	real *hcrita, real *gwl0l, real *aqh, real *bqh, integer *ktold, 
	integer *kbold, integer *nunitd, integer *iunit, integer *nmat, 
	integer *nmatd, integer *nlay, logical *lroot, logical *ltemp, 
	logical *lwdep, logical *lequil, logical *lscreen, logical *qdrain, 
	real *zbotdr, real *basegw, real *rspacing, integer *iposdr, real *
	rkhtop, real *rkhbot, real *rkvtop, real *rkvbot, real *entres, real *
	wetper, real *zintf, real *geofac, logical *linitw, logical *lvarbc, 
	real *xconv, real *tconv, logical *lmeteo, logical *lvapor, integer *
	iver, logical *lprint, logical *lcentrif, logical *lsnow, real *hseep,
	 logical *lflux, logical *lactrsu, integer *ierr)
{
    /* Format strings */
    static char fmt_100[] = "(\002 Date: \002,i3,\002.\002,i2,\002.\002,\002"
	    "    Time: \002,i3,\002:\002,i2,\002:\002,i2)";
    static char fmt_110[] = "(f6.3,i5,f8.3,f8.5)";
    static char fmt_120[] = "(13l6)";

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern integer igetfileversion_(integer *, integer *);
    static real g;
    static integer i__, ii;
    extern /* Subroutine */ int conversion_(char *, char *, real *, real *, 
	    ftnlen, ftnlen);
    static char hed[72];
    static integer iday, mins, isecs;
    static char lunit[5], munit[5], tunit[5];
    extern /* Subroutine */ int getdat_(integer *, integer *, integer *), 
	    gettim_(integer *, integer *, integer *, integer *);
    static integer imonth;
    static logical ldummy;
    static integer ihours;

    /* Fortran I/O blocks */
    static cilist io___1 = { 1, 30, 0, 0, 0 };
    static cilist io___2 = { 1, 30, 0, 0, 0 };
    static cilist io___3 = { 1, 30, 0, "(a)", 0 };
    static cilist io___5 = { 1, 30, 0, 0, 0 };
    static cilist io___6 = { 1, 30, 0, "(a)", 0 };
    static cilist io___8 = { 1, 30, 0, "(a)", 0 };
    static cilist io___10 = { 1, 30, 0, "(a)", 0 };
    static cilist io___12 = { 1, 30, 0, 0, 0 };
    static cilist io___13 = { 1, 30, 0, 0, 0 };
    static cilist io___14 = { 1, 30, 0, 0, 0 };
    static cilist io___15 = { 1, 30, 0, 0, 0 };
    static cilist io___17 = { 1, 30, 0, 0, 0 };
    static cilist io___18 = { 1, 30, 0, 0, 0 };
    static cilist io___19 = { 1, 30, 0, 0, 0 };
    static cilist io___20 = { 1, 30, 0, 0, 0 };
    static cilist io___21 = { 1, 30, 0, 0, 0 };
    static cilist io___22 = { 1, 30, 0, 0, 0 };
    static cilist io___23 = { 1, 30, 0, 0, 0 };
    static cilist io___24 = { 1, 30, 0, 0, 0 };
    static cilist io___25 = { 1, 30, 0, 0, 0 };
    static cilist io___26 = { 1, 30, 0, 0, 0 };
    static cilist io___28 = { 1, 30, 0, 0, 0 };
    static cilist io___29 = { 1, 30, 0, 0, 0 };
    static cilist io___30 = { 1, 30, 0, 0, 0 };
    static cilist io___31 = { 1, 30, 0, 0, 0 };
    static cilist io___32 = { 1, 30, 0, 0, 0 };
    static cilist io___33 = { 1, 30, 0, 0, 0 };
    static cilist io___34 = { 1, 30, 0, 0, 0 };
    static cilist io___35 = { 1, 30, 0, 0, 0 };
    static cilist io___36 = { 1, 30, 0, 0, 0 };
    static cilist io___37 = { 1, 30, 0, 0, 0 };
    static cilist io___38 = { 1, 30, 0, 0, 0 };
    static cilist io___39 = { 1, 30, 0, 0, 0 };
    static cilist io___40 = { 1, 30, 0, 0, 0 };
    static cilist io___41 = { 1, 30, 0, 0, 0 };
    static cilist io___42 = { 1, 30, 0, 0, 0 };
    static cilist io___43 = { 1, 30, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, 0, 0 };
    static cilist io___57 = { 0, 6, 0, 0, 0 };
    static cilist io___58 = { 0, 6, 0, 0, 0 };
    static cilist io___59 = { 0, 6, 0, 0, 0 };
    static cilist io___60 = { 0, 6, 0, 0, 0 };
    static cilist io___62 = { 1, 0, 0, 0, 0 };
    static cilist io___63 = { 1, 0, 0, 0, 0 };
    static cilist io___69 = { 1, 0, 0, fmt_100, 0 };
    static cilist io___70 = { 1, 0, 0, 0, 0 };
    static cilist io___72 = { 1, 50, 0, 0, 0 };
    static cilist io___73 = { 1, 50, 0, 0, 0 };
    static cilist io___74 = { 1, 50, 0, fmt_110, 0 };
    static cilist io___75 = { 1, 50, 0, 0, 0 };
    static cilist io___76 = { 1, 50, 0, 0, 0 };
    static cilist io___77 = { 1, 50, 0, fmt_120, 0 };


    /* Parameter adjustments */
    --iunit;

    /* Function Body */
    *iver = igetfileversion_(&c__30, &c__1);
    i__1 = s_rsle(&io___1);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___2);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsfe(&io___3);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_fio(&c__1, hed, (ftnlen)72);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___5);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsfe(&io___6);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_fio(&c__1, lunit, (ftnlen)5);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsfe(&io___8);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_fio(&c__1, tunit, (ftnlen)5);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsfe(&io___10);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_fio(&c__1, munit, (ftnlen)5);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L901;
    }
    conversion_(lunit, tunit, xconv, tconv, (ftnlen)5, (ftnlen)5);
    i__1 = s_rsle(&io___12);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___13);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*lwat), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*lchem), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*ltemp), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*sinkf), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*lroot), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*shorto), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*lwdep), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*lscreen), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*atmbc), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*lequil), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*iver == 3) {
	i__1 = s_rsle(&io___14);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___15);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lsnow), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&ldummy, (ftnlen)sizeof(logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&ldummy, (ftnlen)sizeof(logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&ldummy, (ftnlen)sizeof(logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    } else if (*iver == 4) {
	i__1 = s_rsle(&io___17);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___18);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lsnow), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&ldummy, (ftnlen)sizeof(logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lmeteo), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lvapor), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lactrsu), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lflux), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    if (*lsnow && ! (*ltemp)) {
	*lsnow = FALSE_;
    }
    i__1 = s_rsle(&io___19);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___20);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*nmat), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*nlay), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*cosalf), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*nmat > *nmatd || *nlay > 10) {
	*ierr = 4;
	return 0;
    }
    i__1 = s_rsle(&io___21);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_rsle(&io___22);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_rsle(&io___23);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*maxit), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*tolth), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*tolh), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_rsle(&io___24);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_rsle(&io___25);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*topinf), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*wlayer), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*kodtop), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*linitw), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L902;
    }
    if (*kodtop == 0) {
	*lvarbc = TRUE_;
    }
    i__1 = s_rsle(&io___26);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L902;
    }
    if (*iver <= 3) {
	ii = 1;
	i__1 = s_rsle(&io___28);
	if (i__1 != 0) {
	    goto L904;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*botinf), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L904;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*qgwlf), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L904;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*freed), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L904;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*seepf), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L904;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*kodbot), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L904;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*qdrain), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L904;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L904;
	}
	ii = 0;
L904:
	if (ii == 1) {
	    *qdrain = FALSE_;
	}
	*hseep = 0.f;
    } else {
	i__1 = s_rsle(&io___29);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*botinf), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*qgwlf), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*freed), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*seepf), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*kodbot), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*qdrain), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*hseep), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
    }
    if (! (*topinf) && *kodtop == -1 || ! (*botinf) && *kodbot == -1 && ! (*
	    qgwlf) && ! (*freed) && ! (*seepf) && ! (*qdrain)) {
	i__1 = s_rsle(&io___30);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = s_rsle(&io___31);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*rtop), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*rbot), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*rroot), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
/* ,hCritS,hCritA */
    } else {
	*rtop = 0.f;
	*rbot = 0.f;
	*rroot = 0.f;
    }
    if (*qgwlf) {
	i__1 = s_rsle(&io___32);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = s_rsle(&io___33);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*gwl0l), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*aqh), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*bqh), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
    }
    if (*qdrain) {
	i__1 = s_rsle(&io___34);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = s_rsle(&io___35);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*iposdr), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = s_rsle(&io___36);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = s_rsle(&io___37);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*zbotdr), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*rspacing), (ftnlen)sizeof(real)
		);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*entres), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
	*zbotdr = -dabs(*zbotdr);
	i__1 = s_rsle(&io___38);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L902;
	}
	if (*iposdr == 1) {
	    i__1 = s_rsle(&io___39);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkhtop), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L902;
	    }
	} else if (*iposdr == 2) {
	    i__1 = s_rsle(&io___40);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*basegw), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkhtop), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*wetper), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L902;
	    }
	} else if (*iposdr == 3) {
	    i__1 = s_rsle(&io___41);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*basegw), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkhtop), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkhbot), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*wetper), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L902;
	    }
	} else if (*iposdr == 4) {
	    i__1 = s_rsle(&io___42);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*basegw), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkvtop), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkvbot), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkhbot), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*wetper), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*zintf), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L902;
	    }
	} else if (*iposdr == 5) {
	    i__1 = s_rsle(&io___43);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*basegw), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkhtop), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkvtop), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkhbot), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*wetper), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*zintf), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*geofac), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L902;
	    }
	}
	*basegw = -dabs(*basegw);
	*zintf = -dabs(*zintf);
    }
/*     Input modifications */
    *rroot = dabs(*rroot);
    *hcrita = 1e10f;
    *hcrita = -dabs(*hcrita);
    if (*topinf) {
	*kodtop = i_sign(&c__3, kodtop);
    }
    if (*botinf) {
	*kodbot = i_sign(&c__3, kodbot);
    }
    if (*atmbc && *kodtop < 0) {
	*hcrits = 0.f;
	*kodtop = -4;
    }
    if (*wlayer) {
	*kodtop = -abs(*kodtop);
    }
    if (*qgwlf) {
	*kodbot = -7;
    }
    if (*freed) {
	*kodbot = -5;
    }
    if (*seepf) {
	*kodbot = -2;
    }
    *ktold = *kodtop;
    *kbold = *kodbot;
    if (*lscreen) {
	s_wsle(&io___44);
	do_lio(&c__9, &c__1, "----------------------------------------------"
		"------", (ftnlen)52);
	e_wsle();
	s_wsle(&io___45);
	do_lio(&c__9, &c__1, "|                                             "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___46);
	do_lio(&c__9, &c__1, "|                    HYDRUS                   "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___47);
	do_lio(&c__9, &c__1, "|                                             "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___48);
	do_lio(&c__9, &c__1, "|   Code for simulating one-dimensional variab"
		"ly   |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___49);
	do_lio(&c__9, &c__1, "|    saturated water flow, heat transport, and"
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___50);
	do_lio(&c__9, &c__1, "|   transport of solutes involved in sequentia"
		"l    |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___51);
	do_lio(&c__9, &c__1, "|         first-order decay reactions         "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___52);
	do_lio(&c__9, &c__1, "|                                             "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___53);
	do_lio(&c__9, &c__1, "|                  version 4.08               "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___54);
	do_lio(&c__9, &c__1, "|                                             "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___55);
	do_lio(&c__9, &c__1, "|           Last modified: January, 2009      "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___56);
	do_lio(&c__9, &c__1, "|                                             "
		"     |", (ftnlen)52);
	e_wsle();
	s_wsle(&io___57);
	do_lio(&c__9, &c__1, "----------------------------------------------"
		"------", (ftnlen)52);
	e_wsle();
	s_wsle(&io___58);
	e_wsle();
	s_wsle(&io___59);
	do_lio(&c__9, &c__1, hed, (ftnlen)72);
	e_wsle();
	s_wsle(&io___60);
	e_wsle();
/*        write(*,*) 'Press Enter to continue' */
/*        read(*,*) */
    }
    ii = 1;
    if (*lprint) {
	ii = *nunitd;
    }
    i__1 = ii;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___62.ciunit = iunit[i__];
	i__2 = s_wsle(&io___62);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, "******* Program HYDRUS", (ftnlen)22);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = e_wsle();
	if (i__2 != 0) {
	    goto L903;
	}
	io___63.ciunit = iunit[i__];
	i__2 = s_wsle(&io___63);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, "******* ", (ftnlen)8);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, hed, (ftnlen)72);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = e_wsle();
	if (i__2 != 0) {
	    goto L903;
	}
	getdat_(&ii, &imonth, &iday);
	gettim_(&ihours, &mins, &isecs, &ii);
	io___69.ciunit = iunit[i__];
	i__2 = s_wsfe(&io___69);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_fio(&c__1, (char *)&iday, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_fio(&c__1, (char *)&imonth, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_fio(&c__1, (char *)&ihours, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_fio(&c__1, (char *)&mins, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_fio(&c__1, (char *)&isecs, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = e_wsfe();
	if (i__2 != 0) {
	    goto L903;
	}
	io___70.ciunit = iunit[i__];
	i__2 = s_wsle(&io___70);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, "Units: L = ", (ftnlen)11);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, lunit, (ftnlen)5);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, ", T = ", (ftnlen)6);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, tunit, (ftnlen)5);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, ", M = ", (ftnlen)6);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = do_lio(&c__9, &c__1, munit, (ftnlen)5);
	if (i__2 != 0) {
	    goto L903;
	}
	i__2 = e_wsle();
	if (i__2 != 0) {
	    goto L903;
	}
/* L11: */
    }
    if (*cosalf <= 1.f) {
	*lcentrif = FALSE_;
    }
    if (*lcentrif) {
	g = *xconv * 9.80665f / *tconv / *tconv;
	*cosalf = *cosalf * *cosalf / g;
    }
    i__1 = s_wsle(&io___72);
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = e_wsle();
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = s_wsle(&io___73);
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_lio(&c__9, &c__1, "CosAlf,MaxIt,TolTh,  TolH", (ftnlen)25);
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = e_wsle();
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = s_wsfe(&io___74);
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*cosalf), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*maxit), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*tolth), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*tolh), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = s_wsle(&io___75);
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = e_wsle();
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = s_wsle(&io___76);
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_lio(&c__9, &c__1, "TopInF,BotInF,AtmBC,SinkF,WLayer,qGWLF,Free"
	    "D,SeepF,lWat,lChem,lTemp,lRoot,lWDep", (ftnlen)79);
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = e_wsle();
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = s_wsfe(&io___77);
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*topinf), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*botinf), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*atmbc), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*sinkf), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*wlayer), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*qgwlf), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*freed), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*seepf), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*lwat), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*lchem), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*ltemp), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*lroot), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = do_fio(&c__1, (char *)&(*lwdep), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L903;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L903;
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
L902:
    *ierr = 2;
    return 0;
/*     Error when writing into an output file */
L903:
    *ierr = 3;
    return 0;
} /* basinf_ */

/* *********************************************************************** */
/* Subroutine */ int conversion_(char *lunit, char *tunit, real *xconv, real *
	tconv, ftnlen lunit_len, ftnlen tunit_len)
{
/*     conversions from m and s to Hydrus units */
    *xconv = 1.f;
    *tconv = 1.f;
    if (s_cmp(lunit, "cm  ", (ftnlen)5, (ftnlen)4) == 0) {
	*xconv = 100.f;
    } else if (s_cmp(lunit, "mm  ", (ftnlen)5, (ftnlen)4) == 0) {
	*xconv = 1e3f;
    }
    if (s_cmp(tunit, "min ", (ftnlen)5, (ftnlen)4) == 0) {
	*tconv = .016666666666666666f;
    } else if (s_cmp(tunit, "hours", (ftnlen)5, (ftnlen)5) == 0) {
	*tconv = 2.7777777777777778e-4f;
    } else if (s_cmp(tunit, "days", (ftnlen)5, (ftnlen)4) == 0) {
	*tconv = 1.1574074074074073e-5f;
    } else if (s_cmp(tunit, "years", (ftnlen)5, (ftnlen)5) == 0) {
	*tconv = 3.1709791983764586e-8f;
    }
    return 0;
} /* conversion_ */

/* *********************************************************************** */
/* Subroutine */ int nodinf_(integer *numnpd, integer *numnp, integer *nobsd, 
	integer *nobs, real *htop, real *hbot, real *x, real *hnew, real *
	hold, integer *matnum, real *htemp, integer *laynum, real *beta, real 
	*ah, real *ak, real *ath, real *conc, real *sorb, real *tempn, real *
	tempo, integer *node, integer *nsd, integer *ns, real *xsurf, logical 
	*lchem, logical *ltemp, logical *lequil, logical *lscreen, logical *
	lbact, real *sorb2, integer *ierr, logical *lprint, logical *lflux, 
	logical *ldualneq)
{
    /* Format strings */
    static char fmt_110[] = "(/\002Nodal point information\002//\002node    "
	    "  x         hOld    MatN LayN  Beta      Ah       AK \002,\002  "
	    "   ATh     Temp    Conc(1...NS)         Sorb(1...NS)\002/)";
    static char fmt_120[] = "(i4,2f11.3,2i5,f8.3,3f9.3,f8.2,10e12.4,10e12.4)";
    static char fmt_130[] = "(/\002 Number of species in the chain : \002,i3)"
	    ;
    static char fmt_140[] = "(///16x,10(15x,a5,i3,\002)\002,7x))";
    static char fmt_200[] = "(/\002         time     \002,10(a29,2x))";
    static char fmt_150[] = "(///16x,10(15x,a5,i3,\002)\002,18x))";
    static char fmt_160[] = "(///16x,10(15x,a5,i3,\002)\002,29x))";
    static char fmt_170[] = "(///16x,10(15x,a5,i3,\002)\002,40x))";
    static char fmt_180[] = "(///16x,10(15x,a5,i3,\002)\002,51x))";
    static char fmt_190[] = "(///16x,10(15x,a5,i3,\002)\002,62x))";
    static char fmt_260[] = "(///16x,10(15x,a5,i3,\002)\002,73x))";
    static char fmt_261[] = "(///14x,10(15x,a5,i3,\002)\002,84x))";
    static char fmt_262[] = "(///14x,10(15x,a5,i3,\002)\002,95x))";
    static char fmt_263[] = "(///14x,10(15x,a5,i3,\002)\002,106x))";
    static char fmt_264[] = "(///14x,10(15x,a5,i3,\002)\002,117x))";
    static char fmt_265[] = "(///14x,10(15x,a5,i3,\002)\002,128x))";
    static char fmt_266[] = "(///14x,10(15x,a5,i3,\002)\002,139x))";
    static char fmt_210[] = "(/\002         time     \002,10(a29,a12,2x))";
    static char fmt_220[] = "(/\002         time     \002,10(a29,2a12,2x))";
    static char fmt_230[] = "(/\002         time     \002,10(a29,3a12,2x))";
    static char fmt_240[] = "(/\002         time     \002,10(a29,4a12,2x))";
    static char fmt_250[] = "(/\002         time     \002,10(a29,5a12,2x))";
    static char fmt_270[] = "(/\002         time     \002,10(a29,6a12,2x))";
    static char fmt_271[] = "(/\002       time     \002,10(a29,7a12,2x))";
    static char fmt_272[] = "(/\002       time     \002,10(a29,8a12,2x))";
    static char fmt_273[] = "(/\002       time     \002,10(a29,9a12,2x))";
    static char fmt_274[] = "(/\002       time     \002,10(a29,10a12,2x))";
    static char fmt_275[] = "(/\002       time     \002,10(a29,11a12,2x))";
    static char fmt_276[] = "(/\002       time     \002,10(a29,12a12,2x))";

    /* System generated locals */
    integer conc_dim1, conc_offset, sorb_dim1, sorb_offset, sorb2_dim1, 
	    sorb2_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    extern integer igetfileversion_(integer *, integer *);
    static real b, c__[5], h__;
    static integer i__, j, l, m, n;
    static real s[5], x1;
    static integer ii;
    static real ax, bx, te, dx, sah, sak;
    static integer nold;
    static real sath;
    static integer iver;
    static char text1[30], text2[30], text3[30];
    static real sbeta;
    static integer nobsa;
    static real sconc[5], shold, ssorb[5], stemp;

    /* Fortran I/O blocks */
    static cilist io___78 = { 0, 6, 0, 0, 0 };
    static cilist io___80 = { 1, 32, 0, 0, 0 };
    static cilist io___83 = { 1, 32, 0, 0, 0 };
    static cilist io___84 = { 1, 32, 0, 0, 0 };
    static cilist io___87 = { 1, 32, 0, 0, 0 };
    static cilist io___97 = { 1, 32, 0, 0, 0 };
    static cilist io___98 = { 1, 32, 0, 0, 0 };
    static cilist io___100 = { 1, 32, 0, 0, 0 };
    static cilist io___102 = { 0, 6, 0, 0, 0 };
    static cilist io___112 = { 1, 50, 0, fmt_110, 0 };
    static cilist io___113 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___114 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___115 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___116 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___117 = { 1, 50, 0, "('end')", 0 };
    static cilist io___118 = { 1, 50, 0, fmt_130, 0 };
    static cilist io___119 = { 1, 32, 0, 0, 0 };
    static cilist io___120 = { 1, 32, 0, 0, 0 };
    static cilist io___125 = { 1, 77, 0, fmt_140, 0 };
    static cilist io___126 = { 1, 77, 0, fmt_200, 0 };
    static cilist io___127 = { 1, 77, 0, fmt_150, 0 };
    static cilist io___128 = { 1, 77, 0, fmt_160, 0 };
    static cilist io___129 = { 1, 77, 0, fmt_170, 0 };
    static cilist io___130 = { 1, 77, 0, fmt_180, 0 };
    static cilist io___131 = { 1, 77, 0, fmt_190, 0 };
    static cilist io___132 = { 1, 77, 0, fmt_260, 0 };
    static cilist io___133 = { 1, 77, 0, fmt_261, 0 };
    static cilist io___134 = { 1, 77, 0, fmt_262, 0 };
    static cilist io___135 = { 1, 77, 0, fmt_263, 0 };
    static cilist io___136 = { 1, 77, 0, fmt_264, 0 };
    static cilist io___137 = { 1, 77, 0, fmt_265, 0 };
    static cilist io___138 = { 1, 77, 0, fmt_266, 0 };
    static cilist io___139 = { 1, 77, 0, fmt_210, 0 };
    static cilist io___140 = { 1, 77, 0, fmt_220, 0 };
    static cilist io___141 = { 1, 77, 0, fmt_230, 0 };
    static cilist io___142 = { 1, 77, 0, fmt_240, 0 };
    static cilist io___143 = { 1, 77, 0, fmt_250, 0 };
    static cilist io___144 = { 1, 77, 0, fmt_270, 0 };
    static cilist io___145 = { 1, 77, 0, fmt_271, 0 };
    static cilist io___146 = { 1, 77, 0, fmt_272, 0 };
    static cilist io___147 = { 1, 77, 0, fmt_273, 0 };
    static cilist io___148 = { 1, 77, 0, fmt_274, 0 };
    static cilist io___149 = { 1, 77, 0, fmt_275, 0 };
    static cilist io___150 = { 1, 77, 0, fmt_276, 0 };


    /* Parameter adjustments */
    --tempo;
    --tempn;
    --ath;
    --ak;
    --ah;
    --beta;
    --laynum;
    --htemp;
    --matnum;
    --hold;
    --hnew;
    --x;
    --node;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    sorb_dim1 = *nsd;
    sorb_offset = 1 + sorb_dim1;
    sorb -= sorb_offset;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    if (*lscreen) {
	s_wsle(&io___78);
	do_lio(&c__9, &c__1, "reading nodal information", (ftnlen)25);
	e_wsle();
    }
    iver = igetfileversion_(&c__32, &c__1);
    i__1 = s_rsle(&io___80);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = s_rsle(&io___83);
	if (i__2 != 0) {
	    goto L901;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L901;
	}
/* L11: */
    }
    i__1 = s_rsle(&io___84);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*numnp), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&ii, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*numnp > *numnpd) {
	*ierr = 3;
	return 0;
    }
    if (*ns > *nsd) {
	*ierr = 5;
	return 0;
    }
/*     Read nodal point information */
    j = *numnp + 1;
L12:
    --j;
    if (! (*lchem) && ! (*ltemp)) {
	i__1 = s_rsle(&io___87);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&x1, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&h__, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&m, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&l, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&b, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&ax, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&bx, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&dx, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	te = 20.f;
    } else if (! (*lchem)) {
	i__1 = s_rsle(&io___97);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&x1, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&h__, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&m, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&l, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&b, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&ax, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&bx, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&dx, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&te, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    } else if (*lequil) {
	i__1 = s_rsle(&io___98);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&x1, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&h__, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&m, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&l, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&b, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&ax, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&bx, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&dx, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&te, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *ns;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__1 = do_lio(&c__4, &c__1, (char *)&c__[ii - 1], (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    } else {
	i__1 = s_rsle(&io___100);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&x1, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&h__, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&m, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&l, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&b, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&ax, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&bx, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&dx, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&te, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *ns;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__1 = do_lio(&c__4, &c__1, (char *)&c__[ii - 1], (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__3 = *ns;
	for (ii = 1; ii <= i__3; ++ii) {
	    i__1 = do_lio(&c__4, &c__1, (char *)&s[ii - 1], (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    n = *numnp - n + 1;
    x[n] = x1;
    hold[n] = h__;
    matnum[n] = m;
    laynum[n] = l;
    beta[n] = b;
    ah[n] = ax;
    ak[n] = bx;
    ath[n] = dx;
    tempo[n] = te;
    i__1 = *ns;
    for (ii = 1; ii <= i__1; ++ii) {
	if (*lchem) {
	    conc[ii + n * conc_dim1] = c__[ii - 1];
	    if (! (*lequil)) {
		sorb[ii + n * sorb_dim1] = s[ii - 1];
	    }
	    if (*lbact || *ldualneq) {
		sorb2[ii + n * sorb2_dim1] = 0.f;
	    }
	}
/* L1: */
    }
    if ((i__1 = j - n) < 0) {
	goto L13;
    } else if (i__1 == 0) {
	goto L18;
    } else {
	goto L14;
    }
L13:
    s_wsle(&io___102);
    do_lio(&c__9, &c__1, "ERROR in NodInf at node =", (ftnlen)25);
    do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
    e_wsle();
    s_stop("", (ftnlen)0);
L14:
    dx = x[nold] - x[n];
    shold = (hold[nold] - hold[n]) / dx;
    sbeta = (beta[nold] - beta[n]) / dx;
    sah = (ah[nold] - ah[n]) / dx;
    sak = (ak[nold] - ak[n]) / dx;
    sath = (ath[nold] - ath[n]) / dx;
    stemp = (tempo[nold] - tempo[n]) / dx;
    if (*lchem) {
	i__1 = *ns;
	for (ii = 1; ii <= i__1; ++ii) {
	    sconc[ii - 1] = (conc[ii + nold * conc_dim1] - conc[ii + n * 
		    conc_dim1]) / dx;
	    ssorb[ii - 1] = (sorb[ii + nold * sorb_dim1] - sorb[ii + n * 
		    sorb_dim1]) / dx;
/* L15: */
	}
    }
    i__1 = n + 1;
    for (i__ = nold - 1; i__ >= i__1; --i__) {
	dx = x[nold] - x[i__];
	hold[i__] = hold[nold] - shold * dx;
	beta[i__] = beta[nold] - sbeta * dx;
	ah[i__] = ah[nold] - sah * dx;
	ak[i__] = ak[nold] - sak * dx;
	ath[i__] = ath[nold] - sath * dx;
	tempo[i__] = tempo[nold] - stemp * dx;
	if (*lchem) {
	    i__2 = *ns;
	    for (ii = 1; ii <= i__2; ++ii) {
		conc[ii + i__ * conc_dim1] = conc[ii + nold * conc_dim1] - 
			sconc[ii - 1] * dx;
		sorb[ii + i__ * sorb_dim1] = sorb[ii + nold * sorb_dim1] - 
			ssorb[ii - 1] * dx;
		if (*lbact || *ldualneq) {
		    sorb2[ii + i__ * sorb2_dim1] = sorb2[ii + nold * 
			    sorb2_dim1] - ssorb[ii - 1] * dx;
		}
/* L16: */
	    }
	}
	matnum[i__] = matnum[i__ + 1];
	laynum[i__] = laynum[i__ + 1];
/* L17: */
    }
    j = n;
L18:
    nold = n;
    if (j > 1) {
	goto L12;
    }
    sbeta = 0.f;
    if (beta[*numnp] > 0.f) {
	sbeta = beta[*numnp] * (x[*numnp] - x[*numnp - 1]) / 2.f;
    }
    i__1 = *numnp - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (beta[i__] > 0.f) {
	    sbeta += beta[i__] * (x[i__ + 1] - x[i__ - 1]) / 2.f;
	}
/* L19: */
    }
    i__1 = *numnp;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (sbeta > 0.f) {
	    beta[i__] /= sbeta;
	} else {
	    beta[i__] = 0.f;
	}
/* L20: */
    }
    *xsurf = x[*numnp];
/*     Print nodal information */
    i__1 = s_wsfe(&io___112);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    for (n = *numnp; n >= 1; --n) {
	if (! (*lchem) && ! (*ltemp)) {
	    i__1 = s_wsfe(&io___113);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__2 = *numnp - n + 1;
	    i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&x[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&hold[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&matnum[n], (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&laynum[n], (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&beta[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ah[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ak[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ath[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	} else if (! (*lchem)) {
	    i__1 = s_wsfe(&io___114);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__2 = *numnp - n + 1;
	    i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&x[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&hold[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&matnum[n], (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&laynum[n], (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&beta[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ah[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ak[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ath[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&tempo[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	} else if (*lequil) {
	    i__1 = s_wsfe(&io___115);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__2 = *numnp - n + 1;
	    i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&x[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&hold[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&matnum[n], (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&laynum[n], (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&beta[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ah[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ak[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ath[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&tempo[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__3 = *ns;
	    for (ii = 1; ii <= i__3; ++ii) {
		i__1 = do_fio(&c__1, (char *)&conc[ii + n * conc_dim1], (
			ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L902;
		}
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	} else {
	    i__1 = s_wsfe(&io___116);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__2 = *numnp - n + 1;
	    i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&x[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&hold[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&matnum[n], (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&laynum[n], (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&beta[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ah[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ak[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&ath[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&tempo[n], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__3 = *ns;
	    for (ii = 1; ii <= i__3; ++ii) {
		i__1 = do_fio(&c__1, (char *)&conc[ii + n * conc_dim1], (
			ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L902;
		}
	    }
	    i__4 = *ns;
	    for (ii = 1; ii <= i__4; ++ii) {
		i__1 = do_fio(&c__1, (char *)&sorb[ii + n * sorb_dim1], (
			ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L902;
		}
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	}
	hnew[n] = hold[n];
	htemp[n] = hold[n];
	tempn[n] = tempo[n];
/* L21: */
    }
    i__1 = s_wsfe(&io___117);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    *hbot = hnew[1];
    *htop = hnew[*numnp];
    i__1 = s_wsfe(&io___118);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_rsle(&io___119);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*nobs), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*nobs > *nobsd) {
	*ierr = 4;
	return 0;
    }
    if (*nobs > 0) {
	i__1 = s_rsle(&io___120);
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *nobs;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = do_lio(&c__3, &c__1, (char *)&node[i__], (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = *nobs;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    node[i__] = *numnp - node[i__] + 1;
/* L22: */
	}
	if (*lprint) {
	    nobsa = min(10,*nobs);
	    s_copy(text3, "Node(", (ftnlen)30, (ftnlen)5);
	    s_copy(text1, "    h        theta    Temp   ", (ftnlen)30, (
		    ftnlen)29);
	    if (*lflux) {
		s_copy(text1, "    h        theta    Flux   ", (ftnlen)30, (
			ftnlen)29);
	    }
	    s_copy(text2, "   Conc     ", (ftnlen)30, (ftnlen)12);
	    if (! (*lchem)) {
		i__1 = s_wsfe(&io___125);
		if (i__1 != 0) {
		    goto L902;
		}
		i__2 = nobsa;
		for (j = 1; j <= i__2; ++j) {
		    i__1 = do_fio(&c__1, text3, (ftnlen)30);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = *numnp - node[j] + 1;
		    i__1 = do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(
			    integer));
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L902;
		}
		i__1 = s_wsfe(&io___126);
		if (i__1 != 0) {
		    goto L902;
		}
		i__3 = nobsa;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__1 = do_fio(&c__1, text1, (ftnlen)30);
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L902;
		}
	    } else {
		if (*ns == 1) {
		    i__1 = s_wsfe(&io___127);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (j = 1; j <= i__3; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 2) {
		    i__1 = s_wsfe(&io___128);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (j = 1; j <= i__2; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 3) {
		    i__1 = s_wsfe(&io___129);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (j = 1; j <= i__3; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 4) {
		    i__1 = s_wsfe(&io___130);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (j = 1; j <= i__2; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 5) {
		    i__1 = s_wsfe(&io___131);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (j = 1; j <= i__3; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 6) {
		    i__1 = s_wsfe(&io___132);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (j = 1; j <= i__2; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 7) {
		    i__1 = s_wsfe(&io___133);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (j = 1; j <= i__3; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 8) {
		    i__1 = s_wsfe(&io___134);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (j = 1; j <= i__2; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 9) {
		    i__1 = s_wsfe(&io___135);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (j = 1; j <= i__3; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 10) {
		    i__1 = s_wsfe(&io___136);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (j = 1; j <= i__2; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 11) {
		    i__1 = s_wsfe(&io___137);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (j = 1; j <= i__3; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 12) {
		    i__1 = s_wsfe(&io___138);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (j = 1; j <= i__2; ++j) {
			i__1 = do_fio(&c__1, text3, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *numnp - node[j] + 1;
			i__1 = do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L902;
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 1) {
		    i__1 = s_wsfe(&io___139);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 2) {
		    i__1 = s_wsfe(&io___140);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 3) {
		    i__1 = s_wsfe(&io___141);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 4) {
		    i__1 = s_wsfe(&io___142);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 5) {
		    i__1 = s_wsfe(&io___143);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 6) {
		    i__1 = s_wsfe(&io___144);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 7) {
		    i__1 = s_wsfe(&io___145);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 8) {
		    i__1 = s_wsfe(&io___146);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 9) {
		    i__1 = s_wsfe(&io___147);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 10) {
		    i__1 = s_wsfe(&io___148);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 11) {
		    i__1 = s_wsfe(&io___149);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__3 = nobsa;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
		if (*ns == 12) {
		    i__1 = s_wsfe(&io___150);
		    if (i__1 != 0) {
			goto L902;
		    }
		    i__2 = nobsa;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, text1, (ftnlen)30);
			if (i__1 != 0) {
			    goto L902;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, text2, (ftnlen)30);
			    if (i__1 != 0) {
				goto L902;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L902;
		    }
		}
	    }
	}
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
/*     Error when writing into an output file */
L902:
    *ierr = 2;
    return 0;
/* L301: */
/* L302: */
} /* nodinf_ */

/* *********************************************************************** */
/* Subroutine */ int initw_(integer *numnp, integer *nmat, integer *matnum, 
	integer *kappa, real *hnew, real *hold, real *htemp, real *pard, real 
	*parw, integer *imodel, real *htop, real *hbot, integer *idualpor, 
	real *thnewim, integer *ierr)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static real thmobile;
    static integer i__, m;
    extern doublereal fh_(integer *, real *, real *);
    static real qe, thtotal;

    /* Parameter adjustments */
    --thnewim;
    --htemp;
    --hold;
    --hnew;
    --kappa;
    --matnum;
    parw -= 12;
    pard -= 12;

    /* Function Body */
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	thtotal = hnew[i__];
	if (*idualpor > 0) {
	    hnew[i__] = hnew[i__] * pard[m * 11 + 2] / (pard[m * 11 + 2] + 
		    pard[m * 11 + 8]);
	    thmobile = hnew[i__];
	    thnewim[i__] = thtotal - thmobile;
	}
	if (kappa[i__] == -1) {
/* Computing MIN */
	    r__1 = (hnew[i__] - pard[m * 11 + 1]) / (pard[m * 11 + 2] - pard[
		    m * 11 + 1]);
	    qe = dmin(r__1,1.f);
	    if (qe < 0.f) {
		goto L901;
	    }
	    hnew[i__] = fh_(imodel, &qe, &pard[m * 11 + 1]);
	} else {
/* Computing MIN */
	    r__1 = (hnew[i__] - parw[m * 11 + 1]) / (parw[m * 11 + 2] - parw[
		    m * 11 + 1]);
	    qe = dmin(r__1,1.f);
	    if (qe < 0.f) {
		goto L901;
	    }
	    hnew[i__] = fh_(imodel, &qe, &parw[m * 11 + 1]);
	}
	hold[i__] = hnew[i__];
	htemp[i__] = hnew[i__];
/* L11: */
    }
    *hbot = hnew[1];
    *htop = hnew[*numnp];
    return 0;
L901:
    *ierr = 1;
    return 0;
} /* initw_ */

/* *********************************************************************** */
/* Subroutine */ int initdualpor_(integer *numnp, integer *nmat, integer *
	matnum, real *par, real *theta, integer *idualpor, real *thnewim, 
	real *tholdim, real *sinkim, real *hnew, real *strans, logical *
	linitw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m;
    extern doublereal fq_(integer *, real *, real *);
    static real se;

    /* Parameter adjustments */
    --strans;
    --hnew;
    --sinkim;
    --tholdim;
    --thnewim;
    --theta;
    --matnum;
    par -= 12;

    /* Function Body */
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	if (*idualpor == 0) {
	    thnewim[i__] = 0.f;
	}
	if (! (*linitw)) {
	    if (*idualpor == 1) {
		se = (theta[i__] - par[m * 11 + 1]) / (par[m * 11 + 2] - par[
			m * 11 + 1]);
		thnewim[i__] = par[m * 11 + 7] + se * (par[m * 11 + 8] - par[
			m * 11 + 7]);
	    } else if (*idualpor == 2) {
		thnewim[i__] = fq_(&c__0, &hnew[i__], &par[m * 11 + 7]);
	    }
	}
	tholdim[i__] = thnewim[i__];
	sinkim[i__] = 0.f;
	strans[i__] = 0.f;
/* L11: */
    }
    return 0;
} /* initdualpor_ */

/* *********************************************************************** */
/* Subroutine */ int matin_(integer *nmat, real *pard, real *parw, real *
	htab1, real *htabn, logical *lscreen, integer *ierr, integer *numnp, 
	real *ah, integer *ihyst, real *ahw, real *athw, real *akw, integer *
	matnum, real *hnew, integer *kappa, real *aths, real *thrr, real *
	conr, real *aks, integer *kappao, integer *imodel, real *xconv, 
	logical *ltable, integer *ikappa, integer *ntabmod, integer *idualpor)
{
    /* Format strings */
    static char fmt_110[] = "(//\002MatNum, Param. array:\002//\002   Mat   "
	    "  Qr     Qs        \002,\002Alfa         n          Ks      l"
	    "\002/)";
    static char fmt_111[] = "(//\002MatNum, Param. array:\002//\002   Mat   "
	    "  Qr     Qs        \002,\002Alfa         n          Ks      l   "
	    "    Qm     Qa     Qk       Kk\002/)";
    static char fmt_112[] = "(//\002MatNum, Param. array:\002//\002   Mat   "
	    "  Qr     Qs        \002,\002Alfa         n          Ks      l   "
	    "    Qm     QsW  AlfaW      KsW\002/)";
    static char fmt_113[] = "(//\002MatNum, Param. array:\002//\002   Mat   "
	    "  Qr     Qs        \002,\002Alfa         n          Ks      l   "
	    "    W2       Alfa2        n2\002/)";
    static char fmt_115[] = "(//\002MatNum, Param. array:\002//\002   Mat   "
	    "  Qr     Qs        \002,\002Alfa         n          Ks      l   "
	    "   QrIm   QsIm  Omega\002/)";
    static char fmt_116[] = "(//\002MatNum, Param. array:\002//\002   Mat   "
	    "  Qr     Qs        \002,\002Alfa         n          Ks      l   "
	    "   QrIm   QsIm   Alfa2        n2  Omega\002/)";
    static char fmt_114[] = "(//\002MatNum, Param. array:\002//\002   Mat   "
	    "  Qr     Qs        Ks\002/)";
    static char fmt_121[] = "(i5,2x,2f7.3,3e12.3,2f7.3,3e12.3)";
    static char fmt_120[] = "(i5,2x,2f7.3,3e12.3,4f7.3,2e12.3)";
    static char fmt_117[] = "(//\002   Saturation   Air-Water Interfacial Ar"
	    "ea\002/)";
    static char fmt_122[] = "(2f10.5)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static real g;
    static integer i__, m;
    static real ae[100];
    extern doublereal fh_(integer *, real *, real *);
    static real row;
    static integer itab, npar;
    static real temp, sint, sigma;
    extern /* Subroutine */ int qromb_(real *, real *, real *, integer *, 
	    real *);
    static real rhentry;

    /* Fortran I/O blocks */
    static cilist io___159 = { 0, 6, 0, 0, 0 };
    static cilist io___160 = { 1, 30, 0, 0, 0 };
    static cilist io___161 = { 1, 30, 0, 0, 0 };
    static cilist io___162 = { 1, 30, 0, 0, 0 };
    static cilist io___163 = { 1, 30, 0, 0, 0 };
    static cilist io___164 = { 0, 6, 0, 0, 0 };
    static cilist io___165 = { 0, 6, 0, 0, 0 };
    static cilist io___166 = { 0, 5, 0, 0, 0 };
    static cilist io___167 = { 1, 30, 0, 0, 0 };
    static cilist io___168 = { 1, 30, 0, 0, 0 };
    static cilist io___170 = { 1, 50, 0, fmt_110, 0 };
    static cilist io___171 = { 1, 50, 0, fmt_111, 0 };
    static cilist io___172 = { 1, 50, 0, fmt_112, 0 };
    static cilist io___173 = { 1, 50, 0, fmt_113, 0 };
    static cilist io___174 = { 1, 50, 0, fmt_115, 0 };
    static cilist io___175 = { 1, 50, 0, fmt_116, 0 };
    static cilist io___176 = { 1, 50, 0, fmt_114, 0 };
    static cilist io___177 = { 1, 30, 0, 0, 0 };
    static cilist io___181 = { 1, 30, 0, 0, 0 };
    static cilist io___182 = { 1, 50, 0, fmt_121, 0 };
    static cilist io___183 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___184 = { 1, 30, 0, 0, 0 };
    static cilist io___185 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___192 = { 0, 50, 0, fmt_117, 0 };
    static cilist io___194 = { 0, 50, 0, fmt_122, 0 };


    /* Parameter adjustments */
    --akw;
    --athw;
    --ahw;
    parw -= 12;
    pard -= 12;
    --kappao;
    --aks;
    --conr;
    --thrr;
    --aths;
    --kappa;
    --hnew;
    --matnum;
    --ah;

    /* Function Body */
    if (*lscreen) {
	s_wsle(&io___159);
	do_lio(&c__9, &c__1, "reading material information", (ftnlen)28);
	e_wsle();
    }
    i__1 = s_rsle(&io___160);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___161);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*htab1), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*htabn), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___162);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___163);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*imodel), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ihyst), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
/*     iModel = 0: van Genuchten */
/*              1: modified van Genuchten (Vogel and Cislerova) */
/*              2: Brooks and Corey */
/*              3: van Genuchte with air entry value of 2 cm */
/*              4: log-normal (Kosugi) */
/*              5: dual-porosity function (Durner) */
/*              6: dual-porosity system with transfer proportional to water content */
/*              7: dual-porosity system with transfer proportional to pressure head */
/*              8: dual-permeability system, not handled by this program */
/*              9: nTabMod: general tables (nTabMod) */
    if (*imodel == 8) {
	s_wsle(&io___164);
	do_lio(&c__9, &c__1, "Dual-permeability models are implemented in di"
		"fferent code !!", (ftnlen)61);
	e_wsle();
	s_wsle(&io___165);
	do_lio(&c__9, &c__1, "Press Enter to continue", (ftnlen)23);
	e_wsle();
	s_rsle(&io___166);
	e_rsle();
	s_stop("", (ftnlen)0);
    }
    if (*imodel < *ntabmod) {
/* Computing MIN */
	r__1 = dabs(*htab1), r__2 = dabs(*htabn);
	*htab1 = -dmin(r__1,r__2);
/* Computing MAX */
	r__1 = dabs(*htab1), r__2 = dabs(*htabn);
	*htabn = -dmax(r__1,r__2);
	*ltable = TRUE_;
	if (*htab1 > -1e-5f && *htabn > -1e-5f || *htab1 == *htabn) {
	    *ltable = FALSE_;
	    *htab1 = *xconv * -1e-4f;
	    *htabn = *xconv * -100.f;
	}
    } else {
	*ltable = TRUE_;
    }
    if (*ihyst > 0) {
	i__1 = s_rsle(&io___167);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___168);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*ikappa), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    } else {
	*ikappa = -1;
    }
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	kappa[i__] = *ikappa;
	kappao[i__] = *ikappa;
/* L11: */
    }
    if (*imodel == 2 || *imodel == 4 || *imodel == 0 && *ihyst == 0) {
	i__1 = s_wsfe(&io___170);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    } else if (*imodel == 1 || *imodel == 3) {
	i__1 = s_wsfe(&io___171);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    } else if (*imodel == 0 && *ihyst > 0) {
	i__1 = s_wsfe(&io___172);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    } else if (*imodel == 5) {
	i__1 = s_wsfe(&io___173);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    } else if (*imodel == 6) {
	i__1 = s_wsfe(&io___174);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    } else if (*imodel == 7) {
	i__1 = s_wsfe(&io___175);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    } else if (*imodel == *ntabmod) {
	i__1 = s_wsfe(&io___176);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    }
    i__1 = s_rsle(&io___177);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    rhentry = *xconv * .02f;
    if (*imodel == 0 || *imodel == 2 || *imodel == 3 || *imodel == 4) {
	npar = 6;
    } else if (*imodel == 1) {
	npar = 10;
    } else if (*imodel == 5) {
	npar = 9;
    } else if (*imodel == 6) {
	npar = 9;
	*imodel = 0;
	*idualpor = 1;
    } else if (*imodel == 7) {
	npar = 11;
	*imodel = 0;
	*idualpor = 2;
    } else if (*imodel == *ntabmod) {
	npar = 3;
    } else {
	npar = 6;
    }
    i__1 = *nmat;
    for (m = 1; m <= i__1; ++m) {
	if (*ihyst == 0) {
	    i__2 = s_rsle(&io___181);
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__3 = npar;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__2 = do_lio(&c__4, &c__1, (char *)&pard[i__ + m * 11], (
			ftnlen)sizeof(real));
		if (i__2 != 0) {
		    goto L901;
		}
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    if (*imodel == 1) {
/* Computing MAX */
		r__1 = pard[m * 11 + 7], r__2 = pard[m * 11 + 2];
		pard[m * 11 + 7] = dmax(r__1,r__2);
/* Computing MIN */
		r__1 = pard[m * 11 + 8], r__2 = pard[m * 11 + 1];
		pard[m * 11 + 8] = dmin(r__1,r__2);
/* Computing MIN */
		r__1 = pard[m * 11 + 9], r__2 = pard[m * 11 + 2];
		pard[m * 11 + 9] = dmin(r__1,r__2);
/* Computing MIN */
		r__1 = pard[m * 11 + 10], r__2 = pard[m * 11 + 5];
		pard[m * 11 + 10] = dmin(r__1,r__2);
	    } else if (*imodel == 3) {
		d__2 = (doublereal) (pard[m * 11 + 3] * rhentry);
		d__3 = (doublereal) pard[m * 11 + 4];
		d__1 = (doublereal) (pow_dd(&d__2, &d__3) + 1.f);
		d__4 = (doublereal) (1.f - 1.f / pard[m * 11 + 4]);
		pard[m * 11 + 7] = pard[m * 11 + 1] + (pard[m * 11 + 2] - 
			pard[m * 11 + 1]) * pow_dd(&d__1, &d__4);
	    }
	    if (*imodel == *ntabmod) {
		i__2 = s_wsfe(&io___182);
		if (i__2 != 0) {
		    goto L902;
		}
		i__2 = do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
		if (i__2 != 0) {
		    goto L902;
		}
		i__3 = npar;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__2 = do_fio(&c__1, (char *)&pard[i__ + m * 11], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L902;
		    }
		}
		i__2 = e_wsfe();
		if (i__2 != 0) {
		    goto L902;
		}
	    } else {
		i__2 = s_wsfe(&io___183);
		if (i__2 != 0) {
		    goto L902;
		}
		i__2 = do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
		if (i__2 != 0) {
		    goto L902;
		}
		i__3 = npar;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__2 = do_fio(&c__1, (char *)&pard[i__ + m * 11], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L902;
		    }
		}
		i__2 = e_wsfe();
		if (i__2 != 0) {
		    goto L902;
		}
	    }
	} else {
/* Hysteresis */
	    i__2 = s_rsle(&io___184);
	    if (i__2 != 0) {
		goto L901;
	    }
	    for (i__ = 1; i__ <= 7; ++i__) {
		i__2 = do_lio(&c__4, &c__1, (char *)&pard[i__ + m * 11], (
			ftnlen)sizeof(real));
		if (i__2 != 0) {
		    goto L901;
		}
	    }
	    i__2 = do_lio(&c__4, &c__1, (char *)&parw[m * 11 + 2], (ftnlen)
		    sizeof(real));
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = do_lio(&c__4, &c__1, (char *)&parw[m * 11 + 3], (ftnlen)
		    sizeof(real));
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = do_lio(&c__4, &c__1, (char *)&parw[m * 11 + 5], (ftnlen)
		    sizeof(real));
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
/* Computing MAX */
	    r__1 = pard[m * 11 + 7], r__2 = pard[m * 11 + 2];
	    pard[m * 11 + 7] = dmax(r__1,r__2);
	    i__2 = s_wsfe(&io___185);
	    if (i__2 != 0) {
		goto L902;
	    }
	    i__2 = do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L902;
	    }
	    for (i__ = 1; i__ <= 7; ++i__) {
		i__2 = do_fio(&c__1, (char *)&pard[i__ + m * 11], (ftnlen)
			sizeof(real));
		if (i__2 != 0) {
		    goto L902;
		}
	    }
	    i__2 = do_fio(&c__1, (char *)&parw[m * 11 + 2], (ftnlen)sizeof(
		    real));
	    if (i__2 != 0) {
		goto L902;
	    }
	    i__2 = do_fio(&c__1, (char *)&parw[m * 11 + 3], (ftnlen)sizeof(
		    real));
	    if (i__2 != 0) {
		goto L902;
	    }
	    i__2 = do_fio(&c__1, (char *)&parw[m * 11 + 5], (ftnlen)sizeof(
		    real));
	    if (i__2 != 0) {
		goto L902;
	    }
	    i__2 = e_wsfe();
	    if (i__2 != 0) {
		goto L902;
	    }
	    parw[m * 11 + 1] = pard[m * 11 + 1];
	    parw[m * 11 + 4] = pard[m * 11 + 4];
	    ahw[m] = pard[m * 11 + 3] / parw[m * 11 + 3];
	    athw[m] = (parw[m * 11 + 2] - parw[m * 11 + 1]) / (pard[m * 11 + 
		    2] - pard[m * 11 + 1]);
	    akw[m] = 1.f;
	    if (*ihyst == 2) {
		akw[m] = parw[m * 11 + 5] / pard[m * 11 + 5];
	    }
	    parw[m * 11 + 7] = parw[m * 11 + 1] + athw[m] * (pard[m * 11 + 7] 
		    - pard[m * 11 + 1]);
	    pard[m * 11 + 8] = pard[m * 11 + 1];
	    pard[m * 11 + 9] = pard[m * 11 + 2];
	    pard[m * 11 + 10] = pard[m * 11 + 5];
	    parw[m * 11 + 8] = parw[m * 11 + 1];
	    parw[m * 11 + 9] = parw[m * 11 + 2];
	    parw[m * 11 + 10] = parw[m * 11 + 5];
	    parw[m * 11 + 6] = pard[m * 11 + 6];
	}
/* L12: */
    }
/*     Hysteresis Update for Initial Pressure Head Distributions */
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	if (*imodel < *ntabmod) {
/* Computing MAX */
	    r__1 = hnew[i__], r__2 = ah[i__] * fh_(imodel, &c_b678, &pard[m * 
		    11 + 1]);
	    hnew[i__] = dmax(r__1,r__2);
	}
	aths[i__] = 1.f;
	aks[i__] = 1.f;
	thrr[i__] = pard[m * 11 + 1];
	conr[i__] = 0.f;
/* L13: */
    }
    return 0;
/*     Calculate the air-water interfacial area */
    itab = 101;
    ae[0] = 0.f;
    g = *xconv * 9.81f;
/* Gravitational acceleration from [m/s2] t */
    temp = 20.f;
/* Computing 2nd power */
    r__1 = temp;
    sigma = 75.6f - temp * .1425f - r__1 * r__1 * 2.38e-4f;
/* Surface tension */
/* Computing 2nd power */
    r__1 = temp - 4.f;
/* Computing 3rd power */
    r__2 = temp - 4.f;
    row = 1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 3.79e-8f;
/* Density of soil */
    row = row * 1e6f / *xconv / *xconv / *xconv;
/* to [g/L3] */
    m = 1;
    s_wsfe(&io___192);
    e_wsfe();
    i__1 = itab;
    for (i__ = 2; i__ <= i__1; ++i__) {
	r__1 = 1.f - (i__ - 1) * .01f;
	qromb_(&c_b682, &r__1, &sint, imodel, &pard[12]);
	ae[i__ - 1] = sint * pard[m * 11 + 2] * row * g / sigma;
	s_wsfe(&io___194);
	r__1 = 1.f - (i__ - 1) * .01f;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ae[i__ - 1], (ftnlen)sizeof(real));
	e_wsfe();
/* L14: */
    }
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
/*     Error when writing into an output file */
L902:
    *ierr = 2;
    return 0;
} /* matin_ */

/* *********************************************************************** */
/* Subroutine */ int hysterin_(integer *numnp, integer *nmat, real *hold, 
	integer *matnum, real *pard, real *parw, real *thnew, real *thold, 
	integer *kappa, real *aths, real *thrr, real *cono, real *conr, real *
	aks, integer *kappao, real *ah, real *ak, integer *ihyst, integer *
	imodel, char *cdatapath, ftnlen cdatapath_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    olist o__1;

    /* Local variables */
    extern integer len_trim__(char *, ftnlen);
    static integer i__, n;
    static char cfilename[260];
    static integer ilengthpath;
    extern /* Subroutine */ int hyster_(integer *, integer *, real *, integer 
	    *, real *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, integer *, 
	    integer *, real *);

    /* Fortran I/O blocks */
    static cilist io___197 = { 1, 35, 0, 0, 0 };
    static cilist io___198 = { 1, 35, 0, 0, 0 };
    static cilist io___200 = { 1, 35, 0, 0, 0 };


    /* Parameter adjustments */
    --ak;
    --ah;
    --kappao;
    --aks;
    --conr;
    --cono;
    --thrr;
    --aths;
    --kappa;
    --thold;
    --thnew;
    --matnum;
    --hold;
    parw -= 12;
    pard -= 12;

    /* Function Body */
    ilengthpath = len_trim__(cdatapath, (ftnlen)260);
/* Writing concatenation */
    i__1[0] = ilengthpath, a__1[0] = cdatapath;
    i__1[1] = 13, a__1[1] = "Hysteresis.in";
    s_cat(cfilename, a__1, i__1, &c__2, (ftnlen)260);
    o__1.oerr = 1;
    o__1.ounit = 35;
    o__1.ofnmlen = 260;
    o__1.ofnm = cfilename;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__2 = f_open(&o__1);
    if (i__2 != 0) {
	goto L901;
    }
    i__2 = s_rsle(&io___197);
    if (i__2 != 0) {
	goto L901;
    }
    i__2 = e_rsle();
    if (i__2 != 0) {
	goto L901;
    }
/* Input file "Hyster.in" */
    i__2 = s_rsle(&io___198);
    if (i__2 != 0) {
	goto L901;
    }
    i__2 = e_rsle();
    if (i__2 != 0) {
	goto L901;
    }
/* n,Theta,Kappa */
    i__2 = *numnp;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__3 = s_rsle(&io___200);
	if (i__3 != 0) {
	    goto L901;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	if (i__3 != 0) {
	    goto L901;
	}
	i__3 = do_lio(&c__4, &c__1, (char *)&thold[i__], (ftnlen)sizeof(real))
		;
	if (i__3 != 0) {
	    goto L901;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&kappao[i__], (ftnlen)sizeof(
		integer));
	if (i__3 != 0) {
	    goto L901;
	}
	i__3 = e_rsle();
	if (i__3 != 0) {
	    goto L901;
	}
	thnew[i__] = thold[i__];
	kappa[i__] = kappao[i__];
/* L11: */
    }
    hyster_(numnp, nmat, &hold[1], &matnum[1], &pard[12], &parw[12], &thnew[1]
	    , &thold[1], &kappa[1], &aths[1], &thrr[1], &cono[1], &conr[1], &
	    aks[1], &kappao[1], &ah[1], &ak[1], ihyst, imodel, &c_b699);
L901:
    return 0;
} /* hysterin_ */

/* *********************************************************************** */
/* Subroutine */ int genmat_(integer *ntab, integer *ntabd, integer *nmat, 
	real *thr, real *ths, real *hsat, real *par, real *htab, real *contab,
	 real *captab, real *consat, real *thetab, integer *imodel, logical *
	lscreen, integer *ntabmod, real *consmax, real *xconv, real *tconv, 
	integer *ierr)
{
    /* Format strings */
    static char fmt_110[] = "(/7x,\002Table of Hydraulic Properties which ar"
	    "e interpolated in simulation\002/7x,65(\002=\002)/)";
    static char fmt_120[] = "(\002  theta         h        log h        C   "
	    "          K\002,\002        log K          S          Kv\002)";
    static char fmt_111[] = "(/7x,\002Hydraulic Properties which are interpo"
	    "lated from input tables in simulation\002/7x,75(\002=\002)/)";
    static char fmt_130[] = "(f8.4,e12.3,e12.4,e12.4,e12.4,e12.4,f10.4,e12.4)"
	    ;
    static char fmt_140[] = "(\002end\002)";

    /* System generated locals */
    integer htab_dim1, htab_offset, contab_dim1, contab_offset, captab_dim1, 
	    captab_offset, thetab_dim1, thetab_offset, i__1, i__2, i__3;
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, m;
    extern doublereal fc_(integer *, real *, real *), fh_(integer *, real *, 
	    real *), fk_(integer *, real *, real *);
    static real qe;
    extern doublereal fq_(integer *, real *, real *), fs_(integer *, real *, 
	    real *);
    static real a10h, a10k, alh, dlh;
    static integer icap;
    static real conv, htab1, htabn;
    extern doublereal convh_(real *, real *, real *, real *, real *);

    /* Fortran I/O blocks */
    static cilist io___202 = { 0, 6, 0, 0, 0 };
    static cilist io___203 = { 1, 50, 0, fmt_110, 0 };
    static cilist io___204 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___210 = { 1, 36, 0, 0, 0 };
    static cilist io___211 = { 1, 36, 0, 0, 0 };
    static cilist io___213 = { 1, 50, 0, fmt_111, 0 };
    static cilist io___214 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___216 = { 1, 50, 0, 0, 0 };
    static cilist io___221 = { 1, 50, 0, fmt_130, 0 };
    static cilist io___222 = { 0, 36, 0, 0, 0 };
    static cilist io___223 = { 0, 36, 0, 0, 0 };
    static cilist io___224 = { 0, 36, 0, 0, 0 };
    static cilist io___225 = { 0, 36, 0, 0, 0 };
    static cilist io___226 = { 0, 36, 0, 0, 0 };
    static cilist io___227 = { 1, 50, 0, 0, 0 };
    static cilist io___228 = { 1, 50, 0, fmt_130, 0 };
    static cilist io___229 = { 1, 50, 0, fmt_140, 0 };


    /* Parameter adjustments */
    thetab_dim1 = *ntabd;
    thetab_offset = 1 + thetab_dim1;
    thetab -= thetab_offset;
    --consat;
    captab_dim1 = *ntabd;
    captab_offset = 1 + captab_dim1;
    captab -= captab_offset;
    contab_dim1 = *ntabd;
    contab_offset = 1 + contab_dim1;
    contab -= contab_offset;
    htab_dim1 = *ntabd;
    htab_offset = 1 + htab_dim1;
    htab -= htab_offset;
    par -= 12;
    --hsat;
    --ths;
    --thr;
    --ntab;

    /* Function Body */
    if (*lscreen) {
	s_wsle(&io___202);
	do_lio(&c__9, &c__1, "generating materials", (ftnlen)20);
	e_wsle();
    }
    if (*imodel < *ntabmod) {
	i__1 = s_wsfe(&io___203);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_wsfe(&io___204);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	htab1 = htab[htab_dim1 + 1];
	htabn = htab[ntab[1] + htab_dim1];
	r__1 = -htabn;
	r__2 = -htab1;
	dlh = (r_lg10(&r__1) - r_lg10(&r__2)) / (ntab[1] - 1);
	i__1 = ntab[1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__1 = -htab1;
	    alh = r_lg10(&r__1) + (i__ - 1) * dlh;
	    d__1 = (doublereal) alh;
	    htab[i__ + htab_dim1] = -pow_dd(&c_b710, &d__1);
/* L11: */
	}
    } else {
	i__1 = s_rsle(&io___210);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___211);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&icap, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
/* (=1; input capacity; =0: do not input c */
	i__1 = s_wsfe(&io___213);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_wsfe(&io___214);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    *consmax = 0.f;
/*      ConSMax=1e+30 */
    i__1 = *nmat;
    for (m = 1; m <= i__1; ++m) {
	if (*imodel < *ntabmod) {
	    hsat[m] = fh_(imodel, &c_b682, &par[m * 11 + 1]);
	    consat[m] = par[m * 11 + 5];
	    if (consat[m] > *consmax) {
		*consmax = consat[m];
	    }
/*          if(ConSat(M).lt.ConSMax) ConSMax=ConSat(M) */
	    thr[m] = par[m * 11 + 1];
	    ths[m] = par[m * 11 + 2];
	    i__2 = s_wsle(&io___216);
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = e_wsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = ntab[1];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		contab[i__ + m * contab_dim1] = fk_(imodel, &htab[i__ + 
			htab_dim1], &par[m * 11 + 1]);
		captab[i__ + m * captab_dim1] = fc_(imodel, &htab[i__ + 
			htab_dim1], &par[m * 11 + 1]);
		thetab[i__ + m * thetab_dim1] = fq_(imodel, &htab[i__ + 
			htab_dim1], &par[m * 11 + 1]);
		qe = fs_(imodel, &htab[i__ + htab_dim1], &par[m * 11 + 1]);
		conv = convh_(&htab[i__ + htab_dim1], &thetab[i__ + m * 
			thetab_dim1], &ths[m], xconv, tconv);
/* Computing MAX */
		r__2 = -htab[i__ + htab_dim1];
		r__1 = dmax(r__2,1e-30f);
		a10h = r_lg10(&r__1);
		a10k = r_lg10(&contab[i__ + m * contab_dim1]);
		i__3 = s_wsfe(&io___221);
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&thetab[i__ + m * thetab_dim1], (
			ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&htab[i__ + htab_dim1], (ftnlen)
			sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&a10h, (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&captab[i__ + m * captab_dim1], (
			ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&contab[i__ + m * contab_dim1], (
			ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&a10k, (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&qe, (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&conv, (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = e_wsfe();
		if (i__3 != 0) {
		    goto L901;
		}
/* L12: */
	    }
	} else if (*imodel == *ntabmod) {
/* Table */
	    s_rsle(&io___222);
	    e_rsle();
	    s_rsle(&io___223);
	    do_lio(&c__3, &c__1, (char *)&ntab[m], (ftnlen)sizeof(integer));
	    e_rsle();
	    s_rsle(&io___224);
	    e_rsle();
	    hsat[m] = 0.f;
	    consat[m] = par[m * 11 + 3];
	    thr[m] = par[m * 11 + 1];
	    ths[m] = par[m * 11 + 2];
	    i__2 = ntab[m];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (icap == 1) {
		    s_rsle(&io___225);
		    do_lio(&c__4, &c__1, (char *)&thetab[i__ + m * 
			    thetab_dim1], (ftnlen)sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&htab[i__ + m * htab_dim1], (
			    ftnlen)sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&contab[i__ + m * 
			    contab_dim1], (ftnlen)sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&captab[i__ + m * 
			    captab_dim1], (ftnlen)sizeof(real));
		    e_rsle();
		} else {
		    s_rsle(&io___226);
		    do_lio(&c__4, &c__1, (char *)&thetab[i__ + m * 
			    thetab_dim1], (ftnlen)sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&htab[i__ + m * htab_dim1], (
			    ftnlen)sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&contab[i__ + m * 
			    contab_dim1], (ftnlen)sizeof(real));
		    e_rsle();
		}
/* L15: */
	    }
	    i__2 = s_wsle(&io___227);
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = e_wsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = ntab[m];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (icap == 0) {
		    if (i__ == 1) {
			captab[i__ + m * captab_dim1] = (thetab[m * 
				thetab_dim1 + 2] - thetab[m * thetab_dim1 + 1]
				) / (htab[m * htab_dim1 + 2] - htab[m * 
				htab_dim1 + 1]);
		    } else if (i__ == ntab[m]) {
			captab[i__ + m * captab_dim1] = (thetab[ntab[m] + m * 
				thetab_dim1] - thetab[ntab[m] - 1 + m * 
				thetab_dim1]) / (htab[ntab[m] + m * htab_dim1]
				 - htab[ntab[m] - 1 + m * htab_dim1]);
		    } else {
			captab[i__ + m * captab_dim1] = (thetab[i__ + 1 + m * 
				thetab_dim1] - thetab[i__ - 1 + m * 
				thetab_dim1]) / (htab[i__ + 1 + m * htab_dim1]
				 - htab[i__ - 1 + m * htab_dim1]);
		    }
		}
		qe = (thetab[i__ + m * thetab_dim1] - par[m * 11 + 1]) / (par[
			m * 11 + 2] - par[m * 11 + 1]);
/* Computing MAX */
		r__2 = -htab[i__ + m * htab_dim1];
		r__1 = dmax(r__2,1e-30f);
		a10h = r_lg10(&r__1);
		a10k = r_lg10(&contab[i__ + m * contab_dim1]);
		i__3 = s_wsfe(&io___228);
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&thetab[i__ + m * thetab_dim1], (
			ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&htab[i__ + m * htab_dim1], (
			ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&a10h, (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&captab[i__ + m * captab_dim1], (
			ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&contab[i__ + m * contab_dim1], (
			ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&a10k, (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = do_fio(&c__1, (char *)&qe, (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
		i__3 = e_wsfe();
		if (i__3 != 0) {
		    goto L901;
		}
/* L16: */
	    }
	}
	i__2 = s_wsfe(&io___229);
	if (i__2 != 0) {
	    goto L901;
	}
	i__2 = e_wsfe();
	if (i__2 != 0) {
	    goto L901;
	}
/* L13: */
    }
    return 0;
/*     Error when writing into an output file */
L901:
    *ierr = 1;
    return 0;
} /* genmat_ */

/* *********************************************************************** */
/* Subroutine */ int tmin_(doublereal *tinit, doublereal *tmax, doublereal *
	tatm, doublereal *told, real *dt, real *dtmax, real *dmul, real *
	dmul2, real *dtmin, doublereal *tprint, doublereal *t, real *dtopt, 
	logical *topinf, logical *botinf, logical *lscreen, integer *itmin, 
	integer *itmax, integer *maxal, real *hcrits, integer *npd, logical *
	atmbc, integer *iver, logical *lprintd, integer *nprstep, doublereal *
	tprintint, logical *lenter, logical *ldayvar, logical *lsinprec, 
	logical *llai, real *rextinct, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern integer igetfileversion_(integer *, integer *);
    static integer i__, mpl, ivera;

    /* Fortran I/O blocks */
    static cilist io___230 = { 0, 6, 0, 0, 0 };
    static cilist io___231 = { 1, 30, 0, 0, 0 };
    static cilist io___232 = { 1, 30, 0, 0, 0 };
    static cilist io___233 = { 1, 30, 0, 0, 0 };
    static cilist io___235 = { 1, 30, 0, 0, 0 };
    static cilist io___236 = { 1, 30, 0, 0, 0 };
    static cilist io___237 = { 1, 30, 0, 0, 0 };
    static cilist io___238 = { 1, 30, 0, 0, 0 };
    static cilist io___239 = { 1, 30, 0, 0, 0 };
    static cilist io___240 = { 1, 30, 0, 0, 0 };
    static cilist io___243 = { 1, 31, 0, 0, 0 };
    static cilist io___244 = { 1, 31, 0, 0, 0 };
    static cilist io___245 = { 1, 31, 0, 0, 0 };
    static cilist io___246 = { 1, 31, 0, 0, 0 };
    static cilist io___247 = { 1, 31, 0, 0, 0 };
    static cilist io___248 = { 1, 31, 0, 0, 0 };
    static cilist io___249 = { 1, 31, 0, 0, 0 };
    static cilist io___250 = { 1, 31, 0, 0, 0 };
    static cilist io___251 = { 1, 31, 0, 0, 0 };
    static cilist io___252 = { 1, 31, 0, 0, 0 };


    /* Parameter adjustments */
    --tprint;

    /* Function Body */
    if (*lscreen) {
	s_wsle(&io___230);
	do_lio(&c__9, &c__1, "reading time information", (ftnlen)24);
	e_wsle();
    }
    i__1 = s_rsle(&io___231);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___232);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___233);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*dt), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*dtmin), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*dtmax), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*dmul), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*dmul2), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*itmin), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*itmax), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&mpl, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (mpl > *npd) {
	*ierr = 2;
	return 0;
    }
    i__1 = s_rsle(&io___235);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___236);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__5, &c__1, (char *)&(*tinit), (ftnlen)sizeof(doublereal))
	    ;
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__5, &c__1, (char *)&(*tmax), (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*iver > 2) {
	i__1 = s_rsle(&io___237);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___238);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lprintd), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*nprstep), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&(*tprintint), (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lenter), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    i__1 = s_rsle(&io___239);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___240);
    if (i__1 != 0) {
	goto L901;
    }
    i__2 = mpl;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_lio(&c__5, &c__1, (char *)&tprint[i__], (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L901;
	}
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    *dtopt = *dt;
    if (*topinf || *botinf || *atmbc) {
	ivera = igetfileversion_(&c__31, &c__1);
	i__1 = s_rsle(&io___243);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___244);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___245);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*maxal), (ftnlen)sizeof(integer)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	if (ivera == 4) {
	    i__1 = s_rsle(&io___246);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_rsle(&io___247);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__8, &c__1, (char *)&(*ldayvar), (ftnlen)sizeof(
		    logical));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__8, &c__1, (char *)&(*lsinprec), (ftnlen)sizeof(
		    logical));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__8, &c__1, (char *)&(*llai), (ftnlen)sizeof(
		    logical));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    if (*llai) {
		i__1 = s_rsle(&io___248);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_rsle();
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = s_rsle(&io___249);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*rextinct), (ftnlen)
			sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_rsle();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	}
	i__1 = s_rsle(&io___250);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___251);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*hcrits), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___252);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    } else {
	*tatm = *tmax;
    }
    tprint[mpl + 1] = *tmax;
    *told = *tinit;
    *t = *tinit + *dt;
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
} /* tmin_ */

/* *********************************************************************** */
/* Subroutine */ int meteoin_(real *latitude, real *altitude, real *
	shortwaverada, real *shortwaveradb, real *longwaverada, real *
	longwaveradb, real *longwaverada1, real *longwaveradb1, real *
	windheight, real *tempheight, integer *icrop, integer *ilai, real *
	cropheight, real *albedo, real *lai, real *xroot, integer *iinterc, 
	real *ainterc, integer *igrowth, real *rgrowth, real *rextinct, 
	integer *iradiation, logical *lenbal, logical *lprint, integer *
	isunsh, integer *irelhum, real *cloudf_ac__, real *cloudf_bc__, 
	logical *lhargr, logical *lmetdaily, real *xconv, integer *ierr)
{
    /* Format strings */
    static char fmt_110[] = "(/\002     Time      Short.Rad.    Long.Rad.   "
	    " Radiation    Sensible      Latent     Heat Flux     Balance\002/"
	    "\002      [d]       [MJ/m2/d]    [MJ/m2/d]    [MJ/m2/d]    [MJ/m"
	    "2/d]   [MJ/m2/d]    [MJ/m2/d]    [MJ/m2/d]\002/)";
    static char fmt_120[] = "(/\002   Time       ET        Evap    Transp   "
	    "  Rns        Rnl    RadTerm  AeroTerm     Prec     Interc   ExIn"
	    "terc\002/\002    [d]     [mm/d]    [mm/d]    [mm/d] [MJ/m2/d]  ["
	    "MJ/m2/d]   [mm/d]    [mm/d]    [mm/d]    [mm/d]    [mm/d]\002/)";
    static char fmt_130[] = "(/\002     Time      Short.Rad.    Long.Rad.   "
	    " Radiation    Sensible      Latent     Heat Flux     Balance    "
	    "    AirTemp\r      AirRH      SolarRad\002/\002      [d]       ["
	    "MJ/m2/d]    [MJ/m2/d]    [MJ/m2/d]    [MJ/m2/d]   [MJ/m2/d]    ["
	    "MJ/m2/d]    [MJ/m2/d]         [C]\r         [%]       [MJ/m2/d"
	    "]\002/)";
    static char fmt_140[] = "(/\002   Time       ET        Evap    Transp   "
	    "  Rns        Rnl    RadTerm  AeroTerm     Prec     Interc   ExIn"
	    "terc   AirTemp\r   AirRH   SolarRad\002/\002    [d]     [mm/d]  "
	    "  [mm/d]    [mm/d] [MJ/m2/d]  [MJ/m2/d]   [mm/d]    [mm/d]    [m"
	    "m/d]    [mm/d]    [mm/d]      [C]\r     [%]    [MJ/m2/d]\002/)";

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern integer igetfileversion_(integer *, integer *);
    static real maxalmet;
    static integer i__, j, iverm;

    /* Fortran I/O blocks */
    static cilist io___254 = { 1, 33, 0, 0, 0 };
    static cilist io___255 = { 1, 33, 0, 0, 0 };
    static cilist io___256 = { 1, 33, 0, 0, 0 };
    static cilist io___258 = { 1, 33, 0, 0, 0 };
    static cilist io___259 = { 1, 33, 0, 0, 0 };
    static cilist io___260 = { 1, 33, 0, 0, 0 };
    static cilist io___261 = { 1, 33, 0, 0, 0 };
    static cilist io___262 = { 1, 33, 0, 0, 0 };
    static cilist io___263 = { 1, 33, 0, 0, 0 };
    static cilist io___264 = { 1, 33, 0, 0, 0 };
    static cilist io___265 = { 1, 33, 0, 0, 0 };
    static cilist io___266 = { 1, 33, 0, 0, 0 };
    static cilist io___267 = { 1, 33, 0, 0, 0 };
    static cilist io___268 = { 1, 33, 0, 0, 0 };
    static cilist io___269 = { 1, 33, 0, 0, 0 };
    static cilist io___270 = { 1, 33, 0, 0, 0 };
    static cilist io___271 = { 1, 33, 0, 0, 0 };
    static cilist io___272 = { 1, 33, 0, 0, 0 };
    static cilist io___273 = { 1, 33, 0, 0, 0 };
    static cilist io___274 = { 1, 33, 0, 0, 0 };
    static cilist io___275 = { 1, 33, 0, 0, 0 };
    static cilist io___276 = { 1, 33, 0, 0, 0 };
    static cilist io___277 = { 1, 33, 0, 0, 0 };
    static cilist io___278 = { 1, 33, 0, 0, 0 };
    static cilist io___279 = { 1, 33, 0, 0, 0 };
    static cilist io___280 = { 1, 33, 0, 0, 0 };
    static cilist io___281 = { 1, 33, 0, 0, 0 };
    static cilist io___282 = { 0, 6, 0, 0, 0 };
    static cilist io___283 = { 0, 6, 0, 0, 0 };
    static cilist io___284 = { 0, 5, 0, 0, 0 };
    static cilist io___285 = { 1, 33, 0, 0, 0 };
    static cilist io___287 = { 1, 33, 0, 0, 0 };
    static cilist io___289 = { 1, 33, 0, 0, 0 };
    static cilist io___290 = { 1, 33, 0, 0, 0 };
    static cilist io___291 = { 1, 33, 0, 0, 0 };
    static cilist io___292 = { 1, 33, 0, 0, 0 };
    static cilist io___293 = { 1, 33, 0, 0, 0 };
    static cilist io___294 = { 1, 33, 0, 0, 0 };
    static cilist io___295 = { 1, 33, 0, 0, 0 };
    static cilist io___296 = { 0, 43, 0, fmt_110, 0 };
    static cilist io___297 = { 0, 43, 0, fmt_120, 0 };
    static cilist io___298 = { 0, 43, 0, fmt_130, 0 };
    static cilist io___299 = { 0, 43, 0, fmt_140, 0 };


/*      Latitude             ! Latitude of the location, [degree, N=+, S=-] */
/*      Altitude             ! Altitude of the location above mean see level [m] */
/*      ShortWaveRadA=0.25   ! fraction of extraterrestrial readiation on */
/*                           ! overcast days, first Angstrom coefficient, a_s, eq.55 */
/*      ShortWaveRadB=0.50   ! Input, fraction of extraterrestrial readiation on */
/*                           !   overcast days, second Angstrom coefficient, b_s, eq.55 */
/*      LAI                  ! Leaf area index [-] */
/*      LongWaveRadA=0.90    ! cloudiness factor, a_c, eq. 59 */
/*      LongWaveRadB=0.10    ! cloudiness factor, b_c, eq. 59 */
/*      LongWaveRadA1=0.34   ! emissivity correlation coefficient, a_l, eq. 60 */
/*      LongWaveRadB1=-0.139 ! emissivity correlation coefficient, b_l, eq. 60 */
/*      WindHeight=200.      ! Measurement hight of wind, [cm] */
/*      TempHeight=190.      ! Measurement hight of temperature, [cm] */
    /* Parameter adjustments */
    rgrowth -= 1001;

    /* Function Body */
    iverm = igetfileversion_(&c__33, &c__1);
    i__1 = s_rsle(&io___254);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___255);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___256);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&maxalmet, (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*iradiation), (ftnlen)sizeof(
	    integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__8, &c__1, (char *)&(*lhargr), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (iverm >= 4) {
	i__1 = s_rsle(&io___258);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___259);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lenbal), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lmetdaily), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    *altitude = 0.f;
    if (*iradiation != 2) {
	i__1 = s_rsle(&io___260);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___261);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*latitude), (ftnlen)sizeof(real)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*altitude), (ftnlen)sizeof(real)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___262);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___263);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*shortwaverada), (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*shortwaveradb), (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___264);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___265);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*longwaverada), (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*longwaveradb), (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___266);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___267);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*longwaverada1), (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*longwaveradb1), (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    i__1 = s_rsle(&io___268);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___269);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*windheight), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*tempheight), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___270);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___271);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*icrop), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*isunsh), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*irelhum), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*iradiation == 1 && *isunsh == 3) {
	i__1 = s_rsle(&io___272);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___273);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*cloudf_ac__), (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*cloudf_bc__), (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    if (*icrop >= 1) {
	i__1 = s_rsle(&io___274);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___275);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*ilai), (ftnlen)sizeof(integer))
		;
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*rextinct), (ftnlen)sizeof(real)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___276);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___277);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*iinterc), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	if (*icrop == 1) {
	    i__1 = s_rsle(&io___278);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_rsle(&io___279);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*cropheight), (ftnlen)
		    sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*albedo), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*lai), (ftnlen)sizeof(real))
		    ;
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*xroot), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    *cropheight = *cropheight * 100.f / *xconv;
/* conversion to cm */
	} else if (*icrop == 2) {
	    i__1 = s_rsle(&io___280);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_rsle(&io___281);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*igrowth), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    if (*igrowth > 1000) {
		s_wsle(&io___282);
		do_lio(&c__9, &c__1, "Number of crop growth data is larger t"
			"han 1000", (ftnlen)46);
		e_wsle();
		s_wsle(&io___283);
		do_lio(&c__9, &c__1, "Press Enter to continue", (ftnlen)23);
		e_wsle();
		s_rsle(&io___284);
		e_rsle();
		s_stop("", (ftnlen)0);
	    }
	    i__1 = s_rsle(&io___285);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = *igrowth;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = s_rsle(&io___287);
		if (i__2 != 0) {
		    goto L901;
		}
		for (j = 1; j <= 5; ++j) {
		    i__2 = do_lio(&c__4, &c__1, (char *)&rgrowth[i__ + j * 
			    1000], (ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L901;
		    }
		}
		i__2 = e_rsle();
		if (i__2 != 0) {
		    goto L901;
		}
		rgrowth[i__ + 2000] = rgrowth[i__ + 2000] * 100.f / *xconv;
/* conversion to cm */
/* L11: */
	    }
	}
	if (*iinterc == 1) {
	    i__1 = s_rsle(&io___289);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_rsle(&io___290);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*ainterc), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
    } else {
	i__1 = s_rsle(&io___291);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___292);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*albedo), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	*ilai = 0.f;
	*iinterc = 0.f;
    }
    i__1 = s_rsle(&io___293);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___294);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___295);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*lprint) {
	if (! (*lmetdaily)) {
	    if (*lenbal) {
		s_wsfe(&io___296);
		e_wsfe();
	    } else {
		s_wsfe(&io___297);
		e_wsfe();
	    }
	} else {
	    if (*lenbal) {
		s_wsfe(&io___298);
		e_wsfe();
	    } else {
		s_wsfe(&io___299);
		e_wsfe();
	    }
	}
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
} /* meteoin_ */

/* *********************************************************************** */
/* Subroutine */ int sinkin_(integer *nmat, logical *lchem, logical *lmosink, 
	logical *lsolred, logical *lsoladd, real *p0, real *poptm, real *p2h, 
	real *p2l, real *p3, real *r2h, real *r2l, real *aosm, real *c50, 
	real *p3c, integer *ns, logical *lmssink, real *crootmax, integer *
	iver, real *omegac, logical *lactrsu, real *omegas, real *spot, real *
	rkm, real *cmin, logical *lomegaw, logical *lscreen, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, imosink, imssink;

    /* Fortran I/O blocks */
    static cilist io___300 = { 0, 6, 0, 0, 0 };
    static cilist io___301 = { 1, 30, 0, 0, 0 };
    static cilist io___302 = { 1, 30, 0, 0, 0 };
    static cilist io___303 = { 1, 30, 0, 0, 0 };
    static cilist io___306 = { 1, 30, 0, 0, 0 };
    static cilist io___307 = { 1, 30, 0, 0, 0 };
    static cilist io___308 = { 1, 30, 0, 0, 0 };
    static cilist io___309 = { 1, 30, 0, 0, 0 };
    static cilist io___310 = { 1, 30, 0, 0, 0 };
    static cilist io___311 = { 1, 30, 0, 0, 0 };
    static cilist io___312 = { 1, 30, 0, 0, 0 };
    static cilist io___313 = { 1, 30, 0, 0, 0 };
    static cilist io___314 = { 1, 30, 0, 0, 0 };
    static cilist io___315 = { 1, 30, 0, 0, 0 };
    static cilist io___316 = { 1, 30, 0, 0, 0 };
    static cilist io___317 = { 1, 30, 0, 0, 0 };
    static cilist io___318 = { 1, 30, 0, 0, 0 };
    static cilist io___320 = { 1, 30, 0, 0, 0 };
    static cilist io___321 = { 1, 30, 0, 0, 0 };


    /* Parameter adjustments */
    --poptm;
    --crootmax;
    --aosm;

    /* Function Body */
    if (*lscreen) {
	s_wsle(&io___300);
	do_lio(&c__9, &c__1, "reading sink information", (ftnlen)24);
	e_wsle();
    }
    i__1 = s_rsle(&io___301);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___302);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*iver <= 2) {
	i__1 = s_rsle(&io___303);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&imosink, (ftnlen)sizeof(integer))
		;
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = do_lio(&c__4, &c__1, (char *)&crootmax[i__], (ftnlen)
		    sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    } else {
	i__1 = s_rsle(&io___306);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&imosink, (ftnlen)sizeof(integer))
		;
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = do_lio(&c__4, &c__1, (char *)&crootmax[i__], (ftnlen)
		    sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*omegac), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    if (imosink == 0) {
	*lmosink = TRUE_;
    } else {
	*lmosink = FALSE_;
    }
    i__1 = s_rsle(&io___307);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*lmosink) {
	i__1 = s_rsle(&io___308);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*p0), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*p2h), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*p2l), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*p3), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*r2h), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*r2l), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___309);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___310);
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *nmat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = do_lio(&c__4, &c__1, (char *)&poptm[i__], (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	*p0 = -dabs(*p0);
	*p2l = -dabs(*p2l);
	*p2h = -dabs(*p2h);
	*p3 = -dabs(*p3);
    } else {
	i__1 = s_rsle(&io___311);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*p0), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*p3), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    if (*lchem) {
	i__1 = s_rsle(&io___312);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___313);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lsolred), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	if (*lsolred) {
	    i__1 = s_rsle(&io___314);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_rsle(&io___315);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__8, &c__1, (char *)&(*lsoladd), (ftnlen)sizeof(
		    logical));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_rsle(&io___316);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    if (*lsoladd) {
		i__1 = s_rsle(&io___317);
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *ns;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__1 = do_lio(&c__4, &c__1, (char *)&aosm[i__], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = e_rsle();
		if (i__1 != 0) {
		    goto L901;
		}
	    } else {
		i__1 = s_rsle(&io___318);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*c50), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*p3c), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *ns;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__1 = do_lio(&c__4, &c__1, (char *)&aosm[i__], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&imssink, (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_rsle();
		if (i__1 != 0) {
		    goto L901;
		}
		if (imssink == 0) {
		    *lmssink = FALSE_;
		} else {
		    *lmssink = TRUE_;
		}
	    }
	}
	if (*ns > 1) {
	    *lactrsu = FALSE_;
	}
/* disable for uptake on last solu */
	if (*lactrsu && *ns == 1) {
/*        if(lActRSU) then              ! Active uptake only for the last solute */
/* disable for uptake on last solu */
	    i__1 = s_rsle(&io___320);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_rsle(&io___321);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*omegas), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*spot), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rkm), (ftnlen)sizeof(real))
		    ;
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*cmin), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__8, &c__1, (char *)&(*lomegaw), (ftnlen)sizeof(
		    logical));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
} /* sinkin_ */

/* *********************************************************************** */
/* Subroutine */ int rootin_(real *trmin, real *trharv, real *xrmin, real *
	xrmax, real *rgr, logical *lscreen, integer *iver, integer *irootin, 
	integer *ngrowth, real *rgrowth, real *trperiod, integer *ierr)
{
    /* Format strings */
    static char fmt_110[] = "(/\002 Root growth information\002/1x,23(\002"
	    "=\002)/\002 tRMin = \002,f10.3,\002 tRHarv = \002,f10.3/\002 xRM"
	    "in = \002,f10.3,\002 xRMax = \002,f10.3/\002 Root growth rate ="
	    " \002,e11.3)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer i__;
    static real rtm;
    static integer irfak;
    static real trmed, xrmed;

    /* Fortran I/O blocks */
    static cilist io___322 = { 0, 6, 0, 0, 0 };
    static cilist io___323 = { 1, 30, 0, 0, 0 };
    static cilist io___324 = { 1, 30, 0, 0, 0 };
    static cilist io___325 = { 1, 30, 0, 0, 0 };
    static cilist io___326 = { 1, 30, 0, 0, 0 };
    static cilist io___327 = { 1, 30, 0, 0, 0 };
    static cilist io___328 = { 0, 6, 0, 0, 0 };
    static cilist io___329 = { 0, 6, 0, 0, 0 };
    static cilist io___330 = { 0, 5, 0, 0, 0 };
    static cilist io___331 = { 1, 30, 0, 0, 0 };
    static cilist io___333 = { 1, 30, 0, 0, 0 };
    static cilist io___334 = { 1, 30, 0, 0, 0 };
    static cilist io___335 = { 1, 30, 0, 0, 0 };
    static cilist io___339 = { 1, 30, 0, 0, 0 };
    static cilist io___341 = { 0, 6, 0, 0, 0 };
    static cilist io___342 = { 1, 50, 0, fmt_110, 0 };


    /* Parameter adjustments */
    rgrowth -= 1001;

    /* Function Body */
    *trperiod = 1e30f;
    if (*lscreen) {
	s_wsle(&io___322);
	do_lio(&c__9, &c__1, "reading of root growth information", (ftnlen)34)
		;
	e_wsle();
    }
    i__1 = s_rsle(&io___323);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*iver == 4) {
	i__1 = s_rsle(&io___324);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___325);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*irootin), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	if (*irootin == 1) {
	    i__1 = s_rsle(&io___326);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_rsle(&io___327);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*ngrowth), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    if (*ngrowth > 1000) {
		s_wsle(&io___328);
		do_lio(&c__9, &c__1, "Number of crop growth data is larger t"
			"han 1000", (ftnlen)46);
		e_wsle();
		s_wsle(&io___329);
		do_lio(&c__9, &c__1, "Press Enter to continue", (ftnlen)23);
		e_wsle();
		s_rsle(&io___330);
		e_rsle();
		s_stop("", (ftnlen)0);
	    }
	    i__1 = s_rsle(&io___331);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = *ngrowth;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = s_rsle(&io___333);
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = do_lio(&c__4, &c__1, (char *)&rgrowth[i__ + 1000], (
			ftnlen)sizeof(real));
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = do_lio(&c__4, &c__1, (char *)&rgrowth[i__ + 5000], (
			ftnlen)sizeof(real));
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = e_rsle();
		if (i__2 != 0) {
		    goto L901;
		}
/* L11: */
	    }
	}
    }
    if (*iver < 4 || *irootin == 2) {
	*irootin = 2;
	i__1 = s_rsle(&io___334);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	if (*iver < 4) {
	    i__1 = s_rsle(&io___335);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&irfak, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*trmin), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&trmed, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*trharv), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*xrmin), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&xrmed, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*xrmax), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	} else {
	    i__1 = s_rsle(&io___339);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&irfak, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*trmin), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&trmed, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*trharv), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*xrmin), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&xrmed, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*xrmax), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*trperiod), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	if (irfak == 1) {
	    trmed = (*trharv + *trmin) / 2.f;
	    xrmed = (*xrmax + *xrmin) / 2.f;
	}
	rtm = trmed - *trmin;
	if (rtm < 1e-20f || xrmed < 1e-10f) {
	    s_wsle(&io___341);
	    do_lio(&c__9, &c__1, "Time(depth) Root Data must be larger then "
		    "Initial Root Growth Time(depth) !!", (ftnlen)76);
	    e_wsle();
	    goto L901;
	}
/* Computing MAX */
	r__1 = 1e-4f, r__2 = *xrmin * (*xrmax - xrmed);
	*rgr = -(1.f / rtm) * log(dmax(r__1,r__2) / (xrmed * (*xrmax - *xrmin)
		));
	i__1 = s_wsfe(&io___342);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_fio(&c__1, (char *)&(*trmin), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_fio(&c__1, (char *)&(*trharv), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_fio(&c__1, (char *)&(*xrmin), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_fio(&c__1, (char *)&(*xrmax), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = do_fio(&c__1, (char *)&(*rgr), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
/*     Error when writing into an output file */
L902:
    *ierr = 2;
    return 0;
} /* rootin_ */

/* *********************************************************************** */
/* Subroutine */ int tempin_(integer *nmat, real *tpar, real *ampl, real *
	tperiod, integer *ktopt, real *ttop, integer *kbott, real *tbot, 
	logical *topinf, logical *botinf, integer *icampbell, integer *iver, 
	real *snowmf, logical *lscreen, integer *ierr)
{
    /* Format strings */
    static char fmt_110[] = "(//\002 Heat transport information\002/1x,26"
	    "(\002=\002)//\002 ample = \002,f10.3//\002   Beta    Qn     Qo  "
	    "      B1         B2         B3\r        Cn         Co         C"
	    "w\002)";
    static char fmt_120[] = "(3f7.3,6e11.3)";

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static real tb, tt;

    /* Fortran I/O blocks */
    static cilist io___343 = { 0, 6, 0, 0, 0 };
    static cilist io___344 = { 1, 30, 0, 0, 0 };
    static cilist io___345 = { 1, 30, 0, 0, 0 };
    static cilist io___347 = { 1, 30, 0, 0, 0 };
    static cilist io___349 = { 1, 30, 0, 0, 0 };
    static cilist io___350 = { 1, 30, 0, 0, 0 };
    static cilist io___351 = { 1, 30, 0, 0, 0 };
    static cilist io___352 = { 1, 50, 0, fmt_110, 0 };
    static cilist io___353 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___354 = { 1, 30, 0, 0, 0 };
    static cilist io___355 = { 1, 30, 0, 0, 0 };


    /* Parameter adjustments */
    tpar -= 11;

    /* Function Body */
    if (*lscreen) {
	s_wsle(&io___343);
	do_lio(&c__9, &c__1, "reading heat transport information", (ftnlen)34)
		;
	e_wsle();
    }
    i__1 = s_rsle(&io___344);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___345);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = *nmat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = s_rsle(&io___347);
	if (i__2 != 0) {
	    goto L901;
	}
	for (j = 1; j <= 9; ++j) {
	    i__2 = do_lio(&c__4, &c__1, (char *)&tpar[j + i__ * 10], (ftnlen)
		    sizeof(real));
	    if (i__2 != 0) {
		goto L901;
	    }
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L901;
	}
/* L11: */
    }
    i__1 = s_rsle(&io___349);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*iver <= 2) {
	i__1 = s_rsle(&io___350);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*ampl), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*tperiod), (ftnlen)sizeof(real))
		;
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	*icampbell = 0;
    } else {
	i__1 = s_rsle(&io___351);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*ampl), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*tperiod), (ftnlen)sizeof(real))
		;
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*icampbell), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*snowmf), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    i__1 = s_wsfe(&io___352);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*ampl), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = *nmat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = s_wsfe(&io___353);
	if (i__2 != 0) {
	    goto L902;
	}
	for (j = 1; j <= 9; ++j) {
	    i__2 = do_fio(&c__1, (char *)&tpar[j + i__ * 10], (ftnlen)sizeof(
		    real));
	    if (i__2 != 0) {
		goto L902;
	    }
	}
	i__2 = e_wsfe();
	if (i__2 != 0) {
	    goto L902;
	}
/* L12: */
    }
    i__1 = s_rsle(&io___354);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___355);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ktopt), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&tt, (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*kbott), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&tb, (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (! (*topinf)) {
	*ttop = tt;
    }
    if (! (*botinf)) {
	*tbot = tb;
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
/*     Error when writing into an output file */
L902:
    *ierr = 2;
    return 0;
} /* tempin_ */

/* *********************************************************************** */
/* Subroutine */ int profil_(integer *n, integer *nmat, real *x, integer *
	matnum, real *xsurf, real *beta, real *ah, real *ak, real *ath, real *
	thr, real *ths, real *cons, real *hs, logical *lscreen, integer *ierr)
{
    /* Format strings */
    static char fmt_110[] = "(//\002    n      depth     THr       THs      "
	    " hs       Ks\002,\002        Ks/KsTop     Beta      Ah        AK"
	    "        ATh\002/)";
    static char fmt_120[] = "(i5,f10.2,2f10.3,f10.1,e12.3,5f10.3)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static integer i__, m;
    static real consn;

    /* Fortran I/O blocks */
    static cilist io___358 = { 0, 6, 0, 0, 0 };
    static cilist io___359 = { 1, 78, 0, fmt_110, 0 };
    static cilist io___363 = { 1, 78, 0, fmt_120, 0 };
    static cilist io___364 = { 1, 78, 0, "('end')", 0 };


    /* Parameter adjustments */
    --ath;
    --ak;
    --ah;
    --beta;
    --matnum;
    --x;
    --hs;
    --cons;
    --ths;
    --thr;

    /* Function Body */
    if (*lscreen) {
	s_wsle(&io___358);
	do_lio(&c__9, &c__1, "printing profile information", (ftnlen)28);
	e_wsle();
    }
    i__1 = s_wsfe(&io___359);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L901;
    }
    consn = cons[matnum[*n]] * ak[*n];
    for (i__ = *n; i__ >= 1; --i__) {
	m = matnum[i__];
	i__1 = s_wsfe(&io___363);
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *n - i__ + 1;
	i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	r__1 = *xsurf - x[i__];
	i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&thr[m], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	r__2 = thr[m] + ath[i__] * (ths[m] - thr[m]);
	i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	r__3 = hs[m] * ah[i__];
	i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	r__4 = cons[m] * ak[i__];
	i__1 = do_fio(&c__1, (char *)&r__4, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	r__5 = cons[m] * ak[i__] / consn;
	i__1 = do_fio(&c__1, (char *)&r__5, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&beta[i__], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&ah[i__], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&ak[i__], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&ath[i__], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
/* L11: */
    }
    i__1 = s_wsfe(&io___364);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L901;
    }
    return 0;
/*     Error when writing into an output file */
L901:
    *ierr = 1;
    return 0;
} /* profil_ */

/* *********************************************************************** */
/*     Read information about solute transport */
/* Subroutine */ int chemin_(logical *lupw, logical *ltdep, integer *nmat, 
	integer *ns, integer *nsd, integer *maxitc, real *chpar, real *tdep, 
	integer *ktopch, real *ctop, integer *kbotch, real *cbot, real *epsi, 
	real *tpulse, real *cumch, real *ctola, real *ctolr, logical *llinear,
	 logical *lequil, logical *lartd, real *pecr, logical *lscreen, real *
	dsurf, real *catm, logical *ltort, logical *lmobim, logical *lbact, 
	logical *lfiltr, integer *imoistdep, real *wdep, integer *nmatd, 
	integer *imodel, real *par, integer *iver, logical *ldualneq, logical 
	*lmassini, logical *leqinit, integer *itort, integer *ierr)
{
    /* Format strings */
    static char fmt_110[] = "(//\002 Solute transport information\002/1x,28"
	    "(\002=\002))";
    static char fmt_120[] = "(/\002 Upstream weighting finite-element metho"
	    "d\002)";
    static char fmt_130[] = "(/\002 Galerkin finite-element method\002)";
    static char fmt_140[] = "(/\002 Artificial dispersion is added when Pecl"
	    "et number is\002,\002 higher than\002,f10.3)";
    static char fmt_150[] = "(//\002 lTDep     lWDep     cTolA     cTolR   M"
	    "axItC\002/l3,6x,i3,e13.3,f10.4,i7///\002 Mat.     Bulk.D.    Dis"
	    "pL    Fraction  Immobile WC\002)";
    static char fmt_160[] = "(i3,f13.4,3f10.4)";
    static char fmt_170[] = "(/\002    Dif.w.      Dif.g.   \002,50(\002-"
	    "\002),\002 (\002,i2,\002.solute)\002)";
    static char fmt_180[] = "(2e12.4/\002 Mat.     KS         Nu         Bet"
	    "a      Henry\r     SinkL1     SinkS1     SinkG1     SinkL1`    S"
	    "inkS1`    SinkG1`\r   SinkL0     SinkS0     SinkG0      Alfa\002)"
	    ;
    static char fmt_190[] = "(i4,14e11.4)";
    static char fmt_200[] = "(/\002 Langmuir nonlinear adsorption isotherm f"
	    "or material \002,i2)";
    static char fmt_210[] = "(/\002 Freundlich nonlinear adsorption isotherm"
	    " for material \002,i2)";
    static char fmt_220[] = "(/\002 No adsorption or linear adsorp. isotherm"
	    " for material \002,i2)";
    static char fmt_222[] = "(/\002 Physical non-equilibrium solute transpor"
	    "t with mobile and imobile water.\002)";
    static char fmt_224[] = "(/\002 Chemical non-equilibrium solute transpor"
	    "t with kinetic and equilibrium sorption sites.\002)";
    static char fmt_230[] = "(/\002 kTopCh      cTop(1...NS)\002/i4,7x,20e10"
	    ".3)";
    static char fmt_240[] = "(/\002 kBotCh      cBot(1...NS)\002/i4,7x,20e10"
	    ".3)";
    static char fmt_250[] = "(/\002 tPulse =   \002,f15.3)";

    /* System generated locals */
    integer chpar_dim1, chpar_offset, wdep_dim1, wdep_offset, i__1, i__2, 
	    i__3;
    real r__1, r__2;

    /* Local variables */
    static integer inonequl, i__, j, m;
    static logical lmoistdep;
    static integer jj;
    extern doublereal fq_(integer *, real *, real *);
    static integer jjj;
    static logical lvar;
    static integer npar2, ibact;

    /* Fortran I/O blocks */
    static cilist io___365 = { 0, 6, 0, 0, 0 };
    static cilist io___366 = { 1, 50, 0, fmt_110, 0 };
    static cilist io___367 = { 1, 30, 0, 0, 0 };
    static cilist io___368 = { 1, 30, 0, 0, 0 };
    static cilist io___369 = { 1, 30, 0, 0, 0 };
    static cilist io___370 = { 1, 30, 0, 0, 0 };
    static cilist io___372 = { 1, 30, 0, 0, 0 };
    static cilist io___373 = { 1, 30, 0, 0, 0 };
    static cilist io___377 = { 1, 50, 0, fmt_120, 0 };
    static cilist io___378 = { 1, 50, 0, fmt_130, 0 };
    static cilist io___379 = { 1, 50, 0, fmt_140, 0 };
    static cilist io___380 = { 1, 50, 0, fmt_150, 0 };
    static cilist io___381 = { 1, 30, 0, 0, 0 };
    static cilist io___383 = { 1, 30, 0, 0, 0 };
    static cilist io___385 = { 1, 50, 0, fmt_160, 0 };
    static cilist io___388 = { 1, 50, 0, fmt_170, 0 };
    static cilist io___389 = { 1, 30, 0, 0, 0 };
    static cilist io___390 = { 1, 30, 0, 0, 0 };
    static cilist io___391 = { 1, 50, 0, fmt_180, 0 };
    static cilist io___392 = { 1, 30, 0, 0, 0 };
    static cilist io___393 = { 1, 30, 0, 0, 0 };
    static cilist io___394 = { 1, 50, 0, fmt_190, 0 };
    static cilist io___395 = { 1, 50, 0, fmt_200, 0 };
    static cilist io___396 = { 1, 50, 0, fmt_210, 0 };
    static cilist io___397 = { 1, 50, 0, fmt_220, 0 };
    static cilist io___398 = { 1, 50, 0, fmt_222, 0 };
    static cilist io___399 = { 1, 50, 0, fmt_224, 0 };
    static cilist io___401 = { 1, 30, 0, 0, 0 };
    static cilist io___402 = { 1, 30, 0, 0, 0 };
    static cilist io___403 = { 1, 30, 0, 0, 0 };
    static cilist io___404 = { 1, 30, 0, 0, 0 };
    static cilist io___405 = { 1, 30, 0, 0, 0 };
    static cilist io___406 = { 1, 30, 0, 0, 0 };
    static cilist io___407 = { 1, 30, 0, 0, 0 };
    static cilist io___408 = { 1, 30, 0, 0, 0 };
    static cilist io___410 = { 1, 30, 0, 0, 0 };
    static cilist io___411 = { 1, 30, 0, 0, 0 };
    static cilist io___412 = { 1, 30, 0, 0, 0 };
    static cilist io___413 = { 1, 30, 0, 0, 0 };
    static cilist io___414 = { 1, 30, 0, 0, 0 };
    static cilist io___415 = { 1, 30, 0, 0, 0 };
    static cilist io___416 = { 1, 30, 0, 0, 0 };
    static cilist io___417 = { 1, 50, 0, fmt_230, 0 };
    static cilist io___418 = { 1, 50, 0, fmt_240, 0 };
    static cilist io___419 = { 1, 30, 0, 0, 0 };
    static cilist io___420 = { 1, 30, 0, 0, 0 };
    static cilist io___421 = { 1, 50, 0, fmt_250, 0 };


    /* Parameter adjustments */
    par -= 12;
    --lmobim;
    --llinear;
    cumch -= 11;
    --cbot;
    --ctop;
    --tdep;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;
    wdep_dim1 = 2 + *nmatd;
    wdep_offset = 1 + wdep_dim1;
    wdep -= wdep_offset;

    /* Function Body */
    if (*lscreen) {
	s_wsle(&io___365);
	do_lio(&c__9, &c__1, "reading solute transport information", (ftnlen)
		36);
	e_wsle();
    }
    i__1 = s_wsfe(&io___366);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_rsle(&io___367);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___368);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*iver <= 2) {
	i__1 = s_rsle(&io___369);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*epsi), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lupw), (ftnlen)sizeof(logical))
		;
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lartd), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*ltdep), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*ctola), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*ctolr), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*maxitc), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*pecr), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*ltort), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    } else {
	i__1 = s_rsle(&io___370);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*epsi), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lupw), (ftnlen)sizeof(logical))
		;
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lartd), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*ltdep), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*ctola), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*ctolr), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*maxitc), (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*pecr), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&(*ns), (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*ltort), (ftnlen)sizeof(logical)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ibact, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lfiltr), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	*lbact = FALSE_;
	if (ibact == 1) {
	    *lbact = TRUE_;
	}
    }
    if (*iver == 4) {
	i__1 = s_rsle(&io___372);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___373);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&inonequl, (ftnlen)sizeof(integer)
		);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&lmoistdep, (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*ldualneq), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*lmassini), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&(*leqinit), (ftnlen)sizeof(
		logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__8, &c__1, (char *)&lvar, (ftnlen)sizeof(logical));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	if (lmoistdep) {
	    *imoistdep = 1;
	}
	if (lvar) {
	    *itort = 1;
	}
    }
    *pecr = dmax(*pecr,.1f);
    if (*lupw) {
	i__1 = s_wsfe(&io___377);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
    } else {
	i__1 = s_wsfe(&io___378);
	if (i__1 != 0) {
	    goto L902;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L902;
	}
	if (*lartd) {
	    i__1 = s_wsfe(&io___379);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*pecr), (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	}
    }
    i__1 = s_wsfe(&io___380);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*ltdep), (ftnlen)sizeof(logical));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*imoistdep), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*ctola), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*ctolr), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*maxitc), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_rsle(&io___381);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    *lequil = TRUE_;
    i__1 = *nmat;
    for (m = 1; m <= i__1; ++m) {
	i__2 = s_rsle(&io___383);
	if (i__2 != 0) {
	    goto L901;
	}
	for (j = 1; j <= 4; ++j) {
	    i__2 = do_lio(&c__4, &c__1, (char *)&chpar[j + m * chpar_dim1], (
		    ftnlen)sizeof(real));
	    if (i__2 != 0) {
		goto L901;
	    }
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L901;
	}
	i__2 = s_wsfe(&io___385);
	if (i__2 != 0) {
	    goto L902;
	}
	i__2 = do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L902;
	}
	for (j = 1; j <= 4; ++j) {
	    i__2 = do_fio(&c__1, (char *)&chpar[j + m * chpar_dim1], (ftnlen)
		    sizeof(real));
	    if (i__2 != 0) {
		goto L902;
	    }
	}
	i__2 = e_wsfe();
	if (i__2 != 0) {
	    goto L902;
	}
	if (chpar[m * chpar_dim1 + 3] < 1.f || chpar[m * chpar_dim1 + 4] > 
		0.f || *lbact) {
	    *lequil = FALSE_;
	}
	lmobim[m] = FALSE_;
	if (! (*lbact) && chpar[m * chpar_dim1 + 4] > 0.f) {
	    lmobim[m] = TRUE_;
	}
	if (! (*lequil) && chpar[m * chpar_dim1 + 1] == 0.f) {
	    goto L903;
	}
/* L11: */
    }
    i__1 = *ns;
    for (jj = 1; jj <= i__1; ++jj) {
	jjj = jj - 1 << 4;
	i__2 = s_wsfe(&io___388);
	if (i__2 != 0) {
	    goto L902;
	}
	i__2 = do_fio(&c__1, (char *)&jj, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L902;
	}
	i__2 = e_wsfe();
	if (i__2 != 0) {
	    goto L902;
	}
	i__2 = s_rsle(&io___389);
	if (i__2 != 0) {
	    goto L901;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L901;
	}
	i__2 = s_rsle(&io___390);
	if (i__2 != 0) {
	    goto L901;
	}
	for (j = 5; j <= 6; ++j) {
	    i__2 = do_lio(&c__4, &c__1, (char *)&chpar[jjj + j + chpar_dim1], 
		    (ftnlen)sizeof(real));
	    if (i__2 != 0) {
		goto L901;
	    }
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L901;
	}
	i__2 = s_wsfe(&io___391);
	if (i__2 != 0) {
	    goto L902;
	}
	for (j = 5; j <= 6; ++j) {
	    i__2 = do_fio(&c__1, (char *)&chpar[jjj + j + chpar_dim1], (
		    ftnlen)sizeof(real));
	    if (i__2 != 0) {
		goto L902;
	    }
	}
	i__2 = e_wsfe();
	if (i__2 != 0) {
	    goto L902;
	}
	i__2 = s_rsle(&io___392);
	if (i__2 != 0) {
	    goto L901;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L901;
	}
	llinear[jj] = TRUE_;
	i__2 = *nmat;
	for (m = 1; m <= i__2; ++m) {
	    chpar[jjj + 5 + m * chpar_dim1] = chpar[jjj + 5 + chpar_dim1];
	    chpar[jjj + 6 + m * chpar_dim1] = chpar[jjj + 6 + chpar_dim1];
	    i__3 = s_rsle(&io___393);
	    if (i__3 != 0) {
		goto L901;
	    }
	    for (j = 7; j <= 20; ++j) {
		i__3 = do_lio(&c__4, &c__1, (char *)&chpar[jjj + j + m * 
			chpar_dim1], (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L901;
		}
	    }
	    i__3 = e_rsle();
	    if (i__3 != 0) {
		goto L901;
	    }
	    i__3 = s_wsfe(&io___394);
	    if (i__3 != 0) {
		goto L902;
	    }
	    i__3 = do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	    if (i__3 != 0) {
		goto L902;
	    }
	    for (j = 7; j <= 20; ++j) {
		i__3 = do_fio(&c__1, (char *)&chpar[jjj + j + m * chpar_dim1],
			 (ftnlen)sizeof(real));
		if (i__3 != 0) {
		    goto L902;
		}
	    }
	    i__3 = e_wsfe();
	    if (i__3 != 0) {
		goto L902;
	    }
	    if ((r__1 = chpar[jjj + 8 + m * chpar_dim1] - 0.f, dabs(r__1)) > 
		    1e-12f) {
		i__3 = s_wsfe(&io___395);
		if (i__3 != 0) {
		    goto L902;
		}
		i__3 = do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
		if (i__3 != 0) {
		    goto L902;
		}
		i__3 = e_wsfe();
		if (i__3 != 0) {
		    goto L902;
		}
	    } else if ((r__1 = chpar[jjj + 9 + m * chpar_dim1] - 1.f, dabs(
		    r__1)) > .001f) {
		i__3 = s_wsfe(&io___396);
		if (i__3 != 0) {
		    goto L902;
		}
		i__3 = do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
		if (i__3 != 0) {
		    goto L902;
		}
		i__3 = e_wsfe();
		if (i__3 != 0) {
		    goto L902;
		}
	    } else {
		i__3 = s_wsfe(&io___397);
		if (i__3 != 0) {
		    goto L902;
		}
		i__3 = do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
		if (i__3 != 0) {
		    goto L902;
		}
		i__3 = e_wsfe();
		if (i__3 != 0) {
		    goto L902;
		}
	    }
	    if (! (*lequil)) {
		if (lmobim[m]) {
		    i__3 = s_wsfe(&io___398);
		    if (i__3 != 0) {
			goto L902;
		    }
		    i__3 = e_wsfe();
		    if (i__3 != 0) {
			goto L902;
		    }
		} else {
		    i__3 = s_wsfe(&io___399);
		    if (i__3 != 0) {
			goto L902;
		    }
		    i__3 = e_wsfe();
		    if (i__3 != 0) {
			goto L902;
		    }
		}
	    }
	    if ((r__1 = chpar[jjj + 8 + m * chpar_dim1] - 0.f, dabs(r__1)) > 
		    1e-12f || (r__2 = chpar[jjj + 9 + m * chpar_dim1] - 1.f, 
		    dabs(r__2)) > .001f) {
		llinear[jj] = FALSE_;
	    }
	    if (*lbact && (chpar[jjj + 18 + m * chpar_dim1] > 0.f || chpar[
		    jjj + 15 + m * chpar_dim1] > 0.f)) {
		llinear[jj] = FALSE_;
	    }
/* L12: */
	}
/* L13: */
    }
    i__1 = (*ns << 4) + 4;
    for (jj = 1; jj <= i__1; ++jj) {
	tdep[jj] = 0.f;
	if (jj <= *ns * 9) {
	    wdep[jj * wdep_dim1 + 1] = 1.f;
	}
	if (jj <= *ns * 9) {
	    wdep[jj * wdep_dim1 + 2] = 0.f;
	}
/* L14: */
    }
    i__1 = *ns;
    for (jj = 1; jj <= i__1; ++jj) {
	for (i__ = 1; i__ <= 10; ++i__) {
	    cumch[i__ + jj * 10] = 0.f;
/* L15: */
	}
	if (*ltdep) {
	    jjj = jj - 1 << 4;
	    if (jj == 1) {
		i__2 = s_rsle(&io___401);
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = e_rsle();
		if (i__2 != 0) {
		    goto L901;
		}
	    }
	    i__2 = s_rsle(&io___402);
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = s_rsle(&io___403);
	    if (i__2 != 0) {
		goto L901;
	    }
	    for (j = 5; j <= 6; ++j) {
		i__2 = do_lio(&c__4, &c__1, (char *)&tdep[jjj + j], (ftnlen)
			sizeof(real));
		if (i__2 != 0) {
		    goto L901;
		}
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = s_rsle(&io___404);
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = s_rsle(&io___405);
	    if (i__2 != 0) {
		goto L901;
	    }
	    for (j = 7; j <= 20; ++j) {
		i__2 = do_lio(&c__4, &c__1, (char *)&tdep[jjj + j], (ftnlen)
			sizeof(real));
		if (i__2 != 0) {
		    goto L901;
		}
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	}
/* L16: */
    }
    i__1 = *ns;
    for (jj = 1; jj <= i__1; ++jj) {
	if (*imoistdep == 1) {
	    if (jj == 1) {
		i__2 = s_rsle(&io___406);
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = e_rsle();
		if (i__2 != 0) {
		    goto L901;
		}
	    }
	    i__2 = s_rsle(&io___407);
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = s_rsle(&io___408);
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&npar2, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    jjj = (jj - 1) * 9;
	    i__2 = s_rsle(&io___410);
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = s_rsle(&io___411);
	    if (i__2 != 0) {
		goto L901;
	    }
	    for (j = 1; j <= 9; ++j) {
		i__2 = do_lio(&c__4, &c__1, (char *)&wdep[(jjj + j) * 
			wdep_dim1 + 1], (ftnlen)sizeof(real));
		if (i__2 != 0) {
		    goto L901;
		}
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = s_rsle(&io___412);
	    if (i__2 != 0) {
		goto L901;
	    }
	    for (j = 1; j <= 9; ++j) {
		i__2 = do_lio(&c__4, &c__1, (char *)&wdep[(jjj + j) * 
			wdep_dim1 + 2], (ftnlen)sizeof(real));
		if (i__2 != 0) {
		    goto L901;
		}
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L901;
	    }
	    i__2 = *nmat;
	    for (m = 1; m <= i__2; ++m) {
		for (j = 1; j <= 9; ++j) {
		    wdep[m + 2 + (jjj + j) * wdep_dim1] = fq_(imodel, &wdep[(
			    jjj + j) * wdep_dim1 + 2], &par[m * 11 + 1]);
/* L17: */
		}
/* L18: */
	    }
	}
/* L19: */
    }
    i__1 = s_rsle(&io___413);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___414);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*ktopch), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__2 = *ns;
    for (jj = 1; jj <= i__2; ++jj) {
	i__1 = do_lio(&c__4, &c__1, (char *)&ctop[jj], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&(*kbotch), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L901;
    }
    i__3 = *ns;
    for (jj = 1; jj <= i__3; ++jj) {
	i__1 = do_lio(&c__4, &c__1, (char *)&cbot[jj], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    if (*ktopch == -2) {
	i__1 = s_rsle(&io___415);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_rsle(&io___416);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*dsurf), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&(*catm), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    i__1 = s_wsfe(&io___417);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*ktopch), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L902;
    }
    i__2 = *ns;
    for (jj = 1; jj <= i__2; ++jj) {
	i__1 = do_fio(&c__1, (char *)&ctop[jj], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_wsfe(&io___418);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*kbotch), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L902;
    }
    i__2 = *ns;
    for (jj = 1; jj <= i__2; ++jj) {
	i__1 = do_fio(&c__1, (char *)&cbot[jj], (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L902;
	}
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = s_rsle(&io___419);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_rsle(&io___420);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&(*tpulse), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = s_wsfe(&io___421);
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = do_fio(&c__1, (char *)&(*tpulse), (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L902;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L902;
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
/*     Error when writing into an output file */
L902:
    *ierr = 2;
    return 0;
/*     Bulk Density is equal to zero */
L903:
    *ierr = 3;
    return 0;
} /* chemin_ */

/* *********************************************************************** */
/* Subroutine */ int opensolutefiles_(integer *ns, char *cdatapath, integer *
	ilengthpath, char *cfilename, integer *ierr, ftnlen cdatapath_len, 
	ftnlen cfilename_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2], i__3;
    olist o__1;

    /* Local variables */
    static integer i__;
    static char ch1[1], ch2[1], cname[12], cname1[13];

    /* Fortran I/O blocks */
    static icilist io___424 = { 0, ch1, 0, "(i1)", 1, 1 };
    static icilist io___426 = { 0, ch1, 0, "(i1)", 1, 1 };
    static icilist io___428 = { 0, ch2, 0, "(i1)", 1, 1 };


    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ <= 9) {
	    s_wsfi(&io___424);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_copy(cname, "solutex.out", (ftnlen)12, (ftnlen)11);
	    *(unsigned char *)&cname[7] = *(unsigned char *)ch1;
/* Writing concatenation */
	    i__2[0] = *ilengthpath, a__1[0] = cdatapath;
	    i__2[1] = 12, a__1[1] = cname;
	    s_cat(cfilename, a__1, i__2, &c__2, (ftnlen)200);
	} else {
	    s_wsfi(&io___426);
	    do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_wsfi(&io___428);
	    i__3 = i__ - 10;
	    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
	    e_wsfi();
	    s_copy(cname1, "solutexx.out", (ftnlen)13, (ftnlen)12);
	    *(unsigned char *)&cname1[7] = *(unsigned char *)ch1;
	    *(unsigned char *)&cname1[8] = *(unsigned char *)ch2;
/* Writing concatenation */
	    i__2[0] = *ilengthpath, a__1[0] = cdatapath;
	    i__2[1] = 13, a__1[1] = cname1;
	    s_cat(cfilename, a__1, i__2, &c__2, (ftnlen)200);
	}
	o__1.oerr = 1;
	o__1.ounit = i__ + 80;
	o__1.ofnmlen = 200;
	o__1.ofnm = cfilename;
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__3 = f_open(&o__1);
	if (i__3 != 0) {
	    goto L901;
	}
/* L11: */
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
} /* opensolutefiles_ */

/* ####################################################################### */

/*     iGetFileVersion - vrati verzi souboru */

/* ####################################################################### */
integer igetfileversion_(integer *fileunit, integer *itext)
{
    /* System generated locals */
    integer ret_val, i__1;
    alist al__1;

    /* Local variables */
    extern integer len_trim__(char *, ftnlen);
    static char cversion[10];
    static shortint iversion;
    extern /* Subroutine */ int emptystr_(char *, ftnlen);
    static integer i__, j;
    static char cline[255];
    static integer ipcplen;
    static char cpcpstr[17];

    /* Fortran I/O blocks */
    static cilist io___433 = { 1, 0, 0, "(a)", 0 };
    static icilist io___437 = { 1, cversion, 0, 0, 10, 1 };
    static cilist io___438 = { 1, 0, 1, 0, 0 };
    static cilist io___440 = { 1, 0, 1, 0, 0 };


    emptystr_(cversion, (ftnlen)10);
    emptystr_(cline, (ftnlen)255);
    ret_val = 0;
    iversion = 0;
    al__1.aerr = 1;
    al__1.aunit = *fileunit;
    i__1 = f_rew(&al__1);
    if (i__1 != 0) {
	goto L1000;
    }
    if (*itext == 1) {
/* Textovy soubor */
	io___433.ciunit = *fileunit;
	i__1 = s_rsfe(&io___433);
	if (i__1 != 0) {
	    goto L1000;
	}
	i__1 = do_fio(&c__1, cline, (ftnlen)255);
	if (i__1 != 0) {
	    goto L1000;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L1000;
	}
	i__ = i_indx(cline, "Pcp_File_Version=", (ftnlen)255, (ftnlen)17);
	if (i__ == 1) {
	    ipcplen = 17;
	    i__1 = len_trim__(cline, (ftnlen)255);
	    for (i__ = ipcplen + 1; i__ <= i__1; ++i__) {
		j = i__ - ipcplen;
		*(unsigned char *)&cversion[j - 1] = *(unsigned char *)&cline[
			i__ - 1];
	    }
	    i__1 = s_rsli(&io___437);
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = do_lio(&c__2, &c__1, (char *)&iversion, (ftnlen)sizeof(
		    shortint));
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1000;
	    }
	    ret_val = iversion;
	} else {
	    al__1.aerr = 1;
	    al__1.aunit = *fileunit;
	    i__1 = f_rew(&al__1);
	    if (i__1 != 0) {
		goto L1000;
	    }
	}
    } else {
/* Binarni soubor */
	io___438.ciunit = *fileunit;
	i__1 = s_rsue(&io___438);
	if (i__1 != 0) {
	    goto L1000;
	}
	i__1 = do_uio(&c__1, cpcpstr, (ftnlen)17);
	if (i__1 != 0) {
	    goto L1000;
	}
	i__1 = e_rsue();
	if (i__1 != 0) {
	    goto L1000;
	}
	i__ = i_indx(cpcpstr, "Pcp_File_Version=", (ftnlen)17, (ftnlen)17);
	if (i__ == 1) {
	    io___440.ciunit = *fileunit;
	    i__1 = s_rsue(&io___440);
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = do_uio(&c__1, (char *)&iversion, (ftnlen)sizeof(shortint));
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = e_rsue();
	    if (i__1 != 0) {
		goto L1000;
	    }
	    ret_val = iversion;
	} else {
	    al__1.aerr = 1;
	    al__1.aunit = *fileunit;
	    i__1 = f_rew(&al__1);
	    if (i__1 != 0) {
		goto L1000;
	    }
	}
    }
L1000:
    return ret_val;
} /* igetfileversion_ */

/* ####################################################################### */

/*     EMPTYSTR - vycisteni stringu */

/* ####################################################################### */
/* Subroutine */ int emptystr_(char *cstring, ftnlen cstring_len)
{
    /* Format strings */
    static char fmt_100[] = "(a1)";

    /* System generated locals */
    integer i__1, i__2;
    icilist ici__1;

    /* Local variables */
    static integer i__, imax;

    /* Fortran I/O blocks */
    static cilist io___443 = { 0, 6, 0, "(1x,a)", 0 };


    imax = i_len(cstring, cstring_len);
    i__1 = imax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici__1.icierr = 1;
	ici__1.icirnum = 1;
	ici__1.icirlen = 1;
	ici__1.iciunit = cstring + (i__ - 1);
	ici__1.icifmt = fmt_100;
	i__2 = s_wsfi(&ici__1);
	if (i__2 != 0) {
	    goto L200;
	}
	i__2 = do_fio(&c__1, " ", (ftnlen)1);
	if (i__2 != 0) {
	    goto L200;
	}
	i__2 = e_wsfi();
	if (i__2 != 0) {
	    goto L200;
	}
    }
    return 0;
L200:
    s_wsfe(&io___443);
    do_fio(&c__1, " Internal error ! ", (ftnlen)18);
    e_wsfe();
    return 0;
} /* emptystr_ */

/* *********************************************************************** */
/* Subroutine */ int init_(real *cosalf, integer *ntab, integer *itcum, 
	integer *tlevel, integer *alevel, integer *plevel, real *hroot, real *
	vroot, integer *iterw, integer *iterc, real *dtmaxc, real *wcumt, 
	real *wcuma, integer *err, logical *lvarbc, integer *nsd, real *croot,
	 real *ccumt, real *ccuma, integer *numnpd, real *sink, real *wc, 
	real *cumq, logical *lmeteo, logical *lbact, logical *lvapor, logical 
	*lenbal, logical *ldayvar, logical *lenter, logical *lfiltr, real *
	tauw, integer *nprstep, integer *ntabmod, integer *idualpor, real *
	dtmaxt, logical *lextrap, logical *lprintd, logical *llai, real *
	rextinct, logical *ldensity, real *excesint, logical *lminstep, 
	logical *lprint, logical *lsnow, real *snowmf, real *snowlayer, real *
	ctemp, integer *isunsh, integer *irelhum, real *xroot, logical *
	lcentrif, real *radius, real *gwl0l, doublereal *tprintint, integer *
	itort, integer *ienhanc, real *hseep, logical *leqinit, logical *
	lsinprec, integer *imoistdep, real *omegac, real *wtransf, logical *
	ldualneq, logical *lflux, integer *icrop, integer *irootin, logical *
	lmassini, logical *lmetdaily, real *sorb2, logical *lactrsu, real *
	omegas, real *spot, logical *lomegaw, real *omegaw, logical *lend, 
	logical *lvaporout, logical *lfluxout)
{
    /* System generated locals */
    integer sorb2_dim1, sorb2_offset, i__1, i__2;

    /* Local variables */
    static integer i__, js, inob;
    static logical lseep;

    /* Parameter adjustments */
    --ntab;
    --ccuma;
    --ccumt;
    --croot;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    --ctemp;
    --wc;
    --sink;
    --cumq;

    /* Function Body */
    *lprint = TRUE_;
    *lminstep = TRUE_;
    *lvarbc = FALSE_;
    *lvaporout = TRUE_;
/* special Nod_inf_v.out file with thermal and iso */
    *lfluxout = FALSE_;
/* special T_Level1.out file with various boundary */
    *cosalf = 1.f;
    ntab[1] = 100;
    *itcum = 0;
    *tlevel = 1;
    *alevel = 1;
    *plevel = 1;
    *hroot = 0.f;
    *vroot = 0.f;
    *iterw = 0;
    *iterc = 0;
    *dtmaxc = 1e30f;
    *dtmaxt = 1e30f;
    *wcumt = 0.f;
    *wcuma = 0.f;
    *wtransf = 0.f;
    *err = 0;
    inob = 1;
    *xroot = 0.f;
    *gwl0l = 0.f;
    *lend = FALSE_;
    i__1 = *nsd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	croot[i__] = 0.f;
	ccumt[i__] = 0.f;
	ccuma[i__] = 0.f;
/* L11: */
    }
    i__1 = *numnpd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sink[i__] = 0.f;
	wc[i__] = 0.f;
	ctemp[i__] = 0.f;
	i__2 = *nsd;
	for (js = 1; js <= i__2; ++js) {
	    sorb2[js + i__ * sorb2_dim1] = 0.f;
/* L12: */
	}
    }
    for (i__ = 1; i__ <= 12; ++i__) {
	cumq[i__] = 0.f;
/* L13: */
    }
/*     New Options - Supported by GUI */
/*     lBact   - Virus transport, ka,kd concept */
/*     lFiltr  - Filtration theory */
/*     lVapor  - Vapor flow */
/*     lEnter  - End the run with pushing Enter key */
/*     nPrStep - Print to the screen and T_Level files at each nPrStep */
/*     lPrintD - Print at a daily (given) interval */
/*     tPrintInt- Print interval */
/*     nTabMod - Kode for the input of the soil hydraulic properties tables */
/*     iDualPor- Dual porosity model; */
/*               = 1: transfer proportional to difference in water contents */
/*               = 2: transfer proportional to difference in pressure heads */
/*     lDualNEq- Both physical and chemical nonequilibrium are considered simultaneously, two-site sorption in
 the mobile zone */
/*               (fraction of equilibrium sites in mobile zone - ChPar(13), Rate of kinetic sorption - ChPar(1
6) */
/*     lMeteo  - Meterorological input to calculate ET */
/*     lLAI    - Distribution of pET based on LAI */
/*     lDayVar - Daily variations in root water uptake and evaporation */
/*     lSinPrec- Sinusoidal distribution of precipitation */
/*     lEnBal  - Evaporation and heat flux is calculated from energy balance */
/*     iSunSh  - =0 Sunshine hours; =1 Cloudeness; =2 Transmission coeff. */
/*     iRelHum - =0 Relative Humidity; =1 Vapor Pressure */
/*     lMetDaily - Daily variations of meteorological variables are generated from daily average, max, and min
 data */
/*     lSnow   - Snow accumulation at the soil surface */
/*     SnowMF  - Amount of snow melted per 1 degree [cm3/cm2/K/d] */
/*     SnowLayer- Thickness of the snow layer */
/*     iTort   - tortuosity in solute transport, =0 for Millington and =1 for Moldrup */
/*     lSeep   - Seepage face initiated by different pressure head */
/*     hSeep   - seepage face with a different bottom pressure */
/*     lEqInit - initial noequilibrium phase is in equilibrium with liquid phase */
/*     lMassIni- initial condition is given in the total concentration [M_solute/M_soil] */
/*     iMoistDep- reaction rates are dependent on the water content (iMoistDep=2) */
/*     OmegaC  - Compensated root water uptake */
/*     lFlux   - Print fluxes in observation nodes instead of temperatures */
/*     lActRSU - Active root solute uptake */
/*     OmegaS  - Compensated root solute active uptake */
/*     SPot    - Potential root solute uptake */
/*     lOmegaW - Reduction of the potential root solute uptake due to reduction of the root water uptake */
    *lbact = FALSE_;
    *lfiltr = FALSE_;
    *lvapor = FALSE_;
    *lenter = TRUE_;
    *lprintd = FALSE_;
    *ldualneq = FALSE_;
    *lmeteo = FALSE_;
    *llai = FALSE_;
    *ldayvar = FALSE_;
    *lsinprec = FALSE_;
    *lenbal = FALSE_;
    *lmetdaily = FALSE_;
    *lsnow = FALSE_;
    lseep = FALSE_;
    *leqinit = FALSE_;
    *lmassini = FALSE_;
    *lflux = FALSE_;
    *lactrsu = FALSE_;
    *lomegaw = FALSE_;
    *nprstep = 1;
    *tprintint = 86400.f;
    *ntabmod = 10;
    *idualpor = 0;
    *isunsh = 0;
    *irelhum = 0;
    *icrop = 0;
    *snowmf = .45f;
/* cm */
    *itort = 0;
    *hseep = 0.f;
    *imoistdep = 0;
    *omegac = 1.f;
    *omegas = 1.f;
    *spot = 0.f;
    *omegaw = 1.f;
    *excesint = 0.f;
    *rextinct = .39f;
    *snowlayer = 0.f;
/*     New Options - Not supported by GUI */
/*     TauW    - Nonequilibrium water flow [Ross and Smettem, 2001] */
/*     lExtrap - Extrapolation of hNew from previous time step */
/*     lDensity- Density dependent flow and transport */
/*     lCentrif- Gravitational acceleration in the cetrifuge. */
/*     Radius  - Distance of the sample from the centre of centrifuge */
/*     iEnhanc - Enhancement factor for vapor flow, =1 OK, =0 no. */
    *lextrap = TRUE_;
    *ldensity = FALSE_;
    *lcentrif = FALSE_;
    *radius = 0.f;
    *irootin = -1;
    *tauw = 0.f;
    *ienhanc = 1;
    return 0;
} /* init_ */

