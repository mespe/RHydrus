/* TIME.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__3 = 3;
static real c_b206 = 365.f;
static real c_b209 = 1.f;
static doublereal c_b225 = 5.253;
static doublereal c_b237 = .14285714285714285;
static doublereal c_b241 = .25;

/* Source file TIME.FOR ||||||||||||||||||||||||||||||||||||||||||||||||| */
/* Subroutine */ int tmcont_(real *dt, real *dtmaxw, real *dtopt, real *dmul, 
	real *dmul2, real *dtmin, integer *iter, doublereal *tprint, 
	doublereal *tatm, doublereal *t, doublereal *tmax, real *dtmaxc, 
	integer *itmin, integer *itmax, logical *lminstep, real *dtinit)
{
    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static doublereal tfix;
    static real dtmax;
    static integer istep;

    if (*lminstep) {
/* Computing MIN */
	r__1 = min(*dtmaxw,*dtmaxc), r__1 = min(r__1,*dtinit);
	dtmax = dmin(r__1,*dtopt);
	*dtopt = dtmax;
	*lminstep = FALSE_;
    } else {
	dtmax = dmin(*dtmaxw,*dtmaxc);
    }
/* Computing MIN */
    d__1 = min(*tprint,*tatm);
    tfix = min(d__1,*tmax);
    if (*iter <= *itmin && tfix - *t >= *dmul * *dtopt) {
/* Computing MIN */
	r__1 = dtmax, r__2 = *dmul * *dtopt;
	*dtopt = dmin(r__1,r__2);
    }
    if (*iter >= *itmax) {
/* Computing MAX */
	r__1 = *dtmin, r__2 = *dmul2 * *dtopt;
	*dtopt = dmax(r__1,r__2);
    }
/* Computing MIN */
    r__1 = *dtopt, r__2 = (real) (tfix - *t);
    *dt = dmin(r__1,r__2);
    istep = 1;
    if (*dt > 0.f) {
	r__1 = (real) (tfix - *t) / *dt;
	istep = r_nint(&r__1);
    }
    if (istep >= 1 && istep <= 10) {
/* Computing MIN */
	r__1 = (real) (tfix - *t) / istep;
	*dt = dmin(r__1,dtmax);
    }
    if (istep == 1) {
	*dt = (real) (tfix - *t);
	if (*dt - dtmax > *dtmin) {
	    *dt /= 2.f;
	}
    }
    if (*dt <= 0.f) {
	*dt = *dtmin / 3.f;
    }
    return 0;
} /* tmcont_ */

/* *********************************************************************** */
doublereal rtime_(shortint *imonth, shortint *iday, shortint *ihours, 
	shortint *imins, shortint *isecs, shortint *i100th)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer noday, nmonth;

    if (*imonth == 1 || *imonth == 3 || *imonth == 5 || *imonth == 7 || *
	    imonth == 8 || *imonth == 10 || *imonth == 12) {
	noday = 31;
    } else if (*imonth == 4 || *imonth == 6 || *imonth == 9 || *imonth == 11) 
	    {
	noday = 30;
    } else if (*imonth == 2) {
	noday = 28;
    }
    nmonth = noday * 24.f * 60.f * 60.f;
    ret_val = nmonth + *iday * 24.f * 60.f * 60.f + *ihours * 60.f * 60.f + *
	    imins * 60.f + *isecs + *i100th / 100.f;
    return ret_val;
} /* rtime_ */

/* *********************************************************************** */
/* Subroutine */ int setbc_(doublereal *tmax, doublereal *tatm, real *rtop, 
	real *rroot, real *rbot, real *hcrita, real *hbot, real *htop, real *
	gwl0l, logical *topinf, logical *botinf, real *ct, real *cbot, 
	integer *ns, real *ttop, real *tbot, real *ampl, logical *ltemp, 
	logical *lchem, integer *kodtop, logical *lvarbc, integer *ierr, 
	logical *lminstep, logical *lmeteo, real *prec, real *rsoil, logical *
	llai, real *rextinct, logical *lcentrif, real *cosalf, real *xconv, 
	real *tconv, integer *imodel, real *htopn, integer *irootin, real *
	xroot, logical *wlayer, logical *llinear, logical *lactrsu, real *
	spot)
{
    /* Format strings */
    static char fmt_101[] = "(a3)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    alist al__1;

    /* Local variables */
    static real g, cb[11], hb, rb;
    static integer jj;
    static real sc, ht, rr, hca, rlai, rpet, hnewt;
    static char theend[3];
    static logical lvargr;
    static real xlimit;
    static integer kodtold;
    static real rtopold;

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 31, 0, fmt_101, 0 };
    static cilist io___9 = { 1, 31, 0, 0, 0 };
    static cilist io___15 = { 1, 31, 0, 0, 0 };
    static cilist io___16 = { 1, 31, 0, 0, 0 };
    static cilist io___19 = { 1, 31, 0, 0, 0 };
    static cilist io___20 = { 1, 31, 0, 0, 0 };
    static cilist io___21 = { 1, 31, 0, 0, 0 };


    /* Parameter adjustments */
    --llinear;
    --cbot;
    --ct;

    /* Function Body */
    kodtold = *kodtop;
    s_rsfe(&io___7);
    do_fio(&c__1, theend, (ftnlen)3);
    e_rsfe();
    if (s_cmp(theend, "end", (ftnlen)3, (ftnlen)3) == 0) {
	*tmax = *tatm;
	return 0;
    } else {
	al__1.aerr = 0;
	al__1.aunit = 31;
	f_back(&al__1);
	if (*irootin != 0) {
	    if (! (*lchem) && ! (*ltemp)) {
		i__1 = s_rsle(&io___9);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*prec), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*rsoil), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rr, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hca, (ftnlen)sizeof(real)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&ht, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_rsle();
		if (i__1 != 0) {
		    goto L901;
		}
	    } else if (*ltemp && ! (*lchem)) {
		i__1 = s_rsle(&io___15);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*prec), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*rsoil), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rr, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hca, (ftnlen)sizeof(real)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&ht, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*ttop), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*tbot), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*ampl), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_rsle();
		if (i__1 != 0) {
		    goto L901;
		}
	    } else {
		i__1 = s_rsle(&io___16);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*prec), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*rsoil), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rr, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hca, (ftnlen)sizeof(real)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&ht, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*ttop), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*tbot), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*ampl), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *ns;
		for (jj = 1; jj <= i__2; ++jj) {
		    i__1 = do_lio(&c__4, &c__1, (char *)&ct[jj], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_lio(&c__4, &c__1, (char *)&cb[jj - 1], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = e_rsle();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	} else {
	    if (! (*lchem) && ! (*ltemp)) {
		i__1 = s_rsle(&io___19);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*prec), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*rsoil), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rr, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hca, (ftnlen)sizeof(real)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&ht, (ftnlen)sizeof(real))
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
	    } else if (*ltemp && ! (*lchem)) {
		i__1 = s_rsle(&io___20);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*prec), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*rsoil), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rr, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hca, (ftnlen)sizeof(real)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&ht, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*ttop), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*tbot), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*ampl), (ftnlen)sizeof(
			real));
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
	    } else {
		i__1 = s_rsle(&io___21);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*prec), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*rsoil), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rr, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hca, (ftnlen)sizeof(real)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&rb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&hb, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&ht, (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*ttop), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*tbot), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&(*ampl), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *ns;
		for (jj = 1; jj <= i__2; ++jj) {
		    i__1 = do_lio(&c__4, &c__1, (char *)&ct[jj], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_lio(&c__4, &c__1, (char *)&cb[jj - 1], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
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
	    }
	}
	if (*llai) {
	    rlai = rr;
	    rpet = *rsoil;
	    sc = 1.f;
	    rr = 0.f;
	    if (rlai > 0.f) {
/* Computing MAX */
		r__1 = 0.f, r__2 = 1.f - exp(-dmax(*rextinct,.1f) * rlai);
		rr = rpet * dmax(r__1,r__2) * sc;
	    }
	    *rsoil = rpet - rr;
	}
    }
/*     Top of the profile */
    if (*topinf) {
	lvargr = FALSE_;
/* Variable gravity field - for Scott Jones */
	if (lvargr) {
	    if (*cosalf != *prec) {
		*lminstep = TRUE_;
	    }
	    *cosalf = *prec;
	    *prec = 0.f;
	}
	rtopold = *rtop;
	*hcrita = -dabs(hca);
	if (*lvarbc) {
	    *rtop = *prec;
	    if ((r__1 = rtopold - *rtop, dabs(r__1)) > dabs(*rtop) * .2f && *
		    rtop < 0.f) {
		*lminstep = TRUE_;
	    }
	    *kodtop = (integer) (*rsoil);
	    *rsoil = 0.f;
	    if (*kodtop == -1 && kodtold == 1 && *prec > 0.f && hnewt > 0.f) {
		hnewt = *xconv * -.01f;
	    }
	} else {
	    if (! (*lmeteo)) {
		*rtop = dabs(*rsoil) - dabs(*prec);
	    }
	    if (*lmeteo) {
		*rtop = -dabs(*prec);
	    }
	    if ((r__1 = rtopold - *rtop, dabs(r__1)) > dabs(*rtop) * .2f && *
		    rtop < 0.f) {
		*lminstep = TRUE_;
	    }
	    if (*rtop > 0.f && rtopold < 0.f && ! (*wlayer)) {
		xlimit = 0.f;
		if (*imodel == 3) {
		    xlimit = *xconv * -.03f;
		}
		if (*kodtop == 4 || *htopn > xlimit) {
		    if (*imodel != 3) {
			xlimit = *xconv * -.01f;
		    }
		    *htopn = xlimit;
		    *kodtop = -4;
		}
	    }
	}
	if (*kodtop == 3 || *lvarbc) {
	    *htop = ht;
	    if ((r__1 = *htop - ht, dabs(r__1)) > dabs(*htop) * .2f) {
		*lminstep = TRUE_;
	    }
	}
	*rroot = dabs(rr);
    }
/*     Bottom of the profile */
    if (*botinf) {
	if (*lcentrif) {
	    g = *xconv * 9.80665f / *tconv / *tconv;
	    *cosalf = hb * hb / g;
	    hb = 0.f;
	    *lminstep = TRUE_;
	}
	if ((r__1 = *rbot - rb, dabs(r__1)) > dabs(*rbot) * .2f) {
	    *lminstep = TRUE_;
	}
	*rbot = rb;
	if ((r__1 = *hbot - hb - *gwl0l, dabs(r__1)) > dabs(*hbot) * .2f) {
	    *lminstep = TRUE_;
	}
	*hbot = hb + *gwl0l;
    }
    if (*lchem) {
	i__1 = *ns;
	for (jj = 1; jj <= i__1; ++jj) {
	    if (! llinear[jj] && ct[jj] > 0.f) {
		*lminstep = TRUE_;
	    }
	    cbot[jj] = cb[jj - 1];
	    if (*lactrsu && *ns == 1) {
		cbot[jj] = 0.f;
		*spot = cb[jj - 1];
	    }
/* L11: */
	}
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
} /* setbc_ */

/* *********************************************************************** */
/* Subroutine */ int setchembc_(real *prec, real *rsoil, integer *ns, real *
	ctop, real *ct, logical *wlayer, real *hnewt, integer *kodtop, 
	integer *ktopch)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer jj;

    /* Parameter adjustments */
    --ct;
    --ctop;

    /* Function Body */
    i__1 = *ns;
    for (jj = 1; jj <= i__1; ++jj) {
	if (*wlayer && *hnewt > 0.f) {
	    ctop[jj] = ctop[jj];
/* this is handled in the main program */
	} else {
	    ctop[jj] = ct[jj];
	    if (abs(*kodtop) == 4 && *ktopch <= 0) {
		if (*prec - *rsoil > 0.f) {
		    ctop[jj] = ct[jj] * *prec / (*prec - *rsoil);
		} else if (*rsoil > 0.f) {
		    ctop[jj] = 0.f;
		}
	    }
	}
/* L11: */
    }
    return 0;
} /* setchembc_ */

/* *********************************************************************** */
/* Subroutine */ int dailyvar_(real *tconv, doublereal *t, real *rroot, real *
	rrootd)
{
    /* Local variables */
    static real pi, tremainder, tday, tperiod;

/*     Temperature, max at 1. p.m. */
/*     Radiation, max at noon, hourly values between 0-6 a.m. and 18-24 p.m. */
/*     represent 1% of daily value, sinusoid in between */
    pi = 3.141592654f;
    tperiod = 1.f;
/* one day  !24.*60.*60.*tConv */
    tday = (real) (*t) / *tconv / 86400;
/*      if(tPeriod.gt.0.) tTopA=tTop+Ampl*sin(2.*PI*sngl(t)/tPeriod-7.*PI/12.) */
/* time in day units */
    tremainder = r_mod(&tday, &tperiod);
    if (tremainder <= .264f || tremainder >= .736f) {
	*rroot = *rrootd * .24f;
    } else {
	*rroot = *rrootd * 2.75f * sin(pi * 2.f * tday / tperiod - pi * 6.f / 
		12.f);
    }
    return 0;
} /* dailyvar_ */

/* *********************************************************************** */
/* Subroutine */ int sinprec_(doublereal *t, doublereal *t1, doublereal *t2, 
	real *rprec, real *rprecd)
{
    /* Local variables */
    static real dt, pi;

/*     Cosinusoidal distribution of precipitation */
    pi = 3.141592654f;
    dt = (real) (*t2 - *t1);
    if (*rprecd > 0.f) {
	*rprec = *rprecd * (cos(pi * 2.f * (real) (*t - *t1) / dt - pi) * 1.f 
		+ 1.f);
    } else {
	*rprec = 0.f;
    }
    return 0;
} /* sinprec_ */

/* *********************************************************************** */
/* Subroutine */ int snow_(real *prec, real *dt, real *temp, real *snowmf, 
	real *snowlayer, real *revap, real *xconv, logical *lminstep, real *
	ctop, real *ct, integer *ns)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static real revapold, snowmelt, q;
    static integer jj;
    static real snowlayero, rtop, snowf, precold;

    /* Parameter adjustments */
    --ct;
    --ctop;

    /* Function Body */
    precold = *prec;
    revapold = *revap;
    q = 1.f;
    if (*snowlayer < *xconv * .001f) {
	if (*temp < -2.f) {
	    q = 1.f;
	} else if (*temp < 2.f) {
	    q = 1.f - (*temp + 2.f) / 4.f;
	} else {
	    q = 0.f;
	}
    }
    rtop = *prec * (1.f - q);
    snowf = *prec * q;
    if (*temp > 0.f && *snowlayer > 0.f) {
	snowmelt = *temp * *snowmf * *dt;
    } else {
	snowmelt = 0.f;
    }
    *snowlayer = *snowlayer + snowf * *dt - snowmelt;
    if (*snowlayer < 0.f) {
	snowmelt += *snowlayer;
	*snowlayer = 0.f;
    } else if (*snowlayer > 0.f && *revap > 0.f) {
	snowlayero = *snowlayer;
	if (*revap * *dt < *snowlayer) {
	    *snowlayer -= *revap * *dt;
	    *revap = 0.f;
	} else {
	    *revap = (*revap * *dt - *snowlayer) / dmax(*dt,1e-8f);
	    *snowlayer = 0.f;
	}
    }
    *prec = rtop + snowmelt / dmax(*dt,1e-8f);
    if ((r__1 = precold - *prec, dabs(r__1)) > dabs(*prec) * .2f && *prec > 
	    0.f) {
	*lminstep = TRUE_;
    }
    if (*snowlayer > *xconv * .001f) {
	i__1 = *ns;
	for (jj = 1; jj <= i__1; ++jj) {
	    if (*snowlayer + *dt * (precold - revapold) > 0.f) {
		ctop[jj] = (*snowlayer * ctop[jj] + *dt * precold * ct[jj]) / 
			(*snowlayer + *dt * (precold - revapold));
	    }
/* L11: */
	}
    }
    return 0;
} /* snow_ */

/* *********************************************************************** */
/* Subroutine */ int meteo_(integer *ikod, logical *lmetdaily, logical *
	ldayvar, doublereal *t, real *dt, doublereal *tinit, doublereal *tmax,
	 doublereal *tatm2, doublereal *tatmn, doublereal *tatm2o, real *
	dtmax, real *rlat, real *ralt, real *shwrada, real *shwradb, real *
	rlwrada, real *rlwradb, real *rlwrada1, real *rlwradb1, real *
	windheight, real *tempheight, integer *icrop, integer *ilai, real *
	rroot, real *xconv, real *tconv, real *rgrowth, integer *ngrowth, 
	integer *iinterc, real *rinterc, real *ainterc, real *excesint, 
	logical *lenbal, real *rextinct, logical *lprint, logical *lhargr, 
	integer *iradiation, integer *isunsh, integer *irelhum, integer *
	imethour, real *cloudf_ac__, real *cloudf_bc__, real *prec, real *
	precc, real *rsoil, real *evapp, real *transp, real *rns, real *rnl, 
	real *radterm, real *aeroterm, real *rst, real *etcomb, real *rad, 
	real *radn, real *rado, real *wind, real *windn, real *windo, real *
	albedo, real *albedon, real *xlai, real *xlain, real *xroot, real *
	xrootn, real *cropheight, real *cropheightn, real *ampl, real *ttop, 
	real *tmaxan, real *tminan, real *tmax1, real *tmaxn, real *tmaxo, 
	real *tmin1, real *tminn, real *tmino, real *tempa, real *tmaxa, real 
	*tmaxao, real *tmina, real *tminao, real *sunhours, real *sunhoursn, 
	real *sunhourso, real *rhmean, real *rhmeann, real *rhmeano, real *
	rhmax, real *rhmaxn, real *rhmaxo, real *rhmin, real *rhminn, real *
	rhmino, real *rh_a__, real *eamean, real *eameann, real *rtop, 
	integer *ierr)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    extern /* Subroutine */ int dailymet_(integer *, doublereal *, real *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, integer *, real *, real *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, integer *, real *,
	     integer *), meteoint_(integer *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, logical *), setmeteo_(real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, integer *, integer *, real *, real *, real *, real *, 
	    doublereal *, real *, real *, real *, integer *, real *, integer *
	    , real *, real *, real *, logical *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, logical *, 
	    integer *, integer *, real *, real *, real *, real *, logical *, 
	    doublereal *, integer *, integer *, integer *), setdaymet_(real *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, integer *, integer *, real *, real *, real *, real *, 
	    doublereal *, real *, real *, real *, integer *, real *, integer *
	    , real *, real *, real *, real *, real *, integer *, real *, real 
	    *, real *, integer *, real *, real *, real *, real *, doublereal *
	    , logical *, real *, real *, real *, real *, real *, real *, real 
	    *, real *, real *, real *, integer *);

    /* Parameter adjustments */
    rgrowth -= 1001;

    /* Function Body */
    *ierr = 0;
    if (*ikod == 1) {
	if (*lmetdaily) {
/* Computing MIN */
	    r__1 = *tconv * 3600.f;
	    *dtmax = dmin(r__1,*dtmax);
/* Maximum time step is less than */
	    *ldayvar = FALSE_;
	    dailymet_(&c__1, t, dt, tinit, tmax, tatm2, tatmn, tmax1, tmaxn, 
		    tmaxo, tmin1, tminn, tmino, tempa, rhmax, rhmaxn, rhmaxo, 
		    rhmin, rhminn, rhmino, rh_a__, irelhum, eamean, eameann, 
		    rad, radn, wind, windn, sunhours, sunhoursn, cropheight, 
		    cropheightn, albedo, albedon, xlai, xlain, xroot, xrootn, 
		    icrop, xconv, ierr);
	    if (*ierr != 0) {
		goto L901;
	    }
	    setdaymet_(rlat, ralt, shwrada, shwradb, rlwrada, rlwradb, 
		    rlwrada1, rlwradb1, windheight, tempheight, icrop, ilai, 
		    cropheight, albedo, xlai, xroot, tatm2, rroot, xconv, 
		    tconv, iinterc, ainterc, ngrowth, &rgrowth[1001], 
		    excesint, rextinct, prec, rsoil, iradiation, wind, rad, 
		    sunhours, isunsh, cloudf_ac__, cloudf_bc__, tempa, rh_a__,
		     t, lenbal, rst, etcomb, evapp, transp, rns, rnl, radterm,
		     aeroterm, rinterc, precc, ierr);
	    if (*ierr == 3) {
		goto L903;
	    }
	} else {
	    setmeteo_(rlat, ralt, shwrada, shwradb, rlwrada, rlwradb, 
		    rlwrada1, rlwradb1, windheight, tempheight, icrop, ilai, 
		    cropheight, albedo, xlai, xroot, tatm2, rroot, xconv, 
		    tconv, iinterc, ainterc, ngrowth, &rgrowth[1001], 
		    excesint, rextinct, lenbal, prec, rsoil, iradiation, 
		    tmaxan, tminan, windn, rhmeann, radn, sunhoursn, lprint, 
		    isunsh, irelhum, ttop, ampl, cloudf_ac__, cloudf_bc__, 
		    lhargr, tinit, &c__1, imethour, ierr);
	    if (*ierr != 0) {
		switch (*ierr) {
		    case 1:  goto L901;
		    case 2:  goto L903;
		}
	    }
	    meteoint_(&c__1, tinit, tatm2o, tatm2, rad, rado, radn, tmaxa, 
		    tmaxao, tmaxan, tmina, tminao, tminan, wind, windo, windn,
		     rhmean, rhmeano, rhmeann, sunhours, sunhourso, sunhoursn,
		     lenbal);
	    if (*ierr != 0) {
		switch (*ierr) {
		    case 1:  goto L902;
		    case 2:  goto L903;
		}
	    }
	}
    } else if (*ikod == 2) {
	if (*lmetdaily) {
	    dailymet_(&c__2, t, dt, tinit, tmax, tatm2, tatmn, tmax1, tmaxn, 
		    tmaxo, tmin1, tminn, tmino, tempa, rhmax, rhmaxn, rhmaxo, 
		    rhmin, rhminn, rhmino, rh_a__, irelhum, eamean, eameann, 
		    rad, radn, wind, windn, sunhours, sunhoursn, cropheight, 
		    cropheightn, albedo, albedon, xlai, xlain, xroot, xrootn, 
		    icrop, xconv, ierr);
	    if (*ierr != 0) {
		goto L901;
	    }
	} else {
	    meteoint_(&c__2, t, tatm2o, tatm2, rad, rado, radn, tmaxa, tmaxao,
		     tmaxan, tmina, tminao, tminan, wind, windo, windn, 
		    rhmean, rhmeano, rhmeann, sunhours, sunhourso, sunhoursn, 
		    lenbal);
	    setmeteo_(rlat, ralt, shwrada, shwradb, rlwrada, rlwradb, 
		    rlwrada1, rlwradb1, windheight, tempheight, icrop, ilai, 
		    cropheight, albedo, xlai, xroot, tatm2, rroot, xconv, 
		    tconv, iinterc, ainterc, ngrowth, &rgrowth[1001], 
		    excesint, rextinct, lenbal, prec, rsoil, iradiation, 
		    tmaxan, tminan, windn, rhmeann, radn, sunhoursn, lprint, 
		    isunsh, irelhum, ttop, ampl, cloudf_ac__, cloudf_bc__, 
		    lhargr, tinit, &c__2, imethour, ierr);
	    if (*ierr != 0) {
		switch (*ierr) {
		    case 1:  goto L901;
		    case 2:  goto L903;
		}
	    }
	}
    } else if (*ikod == 3) {
	if (*lmetdaily) {
	    dailymet_(&c__4, t, dt, tinit, tmax, tatm2, tatmn, tmax1, tmaxn, 
		    tmaxo, tmin1, tminn, tmino, tempa, rhmax, rhmaxn, rhmaxo, 
		    rhmin, rhminn, rhmino, rh_a__, irelhum, eamean, eameann, 
		    rad, radn, wind, windn, sunhours, sunhoursn, cropheight, 
		    cropheightn, albedo, albedon, xlai, xlain, xroot, xrootn, 
		    icrop, xconv, ierr);
	    if (*ierr != 0) {
		goto L901;
	    }
	    setdaymet_(rlat, ralt, shwrada, shwradb, rlwrada, rlwradb, 
		    rlwrada1, rlwradb1, windheight, tempheight, icrop, ilai, 
		    cropheight, albedo, xlai, xroot, tatm2, rroot, xconv, 
		    tconv, iinterc, ainterc, ngrowth, &rgrowth[1001], 
		    excesint, rextinct, prec, rsoil, iradiation, wind, rad, 
		    sunhours, isunsh, cloudf_ac__, cloudf_bc__, tempa, rh_a__,
		     t, lenbal, rst, etcomb, evapp, transp, rns, rnl, radterm,
		     aeroterm, rinterc, precc, ierr);
	    if (*ierr == 3) {
		goto L903;
	    }
	    *rtop = dabs(*rsoil) - dabs(*prec);
	}
	if (*lenbal && ! (*lmetdaily)) {
	    meteoint_(&c__3, t, tatm2o, tatm2, rad, rado, radn, tmaxa, tmaxao,
		     tmaxan, tmina, tminao, tminan, wind, windo, windn, 
		    rhmean, rhmeano, rhmeann, sunhours, sunhourso, sunhoursn, 
		    lenbal);
	}
    }
    return 0;
L901:
    *ierr = 1;
/* 932 */
    return 0;
L902:
    *ierr = 2;
/* 913 */
    return 0;
L903:
    *ierr = 3;
/* 933 */
    return 0;
    return 0;
} /* meteo_ */

/* *********************************************************************** */
/*     Subroutine reading meteorological input data and calculating potential */
/*     evapotranspiration using either Penman-Montheith or Hargreaves equations. */
/* Subroutine */ int setmeteo_(real *latitude, real *altitude, real *
	shortwaverada, real *shortwaveradb, real *longwaverada, real *
	longwaveradb, real *longwaverada1, real *longwaveradb1, real *
	windheight, real *tempheight, integer *icrop, integer *ilai, real *
	cropheight, real *albedo, real *lai, real *xroot, doublereal *tatm, 
	real *rroot, real *xconv, real *tconv, integer *iinterc, real *
	ainterc, integer *ngrowth, real *rgrowth, real *excesint, real *
	rextinct, logical *lenbal, real *prec, real *rsoil, integer *
	iradiation, real *tmax, real *tmin, real *wind_ms__, real *rhmean, 
	real *rad, real *sunhours, logical *lprint, integer *isunsh, integer *
	irelhum, real *taver, real *ampl, real *cloudf_ac__, real *
	cloudf_bc__, logical *lhargr, doublereal *tinit, integer *iinit, 
	integer *imethour, integer *ierr)
{
    /* Format strings */
    static char fmt_101[] = "(a3)";
    static char fmt_120[] = "(f8.2,10f10.3)";

    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4;
    alist al__1;

    /* Local variables */
    static real wind_kmd__, aeroterm;
    static integer i__, j, k;
    extern /* Subroutine */ int radglobal_(real *, real *, real *, real *, 
	    real *, real *, real *), intercept_(real *, integer *, real *, 
	    real *, real *, real *, real *, real *);
    static real ra, rc, sc, pi, rn, tt, xx, yy;
    extern /* Subroutine */ int radlongnet_(real *, real *, real *, real *, 
	    real *, real *, real *), cloudiness_(real *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *);
    static real raa, n_n__, scf, rnl, rns, sum, row;
    extern /* Subroutine */ int aero_(real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, integer *);
    static real sine, sine1, eadew;
    extern /* Subroutine */ int table_(integer *, real *, integer *, integer *
	    , doublereal *, real *, real *, real *, real *);
    static real omega, precc, dayno, evapp, cover, rconv, rad_cs__;
    static char theend[3];
    static real etcomb, cloudf, transp, hourno, ttconv, rad_csh__, ea_tmin__, 
	    ea_tmax__, radterm;
    static doublereal tatmold;
    static real rinterc;
    extern /* Subroutine */ int cropres_(real *, integer *, real *, real *, 
	    integer *, real *, real *);

    /* Fortran I/O blocks */
    static cilist io___52 = { 0, 33, 0, fmt_101, 0 };
    static cilist io___54 = { 1, 33, 0, 0, 0 };
    static cilist io___56 = { 1, 33, 0, 0, 0 };
    static cilist io___90 = { 0, 43, 0, fmt_120, 0 };


/*     DayNo   - Day number */
/*     Rad     - Net radiation flux at the surface [MJ/m2/d] */
/*     TMax    - Maximum temperature [C] */
/*     TMin    - Minimum temperature [C] */
/*     RHMean  - Relative humidity [%] */
/*     Wind_kmd- Average daily wind speed [km/d] */
/*     Albedo  - Albedo [-] */
/*     xRoot   - Rooting depth [L] */
/*     iRadiation - = 0: potential radiation, = 1: solar radiation, = 2: net radiation */
/*     iInterc - interception */
/*               =1: uses LAI and Soil Cover (SC) */
/*               =2; uses only SC */
    /* Parameter adjustments */
    rgrowth -= 1001;

    /* Function Body */
    pi = 3.141592654f;
    raa = 0.f;
    rc = 60.f;
    scf = 0.f;
/*     Conversion to mm/d from L/T */
    rconv = *xconv * .001f;
    ttconv = *tconv * 86400.f;
    tatmold = *tatm;
    s_rsfe(&io___52);
    do_fio(&c__1, theend, (ftnlen)3);
    e_rsfe();
    if (s_cmp(theend, "end", (ftnlen)3, (ftnlen)3) == 0) {
/*        tMax=tAtm */
	return 0;
    } else {
	al__1.aerr = 0;
	al__1.aunit = 33;
	f_back(&al__1);
	if (*icrop == 3) {
	    i__1 = s_rsle(&io___54);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rad), (ftnlen)sizeof(real))
		    ;
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmax), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmin), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rhmean), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&wind_kmd__, (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*sunhours), (ftnlen)sizeof(
		    real));
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
	} else {
	    i__1 = s_rsle(&io___56);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rad), (ftnlen)sizeof(real))
		    ;
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmax), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmin), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rhmean), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&wind_kmd__, (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*sunhours), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
    }
    dayno = (real) (*tatm) / ttconv;
/* Conversion to [d] */
    if (*iinit == 1) {
/* Check time interval of m */
	if (*tatm - *tinit <= ttconv * .9999f) {
/* short interval */
	    *imethour = 1;
	} else {
	    *imethour = 0;
/* daily data */
	}
    }
    if (*icrop == 2) {
	i__ = 1000;
	j = 5;
	table_(ngrowth, &rgrowth[1001], &i__, &j, tatm, cropheight, albedo, 
		lai, xroot);
    }
    *wind_ms__ = wind_kmd__ / 86.4f;
/* Conversion, [m/s] */
    if (*lenbal) {
	return 0;
    }
    cropres_(&rc, ilai, lai, rextinct, icrop, cropheight, &scf);
    *taver = (*tmax + *tmin) / 2.f;
    if (! (*lhargr)) {
	*ampl = (*tmax - *tmin) / 2.f;
	ea_tmax__ = exp(*tmax * 17.27f / (*tmax + 237.3f)) * .6108f;
/* Equati */
	ea_tmin__ = exp(*tmin * 17.27f / (*tmin + 237.3f)) * .6108f;
/* Equati */
	if (*irelhum == 0) {
/*         When average relative humidity is the input */
	    eadew = *rhmean / (50.f / ea_tmin__ + 50.f / ea_tmax__);
/* Equati */
	} else if (*irelhum == 1) {
/*         When average daily vapor pressure is the input */
	    eadew = *rhmean;
	    *rhmean = eadew * (50.f / ea_tmin__ + 50.f / ea_tmax__);
/* Equati */
	}
    }
/*     Net shortwave radiation */
    dayno = r_mod(&dayno, &c_b206);
    if (*iradiation != 2 || *lhargr) {
	radglobal_(&ra, latitude, &dayno, &omega, &xx, &yy, &sc);
	if (*lhargr) {
	    goto L11;
	}
	cloudiness_(&cloudf, isunsh, sunhours, &omega, longwaveradb, 
		longwaverada, &n_n__, &cover, &tt);
	if (*iradiation == 0 || *isunsh == 3) {
	    if (*iradiation == 0) {
		*rad = ra * (*shortwaverada + *shortwaveradb * n_n__);
	    }
	    if (*isunsh == 3) {
		tt = *rad / ra;
/* Computing MAX */
/* Computing MIN */
		r__3 = 1.f, r__4 = 2.33f - tt * 3.33f;
		r__1 = .1f, r__2 = dmin(r__3,r__4);
		cover = dmax(r__1,r__2);
		n_n__ = 1.f - cover;
	    }
	    if (*iradiation == 1 && *isunsh == 3) {
		rad_cs__ = ra * (*shortwaverada + *shortwaveradb * 1.f);
/* solar */
		if (*imethour == 1) {
/* Short term interval, calculate da */
		    sum = 0.f;
		    for (k = 1; k <= 24; ++k) {
			sine1 = xx + yy * cos(pi * 2.f / 24.f * (k - 12.f));
			sum += dmax(sine1,0.f) / 24.f;
/* L10: */
		    }
		    hourno = r_mod(&dayno, &c_b209) * 24.f;
		    sine = xx + yy * cos(pi * 2.f / 24.f * (hourno - 12.f));
/* Hourly variat */
/* Computing MAX */
		    r__1 = sine * rad_cs__ / sum;
		    rad_csh__ = dmax(r__1,0.f);
		    if (rad_csh__ <= 1e-4f) {
			cloudf = *cloudf_ac__ * .6f + *cloudf_bc__;
/* averag */
		    } else if (*rad >= rad_csh__) {
			cloudf = *cloudf_ac__ * 1.f + *cloudf_bc__;
/* Equati */
		    } else {
/* Computing MAX */
			r__1 = *cloudf_ac__ * *rad / rad_csh__ + *cloudf_bc__;
			cloudf = dmax(r__1,.01f);
/* Equat */
		    }
		} else {
/* Daily interval */
		    cloudf = *cloudf_ac__ * *rad / rad_cs__ + *cloudf_bc__;
/* Equati */
		}
	    }
	}
	rns = (1.f - *albedo) * *rad;
/*       Calculate net longwave radiation */
/* Equati */
	radlongnet_(&rnl, tmax, tmin, longwaverada1, longwaveradb1, &eadew, &
		cloudf);
/*       Net radiation */
	rn = rns - rnl;
/* Equati */
    } else {
	rn = *rad;
    }
/*     Calculate aerodynamic and radiation terms of the Penman-Montheith equation */
    aero_(&aeroterm, &radterm, &rn, cropheight, windheight, tempheight, 
	    wind_ms__, altitude, taver, tmax, tmin, &rc, &raa, &ea_tmax__, &
	    ea_tmin__, &eadew, ierr);
    if (*ierr == 3) {
	return 0;
    }
/*     Evapotranspiration */
/* Computing MAX */
    r__1 = 0.f, r__2 = radterm + aeroterm;
    etcomb = dmax(r__1,r__2);
/* Equati */
    row = 1e3f;
/* (ms) water density [kg/m3] */
/* Computing 2nd power */
    r__1 = *taver - 4.f;
/* Computing 3rd power */
    r__2 = *taver - 4.f;
    row = (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 3.79e-8f) * 
	    1e3f;
    etcomb = etcomb / row * 1e3f;
/* conversion from [kg/m2/d] to [mm/d] */
L11:
    if (*lhargr) {
/* Hargreaves Formula */
	etcomb = ra * 9.3839999999999993e-4f * (*taver + 17.8f) * sqrt(*tmax 
		- *tmin);
/* [mm/d], 0.4 */
    }
/*     Potential Evaporation and Transpirations [mm/d] */
    evapp = etcomb * (1.f - scf);
    transp = etcomb * scf;
/*     Calculate interception */
    *prec = *prec / rconv * ttconv;
/* to mm/d */
    precc = *prec;
    intercept_(&rinterc, iinterc, lai, ainterc, &scf, prec, excesint, &transp)
	    ;
/*     conversion from mm/d to L/T */
/* Computing MAX */
    r__1 = transp - rinterc;
    *rroot = dmax(r__1,0.f) * rconv / ttconv;
    *rsoil = evapp * rconv / ttconv;
    *prec = *prec * rconv / ttconv;
    if (dayno == 0.f) {
	dayno = 365.f;
    }
    if (*lprint) {
	s_wsfe(&io___90);
	do_fio(&c__1, (char *)&dayno, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&etcomb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&evapp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&transp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rns, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rnl, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&radterm, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&aeroterm, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&precc, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rinterc, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*excesint), (ftnlen)sizeof(real));
	e_wsfe();
    }
    return 0;
/*     Error when reading from an input file */
L901:
    *ierr = 1;
    return 0;
} /* setmeteo_ */

/* *********************************************************************** */
/*     Calculate aerodynamic and radiation terms of the Penman-Montheith equation */
/* Subroutine */ int aero_(real *aeroterm, real *radterm, real *rn, real *
	cropheight, real *windheight, real *tempheight, real *wind_ms__, real 
	*altitude, real *taver, real *tmax, real *tmin, real *rc, real *raa, 
	real *ea_tmax__, real *ea_tmin__, real *eadew, integer *ierr)
{
    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static real aerotcff, aerdynres, dl0, dlt, patm, dl_dl__, gamma, gm_dl__, 
	    gamma1, lambda, eamean;

    dl0 = *cropheight * .667f;
    if (dl0 >= *windheight || dl0 >= *tempheight) {
	goto L901;
    }
    aerdynres = log((*windheight - dl0) / (*cropheight * .123f)) * log((*
	    tempheight - dl0) / (*cropheight * .0123f)) / .16809999999999997f;
/* Equati */
    if (*wind_ms__ > 0.f) {
	*raa = aerdynres / *wind_ms__;
    }
/* Equati */
    aerotcff = 187340.42880000002f / aerdynres / 1.01f;
/* Equati */
    d__1 = (doublereal) ((293.f - *altitude * .0065f) / 293.f);
    patm = pow_dd(&d__1, &c_b225) * 101.3f;
/* Equati */
    lambda = 2.501f - *taver * .002361f;
/* Equati */
    gamma = patm * .0016286f / lambda;
/* Equati */
    gamma1 = gamma;
    if (*raa > 0.f) {
	gamma1 = gamma * (*rc / *raa + 1);
    }
/* Equati */
/* Computing 2nd power */
    r__1 = *tmax + 237.3f;
/* Computing 2nd power */
    r__2 = *tmin + 237.3f;
    dlt = *ea_tmax__ * 2049.f / (r__1 * r__1) + *ea_tmin__ * 2049.f / (r__2 * 
	    r__2);
/* Equati */
    dl_dl__ = dlt / (dlt + gamma1);
/* Equati */
    gm_dl__ = gamma / (dlt + gamma1);
/* Equati */
    eamean = (*ea_tmax__ + *ea_tmin__) / 2.f;
/* Equati */
    *aeroterm = gm_dl__ * aerotcff / (*taver + 273.f) * *wind_ms__ * (eamean 
	    - *eadew);
/*     Radiation term */
/* Equat */
    *radterm = 0.f;
    if (lambda > 0.f) {
	*radterm = dl_dl__ * *rn / lambda;
    }
/* Equati */
    return 0;
L901:
    *ierr = 3;
    return 0;
} /* aero_ */

/* *********************************************************************** */
/*     Calculate global radiation */
/* Subroutine */ int radglobal_(real *ra, real *latitude, real *dayno, real *
	omega, real *xx, real *yy, real *sc)
{
    /* Local variables */
    static real soldeclin, dr, pi, lat1, lat2;

    pi = 3.141592654f;
    *sc = 118.08f;
/* Solar Cons */
    lat1 = (*latitude - r_int(latitude)) * 1.6666666666666667f + r_int(
	    latitude);
    lat2 = lat1 * pi / 180.f;
    soldeclin = sin(pi * 2.f / 365.f * *dayno - 1.39f) * .4093f;
/* Equation 2 */
    *omega = acos(-tan(soldeclin) * tan(lat2));
/* Equation 2 */
    *xx = sin(soldeclin) * sin(lat2);
/* Equation 1 */
    *yy = cos(soldeclin) * cos(lat2);
/* Equation 1 */
    dr = cos(pi * 2.f / 365.f * *dayno) * .033f + 1.f;
/* Equation 2 */
    *ra = *sc / pi * dr * (*omega * *xx + sin(*omega) * *yy);
/* Equation 1 */
    return 0;
} /* radglobal_ */

/* *********************************************************************** */
/*     Calculate cloudiness factor */
/*     CloudF	- Cloudiness factor [-] */
/*     Cover	  - Cloud cover fraction */
/* Subroutine */ int cloudiness_(real *cloudf, integer *isunsh, real *
	sunhours, real *omega, real *longwaveradb, real *longwaverada, real *
	n_n__, real *cover, real *tt)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static real pi, nn;

    pi = 3.141592654f;
    if (*isunsh == 0) {
/* sunshine hours */
	nn = 24.f / pi * *omega;
/* Equation 2 */
/* Computing MIN */
	r__1 = *sunhours / nn;
	*n_n__ = dmin(r__1,1.f);
/* Equation 5 */
	*cloudf = *longwaveradb + *longwaverada * *n_n__;
/* Equation 5 */
	*cover = 1.f - *n_n__;
    } else if (*isunsh == 1) {
/* cloudiness */
	*cloudf = *sunhours;
	*n_n__ = (*cloudf - *longwaveradb) / *longwaverada;
	*cover = 1.f - *n_n__;
    } else if (*isunsh == 2) {
/* transmission coeff */
	*tt = *sunhours;
/* Computing MAX */
/* Computing MIN */
	r__3 = 1.f, r__4 = 2.33f - *tt * 3.33f;
	r__1 = 1e-4f, r__2 = dmin(r__3,r__4);
	*cover = dmax(r__1,r__2);
	*n_n__ = 1.f - *cover;
	*cloudf = *longwaveradb + *longwaverada * *n_n__;
/* Equation 5 */
    }
    return 0;
} /* cloudiness_ */

/* *********************************************************************** */
/*     Calculate crop canopy resistance, rc */
/*     SunHours- Bright sunshine hours per day [hr] */
/*               alternatively, it can be Tt instead of SunHours */
/*     iLAI    - LAI is calculated from crop height for grass (1), for alfalfa (2), and soil cover (3) */
/*     CropHeight - Crop Height [cm] */
/*     LAI     - Leaf Area Index [-] */
/* Subroutine */ int cropres_(real *rc, integer *ilai, real *lai, real *
	rextinct, integer *icrop, real *cropheight, real *scf)
{
    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static logical lcrop;

    if (*icrop == 0 || *cropheight <= 0.f) {
	*cropheight = .1f;
/* cm */
	lcrop = FALSE_;
    } else {
	lcrop = TRUE_;
    }
    *rc = 0.f;
    if (lcrop) {
	if (*ilai == 1) {
/* only clipped grass */
	    *lai = *cropheight * .24f;
/* Equati */
	} else if (*ilai == 2) {
/* alfalfa */
	    *lai = log(*cropheight) * 1.5f + 5.5f;
/* Equati */
	} else if (*ilai == 3) {
/* SCF */
	    if (*lai < 1.f) {
		*lai = -log(1.f - *lai) / dmax(*rextinct,.1f);
	    } else {
		*lai = 10.f;
	    }
	}
	if (*lai > 0.f) {
	    *rc = 200.f / *lai;
	}
/* Equati */
	if (*lai > 0.f) {
/* Computing MAX */
	    r__1 = 0.f, r__2 = 1.f - exp(-dmax(*rextinct,.1f) * *lai);
	    *scf = dmax(r__1,r__2);
	}
    } else {
	*rc = 0.f;
    }
    return 0;
} /* cropres_ */

/* *********************************************************************** */
/*     Calculate net longwave radiation */
/* Subroutine */ int radlongnet_(real *rnl, real *tmax, real *tmin, real *
	longwaverada1, real *longwaveradb1, real *ea, real *cloudf)
{
    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static real emissivity, sigma;

/* Computing 4th power */
    r__1 = *tmax + 273.16f, r__1 *= r__1;
/* Computing 4th power */
    r__2 = *tmin + 273.16f, r__2 *= r__2;
    sigma = (r__1 * r__1 + r__2 * r__2) * 2.45e-9f;
/* Equati */
    emissivity = *longwaverada1 + *longwaveradb1 * sqrt(*ea);
/* Equati */
    *rnl = sigma * *cloudf * emissivity;
/* Equati */
    return 0;
} /* radlongnet_ */

/* *********************************************************************** */
/*     Calculate interception [mm] */
/* Subroutine */ int intercept_(real *rinterc, integer *iinterc, real *lai, 
	real *ainterc, real *scf, real *prec, real *excesint, real *transp)
{
    /* System generated locals */
    real r__1;

    *rinterc = 0.f;
    if (*iinterc > 0) {
	if (*iinterc == 1 && *lai > 0.f) {
/* Computing MIN */
	    r__1 = *ainterc * *lai * (1.f - 1.f / (*scf * *prec / *ainterc / *
		    lai + 1));
	    *rinterc = dmin(r__1,*prec);
/* Newly intercepted */
	    if (*rinterc + *excesint > *ainterc * *lai) {
		*rinterc = *ainterc * *lai - *excesint;
	    }
/* can not be more than max interception */
	}
/* Computing MAX */
	r__1 = *prec - *rinterc;
	*prec = dmax(r__1,0.f);
	*rinterc = *excesint + *rinterc;
/* Old interception + new interception */
	*excesint = 0.f;
	if (*transp - *rinterc < 0.f) {
	    *excesint = -(*transp) + *rinterc;
	}
    }
    return 0;
} /* intercept_ */

/* *********************************************************************** */
/*     Interpolate meteo variables in time */
/*     i = 1 (initialize), i = 2 (move new into old), i = 3 (interpolate) */
/* Subroutine */ int meteoint_(integer *i__, doublereal *t, doublereal *
	tatm2o, doublereal *tatm2, real *rad, real *rado, real *radn, real *
	tmaxa, real *tmaxao, real *tmaxan, real *tmina, real *tminao, real *
	tminan, real *wind, real *windo, real *windn, real *rhmean, real *
	rhmeano, real *rhmeann, real *sunhours, real *sunhourso, real *
	sunhoursn, logical *lenbal)
{
/* for lEnBal, there is no need to interpolate Rad a */
    if (*i__ == 1) {
	*tatm2o = *t;
	*rado = *radn;
	*rad = *radn;
	*tmaxao = *tmaxan;
	*tmaxa = *tmaxan;
	*tminao = *tminan;
	*tmina = *tminan;
	*windo = *windn;
	*wind = *windn;
	*rhmeano = *rhmeann;
	*rhmean = *rhmeann;
	*sunhourso = *sunhoursn;
	*sunhours = *sunhoursn;
    } else if (*i__ == 2) {
	*tatm2o = *tatm2;
	*rado = *radn;
	*tmaxao = *tmaxan;
	*tminao = *tminan;
	*windo = *windn;
	*rhmeano = *rhmeann;
	*sunhourso = *sunhoursn;
    } else if (*i__ == 3) {
	*rad = *rado + (*radn - *rado) * ((real) (*t) - *tatm2o) / (*tatm2 - *
		tatm2o);
	if (*lenbal) {
	    *rad = *radn;
	}
	*tmaxa = *tmaxao + (*tmaxan - *tmaxao) * ((real) (*t) - *tatm2o) / (*
		tatm2 - *tatm2o);
	*tmina = *tminao + (*tminan - *tminao) * ((real) (*t) - *tatm2o) / (*
		tatm2 - *tatm2o);
	*wind = *windo + (*windn - *windo) * ((real) (*t) - *tatm2o) / (*
		tatm2 - *tatm2o);
	*rhmean = *rhmeano + (*rhmeann - *rhmeano) * ((real) (*t) - *tatm2o) /
		 (*tatm2 - *tatm2o);
	*sunhours = *sunhourso + (*sunhoursn - *sunhourso) * ((real) (*t) - *
		tatm2o) / (*tatm2 - *tatm2o);
	if (*lenbal) {
	    *sunhours = *sunhoursn;
	}
    }
    return 0;
} /* meteoint_ */

/* *********************************************************************** */
/*     Evaluate surface energy balance */
/* Subroutine */ int evapor_(doublereal *t, real *temps, real *tmaxa, real *
	tmina, real *rad, real *htop, real *tempheight, real *windheight, 
	real *wind, real *rhmean, real *heatflux, real *rtop, real *prec, 
	real *tperiod, real *latitude, real *albedo, real *sunhours, real *
	thetat, real *xconv, real *tconv, integer *iradiation, real *rns, 
	real *rnl, real *rn, real *sensflux, real *evap, real *lat, real *
	const__, integer *isunsh, real *r_h__, logical *lmetdaily, real *rst, 
	real *tempa, real *rh_a__, real *shortwaverada, real *shortwaveradb, 
	real *longwaverada, real *longwaveradb, integer *imethour)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static real g, h__, r__;
    extern /* Subroutine */ int radglobal_(real *, real *, real *, real *, 
	    real *, real *, real *);
    static real ca, ra, sc, pi, hr, tt, xx, yy;
    extern /* Subroutine */ int cloudiness_(real *, integer *, real *, real *,
	     real *, real *, real *, real *, real *);
    static real n_n__, rah;
    extern /* Subroutine */ int radlongnet1_(real *, real *, real *, real *, 
	    real *, real *, real *);
    static real roa, r_v__, r_s__, rov, row, radh, sine, xmol, evap1, omega, 
	    tmean, dayno, cover, rconv, rovsa, rovss, cloudf, tkelva, tkelvs, 
	    hourno;
    extern /* Subroutine */ int aerores_(real *, real *, real *, real *, real 
	    *, real *, real *, real *);
    extern doublereal xlatent_(real *);

/*     Rn      - Net radiation [MJ/m2/d] */
/*     Ra      - Potential radiation [MJ/m2/d] */
/*     Rad     - Solar radiation flux at the surface [MJ/m2/d] */
/*     Rns     - Net short wave radiation [MJ/m2/d] */
/*     Rst     - Daily variated net short wave radiation [MJ/m2/d] */
/*     Rnl     - Net long wave radiation [MJ/m2/d] */
/*     Rnlu    - Outgoing long wave radiation [MJ/m2/d] */
/*     Rnld    - Incoming long wave radiation [MJ/m2/d] */
/*     SC      - Solar Constant [MJ/m2/d] (=118.08) */
/*     g       - gravitational acceleration [m/s2] (9.81) */
/*     xMol    - molecular weight of water [kg/mol] (0.018015) */
/*     R       - universal gas constant [J/mol/K] (8.314) */
/*     row     - density of soil water [kg/m3] */
/*     Lat     - mass latent heat of vaporization of water [J/kg, m2/s2] */
/*     xLatent - volumetric latent heat of vaporization of water [J/m3,kg/m/s2] */
/*     TempS   - soil temperature */
/*     TempA   - atmosphere temperature */
/*     rovs    - saturated vapor  density [kg/m3] */
/*     rov     - vapor density [kg/m3] */
/*     Hr      - relative humidity [-] */
/*     uu      - friction velocity */
/*     Wind    - wind speed [m/s] */
/*     SensFlux - sensible flux [W/m2] */
/*     HeatFlux - heat flux fo soil [W/m2] */
/*     Evap    - evaporation flux [kg/s/m2] */
/*     sigma   - Stephan-Boltzmann constant [5.6697e-8 J/s/m2/K4], [4.899e-9 MJ/d/m2/K4] */
/*     Tt      - Transmission coefficient (either given or calculated from potential and */
/*               measured solar radiation) */
/*     r_v		  - Aerodynamic resistance to vapor flow [s/m] */
/*     r_h		  - Aerodynamic resistance to heat flow [s/m] (=r_v) */
/*     r_s		  - Soil surface resistance [s/m] */
/*     CloudF	- Cloudiness factor [-] */
/*     Cover	  - Cloud cover fraction */
    pi = 3.141592654f;
    g = 9.81f;
    xmol = .018015f;
    r__ = 8.314f;
    ca = 1200.f;
/*     Hourly variations of atmospheric temperature */
    if (! (*lmetdaily)) {
	tmean = (*tmaxa + *tmina) / 2.f;
	*tempa = tmean;
	if (*tperiod > 0.f) {
	    *tempa = tmean + (*tmaxa - *tmina) / 2.f * sin(pi * 2.f * (real) (
		    *t) / *tperiod - pi * 7.f / 12.f);
	}
    }
    tkelvs = *temps + 273.15f;
    tkelva = *tempa + 273.15f;
    rovsa = exp(31.3716f - 6014.79f / tkelva - tkelva * .00792495f) * .001f / 
	    tkelva;
    if (*lmetdaily) {
	*rhmean = *rh_a__;
    }
    roa = *rhmean / 100.f * rovsa;
/*     Net shortwave radiation (Cloud cover fraction to calculate Rnl) */
    if (*iradiation != 2) {
	r__1 = (real) (*t);
	dayno = r_mod(&r__1, &c_b206);
	radglobal_(&ra, latitude, &dayno, &omega, &xx, &yy, &sc);
	cloudiness_(&cloudf, isunsh, sunhours, &omega, longwaveradb, 
		longwaverada, &n_n__, &cover, &tt);
	if (*iradiation == 0 || *isunsh == 3) {
	    if (*iradiation == 0) {
		*rad = ra * (*shortwaverada + *shortwaveradb * n_n__);
	    }
	    if (*isunsh == 3) {
		if (*imethour == 1) {
/* short term interval radiation da */
		    hourno = r_mod(&dayno, &c_b209) * 24.f;
		    sine = xx + yy * cos(pi * 2.f / 24.f * (hourno - 12.f));
/* Computing MAX */
		    r__1 = sc * sine;
		    rah = dmax(r__1,0.f);
/* Hourly variations of extraterres */
		    if (rah <= 1e-4f) {
/* night time */
			cover = .6f;
/* average value calculated from Ra */
		    } else {
/* Computing MIN */
			r__1 = *rad / rah;
			tt = dmin(r__1,1.f);
/* Computing MAX */
/* Computing MIN */
			r__3 = 1.f, r__4 = 2.33f - tt * 3.33f;
			r__1 = .1f, r__2 = dmin(r__3,r__4);
			cover = dmax(r__1,r__2);
		    }
		} else {
/* daily interval radiation data */
		    tt = *rad / ra;
/* Computing MAX */
/* Computing MIN */
		    r__3 = 1.f, r__4 = 2.33f - tt * 3.33f;
		    r__1 = .1f, r__2 = dmin(r__3,r__4);
		    cover = dmax(r__1,r__2);
		}
	    }
	}
	if (*iradiation == 0 || *imethour == 0) {
	    if (*isunsh != 2) {
		tt = *rad / ra;
	    }
	    hourno = r_mod(&dayno, &c_b209) * 24.f;
	    sine = xx + yy * cos(pi * 2.f / 24.f * (hourno - 12.f));
/* Hourly */
/* Computing MAX */
	    r__1 = sc * tt * sine;
	    radh = dmax(r__1,0.f);
/* Equati */
	} else {
	    radh = *rad;
	}
/*       Surface Albedo van Bavel and Hillel (1976) */
	*albedo = .25f;
	if (*thetat > .25f) {
	    *albedo = .1f;
	}
	if (*thetat > .1f && *thetat <= .25f) {
	    *albedo = .35f - *thetat;
	}
	if (*lmetdaily) {
	    radh = *rst;
	}
	*rns = (1.f - *albedo) * radh;
/*       Net longwave radiation */
/* Equati */
	radlongnet1_(rnl, tempa, rhmean, &cover, thetat, &tkelva, &tkelvs);
	*rn = *rns + *rnl;
/* Equati */
    } else {
	*rn = *rad;
    }
/*     Calculate aerodynamic resistance to vapor flow */
    aerores_(&r_v__, tempheight, windheight, wind, &tkelvs, &tkelva, &ca, &g);
/*     Soil surface resistance (van de Griend and Owe, 1994) */
    r_s__ = 10.f;
    if (*thetat < .15f) {
	r_s__ = exp((.15f - *thetat) * 35.63f) * 10.f;
    }
/*     Aerodynamic resistances to heat transfer */
    *r_h__ = r_v__;
    *sensflux = ca * (*temps - *tempa) / *r_h__;
    h__ = *htop / *xconv;
/* conversions to m */
    hr = exp(h__ * xmol * g / r__ / tkelvs);
/*      if(hTop.lt.0.999*hCritA.and.TempS.gt.TempS1) Hr=0.0001 */
    rovss = exp(31.3716f - 6014.79f / tkelvs - tkelvs * .00792495f) * .001f / 
	    tkelvs;
    rov = rovss * hr;
/* Computing MAX */
    r__1 = 0.f, r__2 = (rov - roa) / (r_v__ + r_s__);
    *evap = dmax(r__1,r__2);
/* Computing 2nd power */
    r__1 = *temps - 4.f;
/* Computing 3rd power */
    r__2 = *temps - 4.f;
    row = (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 3.79e-8f) * 
	    1e3f;
    *lat = xlatent_(tempa) / row;
/*     Conversion from [MJ/m2/d] to [J/m2/s,W/m2] */
    rconv = 11.574074074074074f;
    *heatflux = *rn * rconv - *sensflux - *lat * *evap;
/*     Conversion to HYDRUS units */
/* [W/m2]=[MJ/m2/d]*conv-[W/m */
    evap1 = *evap / row * *xconv / *tconv;
/* [L/T]=[kg/m2/s][m3/kg]*con */
    *const__ = *lat * row / *xconv * *tconv / *tconv / *tconv / *tconv;
/* [J/kg][kg/m3]*conv= */
    *rtop = -(*prec) + evap1;
/*     Unit conversion from W/m2 [kg/s3] to [kg/T3] */
    *heatflux = -(*heatflux) / *tconv / *tconv / *tconv;
/* Negative is in the pro */
    return 0;
} /* evapor_ */

/* *********************************************************************** */
/*     Net longwave radiation */
/* Subroutine */ int radlongnet1_(real *rnl, real *tempa, real *rhmean, real *
	cover, real *thetat, real *tkelva, real *tkelvs)
{
    /* System generated locals */
    real r__1, r__2, r__3;
    doublereal d__1;

    /* Local variables */
    static real ea, es, rld, rlu, epsi, sigma, epsia, epsis;

    sigma = 4.899e-9f;
    es = exp(*tempa * 17.27f / (*tempa + 237.3f)) * .6108f;
/* Equation */
    ea = *rhmean / 100.f * es;
    d__1 = (doublereal) (ea / *tkelva);
    epsi = pow_dd(&d__1, &c_b237) * 1.24f;
/* Brutsaer */
/* Computing MAX */
/* Computing MIN */
    r__3 = (1.f - *cover * .84f) * epsi + *cover * .84f;
    r__1 = 0.f, r__2 = dmin(r__3,1.f);
    epsia = dmax(r__1,r__2);
/* Computing MIN */
    r__1 = *thetat * .18f + .9f;
    epsis = dmin(r__1,1.f);
/* Computing 4th power */
    r__1 = *tkelvs, r__1 *= r__1;
    rlu = epsis * sigma * (r__1 * r__1);
/* Computing 4th power */
    r__1 = *tkelva, r__1 *= r__1;
    rld = epsia * sigma * (r__1 * r__1);
    *rnl = rld - rlu;
    return 0;
} /* radlongnet1_ */

/* *********************************************************************** */
/*     Aerodynamic resistance to vapor flow (Camillo and Gurney, 1986) */
/* Subroutine */ int aerores_(real *r_v__, real *tempheight, real *windheight,
	 real *wind, real *tkelvs, real *tkelva, real *ca, real *g)
{
    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static real dl, pi, mo, rk, zh, zm, uu, xx, psih, zeta, psim, diff0, 
	    difft, rvmin, rvmax, theight, wheight;

/*     rK      - von Karman constant (=0.41) */
/*     psim    - atmospheric stability factor for momentum */
/*     psih    - atmospheric stability factor for heat */
/*     zm		  - roughness parameter for momentum [m] */
/*     zh		  - roughness parameter for heat transport [m] */
/*     dl		  - displacement level for heat transport [m] (0) */
/*     Mo		  - Monin-Obukhov scaling length */
/*     zeta	  - Unitless height */
    rk = .41f;
    pi = 3.141592654f;
    zm = .001f;
/* roughness parameter for momentum [m] */
    zh = zm;
/* roughness parameter for heat transport */
    dl = 0.f;
/* displacement level for heat transport */
    theight = *tempheight / 100.f;
/* conversions to m */
    wheight = *windheight / 100.f;
    if (*wind > 0.f) {
	if ((r__1 = *tkelvs - *tkelva, dabs(r__1)) < .01f) {
/* The atmosphere is neutral */
	    *r_v__ = 1.f / *wind / rk / rk * log((theight - dl) / zh) * log((
		    wheight - dl) / zm);
	    goto L12;
	}
	psim = 0.f;
/* The atmosphere is not neutral */
	psih = 0.f;
	for (i__ = 1; i__ <= 6; ++i__) {
	    uu = *wind * rk / (log((wheight - dl + zm) / zm) + psim);
	    *r_v__ = 1.f / uu / rk * (log((theight - dl + zh) / zh) + psih);
	    if (i__ == 1) {
		rvmin = *r_v__ * .1f;
		rvmax = *r_v__ * 10.f;
	    } else {
		if (*r_v__ > rvmax) {
		    *r_v__ = rvmax;
		    goto L12;
		} else if (*r_v__ < rvmin) {
		    *r_v__ = rvmin;
		    goto L12;
		}
	    }
	    mo = -(*ca) * *tkelva * uu * uu * uu / rk / *g / (*ca * (*tkelvs 
		    - *tkelva) / *r_v__);
/* Monin- */
	    zeta = (theight - dl) / mo;
/* unitless height */
	    if (zeta < 0.f) {
/* unstable */
/* Computing MAX */
		r__1 = 0.f, r__2 = 1.f - zeta * 16;
		d__1 = (doublereal) dmax(r__1,r__2);
		xx = pow_dd(&d__1, &c_b241);
		psih = log((xx * xx + 1) / 2.f) * -2.f;
		psim = log((xx + 1) / 2.f) * -2.f - log((xx * xx + 1) / 2.f) 
			+ atan(xx) * 2.f - pi / 2.f;
	    } else if (zeta > 0.f) {
/* stable */
		if (zeta < 1.f) {
		    psih = zeta * 5.f;
		    psim = psih;
		} else if (zeta >= 1.f) {
		    psih = 5.f;
		    psim = psih;
		}
	    }
/* L11: */
	}
    } else {
/* no wind */
	diff0 = 2.12e-5f;
/* Computing 2nd power */
	r__1 = *tkelva / 273.15f;
	difft = diff0 * (r__1 * r__1);
	*r_v__ = theight / difft;
    }
L12:
    return 0;
} /* aerores_ */

/* *********************************************************************** */
/*     Adjusts evaporation and heat flux based on the difference between */
/*     potential and actual evaporation */
/* Subroutine */ int updateenergy_(doublereal *t, real *vtop, real *rtop, 
	real *heatfl, real *temps, real *rns, real *rnl, real *rn, real *evap,
	 real *lat, real *sensflux, real *xconv, real *tconv, real *const__, 
	integer *tlevel, integer *nprstep, real *r_h__, real *dz, real *rlamb,
	 integer *itemp, logical *lprint, logical *lmetdaily, real *tempa, 
	real *rh_a__, real *rst)
{
    /* Format strings */
    static char fmt_110[] = "(12e13.5)";
    static char fmt_120[] = "(12e13.5)";

    /* System generated locals */
    real r__1, r__2, r__3;

    /* Local variables */
    static real deltasens, ca, cw, row, rconv, deltae, deltat, balance;

    /* Fortran I/O blocks */
    static cilist io___179 = { 0, 43, 0, fmt_110, 0 };
    static cilist io___180 = { 0, 43, 0, fmt_120, 0 };


/*     DeltaE    - Difference of Latent heat term [kg/s3] */
/*     DeltaT    - Approximate temperature change at surface [C] */
/*     DeltaSens - Difference of Sensible heat term [kg/s3] */
/*     rLamb	    - Thermal conductivity of soil surface [kg/m/s3/K] */
/*     dz        - Distance between 1st and 2nd lattice [L] */
    deltae = 0.f;
    if (*rtop > 0.f && *rtop > *vtop) {
/* Computing 2nd power */
	r__1 = *temps - 4.f;
/* Computing 3rd power */
	r__2 = *temps - 4.f;
	row = (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 3.79e-8f)
		 * 1e3f;
	deltae = (*rtop - dmax(0.f,*vtop)) * *const__;
/* [kg/T3] */
	deltae = deltae * *tconv * *tconv * *tconv;
/* to [kg/s3 */
	if (deltae > 0.f && dabs(*vtop) > 0.f) {
	    ca = 1200.f;
/* [J/m3/K] */
	    cw = 4.18e6f;
	    *rlamb = *rlamb / *xconv * *tconv * *tconv * *tconv;
/* to [kg/m/ */
	    deltat = deltae / (ca / *r_h__ + *rlamb / (*dz / *xconv) - cw * (*
		    vtop * *tconv / *xconv));
	    deltasens = ca * deltat / *r_h__;
/* [kg/s3] */
	    *sensflux += deltasens;
/* [kg/s3] */
	    *heatfl -= (deltae - deltasens) / *tconv / *tconv / *tconv;
/* [kg/T3] */
	    *evap -= (*rtop - dmax(0.f,*vtop)) / *xconv * *tconv * row;
	}
    }
/*     Conversion from [J/m2/s,W/m2] to [MJ/m2/d] */
    rconv = 11.574074074074074f;
    balance = *rn - *sensflux / rconv - *lat * *evap / rconv + *heatfl / 
	    rconv * *tconv * *tconv * *tconv;
    if (*lprint && (r__1 = (real) ((*tlevel + *nprstep - 1) / *nprstep) - (*
	    tlevel + *nprstep - 1) / (real) (*nprstep), dabs(r__1)) < 1e-4f &&
	     *itemp == 1) {
	if (! (*lmetdaily)) {
	    s_wsfe(&io___179);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*rns), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*rnl), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*rn), (ftnlen)sizeof(real));
	    r__1 = *sensflux / rconv;
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    r__2 = *lat * *evap / rconv;
	    do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
	    r__3 = *heatfl / rconv * *tconv * *tconv * *tconv;
	    do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&balance, (ftnlen)sizeof(real));
	    e_wsfe();
	} else {
	    s_wsfe(&io___180);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*rns), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*rnl), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*rn), (ftnlen)sizeof(real));
	    r__1 = *sensflux / rconv;
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    r__2 = *lat * *evap / rconv;
	    do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
	    r__3 = *heatfl / rconv * *tconv * *tconv * *tconv;
	    do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&balance, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*tempa), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*rh_a__), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*rst), (ftnlen)sizeof(real));
	    e_wsfe();
	}
    }
    return 0;
} /* updateenergy_ */

/* *********************************************************************** */
/*     Read Meteo.in and generate daily variations of air temperature */
/*     and relative humidity (April 2008 Masaru Sakai) */
/* Subroutine */ int dailymet_(integer *ikod, doublereal *t, real *dt, 
	doublereal *tinit, doublereal *tend, doublereal *tatm, doublereal *
	tatmn, real *tmax, real *tmaxn, real *tmaxo, real *tmin, real *tminn, 
	real *tmino, real *tempa, real *rhmax, real *rhmaxn, real *rhmaxo, 
	real *rhmin, real *rhminn, real *rhmino, real *rh_a__, integer *
	irelhum, real *eamean, real *eameann, real *rad, real *radn, real *
	wind_ms__, real *wind_msn__, real *sunhours, real *sunhoursn, real *
	cropheight, real *cropheightn, real *albedo, real *albedon, real *lai,
	 real *lain, real *xroot, real *xrootn, integer *icrop, real *xconv, 
	integer *ierr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static real wind_kmd__;
    extern /* Subroutine */ int humidity_(real *, real *, real *, integer *, 
	    real *, real *, real *);
    static real wind_kmdn__;
    extern /* Subroutine */ int dailytemp_(doublereal *, real *, doublereal *,
	     doublereal *, doublereal *, real *, real *, real *, real *, real 
	    *, real *, real *);
    static real rhmean;
    static integer idaily;
    static real rhmeann;

    /* Fortran I/O blocks */
    static cilist io___182 = { 1, 33, 0, 0, 0 };
    static cilist io___185 = { 1, 33, 0, 0, 0 };
    static cilist io___186 = { 1, 33, 0, 0, 0 };
    static cilist io___189 = { 1, 33, 0, 0, 0 };


/*     iDaily    - =1: initial, =2: update, =3: tFinal-1 (last day), =4: 0:00<t<24:00 */
/*     tInit     - Initial time */
/*     tEnd	    - Final time */
/*     tAtm,N    - Current time, next DOY */
/*     TMax,N(O)	- Max temperature at current, next, and previous DOY */
/*     TMin,N(O)	- Min temperature at current, next, and previous DOY */
/*     TempA		  - Air temperature at Current Time Step */
/*     RHMax,N(O)- Max relative humidity at current, next, previous DOY */
/*     RHMin,N(O)- Min relative humidity at current, next, previous DOY */
/*     RH_A		  - Relative humidity at the current time step */
/*     EaMean,N  - Daily averaged actual vapor pressure */
/*     Rad,N     - Radiation data at the current and next DOY */
/*     Wind_ms,N	- Wind speed (m/s) at the current and next DOY */
/*     SunHours,N- Sunshine hours at the current and next DOY */
/*     CropHeight,N - Crop height at the current and next DOY */
/*     Albedo,N  - Albedo at the current and next DOY */
/*     LAI,N		  - LAI at the current and next DOY */
/*     xRoot,N	  - Root depth  at the current and next DOY */
/*     Read meteorological data from Meteo.In */
/*     iDaily=1: Initialy read 1st and 2nd day's data */
/*     iDaily=2: In the middle of calculations, update data and then read the next day's data */
/*     iDaily=3: At the tFinal-1 (last day) or 0:00<t<24:00 (iDaily=4), do not read data */
/*     At the end of calculations, skip this subroutine */
    idaily = *ikod;
    if ((d__1 = *tend - *t, abs(d__1)) <= *dt * .001f) {
	return 0;
    }
    if ((d__1 = *tend - 1.f - *t, abs(d__1)) <= *dt * .001f) {
	idaily = 3;
    }
    if (idaily == 1) {
	if (*icrop == 3) {
	    i__1 = s_rsle(&io___182);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rad), (ftnlen)sizeof(real))
		    ;
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmax), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmin), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&rhmean, (ftnlen)sizeof(real))
		    ;
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&wind_kmd__, (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*sunhours), (ftnlen)sizeof(
		    real));
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
	} else {
	    i__1 = s_rsle(&io___185);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tatm), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*rad), (ftnlen)sizeof(real))
		    ;
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmax), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmin), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&rhmean, (ftnlen)sizeof(real))
		    ;
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&wind_kmd__, (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*sunhours), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	humidity_(tmax, tmin, &rhmean, irelhum, rhmax, rhmin, eamean);
/* Es */
	*wind_ms__ = wind_kmd__ / 86.4f;
/* Conversion to m/s */
    }
/*     Update Data */
    if (idaily == 2 || idaily == 3) {
	*tatm = *tatmn;
	*tmaxo = *tmax;
	*tmax = *tmaxn;
	*tmino = *tmin;
	*tmin = *tminn;
	*rhmaxo = *rhmax;
	*rhmax = *rhmaxn;
	*rhmino = *rhmin;
	*rhmin = *rhminn;
	*eamean = *eameann;
	*rad = *radn;
	*wind_ms__ = *wind_msn__;
	*sunhours = *sunhoursn;
	if (*icrop == 3) {
	    *cropheight = *cropheightn;
	    *albedo = *albedon;
	    *lai = *lain;
	    *xroot = *xrootn;
	}
    }
/*     Read Next DOY's data */
    if (idaily <= 2) {
	if (*icrop == 3) {
	    i__1 = s_rsle(&io___186);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tatmn), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*radn), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmaxn), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tminn), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&rhmeann, (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&wind_kmdn__, (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*sunhoursn), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*cropheightn), (ftnlen)
		    sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*albedon), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*lain), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*xrootn), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	    *cropheightn = *cropheightn * 100.f / *xconv;
/* conversion to cm */
	} else {
	    i__1 = s_rsle(&io___189);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tatmn), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*radn), (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tmaxn), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*tminn), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&rhmeann, (ftnlen)sizeof(real)
		    );
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&wind_kmdn__, (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_lio(&c__4, &c__1, (char *)&(*sunhoursn), (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	humidity_(tmaxn, tminn, &rhmeann, irelhum, rhmaxn, rhminn, eameann);
/* Estimate RHMax,RHmin */
	*wind_msn__ = wind_kmdn__ / 86.4f;
/* Conversion to m/s */
    }
    if (idaily == 2) {
	return 0;
    }
/*     Calculate temperature, TempA, at the current time step */
/* do not calculate temp an */
    dailytemp_(t, dt, tinit, tend, tatm, tmax, tmaxn, tmaxo, tmin, tminn, 
	    tmino, tempa);
/*     Calculate relative humidity, RH_A, at the current time step */
    dailytemp_(t, dt, tinit, tend, tatm, rhmin, rhminn, rhmino, rhmax, rhmaxn,
	     rhmaxo, rh_a__);
    return 0;
/*     Error */
L901:
    *ierr = 1;
    return 0;
} /* dailymet_ */

/* *********************************************************************** */
/*     Calculate temperature (or relative humidity) from daily max and min values */
/* Subroutine */ int dailytemp_(doublereal *t, real *dt, doublereal *tinit, 
	doublereal *tend, doublereal *tatm, real *tmax, real *tmaxn, real *
	tmaxo, real *tmin, real *tminn, real *tmino, real *tempa)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static real pi, tmina, tmaxa, mintime, maxtime;

/*     maxTime	- The time of maximum temperature 13:00 */
/*     minTime	- The time of minimum temperature 1:00 */
    pi = 3.141592654f;
    maxtime = .54166666666666663f;
    mintime = .041666666666666664f;
/*     Calculate maximum temperature, TMaxA, for current time */
    if (*t < *tinit + maxtime) {
/* Initial - 1 */
	tmaxa = *tmax;
    } else if (*t >= *tend - 1.f + maxtime) {
/* 13:00 of l */
	tmaxa = *tmax;
    } else if (*t >= *tatm - 1.f + maxtime && (*t < *tatm || (d__1 = *t - *
	    tatm, abs(d__1)) <= *dt * .001f)) {
/* 13:00-24:0 */
	tmaxa = (*tmaxn - *tmax) * (*t - (*tatm - 1.f + maxtime)) + *tmax;
    } else if ((*t > *tatm - 1.f || (d__1 = *t - (*tatm - 1.f), abs(d__1)) <= 
	    *dt * .001f) && *t < *tatm - 1.f + maxtime) {
/* 0:00-13:00 */
	tmaxa = (*tmax - *tmaxo) * (*t - (*tatm - 2.f + maxtime)) + *tmaxo;
    }
/*     Calculate minimum temperature, TMinA, for current time */
    if (*t < *tinit + mintime) {
/* Initial - 1 */
	tmina = *tmin;
    } else if (*t >= *tend - 1.f + mintime) {
/* 1:00 of la */
	tmina = *tmin;
    } else if (*t >= *tatm - 1.f + mintime && (*t < *tatm || (d__1 = *t - *
	    tatm, abs(d__1)) <= *dt * .001f)) {
/* 1:00-24:00 */
	tmina = (*tminn - *tmin) * (*t - (*tatm - 1.f + mintime)) + *tmin;
    } else if ((*t > *tatm - 1.f || (d__1 = *t - (*tatm - 1.f), abs(d__1)) <= 
	    *dt * .001f) && *t < *tatm - 1.f + mintime) {
/* 0:00-1:00 */
	tmina = (*tmin - *tmino) * (*t - (*tatm - 2.f + mintime)) + *tmino;
    }
/*     Calculate current temperature */
    *tempa = (tmaxa + tmina) / 2.f + (tmaxa - tmina) / 2.f * cos(pi * 2.f * (*
	    t - maxtime));
    return 0;
} /* dailytemp_ */

/* *********************************************************************** */
/*     Calculate max and min relative humidity from daily max and min */
/*     temperatures and average RH */
/* Subroutine */ int humidity_(real *tmax, real *tmin, real *rhmean, integer *
	irelhum, real *rhmax, real *rhmin, real *eamean)
{
    /* Local variables */
    static real es_tmin__, es_tmax__;

/*     iRelHum	- =0: Input RH, =1: Input Vapor pressure for RHMean */
/*     EaMean	- Average Saturation vapor pressure */
/*     Es_TMax	- Maximum Saturation vapor pressure */
/*     Es_TMin	- Minimum Saturation vapor pressure */
/*     RHMax   - Estimated Max RH */
/*     RHMin   - Estimated Min RH */
    es_tmax__ = exp(*tmax * 17.27f / (*tmax + 237.3f)) * .6108f;
/* Equati */
    es_tmin__ = exp(*tmin * 17.27f / (*tmin + 237.3f)) * .6108f;
/* Equati */
    if (*irelhum == 0) {
	*rhmin = *rhmean * 2.f * es_tmin__ / (es_tmax__ + es_tmin__);
	*rhmax = *rhmean * 2.f * es_tmax__ / (es_tmax__ + es_tmin__);
	*eamean = *rhmin / 100.f * es_tmax__;
    } else {
	*eamean = *rhmean;
	*rhmin = *eamean / es_tmax__ * 100.f;
	*rhmax = *eamean / es_tmin__ * 100.f;
    }
    if (*rhmax >= 100.f) {
	if (*irelhum == 1) {
	    *rhmean = *eamean * (50.f / es_tmin__ + 50.f / es_tmax__);
	}
/* Equ */
	*rhmin = *rhmean - (100.f - *rhmean);
	*rhmax = 100.f;
    }
    return 0;
} /* humidity_ */

/* *********************************************************************** */
/* Subroutine */ int setdaymet_(real *latitude, real *altitude, real *
	shortwaverada, real *shortwaveradb, real *longwaverada, real *
	longwaveradb, real *longwaverada1, real *longwaveradb1, real *
	windheight, real *tempheight, integer *icrop, integer *ilai, real *
	cropheight, real *albedo, real *lai, real *xroot, doublereal *tatm, 
	real *rroot, real *xconv, real *tconv, integer *iinterc, real *
	ainterc, integer *ngrowth, real *rgrowth, real *excesint, real *
	rextinct, real *prec, real *rsoil, integer *iradiation, real *
	wind_ms__, real *rad, real *sunhours, integer *isunsh, real *
	cloudf_ac__, real *cloudf_bc__, real *tempa, real *rh_a__, doublereal 
	*t, logical *lenbal, real *rst, real *etcomb, real *evapp, real *
	transp, real *rns, real *rnl, real *radterm, real *aeroterm, real *
	rinterc, real *precc, integer *ierr)
{
    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, k;
    extern /* Subroutine */ int radglobal_(real *, real *, real *, real *, 
	    real *, real *, real *);
    static real t2;
    extern /* Subroutine */ int intercept_(real *, integer *, real *, real *, 
	    real *, real *, real *, real *);
    static real ea, ra, rc, sc, es, pi, rn, tt, xx, yy;
    extern /* Subroutine */ int radlongnet_(real *, real *, real *, real *, 
	    real *, real *, real *), cloudiness_(real *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *);
    static real raa, n_n__, scf, sum, row;
    extern /* Subroutine */ int aero_(real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, integer *);
    static real sine, sine1;
    extern /* Subroutine */ int table_(integer *, real *, integer *, integer *
	    , doublereal *, real *, real *, real *, real *);
    static real omega, dayno, cover, rconv, rad_cs__, cloudf, ttconv;
    extern /* Subroutine */ int cropres_(real *, integer *, real *, real *, 
	    integer *, real *, real *);

/*     tAtm	- Current DOY */
/*     TempA - Air tempeature [C] */
/*     RH_A  - Air relative humidity [%] */
/*     Ea    - Daily averaged vapor pressure [kPa] */
/*     Es		- Saturation vapor pressure at current time step [kPa] */
/*     Rad   - Incoming solar radiation [MJ/m2/d] */
/*     Rst   - Solar radiation at current time step [MJ/m2/d] */
    /* Parameter adjustments */
    rgrowth -= 1001;

    /* Function Body */
    pi = 3.141592654f;
    raa = 0.f;
    rc = 60.f;
    scf = 0.f;
/*     Conversion to mm/d from L/T */
    rconv = *xconv * .001f;
    ttconv = *tconv * 86400.f;
    if (*icrop == 2) {
	i__ = 1000;
	j = 5;
	table_(ngrowth, &rgrowth[1001], &i__, &j, t, cropheight, albedo, lai, 
		xroot);
    }
    cropres_(&rc, ilai, lai, rextinct, icrop, cropheight, &scf);
/*     Calculate vapor pressure */
    es = exp(*tempa * 17.27f / (*tempa + 237.3f)) * .6108f;
/* Equation */
    ea = es * *rh_a__ / 100.f;
/*     Daily Variated Net shortwave radiation, Rst */
    r__1 = (real) (*tatm);
    dayno = r_mod(&r__1, &c_b206);
    if (*iradiation != 2) {
	radglobal_(&ra, latitude, &dayno, &omega, &xx, &yy, &sc);
/*       Calculate daily average incoming shortwave radiation */
	if (*iradiation == 0) {
	    cloudiness_(&cloudf, isunsh, sunhours, &omega, longwaveradb, 
		    longwaverada, &n_n__, &cover, &tt);
	    *rad = ra * (*shortwaverada + *shortwaveradb * n_n__);
/* Equation */
	}
/*       Average of Rst of every hour should be daily averaged Rs value */
	sum = 0.f;
	for (k = 1; k <= 24; ++k) {
	    sine1 = xx + yy * cos(pi * 2.f / 24.f * (k - 12.f));
	    sum += dmax(sine1,0.f) / 24.f;
/* L10: */
	}
	t2 = (*t - dayno) * 24.f;
	sine = xx + yy * cos(pi * 2.f / 24.f * (t2 - 12.f));
/* Computing MAX */
	r__1 = sine * *rad / sum;
	*rst = dmax(r__1,0.f);
	*rns = (1.f - *albedo) * *rst;
/*       Calculate Cloudiness Factor for net long wave radiation */
/* Equati */
	cloudiness_(&cloudf, isunsh, sunhours, &omega, longwaveradb, 
		longwaverada, &n_n__, &cover, &tt);
	if (*iradiation == 1 && *isunsh == 3) {
/* measured solar */
	    rad_cs__ = ra * (*shortwaverada + *shortwaveradb * 1.f);
/* solar radiatio */
	    cloudf = *cloudf_ac__ * *rad / rad_cs__ + *cloudf_bc__;
/* Equati */
	}
/*       Net longwave radiation */
	radlongnet_(rnl, tempa, tempa, longwaverada1, longwaveradb1, &ea, &
		cloudf);
/*       Net radiation */
	rn = *rns - *rnl;
/* Equati */
    } else {
/* Use Measured Net Radiation data */
	rn = *rad;
    }
    if (*lenbal) {
	return 0;
    }
/*     Calculate aerodynamic and radiation terms of the Penman-Montheith equation */
    aero_(aeroterm, radterm, &rn, cropheight, windheight, tempheight, 
	    wind_ms__, altitude, tempa, tempa, tempa, &rc, &raa, &es, &es, &
	    ea, ierr);
    if (*ierr == 3) {
	return 0;
    }
/*     Evapotranspiration */
/* Computing MAX */
    r__1 = 0.f, r__2 = *radterm + *aeroterm;
    *etcomb = dmax(r__1,r__2);
/* Equati */
    row = 1e3f;
/* (ms) water density [kg/m3] */
/* Computing 2nd power */
    r__1 = *tempa - 4.f;
/* Computing 3rd power */
    r__2 = *tempa - 4.f;
    row = (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 3.79e-8f) * 
	    1e3f;
    *etcomb = *etcomb / row * 1e3f;
/*     Potential Evaporation and Transpirations [mm/d] */
/* (ms) conversion [kg/m2/d] to [mm/d] */
    *evapp = *etcomb * (1.f - scf);
    *transp = *etcomb * scf;
/*     Calculate interception */
    *prec = *prec / rconv * ttconv;
/* to mm/d */
    *precc = *prec;
    intercept_(rinterc, iinterc, lai, ainterc, &scf, prec, excesint, transp);
/*     conversion from mm/d to L/T */
/* Computing MAX */
    r__1 = *transp - *rinterc;
    *rroot = dmax(r__1,0.f) * rconv / ttconv;
    *rsoil = *evapp * rconv / ttconv;
    *prec = *prec * rconv / ttconv;
    return 0;
/* 120   format(f8.2,14f10.3) */
} /* setdaymet_ */

/* *********************************************************************** */
/* 	    Output for Penman-Monteith with Daily variated meteorological information */
/* Subroutine */ int daymeteoout_(doublereal *t, real *etcomb, real *evapp, 
	real *transp, real *rns, real *rnl, real *radterm, real *aeroterm, 
	real *precc, real *rinterc, real *excesint, real *tempa, real *rh_a__,
	 real *rst, logical *lprint)
{
    /* Format strings */
    static char fmt_120[] = "(f8.2,14f10.3)";

    /* Fortran I/O blocks */
    static cilist io___225 = { 0, 43, 0, fmt_120, 0 };


    if (*lprint) {
	s_wsfe(&io___225);
	do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*etcomb), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*evapp), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*transp), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*rns), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*rnl), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*radterm), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*aeroterm), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*precc), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*rinterc), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*excesint), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*tempa), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*rh_a__), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*rst), (ftnlen)sizeof(real));
	e_wsfe();
    }
    return 0;
} /* daymeteoout_ */

/* *********************************************************************** */
/* Subroutine */ int table_(integer *n1, real *rtable, integer *n, integer *m,
	 doublereal *t, real *a, real *b, real *c__, real *d__)
{
    /* System generated locals */
    integer rtable_dim1, rtable_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static real y[4], dt;

    /* Parameter adjustments */
    rtable_dim1 = *n;
    rtable_offset = 1 + rtable_dim1;
    rtable -= rtable_offset;

    /* Function Body */
    if (*t <= rtable[rtable_dim1 + 1]) {
	i__1 = *m - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    y[i__ - 1] = rtable[(i__ + 1) * rtable_dim1 + 1];
/* L11: */
	}
    } else if (*t >= rtable[*n1 + rtable_dim1]) {
	i__1 = *m - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    y[i__ - 1] = rtable[*n1 + (i__ + 1) * rtable_dim1];
/* L12: */
	}
    } else {
	i__1 = *n1;
	for (j = 2; j <= i__1; ++j) {
	    if (*t > rtable[j - 1 + rtable_dim1] && *t <= rtable[j + 
		    rtable_dim1]) {
		dt = *t - rtable[j - 1 + rtable_dim1];
		i__2 = *m - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__ - 1] = rtable[j - 1 + (i__ + 1) * rtable_dim1] + (
			    rtable[j + (i__ + 1) * rtable_dim1] - rtable[j - 
			    1 + (i__ + 1) * rtable_dim1]) * dt / (rtable[j + 
			    rtable_dim1] - rtable[j - 1 + rtable_dim1]);
/* L13: */
		}
	    }
/* L14: */
	}
    }
    *a = y[0];
    *b = y[1];
    *c__ = y[2];
    *d__ = y[3];
    return 0;
} /* table_ */

