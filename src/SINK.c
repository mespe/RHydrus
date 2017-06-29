/* SINK.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;

/* Source file SINK.FOR ||||||||||||||||||||||||||||||||||||||||||||||||| */
/* Subroutine */ int setsnk_(integer *n, integer *nmat, integer *matnum, real 
	*x, real *hroot, real *vroot, real *sink, real *tpot, real *hnew, 
	logical *lmosink, logical *lsolred, logical *lsoladd, real *p0, real *
	poptm, real *p2h, real *p2l, real *p3, real *r2h, real *r2l, real *
	aosm, real *c50, real *p3c, real *beta, logical *lchem, integer *ns, 
	integer *nsd, real *conc, real *croot, logical *lmssink, real *thnew, 
	real *pard, real *dt, real *omegac, integer *imodel, real *con, 
	logical *lomegaw, real *omegaw, real *rbot)
{
    /* System generated locals */
    integer conc_dim1, conc_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static integer i__, j, m, ii;
    extern doublereal fq_(integer *, real *, real *);
    static real hr, dxm, sum1, sum2, alfa, cred, hred;
    static integer iter;
    static real pmin;
    extern doublereal falfa_(logical *, real *, real *, real *, real *, real *
	    , real *, real *, real *, real *);
    static real salfa, omega, aroot;
    static integer istep, nstep;
    extern doublereal fsalfa_(logical *, real *, real *, real *);
    static logical lhanks;
    static real compen, xconst, thlimit;

    /* Parameter adjustments */
    --con;
    --thnew;
    --beta;
    --hnew;
    --sink;
    --x;
    --matnum;
    pard -= 12;
    --poptm;
    --croot;
    --aosm;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    compen = 1.f;
    nstep = 1;
    if (*omegac < 1.f) {
	nstep = 2;
    }
    omega = 0.f;
    *vroot = 0.f;
    *hroot = 0.f;
    aroot = 0.f;
    i__1 = *ns;
    for (ii = 1; ii <= i__1; ++ii) {
	croot[ii] = 0.f;
/* L11: */
    }
    lhanks = FALSE_;
    if (lhanks) {
	hr = *p3;
	for (iter = 1; iter <= 5; ++iter) {
	    sum1 = 0.f;
	    sum2 = 0.f;
	    xconst = 1.f;
/* Penalty for distance */
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		if (beta[i__] > 0.f) {
		    if (i__ == *n) {
			dxm = (x[i__] - x[i__ - 1]) / 2.f;
		    } else {
			dxm = (x[i__ + 1] - x[i__ - 1]) / 2.f;
		    }
		    if (hnew[i__] > hr - xconst * x[i__]) {
			sum1 += con[i__] * beta[i__] * (hnew[i__] + xconst * 
				x[i__]) * dxm;
			sum2 += con[i__] * beta[i__] * dxm;
		    }
		}
/* L21: */
	    }
/* Computing MAX */
	    r__1 = (sum1 - *tpot) / sum2;
	    hr = dmax(r__1,*p3);
/* L22: */
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (beta[i__] > 0.f) {
		if (i__ == *n) {
		    dxm = (x[i__] - x[i__ - 1]) / 2.f;
		} else {
		    dxm = (x[i__ + 1] - x[i__ - 1]) / 2.f;
		}
/* Computing MAX */
		r__1 = -con[i__] * beta[i__] * (hr - xconst * x[i__] - hnew[
			i__]);
		sink[i__] = dmax(r__1,0.f);
/*            Sink(i)=-Con(i)*Beta(i)*(hR-xConst*x(i)-hNew(i)) */
		*vroot += sink[i__] * dxm;
		*hroot += hnew[i__] * dxm;
		aroot += dxm;
	    }
/* L23: */
	}
	if (aroot > .001f) {
	    *hroot /= aroot;
	}
	return 0;
    }
    i__1 = nstep;
    for (istep = 1; istep <= i__1; ++istep) {
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    if (beta[i__] > 0.f) {
		if (i__ == *n) {
		    dxm = (x[i__] - x[i__ - 1]) / 2.f;
		} else {
		    dxm = (x[i__ + 1] - x[i__ - 1]) / 2.f;
		}
		m = matnum[i__];
		hred = hnew[i__];
		salfa = 1.f;
		if (*lchem && *lsolred) {
		    cred = 0.f;
		    i__3 = *ns;
		    for (j = 1; j <= i__3; ++j) {
			cred += aosm[j] * conc[j + i__ * conc_dim1];
/* L15: */
		    }
		    if (*lsoladd) {
			hred += cred;
		    } else {
			salfa = fsalfa_(lmssink, &cred, c50, p3c);
		    }
		}
		alfa = falfa_(lmosink, tpot, &hred, p0, &poptm[m], p2h, p2l, 
			p3, r2h, r2l);
		if (istep != nstep) {
		    omega += alfa * salfa * beta[i__] * dxm;
		    goto L13;
		} else {
		    compen = 1.f;
		    if (omega < *omegac && omega > 0.f) {
			compen = *omegac;
		    }
		    if (omega >= *omegac) {
			compen = omega;
		    }
		}
		sink[i__] = alfa * salfa * beta[i__] * *tpot / compen;
		if (thnew[i__] - 2.5e-4f < pard[matnum[i__] * 11 + 1]) {
		    sink[i__] = 0.f;
		}
		if (*lmosink) {
		    pmin = *p3;
		}
		if (! (*lmosink)) {
		    pmin = *p0 * 10.f;
		}
		thlimit = fq_(imodel, &pmin, &pard[matnum[i__] * 11 + 1]);
/*            Sink(i)=min(Sink(i),0.5*(ThNew(i)-ParD(1,MatNum(i)))/dt) */
/* Computing MIN */
/* Computing MAX */
		r__3 = 0.f, r__4 = (thnew[i__] - thlimit) * .5f / *dt;
		r__1 = sink[i__], r__2 = dmax(r__3,r__4);
		sink[i__] = dmin(r__1,r__2);
		*vroot += sink[i__] * dxm;
		*hroot += hnew[i__] * dxm;
		i__3 = *ns;
		for (ii = 1; ii <= i__3; ++ii) {
		    if (*lchem) {
			croot[ii] += conc[ii + i__ * conc_dim1] * dxm;
		    }
/* L12: */
		}
		aroot += dxm;
	    } else {
		sink[i__] = 0.f;
	    }
	    if (beta[i__] < 0.f) {
/* Eddy Woehling's modification */
		if (i__ == *n) {
		    dxm = (x[i__] - x[i__ - 1]) / 2.f;
		} else {
		    dxm = (x[i__ + 1] - x[i__ - 1]) / 2.f;
		}
		sink[i__] = beta[i__] * *rbot;
/* Computing MAX */
		r__1 = sink[i__], r__2 = (thnew[i__] - pard[matnum[i__] * 11 
			+ 2]) * .5f / *dt;
		sink[i__] = dmax(r__1,r__2);
	    }
L13:
	    ;
	}
/* L16: */
    }
    if (aroot > .001f) {
	*hroot /= aroot;
	i__1 = *ns;
	for (ii = 1; ii <= i__1; ++ii) {
	    croot[ii] /= aroot;
/* L14: */
	}
    }
    if (*lomegaw && *tpot > 0.f) {
	*omegaw = *vroot / *tpot;
    }
    return 0;
} /* setsnk_ */

/* *********************************************************************** */
/*     Subroutine calculating root solute uptake with and without compensation */
/* Subroutine */ int setssnk_(integer *js, integer *ns, integer *n, 
	doublereal *t, real *x, real *beta, real *sink, real *sinks, integer *
	nsd, real *conc, real *omegaw, real *crootmax, logical *lactrsu, real 
	*omegas, real *spot, real *rkm, real *cmin)
{
    /* Format strings */
    static char fmt_100[] = "(3x,e14.7,1x,4e12.4)";

    /* System generated locals */
    integer conc_dim1, conc_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static real auptakea, spuptake;
    static integer i__;
    static real sauptakea, sauptakep, cc, sauptakean, dxm, omega;
    static logical llast;
    static integer istep, nstep;
    static real omega1, compen;

    /* Fortran I/O blocks */
    static cilist io___37 = { 0, 78, 0, fmt_100, 0 };


/*     Inputs: */
/*     SPot      - potential root solute uptake */
/*     OmegaS    - solute stress index */
/*     rKM       - Michaelis-Menten constant */
/*     lActRSU   - consider active root solute uptake */
/*     cRootMax  - maximum concentration for the passive solute uptake */
/*     From Water Flow */
/*     Sink(i)   - Root water uptake */
/*     OmegaW    - ratio of actual and potential transpiration */
/*     SPUptake  - passive root solute uptake (step 1) */
/*     SAUptakeP - potential active solute uptake (step 1) */
/*     SAUptakeA - uncompensated actual active solute uptake (step 2) */
/*     SAUptakeA - compensated actual active solute uptake (step 3) */
/*     SinkS(i)  - local active solute uptake */
/*     Initialization */
    /* Parameter adjustments */
    --sinks;
    --sink;
    --beta;
    --x;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    compen = 1.f;
    nstep = 1;
    if (*lactrsu) {
	nstep = 2;
    }
    if (*lactrsu && *omegas < 1.f) {
	nstep = 3;
    }
/*     step 1: Passive uptake */
/*     step 2: Active uptake without compensation */
/*     step 3: Active uptake with compensation */
    llast = FALSE_;
/* Active uptake only for the last so */
    if (llast && *js < *ns) {
	nstep = 1;
    }
    omega = 0.f;
    spuptake = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sinks[i__] = 0.f;
/* L10: */
    }
    i__1 = nstep;
    for (istep = 1; istep <= i__1; ++istep) {
	sauptakea = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (beta[i__] > 0.f) {
		if (i__ == *n) {
		    dxm = (x[i__] - x[i__ - 1]) / 2.f;
		} else if (i__ == 1) {
		    dxm = (x[i__] - x[i__ + 1]) / 2.f;
		} else {
		    dxm = (x[i__ + 1] - x[i__ - 1]) / 2.f;
		}
/* Computing MAX */
		r__1 = conc[*js + i__ * conc_dim1] - *cmin;
		cc = dmax(r__1,0.f);
		if (istep == 1) {
/* Computing MAX */
/* Computing MIN */
		    r__2 = conc[*js + i__ * conc_dim1];
		    r__1 = dmin(r__2,*crootmax);
		    sinks[i__] = sink[i__] * dmax(r__1,0.f);
		    spuptake += sinks[i__] * dxm;
/*             This is needed only for the last node, but that node may not have beta */
/* Computing MAX */
		    r__1 = *spot * *omegaw - spuptake;
		    sauptakep = dmax(r__1,0.f);
		} else if (istep == 2) {
		    auptakea = cc / (*rkm + cc) * beta[i__] * sauptakep;
		    omega += auptakea * dxm;
		    if (nstep == 2) {
			sinks[i__] += auptakea;
		    }
/*             This is needed only for the last node, but that node may not have beta */
		    sauptakea = omega;
		    sauptakean = omega;
		    if (sauptakep > 0.f) {
			omega1 = omega / sauptakep;
		    }
		} else if (istep == 3) {
/*             This is needed only for the first node, but that node may not have beta */
		    if (omega1 < *omegas && omega1 > 0.f) {
			compen = *omegas;
		    }
		    if (omega1 >= *omegas) {
			compen = omega1;
		    }
		    if (compen > 0.f) {
			auptakea = cc / (*rkm + cc) * beta[i__] * sauptakep / 
				compen;
		    }
		    sinks[i__] += auptakea;
		    sauptakea += auptakea * dxm;
		}
	    } else {
		sinks[i__] = 0.f;
	    }
/* L11: */
	}
	if (istep == nstep && *js == *ns) {
	    s_wsfe(&io___37);
	    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&spuptake, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sauptakep, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sauptakea, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sauptakean, (ftnlen)sizeof(real));
	    e_wsfe();
	}
/* the */
/* L12: */
    }
    return 0;
} /* setssnk_ */

/* *********************************************************************** */
doublereal fsalfa_(logical *lmode, real *cred, real *c50, real *p3c)
{
    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1, d__2;

    if (*lmode) {
	ret_val = 0.f;
	if (dabs(*c50) > 0.f) {
	    d__1 = (doublereal) (*cred / *c50);
	    d__2 = (doublereal) (*p3c);
	    ret_val = 1.f / (pow_dd(&d__1, &d__2) + 1.f);
	}
    } else {
	if (*cred <= *c50) {
	    ret_val = 1.f;
	} else {
/* Computing MAX */
	    r__1 = 0.f, r__2 = 1.f - (*cred - *c50) * *p3c * .01f;
	    ret_val = dmax(r__1,r__2);
	}
    }
    return ret_val;
} /* fsalfa_ */

/* *********************************************************************** */
doublereal falfa_(logical *lmosink, real *tpot, real *h__, real *p0, real *p1,
	 real *p2h, real *p2l, real *p3, real *r2h, real *r2l)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Local variables */
    static real p2;

    if (*lmosink) {
	if (*tpot < *r2l) {
	    p2 = *p2l;
	}
	if (*tpot > *r2h) {
	    p2 = *p2h;
	}
	if (*tpot >= *r2l && *tpot <= *r2h) {
	    p2 = *p2h + (*r2h - *tpot) / (*r2h - *r2l) * (*p2l - *p2h);
	}
	ret_val = 0.f;
	if (*h__ > *p3 && *h__ < p2) {
	    ret_val = (*h__ - *p3) / (p2 - *p3);
	}
	if (*h__ >= p2 && *h__ <= *p1) {
	    ret_val = 1.f;
	}
	if (*h__ > *p1 && *h__ < *p0 && *p0 - *p1 > 0.f) {
	    ret_val = (*h__ - *p0) / (*p1 - *p0);
	}
/*       Uptake even at full saturation, when both P1 and P2 are equal to zero */
	if (*h__ >= p2 && *p1 == 0.f && *p0 == 0.f) {
	    ret_val = 1.f;
	}
    } else {
	d__1 = (doublereal) (*h__ / *p0);
	d__2 = (doublereal) (*p3);
	ret_val = 1.f / (pow_dd(&d__1, &d__2) + 1.f);
    }
    return ret_val;
} /* falfa_ */

/* *********************************************************************** */
/* Subroutine */ int setrg_(integer *numnp, real *x, real *beta, doublereal *
	t, real *trmin, real *trharv, real *xrmin, real *xrmax, real *rgr, 
	real *xroot, logical *lroot, integer *irootin, integer *ngrowth, real 
	*rgrowth, real *trperiod)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__, j;
    static real tt, xr;
    extern /* Subroutine */ int table_(integer *, real *, integer *, integer *
	    , doublereal *, real *, real *, real *, real *);
    static real sbeta, troot, rdummy;

    /* Parameter adjustments */
    --beta;
    --x;
    rgrowth -= 1001;

    /* Function Body */
    if (*lroot && *irootin == 1) {
	i__ = 1000;
	j = 5;
	table_(ngrowth, &rgrowth[1001], &i__, &j, t, &rdummy, &rdummy, &
		rdummy, xroot);
    }
    if (*lroot && *irootin == 2) {
	r__1 = (real) (*t);
	troot = r_mod(&r__1, trperiod);
	if (troot < *trmin || troot > *trharv) {
	    i__1 = *numnp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		beta[i__] = 0.f;
/* L11: */
	    }
	    return 0;
	}
	xr = *xrmax;
	if (*xrmin <= .001f) {
	    *xrmin = .001f;
	}
	tt = troot - *trmin;
	xr = *xrmax * *xrmin / (*xrmin + (*xrmax - *xrmin) * exp(-(*rgr) * tt)
		);
    } else {
	xr = *xroot;
    }
    sbeta = 0.f;
    i__1 = *numnp - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__] < x[*numnp] - xr) {
	    beta[i__] = 0.f;
	} else if (x[i__] < x[*numnp] - xr * .2f) {
	    beta[i__] = 2.08333f / xr * (1 - (x[*numnp] - x[i__]) / xr);
	} else {
	    beta[i__] = 1.66667f / xr;
	}
	if (i__ != *numnp) {
	    sbeta += beta[i__] * (x[i__ + 1] - x[i__ - 1]) / 2.f;
	} else {
	    sbeta += beta[i__] * (x[i__] - x[i__ - 1]) / 2.f;
	}
/* L12: */
    }
    if (sbeta < 1e-4f) {
	beta[*numnp - 1] = 1.f / ((x[*numnp] - x[*numnp - 2]) / 2.f);
    } else {
	i__1 = *numnp - 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    beta[i__] /= sbeta;
/* L13: */
	}
    }
    return 0;
} /* setrg_ */

/* |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| */
/*     Jasper Vrugt function */
/* Subroutine */ int rootinn_(integer *numnp, real *beta, real *z__)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__;
    static real r1, r2, z0, za, zm, zmax, sbeta, rroota;

    /* Fortran I/O blocks */
    static cilist io___46 = { 0, 44, 0, 0, 0 };
    static cilist io___47 = { 0, 44, 0, 0, 0 };


/*     read input */
    /* Parameter adjustments */
    --z__;
    --beta;

    /* Function Body */
    s_rsle(&io___46);
    e_rsle();
    s_rsle(&io___47);
    do_lio(&c__4, &c__1, (char *)&zm, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&z0, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&za, (ftnlen)sizeof(real));
    e_rsle();
/*     coordinate of the surface */
    zmax = z__[*numnp];
/*     calculate non-normalized uptake intensity */
    r1 = 0.f;
    r2 = 0.f;
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dabs(zm) > 1e-5f) {
	    r1 = (zm - (zmax - z__[i__])) / zm;
	    rroota = za;
	    if (zmax - z__[i__] > z0) {
		rroota = 1.f;
	    }
	    r2 = rroota / zm * (r__1 = z0 - (zmax - z__[i__]), dabs(r__1));
	    r2 = exp(-r2);
	}
/* Computing MAX */
	r__1 = r1 * r2;
	beta[i__] = dmax(r__1,0.f);
/* L11: */
    }
/*     normalize uptake intensity */
    sbeta = beta[*numnp] * (z__[*numnp] - z__[*numnp - 1]) / 2.f;
    i__1 = *numnp - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	sbeta += beta[i__] * (z__[i__ + 1] - z__[i__ - 1]) / 2.f;
/* L12: */
    }
    i__1 = *numnp;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (sbeta > 0.f) {
	    beta[i__] /= sbeta;
	} else {
	    beta[i__] = 0.f;
	}
/* L13: */
    }
    return 0;
} /* rootinn_ */

