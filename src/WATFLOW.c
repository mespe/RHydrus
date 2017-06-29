/* WATFLOW.f -- translated by f2c (version 12.02.01).
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

static real c_b3 = .9999f;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static doublereal c_b20 = .5;
static doublereal c_b65 = 2.3333333333333335;

/* Source file WATFLOW.FOR |||||||||||||||||||||||||||||||||||||||||||||| */
/* Subroutine */ int watflow_(integer *numnp, integer *ntab, integer *ntabd, 
	integer *nmat, real *htab, real *contab, real *captab, real *hnew, 
	real *hold, integer *matnum, real *pard, real *parw, real *con, real *
	cap, real *consat, real *ah, real *ak, real *ath, real *hsat, real *
	htemp, integer *kodtop, integer *kodbot, real *rtop, real *rbot, real 
	*cosalf, doublereal *t, real *dt, real *x, real *sink, doublereal *p, 
	doublereal *r__, doublereal *s, logical *freed, logical *seepf, 
	logical *qgwlf, real *aqh, real *bqh, real *gwl0l, real *htop, real *
	hbot, real *hcrita, real *hcrits, logical *wlayer, integer *iter, 
	integer *itcum, logical *topinf, integer *ktold, integer *kbold, real 
	*tolth, real *tolh, integer *maxit, real *dtmin, doublereal *told, 
	real *dtopt, logical *convgf, real *thetab, real *thnew, real *thold, 
	real *thr, real *ths, logical *lwtdep, real *tempn, integer *kappa, 
	integer *kappao, real *aths, real *thrr, real *cono, real *conr, real 
	*aks, real *ahw, real *athw, real *akw, integer *ihyst, integer *
	imodel, logical *qdrain, real *zbotdr, real *basegw, real *rspacing, 
	integer *iposdr, real *rkhtop, real *rkhbot, real *rkvtop, real *
	rkvbot, real *entres, real *wetper, real *zintf, real *geofac, 
	logical *ltable, logical *lvapor, real *xconv, real *tconv, real *
	conlt, real *convt, real *convh, real *tauw, real *theq, real *thvnew,
	 real *thvold, integer *ntabmod, integer *idualpor, real *thnewim, 
	real *tholdim, real *sinkim, real *vtop, real *tempo, integer *itemp, 
	real *wtransf, logical *ldensity, real *conc, integer *nsd, integer *
	ienhanc, logical *lcentrif, real *radius, real *hseep)
{
    /* System generated locals */
    integer contab_dim1, contab_offset, captab_dim1, captab_offset, htab_dim1,
	     htab_offset, thetab_dim1, thetab_offset, conc_dim1, conc_offset, 
	    i__1;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static integer i__, m;
    extern doublereal fh_(integer *, real *, real *);
    static doublereal pb, rb, sb;
    static real th;
    static doublereal pt, rt, st;
    static real rate;
    static doublereal rmin;
    static real rmax, epsh;
    extern /* Subroutine */ int hyst_(integer *, integer *, real *, real *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    integer *, integer *), reset_(integer *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, logical *, real *, real *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, logical *, real *, real *, real *, real *
	    , real *, real *, logical *, real *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, logical *, logical *, real *, real *, real *, 
	    real *, real *, real *, integer *, real *, logical *, real *, 
	    integer *, logical *, real *), shift_(integer *, integer *, real *
	    , real *, real *, real *, real *, real *, logical *, real *, real 
	    *, real *, logical *, integer *, logical *, real *, real *, real *
	    , real *, logical *, logical *, real *, real *, real *, real *, 
	    integer *, real *, logical *, real *, integer *, logical *, real *
	    , real *), gauss_(integer *, integer *, integer *, real *, real *,
	     real *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static real epsth;
    extern /* Subroutine */ int setmat_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, logical *, real *, integer *, 
	    real *, integer *, real *, real *, real *, real *, real *, real *,
	     real *, integer *, logical *, logical *, real *, real *, real *, 
	    real *, real *, real *, integer *, real *, logical *, real *, 
	    integer *, integer *);
    static logical itcrit;
    extern /* Subroutine */ int hyster_(integer *, integer *, real *, integer 
	    *, real *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, integer *, 
	    integer *, real *), dualpor_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *, integer *, real *, real *, real *, real *);

    /* Fortran I/O blocks */
    static cilist io___16 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    --tempo;
    --sinkim;
    --tholdim;
    --thnewim;
    --thvold;
    --thvnew;
    --theq;
    --convh;
    --convt;
    --conlt;
    --aks;
    --conr;
    --cono;
    --thrr;
    --aths;
    --kappao;
    --kappa;
    --tempn;
    --thold;
    --thnew;
    --s;
    --r__;
    --p;
    --sink;
    --x;
    --htemp;
    --ath;
    --ak;
    --ah;
    --cap;
    --con;
    --matnum;
    --hold;
    --hnew;
    --akw;
    --athw;
    --ahw;
    --ths;
    --thr;
    thetab_dim1 = *ntabd;
    thetab_offset = 1 + thetab_dim1;
    thetab -= thetab_offset;
    --hsat;
    --consat;
    parw -= 12;
    pard -= 12;
    captab_dim1 = *ntabd;
    captab_offset = 1 + captab_dim1;
    captab -= captab_offset;
    contab_dim1 = *ntabd;
    contab_offset = 1 + contab_dim1;
    contab -= contab_offset;
    htab_dim1 = *ntabd;
    htab_offset = 1 + htab_dim1;
    htab -= htab_offset;
    --ntab;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    rmax = 1e10f;
    rmin = 1e-100;
/*     Nonequilibrium transport [Ross and Smettem, 2001] */
    rate = 1.f;
    if (*tauw > 0.f) {
/* Computing MIN */
/* Computing MAX */
	r__3 = 1e-6f, r__4 = 1.f - exp(-(*dt) / *tauw);
	r__1 = 1.f, r__2 = dmax(r__3,r__4);
	rate = dmin(r__1,r__2);
    }
/*     Dual porosity mass transfer */
    if (*idualpor > 0) {
	dualpor_(numnp, nmat, &matnum[1], idualpor, &thold[1], &ths[1], &thr[
		1], &thnewim[1], &tholdim[1], &pard[12], &sinkim[1], dt, 
		imodel, &hnew[1], hcrita, &x[1], wtransf);
    }
L11:
    *iter = 0;
    *convgf = TRUE_;
/*     End of ponding (works for both BC and VG) */
    if (*wlayer && hnew[*numnp] > 0.f && hnew[*numnp] < *xconv * 5e-5f && *
	    rtop >= 0.f) {
	hnew[*numnp] = fh_(imodel, &c_b3, &pard[matnum[*numnp] * 11 + 1]);
	hold[*numnp] = hnew[*numnp];
	htemp[*numnp] = hnew[*numnp];
    }
L12:
/*     Generate terms of matrix equation and solve by Gauss elimination */
    if (*ihyst != 3) {
	setmat_(numnp, &ntab[1], ntabd, nmat, &htab[htab_offset], &contab[
		contab_offset], &captab[captab_offset], &hnew[1], &matnum[1], 
		&pard[12], &con[1], &cap[1], &consat[1], &ah[1], &ak[1], &ath[
		1], &hsat[1], &htemp[1], &thetab[thetab_offset], &theq[1], &
		thr[1], &ths[1], lwtdep, &tempn[1], iter, &cono[1], &kappa[1],
		 &aths[1], &thrr[1], &conr[1], &aks[1], &ahw[1], &athw[1], &
		akw[1], imodel, ltable, lvapor, &thvnew[1], &conlt[1], &convt[
		1], &convh[1], xconv, tconv, ntabmod, hcrita, ldensity, &conc[
		conc_offset], nsd, ienhanc);
	if (*iter == 2 && *ihyst > 0) {
	    hyster_(numnp, nmat, &hold[1], &matnum[1], &pard[12], &parw[12], &
		    thnew[1], &thold[1], &kappa[1], &aths[1], &thrr[1], &cono[
		    1], &conr[1], &aks[1], &kappao[1], &ah[1], &ak[1], ihyst, 
		    imodel, tolth);
	}
    } else {
/* Bob Lenhard hysteresis */
	hyst_(numnp, nmat, &pard[12], &parw[12], &matnum[1], &kappa[1], &hnew[
		1], &hold[1], &theq[1], &con[1], &cap[1], &c__0, &c__2);
    }
    reset_(numnp, rtop, rbot, cosalf, dt, &x[1], &hold[1], &con[1], &cap[1], 
	    wlayer, &hnew[1], &sink[1], &p[1], &r__[1], &s[1], &pb, &rb, &sb, 
	    &pt, &rt, &st, freed, qgwlf, aqh, bqh, gwl0l, &thnew[1], &thold[1]
	    , vtop, qdrain, zbotdr, basegw, rspacing, iposdr, rkhtop, rkhbot, 
	    rkvtop, rkvbot, entres, wetper, zintf, geofac, &theq[1], &rate, 
	    lvapor, lwtdep, &tempn[1], &conlt[1], &convt[1], &convh[1], &
	    thvnew[1], &thvold[1], idualpor, &sinkim[1], ldensity, &conc[
	    conc_offset], nsd, lcentrif, radius);
    shift_(numnp, kodtop, rtop, rbot, htop, hbot, hcrita, cosalf, wlayer, &
	    con[1], &hnew[1], &x[1], topinf, kodbot, seepf, &thnew[1], &thold[
	    1], &sink[1], dt, lvapor, lwtdep, &conlt[1], &tempn[1], &convh[1],
	     &convt[1], idualpor, &sinkim[1], ldensity, &conc[conc_offset], 
	    nsd, lcentrif, radius, hseep);
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	htemp[i__] = hnew[i__];
/* L13: */
    }
    gauss_(numnp, kodtop, kodbot, htop, hbot, &hnew[1], &p[1], &r__[1], &s[1],
	     &pb, &rb, &sb, &pt, &rt, &st, &rmin);
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((r__1 = hnew[i__], dabs(r__1)) > rmax) {
	    hnew[i__] = r_sign(&rmax, &hnew[i__]);
	}
	if (abs(*kodtop) == 4 && hnew[i__] < *hcrita && i__ == *numnp) {
	    hnew[i__] = *hcrita;
	}
	if (abs(*kodtop) == 4 && hnew[i__] < *hcrita && i__ > *numnp * 9 / 10 
		&& sink[i__] <= 0.f) {
	    hnew[i__] = *hcrita;
	}
/* L17: */
    }
    ++(*iter);
    ++(*itcum);
/*     Test for convergence */
    itcrit = TRUE_;
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	epsth = 0.f;
	epsh = 0.f;
	if (htemp[i__] < hsat[m] && hnew[i__] < hsat[m]) {
/* .and.TauW.eq.0. */
	    th = thnew[i__] + cap[i__] * (hnew[i__] - htemp[i__]) / (ths[m] - 
		    thr[m]) / ath[i__] * rate;
	    epsth = (r__1 = thnew[i__] - th, dabs(r__1));
/*          if(TauW.gt.0) EpsH=abs(hNew(i)-hTemp(i))-abs(0.05*hNew(i)) */
	} else {
	    epsh = (r__1 = hnew[i__] - htemp[i__], dabs(r__1));
	}
	if (epsth > *tolth || epsh > *tolh || (r__1 = hnew[i__], dabs(r__1)) 
		> rmax * .999f) {
	    itcrit = FALSE_;
	    if ((r__1 = hnew[i__], dabs(r__1)) > rmax * .999f) {
		*iter = *maxit;
	    }
	    goto L15;
	}
/* L14: */
    }
L15:
    if (! itcrit || (*iter <= 1 || *iter <= 2 && *ihyst > 0)) {
	if (*iter < *maxit) {
	    goto L12;
	} else if (*dt <= *dtmin) {
	    *convgf = FALSE_;
	    s_wsle(&io___16);
	    do_lio(&c__9, &c__1, " The numerical solution has not converged "
		    "! ", (ftnlen)44);
	    e_wsle();
	    return 0;
	} else {
	    i__1 = *numnp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*ihyst > 0) {
		    kappa[i__] = kappao[i__];
		}
		hnew[i__] = hold[i__];
		htemp[i__] = hold[i__];
		tempn[i__] = tempo[i__];
/* L16: */
	    }
	    *kodtop = *ktold;
	    *kodbot = *kbold;
/* Computing MAX */
	    r__1 = *dt / 3;
	    *dt = dmax(r__1,*dtmin);
	    *dtopt = *dt;
	    *t = *told + *dt;
	    if (*tauw > 0.f) {
/* Computing MIN */
/* Computing MAX */
		r__3 = 1e-6f, r__4 = 1.f - exp(-(*dt) / *tauw);
		r__1 = 1.f, r__2 = dmax(r__3,r__4);
		rate = dmin(r__1,r__2);
	    }
	    *itemp = 0;
	    goto L11;
	}
    }
    if (itcrit) {
	i__1 = *numnp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    thnew[i__] += cap[i__] * (hnew[i__] - htemp[i__]) * rate;
/* L18: */
	}
    }
    if (*wlayer) {
	if (hnew[*numnp] > *hcrits) {
	    *kodtop = 4;
	    *htop = *hcrits;
	}
    }
    if (*ihyst == 3) {
	hyst_(numnp, nmat, &pard[12], &parw[12], &matnum[1], &kappa[1], &hnew[
		1], &hold[1], &thnew[1], &con[1], &cap[1], &c__0, &c__3);
    }
    return 0;
} /* watflow_ */

/* *********************************************************************** */
/* Subroutine */ int reset_(integer *n, real *rtop, real *rbot, real *cosalf, 
	real *dt, real *x, real *hold, real *con, real *cap, logical *wlayer, 
	real *hnew, real *sink, doublereal *p, doublereal *r__, doublereal *s,
	 doublereal *pb, doublereal *rb, doublereal *sb, doublereal *pt, 
	doublereal *rt, doublereal *st, logical *freed, logical *qgwlf, real *
	aqh, real *bqh, real *gwl0l, real *thnew, real *thold, real *vtop, 
	logical *qdrain, real *zbotdr, real *basegw, real *rspacing, integer *
	iposdr, real *rkhtop, real *rkhbot, real *rkvtop, real *rkvbot, real *
	entres, real *wetper, real *zintf, real *geofac, real *theq, real *
	rate, logical *lvapor, logical *lwtdep, real *temp, real *conlt, real 
	*convt, real *convh, real *thvnew, real *thvold, integer *idualpor, 
	real *sinkim, logical *ldensity, real *conc, integer *nsd, logical *
	lcentrif, real *radius)
{
    /* System generated locals */
    integer conc_dim1, conc_offset, i__1;
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static doublereal b;
    static integer i__;
    static doublereal a2, a3, f2;
    static real dx, fre, dxb;
    extern doublereal fqh_(real *, real *, real *);
    static real dxa;
    extern doublereal fro_(integer *, real *);
    static real cona, conb, grav;
    static logical lgeom;
    static real conta, contb;
    extern doublereal fqdrain_(real *, real *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --sinkim;
    --thvold;
    --thvnew;
    --convh;
    --convt;
    --conlt;
    --temp;
    --theq;
    --thold;
    --thnew;
    --s;
    --r__;
    --p;
    --sink;
    --hnew;
    --cap;
    --con;
    --hold;
    --x;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    lgeom = FALSE_;
/* Arithmetic average (false), geometric average (t */
    fre = 1.f;
    grav = *cosalf;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	thnew[i__] = thold[i__] + (theq[i__] - thold[i__]) * *rate;
/* L10: */
    }
/*     Finite differences */
/*     Bottom BC */
    dxb = x[2] - x[1];
    dx = dxb / 2.f;
    conb = (con[1] + con[2]) / 2.f;
/* Arithmetic average */
    if (lgeom) {
	d__1 = (doublereal) (con[1] * con[2]);
	conb = pow_dd(&d__1, &c_b20);
    }
/* Geometric average */
    if (*lcentrif) {
	grav = *cosalf * (*radius + (r__1 = (x[1] + x[2]) / 2.f, dabs(r__1)));
    }
    b = conb * grav;
    if (*ldensity) {
	r__1 = (conc[conc_dim1 + 1] + conc[(conc_dim1 << 1) + 1]) / 2.f;
	fre = fro_(&c__1, &r__1);
	b = conb * grav * fre;
	fre = fro_(&c__1, &conc[conc_dim1 + 1]);
    }
    if (*lvapor) {
	conb += (convh[1] + convh[2]) / 2.f;
    }
    s[1] = -conb / dxb;
    if (*freed) {
	*rbot = -conb * grav * fre;
    }
    f2 = cap[1] * dx / *dt * fre * *rate;
    *rb = conb / dxb + f2;
    *sb = -conb / dxb;
    if (*qgwlf) {
	r__1 = hnew[1] - *gwl0l;
	*rbot = fqh_(&r__1, aqh, bqh);
    }
    if (*qdrain) {
	r__1 = x[1] + hnew[1];
	*rbot = fqdrain_(&r__1, zbotdr, basegw, rspacing, iposdr, rkhtop, 
		rkhbot, rkvtop, rkvbot, entres, wetper, zintf, geofac);
    }
    *pb = b - sink[1] * dx + f2 * hnew[1] - (thnew[1] - thold[1]) * dx / *dt *
	     fre + *rbot;
    if (*idualpor > 0) {
	*pb -= sinkim[1] * dx;
    }
    if (*lvapor || *lwtdep) {
	contb = 0.f;
	if (*lvapor) {
	    contb += (convt[1] + convt[2]) / 2.f;
	}
	if (*lwtdep) {
	    contb += (conlt[1] + conlt[2]) / 2.f;
	}
	*pb = *pb + contb * (temp[2] - temp[1]) / dxb - (thvnew[1] - thvold[1]
		) * dx / *dt;
    }
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dxa = x[i__] - x[i__ - 1];
	dxb = x[i__ + 1] - x[i__];
	dx = (dxa + dxb) / 2.f;
	cona = (con[i__] + con[i__ - 1]) / 2.f;
	conb = (con[i__] + con[i__ + 1]) / 2.f;
	if (lgeom) {
	    d__1 = (doublereal) (con[i__] * con[i__ - 1]);
	    cona = pow_dd(&d__1, &c_b20);
	}
	if (lgeom) {
	    d__1 = (doublereal) (con[i__] * con[i__ + 1]);
	    conb = pow_dd(&d__1, &c_b20);
	}
	if (*lcentrif) {
	    grav = *cosalf * (*radius + (r__1 = x[i__], dabs(r__1)));
	}
	b = (cona - conb) * grav;
	if (*ldensity) {
	    r__1 = (conc[i__ * conc_dim1 + 1] + conc[(i__ - 1) * conc_dim1 + 
		    1]) / 2.f;
	    r__2 = (conc[i__ * conc_dim1 + 1] + conc[(i__ + 1) * conc_dim1 + 
		    1]) / 2.f;
	    b = (cona * fro_(&c__1, &r__1) - conb * fro_(&c__1, &r__2)) * 
		    grav;
	    fre = fro_(&c__1, &conc[i__ * conc_dim1 + 1]);
	}
	if (*lcentrif) {
	    b += *cosalf * con[i__] * dx;
	}
	if (*lvapor) {
	    cona += (convh[i__] + convh[i__ - 1]) / 2.f;
	    conb += (convh[i__] + convh[i__ + 1]) / 2.f;
	}
	a2 = cona / dxa + conb / dxb;
	a3 = -conb / dxb;
	f2 = cap[i__] * dx / *dt * fre * *rate;
	r__[i__] = a2 + f2;
	p[i__] = f2 * hnew[i__] - (thnew[i__] - thold[i__]) * dx / *dt * fre 
		- b - sink[i__] * dx;
	if (*idualpor > 0) {
	    p[i__] -= sinkim[i__] * dx;
	}
	s[i__] = a3;
	if (*lvapor || *lwtdep) {
	    conta = 0.f;
	    if (*lvapor) {
		conta += (convt[i__] + convt[i__ - 1]) / 2.f;
	    }
	    if (*lwtdep) {
		conta += (conlt[i__] + conlt[i__ - 1]) / 2.f;
	    }
	    contb = 0.f;
	    if (*lvapor) {
		contb += (convt[i__] + convt[i__ + 1]) / 2.f;
	    }
	    if (*lwtdep) {
		contb += (conlt[i__] + conlt[i__ + 1]) / 2.f;
	    }
	    p[i__] = p[i__] + contb * (temp[i__ + 1] - temp[i__]) / dxb - 
		    conta * (temp[i__] - temp[i__ - 1]) / dxa - (thvnew[i__] 
		    - thvold[i__]) * dx / *dt;
	}
/* L11: */
    }
/*     Top BC */
    dxa = x[*n] - x[*n - 1];
    dx = dxa / 2.f;
    cona = (con[*n] + con[*n - 1]) / 2.f;
    if (lgeom) {
	d__1 = (doublereal) (con[*n] * con[*n - 1]);
	cona = pow_dd(&d__1, &c_b20);
    }
    if (*lcentrif) {
	grav = *cosalf * (*radius + (r__1 = (x[*n] + x[*n - 1]) / 2.f, dabs(
		r__1)));
    }
    b = cona * grav;
    if (*lvapor) {
	cona += (convh[*n] + convh[*n - 1]) / 2.f;
    }
    if (*ldensity) {
	r__1 = (conc[*n * conc_dim1 + 1] + conc[(*n - 1) * conc_dim1 + 1]) / 
		2.f;
	b = cona * grav * fro_(&c__1, &r__1);
	fre = fro_(&c__1, &conc[*n * conc_dim1 + 1]);
    }
    f2 = cap[*n] * dx / *dt * fre * *rate;
    *rt = cona / dxa + f2;
    *st = -cona / dxa;
    *pt = f2 * hnew[*n] - (thnew[*n] - thold[*n]) * dx / *dt * fre - sink[*n] 
	    * dx - b;
    if (*idualpor > 0) {
	*pt -= sinkim[*n] * dx;
    }
    if (*lvapor || *lwtdep) {
	conta = 0.f;
	if (*lvapor) {
	    conta += (convt[*n] + convt[*n - 1]) / 2.f;
	}
	if (*lwtdep) {
	    conta += (conlt[*n] + conlt[*n - 1]) / 2.f;
	}
	*pt = *pt - conta * (temp[*n] - temp[*n - 1]) / dxa - (thvnew[*n] - 
		thvold[*n]) * dx / *dt;
    }
    *vtop = -((real) (*st)) * hnew[*n - 1] - (real) (*rt) * hnew[*n] + (real) 
	    (*pt);
    *pt -= *rtop;
    if (*wlayer) {
	if (hnew[*n] > 0.f) {
	    *rt += 1.f / *dt;
/* Computing MAX */
	    r__1 = hold[*n];
	    *pt += dmax(r__1,0.f) / *dt;
	} else {
/* Computing MAX */
	    r__1 = hold[*n];
	    *pt += dmax(r__1,0.f) / *dt;
	}
    }
    return 0;
} /* reset_ */

/* *********************************************************************** */
/* Subroutine */ int gauss_(integer *n, integer *kodtop, integer *kodbot, 
	real *htop, real *hbot, real *hnew, doublereal *p, doublereal *r__, 
	doublereal *s, doublereal *pb, doublereal *rb, doublereal *sb, 
	doublereal *pt, doublereal *rt, doublereal *st, doublereal *rmin)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;

/*     Forward */
    /* Parameter adjustments */
    --s;
    --r__;
    --p;
    --hnew;

    /* Function Body */
    if (*kodbot >= 0) {
	p[2] -= s[1] * *hbot;
    } else {
	if (abs(*rb) < *rmin) {
	    *rb = *rmin;
	}
	p[2] -= *pb * s[1] / *rb;
	r__[2] -= *sb * s[1] / *rb;
    }
    i__1 = *n - 1;
    for (i__ = 3; i__ <= i__1; ++i__) {
	if ((d__1 = r__[i__ - 1], abs(d__1)) < *rmin) {
	    r__[i__ - 1] = *rmin;
	}
	p[i__] -= p[i__ - 1] * s[i__ - 1] / r__[i__ - 1];
	r__[i__] -= s[i__ - 1] * s[i__ - 1] / r__[i__ - 1];
/* L11: */
    }
    if (*kodtop > 0) {
	p[*n - 1] -= s[*n - 1] * *htop;
    } else {
	if ((d__1 = r__[*n - 1], abs(d__1)) < *rmin) {
	    r__[*n - 1] = *rmin;
	}
	p[*n] = *pt - p[*n - 1] * *st / r__[*n - 1];
	r__[*n] = *rt - s[*n - 1] * *st / r__[*n - 1];
    }
/*     Back */
    if ((d__1 = r__[*n - 1], abs(d__1)) < *rmin) {
	r__[*n - 1] = *rmin;
    }
    if (*kodtop > 0) {
	hnew[*n] = *htop;
	hnew[*n - 1] = (real) (p[*n - 1] / r__[*n - 1]);
    } else {
	hnew[*n] = (real) (p[*n] / r__[*n]);
	hnew[*n - 1] = (real) ((p[*n - 1] - s[*n - 1] * hnew[*n]) / r__[*n - 
		1]);
    }
    for (i__ = *n - 2; i__ >= 2; --i__) {
	if ((d__1 = r__[i__], abs(d__1)) < *rmin) {
	    r__[i__] = *rmin;
	}
	hnew[i__] = (real) ((p[i__] - s[i__] * hnew[i__ + 1]) / r__[i__]);
/* L12: */
    }
    if (*kodbot >= 0) {
	hnew[1] = *hbot;
    } else {
	if (abs(*rb) < *rmin) {
	    *rb = *rmin;
	}
	hnew[1] = (real) ((*pb - *sb * hnew[2]) / *rb);
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L13: */
    }
    return 0;
} /* gauss_ */

/* *********************************************************************** */
/* Subroutine */ int shift_(integer *n, integer *kodtop, real *rtop, real *
	rbot, real *htop, real *hbot, real *hcrita, real *cosalf, logical *
	wlayer, real *con, real *hnew, real *x, logical *topinf, integer *
	kodbot, logical *seepf, real *thnew, real *thold, real *sink, real *
	dt, logical *lvapor, logical *lwtdep, real *conlt, real *temp, real *
	convh, real *convt, integer *idualpor, real *sinkim, logical *
	ldensity, real *conc, integer *nsd, logical *lcentrif, real *radius, 
	real *hseep)
{
    /* System generated locals */
    integer conc_dim1, conc_offset;
    real r__1;

    /* Local variables */
    static integer m;
    static real dx, fre;
    extern doublereal fro_(integer *, real *);
    static real grav, vbot, vtop;

    /* Parameter adjustments */
    --sinkim;
    --convt;
    --convh;
    --temp;
    --conlt;
    --sink;
    --thold;
    --thnew;
    --x;
    --hnew;
    --con;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    fre = 1.f;
    grav = *cosalf;
/*     Seepage face at the bottom */
    if (*seepf) {
	dx = x[2] - x[1];
	if (*ldensity) {
	    fre = fro_(&c__1, &conc[conc_dim1 + 1]);
	}
	if (*lcentrif) {
	    grav = *cosalf * (*radius + (r__1 = (x[2] + x[1]) / 2.f, dabs(
		    r__1)));
	}
	vbot = -(con[1] + con[2]) / 2.f * ((hnew[2] - hnew[1]) / dx + grav * 
		fre) - dx / 2.f * fre * ((thnew[1] - thold[1]) / *dt + sink[1]
		);
	if (*kodbot >= 0) {
	    if (vbot > 0.f) {
		*kodbot = -2;
		*rbot = 0.f;
	    }
	} else {
	    if (hnew[1] >= *hseep) {
		*kodbot = 2;
		*hbot = *hseep;
	    }
	}
    }
/*     Atmospheric boundary condition */
    if (*topinf && (abs(*kodtop) == 4 || abs(*kodtop) == 1 && *rtop > 0.f)) {
	if (*kodtop > 0) {
	    m = *n - 1;
	    dx = x[*n] - x[m];
	    if (*ldensity) {
		fre = fro_(&c__1, &conc[*n * conc_dim1 + 1]);
	    }
	    if (*lcentrif) {
		grav = *cosalf * (*radius + (r__1 = (x[*n] + x[m]) / 2.f, 
			dabs(r__1)));
	    }
	    vtop = -(con[*n] + con[m]) / 2.f * ((hnew[*n] - hnew[m]) / dx + 
		    grav * fre) - (thnew[*n] - thold[*n]) * fre * dx / 2.f / *
		    dt - sink[*n] * dx / 2.f;
	    if (*idualpor > 0) {
		vtop -= sinkim[*n] * dx / 2.f;
	    }
	    if (*lwtdep) {
		vtop -= (conlt[*n] + conlt[m]) / 2.f * (temp[*n] - temp[m]) / 
			dx;
	    }
	    if (*lvapor) {
		vtop = vtop - (convh[*n] + convh[m]) / 2.f * (hnew[*n] - hnew[
			m]) / dx - (convt[*n] + convt[m]) / 2.f * (temp[*n] - 
			temp[m]) / dx;
	    }
	    if (dabs(vtop) > dabs(*rtop) || vtop * *rtop <= 0.f) {
		if (abs(*kodtop) == 4) {
		    *kodtop = -4;
		}
	    }
	    if (*kodtop == 4 && hnew[*n] <= *hcrita * .99f && *rtop < 0.f) {
		*kodtop = -4;
	    }
	} else {
	    if (! (*wlayer)) {
		if (hnew[*n] > 0.f) {
		    if (abs(*kodtop) == 4) {
			*kodtop = 4;
		    }
		    if (abs(*kodtop) == 1) {
			*kodtop = 1;
		    }
		    *htop = 0.f;
		}
	    }
	    if (hnew[*n] <= *hcrita) {
		if (abs(*kodtop) == 4) {
		    *kodtop = 4;
		}
		if (abs(*kodtop) == 1) {
		    *kodtop = 1;
		}
		*htop = *hcrita;
	    }
	}
    }
    return 0;
} /* shift_ */

/* *********************************************************************** */
/* Subroutine */ int setmat_(integer *numnp, integer *ntab, integer *ntabd, 
	integer *nmat, real *htab, real *contab, real *captab, real *hnew, 
	integer *matnum, real *pard, real *con, real *cap, real *consat, real 
	*ah, real *ak, real *ath, real *hsat, real *htemp, real *thetab, real 
	*theta, real *thr, real *ths, logical *lwtdep, real *temp, integer *
	iter, real *cono, integer *kappa, real *aths, real *thrr, real *conr, 
	real *aks, real *ahw, real *athw, real *akw, integer *imodel, logical 
	*ltable, logical *lvapor, real *thetav, real *conlt, real *convt, 
	real *convh, real *xconv, real *tconv, integer *ntabmod, real *hcrita,
	 logical *ldensity, real *conc, integer *nsd, integer *ienhanc)
{
    /* System generated locals */
    integer htab_dim1, htab_offset, contab_dim1, contab_offset, captab_dim1, 
	    captab_offset, thetab_dim1, thetab_offset, conc_dim1, conc_offset,
	     i__1, i__2;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    extern /* Subroutine */ int convapor_(integer *, integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, logical *, real *, integer *);
    static integer i__, j, m;
    extern doublereal fc_(integer *, real *, real *);
    static real dh;
    extern doublereal fk_(integer *, real *, real *);
    static real at, bt;
    extern doublereal fq_(integer *, real *, real *);
    static integer it;
    static real hi1, hi2, dlh, him;
    extern doublereal fro_(integer *, real *);
    static real alh1, capi, coni, thei;
    extern /* Subroutine */ int vaporcontent_(integer *, integer *, real *, 
	    real *, real *, real *, integer *, real *, real *);
    static real tempr;

    /* Parameter adjustments */
    --convh;
    --convt;
    --conlt;
    --thetav;
    --aks;
    --conr;
    --thrr;
    --aths;
    --kappa;
    --cono;
    --temp;
    --theta;
    --htemp;
    --ath;
    --ak;
    --ah;
    --cap;
    --con;
    --matnum;
    --hnew;
    --akw;
    --athw;
    --ahw;
    --ths;
    --thr;
    thetab_dim1 = *ntabd;
    thetab_offset = 1 + thetab_dim1;
    thetab -= thetab_offset;
    --hsat;
    --consat;
    pard -= 12;
    captab_dim1 = *ntabd;
    captab_offset = 1 + captab_dim1;
    captab -= captab_offset;
    contab_dim1 = *ntabd;
    contab_offset = 1 + contab_dim1;
    contab -= contab_offset;
    htab_dim1 = *ntabd;
    htab_offset = 1 + htab_dim1;
    htab -= htab_offset;
    --ntab;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    if (*imodel < *ntabmod) {
	r__1 = -htab[htab_dim1 + 1];
	alh1 = r_lg10(&r__1);
	r__1 = -htab[ntab[1] + htab_dim1];
	dlh = (r_lg10(&r__1) - alh1) / (ntab[1] - 1);
    }
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	at = 1.f;
	bt = 1.f;
	if (*lwtdep) {
/* Temperature dependence */
	    tempr = 20.f;
/* Computing 2nd power */
	    r__1 = temp[i__];
/* Computing 2nd power */
	    r__2 = tempr;
	    at = (75.6f - temp[i__] * .1425f - r__1 * r__1 * 2.38e-4f) / (
		    75.6f - tempr * .1425f - r__2 * r__2 * 2.38e-4f);
/* Su */
/* Computing 2nd power */
	    r__1 = temp[i__] - 4.f;
/* Computing 3rd power */
	    r__2 = temp[i__] - 4.f;
/* Computing 2nd power */
	    r__3 = tempr - 4.f;
/* Computing 3rd power */
	    r__4 = tempr - 4.f;
	    bt = (1.787f - tempr * .007f) / (tempr * .03225f + 1.f) / ((
		    1.787f - temp[i__] * .007f) / (temp[i__] * .03225f + 1.f))
		     * (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 
		    3.79e-8f) / (1.f - r__3 * r__3 * 7.37e-6f + r__4 * (r__4 *
		     r__4) * 3.79e-8f);
/* Dy */
/* De */
	}
	if (*ldensity) {
	    bt = bt * fro_(&c__1, &conc[i__ * conc_dim1 + 1]) / fro_(&c__2, &
		    conc[i__ * conc_dim1 + 1]);
	}
/* Bu */
	m = matnum[i__];
	if (kappa[i__] == -1) {
/* Computing MIN */
	    r__1 = hsat[m], r__2 = htemp[i__] / ah[i__] / at;
	    hi1 = dmin(r__1,r__2);
/* Computing MIN */
	    r__1 = hsat[m], r__2 = hnew[i__] / ah[i__] / at;
	    hi2 = dmin(r__1,r__2);
	} else if (kappa[i__] == 1) {
/* Computing MIN */
	    r__1 = hsat[m], r__2 = htemp[i__] / ah[i__] / ahw[m] / at;
	    hi1 = dmin(r__1,r__2);
/* Computing MIN */
	    r__1 = hsat[m], r__2 = hnew[i__] / ah[i__] / ahw[m] / at;
	    hi2 = dmin(r__1,r__2);
	}
	him = hi1 * .1f + hi2 * .9f;
	if (*imodel < *ntabmod) {
/* Conductivity */
	    if (hi1 >= hsat[m] && hi2 >= hsat[m]) {
		coni = consat[m];
	    } else if (him > htab[ntab[1] + htab_dim1] && him <= htab[
		    htab_dim1 + 1] && *ltable) {
		r__1 = -him;
		it = (integer) ((r_lg10(&r__1) - alh1) / dlh) + 1;
		dh = (him - htab[it + htab_dim1]) / (htab[it + 1 + htab_dim1] 
			- htab[it + htab_dim1]);
		coni = contab[it + m * contab_dim1] + (contab[it + 1 + m * 
			contab_dim1] - contab[it + m * contab_dim1]) * dh;
	    } else {
		coni = fk_(imodel, &him, &pard[m * 11 + 1]);
	    }
	} else if (*imodel == *ntabmod) {
/* Tables */
	    if (hi1 >= hsat[m] && hi2 >= hsat[m]) {
		coni = consat[m];
	    } else if (him >= htab[ntab[m] + m * htab_dim1] && him <= htab[m *
		     htab_dim1 + 1]) {
		it = 1;
		i__2 = ntab[m] - 1;
		for (j = 1; j <= i__2; ++j) {
		    if (him >= htab[j + 1 + m * htab_dim1] && him < htab[j + 
			    m * htab_dim1]) {
			it = j;
		    }
/* L12: */
		}
		dh = (him - htab[it + m * htab_dim1]) / (htab[it + 1 + m * 
			htab_dim1] - htab[it + m * htab_dim1]);
		coni = contab[it + m * contab_dim1] + (contab[it + 1 + m * 
			contab_dim1] - contab[it + m * contab_dim1]) * dh;
	    } else {
		if (him > htab[m * htab_dim1 + 1]) {
		    coni = contab[m * contab_dim1 + 1] + (consat[m] - contab[
			    m * contab_dim1 + 1]) * (htab[m * htab_dim1 + 1] 
			    - hi2) / htab[m * htab_dim1 + 1];
		} else if (him < htab[ntab[m] + m * htab_dim1]) {
		    r__1 = -him;
		    r__2 = -htab[ntab[m] + m * htab_dim1];
		    r__3 = -htab[ntab[m] + m * htab_dim1];
		    coni = contab[ntab[m] + m * contab_dim1] - contab[ntab[m] 
			    + m * contab_dim1] * (r_lg10(&r__1) - r_lg10(&
			    r__2)) / (10.f - r_lg10(&r__3));
		}
	    }
	}
	if (*imodel < *ntabmod) {
/* Capacity and water content */
	    if (him >= hsat[m]) {
		capi = 0.f;
		thei = ths[m];
	    } else if (him >= htab[ntab[1] + htab_dim1] && him <= htab[
		    htab_dim1 + 1] && *ltable) {
		r__1 = -him;
		it = (integer) ((r_lg10(&r__1) - alh1) / dlh) + 1;
		dh = (him - htab[it + htab_dim1]) / (htab[it + 1 + htab_dim1] 
			- htab[it + htab_dim1]);
		capi = captab[it + m * captab_dim1] + (captab[it + 1 + m * 
			captab_dim1] - captab[it + m * captab_dim1]) * dh;
		thei = thetab[it + m * thetab_dim1] + (thetab[it + 1 + m * 
			thetab_dim1] - thetab[it + m * thetab_dim1]) * dh;
	    } else {
		capi = fc_(imodel, &him, &pard[m * 11 + 1]);
		thei = fq_(imodel, &him, &pard[m * 11 + 1]);
	    }
	} else if (*imodel == *ntabmod) {
/* Tables */
	    if (hi2 >= hsat[m]) {
		capi = 0.f;
		thei = ths[m];
	    } else if (hi2 >= htab[ntab[m] + m * htab_dim1] && hi2 <= htab[m *
		     htab_dim1 + 1]) {
		it = 1;
		i__2 = ntab[m] - 1;
		for (j = 1; j <= i__2; ++j) {
		    if (hi2 >= htab[j + 1 + m * htab_dim1] && hi2 <= htab[j + 
			    m * htab_dim1]) {
			it = j;
		    }
/* L13: */
		}
		dh = (hi2 - htab[it + m * htab_dim1]) / (htab[it + 1 + m * 
			htab_dim1] - htab[it + m * htab_dim1]);
		capi = captab[it + m * captab_dim1] + (captab[it + 1 + m * 
			captab_dim1] - captab[it + m * captab_dim1]) * dh;
		thei = thetab[it + m * thetab_dim1] + (thetab[it + 1 + m * 
			thetab_dim1] - thetab[it + m * thetab_dim1]) * dh;
	    } else {
		if (hi2 > htab[m * htab_dim1 + 1]) {
		    capi = captab[m * captab_dim1 + 1] * hi2 / htab[m * 
			    htab_dim1 + 1];
		    thei = thetab[m * thetab_dim1 + 1] + (ths[m] - thetab[m * 
			    thetab_dim1 + 1]) * (htab[m * htab_dim1 + 1] - 
			    hi2) / htab[m * htab_dim1 + 1];
		} else if (hi2 < htab[ntab[m] + m * htab_dim1]) {
		    r__1 = -hi2;
		    r__2 = -htab[ntab[m] + m * htab_dim1];
		    r__3 = -htab[ntab[m] + m * htab_dim1];
		    capi = captab[ntab[m] + m * captab_dim1] - captab[ntab[m] 
			    + m * captab_dim1] * (r_lg10(&r__1) - r_lg10(&
			    r__2)) / (6.f - r_lg10(&r__3));
		    thei = thr[m] + (thetab[ntab[m] + m * thetab_dim1] - thr[
			    m]) * (hi2 + 1e6f) / (htab[ntab[m] + m * 
			    htab_dim1] + 1e6f);
		}
	    }
	}
	if (kappa[i__] == -1) {
/* Drying */
	    con[i__] = coni * ak[i__] * bt * aks[i__];
	    cap[i__] = capi * ath[i__] * aths[i__] / ah[i__] / at;
	    theta[i__] = thr[m] + (thei - thr[m]) * ath[i__] * aths[i__];
	} else {
/* Wetting */
	    con[i__] = conr[i__] + coni * ak[i__] * bt * aks[i__] * akw[m];
	    cap[i__] = capi * ath[i__] * aths[i__] * athw[m] / ah[i__] / ahw[
		    m] / at;
	    theta[i__] = thrr[i__] + athw[m] * ath[i__] * aths[i__] * (thei - 
		    thr[m]);
	}
	if (*iter == 0) {
	    cono[i__] = con[i__];
	}
/* L11: */
    }
    if (*lvapor || *lwtdep) {
/* Vapor flow */
	convapor_(numnp, nmat, &matnum[1], &hnew[1], &temp[1], &con[1], &
		theta[1], &ths[1], &conlt[1], &convt[1], &convh[1], xconv, 
		tconv, lvapor, hcrita, ienhanc);
	if (*lvapor) {
	    vaporcontent_(numnp, nmat, &theta[1], &thetav[1], &temp[1], &hnew[
		    1], &matnum[1], &ths[1], xconv);
	}
    }
    return 0;
} /* setmat_ */

/* *********************************************************************** */
doublereal fqh_(real *gwl, real *aqh, real *bqh)
{
    /* System generated locals */
    real ret_val;

    ret_val = *aqh * exp(*bqh * dabs(*gwl));
    return ret_val;
} /* fqh_ */

/* *********************************************************************** */
doublereal fqdrain_(real *gwl, real *zbotdr, real *basegw, real *rspacing, 
	integer *iposdr, real *khtop, real *khbot, real *kvtop, real *kvbot, 
	real *entres, real *wetper, real *zintf, real *geofac)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static integer i__;
    static real x, dh, fx, eqd, dbot, rrad, rhor, rver, zimp, totres;

/*    ------------------------------------------------------------------- */
/*     Purpose: determines the drainage flux */
/*     Based on the SWAP model by van Dam et al. 1997 */
/*     iPosDr   kod for the position of the drain.......................I */
/*              =1: Homogeneous profile, drain on top of impervious layer */
/*              =2: Homogeneous profile, drain above impervious layer */
/*              =3: Heterogeneous profile, drain at interface between */
/*                                         both soil layers */
/*              =4: Heterogeneous profile, drain in bottom layer */
/*              =5: Heterogeneous profile, drain in top layer */
/*     Input    Calculation: GWL,Pond */
/*              Common: zBotDr,rSpacing,Entres */
/*              iPosDr=1: KhTop */
/*              iPosDr=2: BaseGW,KhTop,WetPer */
/*              iPosDr=3: BaseGW,KhTop,KhBot,WetPer */
/*              iPosDr=4: BaseGW,KvTop,KvBot,KhBot,WetPer,zInTF */
/*              iPosDr=5: BaseGW,KhTop,KvTop,KhBot,WetPer,zInTF,GeoFac */
/*     GWL      ground water level...................................I(s) */
/*     zBotDr   coordinate of the bottom of the drainage................I */
/*     Pond     Ponding (cm).........................................I(s) */
/*     BaseGW   coordinate of the impervious layer cm (2,3,4,5).........I */
/*     rSpacing drain spacing (1,2,3,4,5)...............................I */
/*     KhTop    horizontal saturated hydraulic conductivity above drain */
/*              (cm/d) (1,2,3,5)........................................I */
/*     KhBot    horizontal saturated hydraulic conductivity below drain */
/*              (cm/d) (3,4,5)..........................................I */
/*     KvTop    vertical saturated hydraulic conductivity above drain */
/*              (cm/d) (4,5)............................................I */
/*     KvBot    vertical saturated hydraulic conductivity below drain */
/*              (cm/d) (4)..............................................I */
/*     Entres   entrance resistance into the drain and/or ditches (d) */
/*              (1,2,3,4,5).............................................I */
/*     WetPer   wet perimeter (cm) (2,3,4,5)............................I */
/*     zInTF    level of the transition between the upper and lower */
/*              soil layer (cm) (4,5)...................................I */
/*     GeoFac   geometry factor (5).....................................I */
/*     FqDrain  drainage flux (cm/d)....................................O */
/*     dh       hydraulic difference between drain and the middle of */
/*              the spacing */
/*     zImp     adjusted coordinate of the impervious layer cm (2,3,4,5) */
/*     dBot     depth to the impervious layer below the drain */
/*     EqD      equilvalent depth (cm) */
/*     TotRes   total drainage resistance */
/*     RVer     vertical drainage resistance */
/*     RHor     horizontal drainage resistance */
/*     RRad     radial drainage resistance */
/*     x        typical length variable */
/*     ------------------------------------------------------------------ */
/*     global variables */
/*     local variables */
/*     drainage flux calculated according to Hooghoudt or Ernst */
    dh = *gwl - *zbotdr;
/*     contributing layer below drains limited to 1/4 L */
    if (*iposdr > 1) {
/* Computing MAX */
	r__1 = *basegw, r__2 = *zbotdr - *rspacing * .25f;
	zimp = dmax(r__1,r__2);
	dbot = *zbotdr - zimp;
	if (dbot < 0.f) {
	    s_stop("Error - Bocodrb: dBot negative", (ftnlen)30);
	}
    }
/*     no infiltration allowed */
    if (dh < 1e-10f) {
	ret_val = 0.f;
	return ret_val;
    }
/*     case 1: homogeneous, on top of impervious layer */
    if (*iposdr == 1) {
/*     calculation of drainage resistance and drainage flux */
	totres = *rspacing * *rspacing / (*khtop * 4.f * dabs(dh)) + *entres;
/*     case 2,3: in homogeneous profile or at interface of 2 layers */
    } else if (*iposdr == 2 || *iposdr == 3) {
/*       calculation of equivalent depth */
	x = dbot * 6.2831799999999998f / *rspacing;
	if (x > .5f) {
	    fx = 0.f;
	    for (i__ = 1; i__ <= 5; i__ += 2) {
		fx += exp(i__ * -2.f * x) * 4.f / (i__ * (1.f - exp(i__ * 
			-2.f * x)));
/* L10: */
	    }
	    eqd = *rspacing * 3.14159f / 8.f / (log(*rspacing / *wetper) + fx)
		    ;
	} else {
	    if (x < 1e-6f) {
		eqd = dbot;
	    } else {
		fx = 9.8695877280999991f / (x * 4.f) + log(x / 
			6.2831799999999998f);
		eqd = *rspacing * 3.14159f / 8.f / (log(*rspacing / *wetper) 
			+ fx);
	    }
	}
	if (eqd > dbot) {
	    eqd = dbot;
	}
/*       calculation of drainage resistance & drainage flux */
	if (*iposdr == 2) {
	    totres = *rspacing * *rspacing / (*khtop * 8.f * eqd + *khtop * 4 
		    * dabs(dh)) + *entres;
	} else if (*iposdr == 3) {
	    totres = *rspacing * *rspacing / (*khbot * 8.f * eqd + *khtop * 4 
		    * dabs(dh)) + *entres;
	}
/*     case 4: drain in bottom layer */
    } else if (*iposdr == 4) {
	if (*zbotdr > *zintf) {
	    s_stop("Error - check zInTF and zBotDr", (ftnlen)30);
	}
/* Computing MAX */
	r__1 = *gwl - *zintf;
	rver = dmax(r__1,0.f) / *kvtop + (dmin(*zintf,*gwl) - *zbotdr) / *
		kvbot;
	rhor = *rspacing * *rspacing / (*khbot * 8.f * dbot);
	rrad = *rspacing / (sqrt(*khbot * *kvbot) * 3.14159f) * log(dbot / *
		wetper);
	totres = rver + rhor + rrad + *entres;
/*     case 5: drain in top layer */
    } else if (*iposdr == 5) {
	if (*zbotdr < *zintf) {
	    s_stop("Error - check zInTF and zBotDr", (ftnlen)30);
	}
	rver = (*gwl - *zbotdr) / *kvtop;
	rhor = *rspacing * *rspacing / (*khtop * 8.f * (*zbotdr - *zintf) + *
		khbot * 8.f * (*zintf - zimp));
	rrad = *rspacing / (sqrt(*khtop * *kvtop) * 3.14159f) * log(*geofac * 
		(*zbotdr - *zintf) / *wetper);
	totres = rver + rhor + rrad + *entres;
    }
    ret_val = -dh / totres;
    return ret_val;
} /* fqdrain_ */

/* *********************************************************************** */
/*     To calculate the velocities */
/* Subroutine */ int veloc_(integer *n, real *hnew, real *con, real *x, real *
	cosalf, real *v, real *thnew, real *thold, real *sink, real *dt, 
	logical *lvapor, logical *lwtdep, real *conlt, real *convt, real *
	convh, real *temp, real *vv, real *thvnew, real *thvold, logical *
	ldensity, real *conc, integer *nsd, logical *lcentrif, real *radius)
{
    /* System generated locals */
    integer conc_dim1, conc_offset, i__1;
    real r__1;

    /* Local variables */
    static integer i__, m;
    static real va, vb, dx1, fre, dxa, dxb;
    extern doublereal fro_(integer *, real *);
    static real dxn, vta, vtb, vva, vvb, grav;

    /* Parameter adjustments */
    --thvold;
    --thvnew;
    --vv;
    --temp;
    --convh;
    --convt;
    --conlt;
    --sink;
    --thold;
    --thnew;
    --v;
    --x;
    --con;
    --hnew;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    fre = 1.f;
    grav = *cosalf;
    m = *n - 1;
    dxn = x[*n] - x[m];
    if (*ldensity) {
	fre = (fro_(&c__1, &conc[*n * conc_dim1 + 1]) + fro_(&c__1, &conc[m * 
		conc_dim1 + 1])) / 2.f;
    }
    if (*lcentrif) {
	grav = *cosalf * (*radius + (r__1 = (x[*n] + x[m]) / 2.f, dabs(r__1)))
		;
    }
    v[*n] = -(con[*n] + con[m]) / 2.f * ((hnew[*n] - hnew[m]) / dxn + fre * 
	    grav) - dxn / 2.f * (fre * (thnew[*n] - thold[*n]) / *dt + sink[*
	    n]);
    if (*lwtdep) {
	v[*n] -= (conlt[*n] + conlt[m]) / 2.f * (temp[*n] - temp[m]) / dxn;
    }
    vv[*n] = 0.f;
    if (*lvapor) {
	vv[*n] = -(convh[*n] + convh[m]) / 2.f * (hnew[*n] - hnew[m]) / dxn - 
		(convt[*n] + convt[m]) / 2.f * (temp[*n] - temp[m]) / dxn - 
		dxn / 2.f * (thvnew[*n] - thvold[*n]) / *dt;
    }
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dxa = x[i__ + 1] - x[i__];
	dxb = x[i__] - x[i__ - 1];
	if (*ldensity) {
	    fre = (fro_(&c__1, &conc[i__ * conc_dim1 + 1]) + fro_(&c__1, &
		    conc[(i__ + 1) * conc_dim1 + 1])) / 2.f;
	}
	if (*lcentrif) {
	    grav = *cosalf * (*radius + (r__1 = (x[i__ + 1] + x[i__]) / 2.f, 
		    dabs(r__1)));
	}
	va = -(con[i__] + con[i__ + 1]) / 2.f * ((hnew[i__ + 1] - hnew[i__]) /
		 dxa + fre * grav);
	if (*ldensity) {
	    fre = (fro_(&c__1, &conc[i__ * conc_dim1 + 1]) + fro_(&c__1, &
		    conc[(i__ - 1) * conc_dim1 + 1])) / 2.f;
	}
	if (*lcentrif) {
	    grav = *cosalf * (*radius + (r__1 = (x[i__] + x[i__ - 1]) / 2.f, 
		    dabs(r__1)));
	}
	vb = -(con[i__] + con[i__ - 1]) / 2.f * ((hnew[i__] - hnew[i__ - 1]) /
		 dxb + fre * grav);
	v[i__] = (va * dxb + vb * dxa) / (dxa + dxb);
	if (*lwtdep) {
	    vta = -(conlt[i__] + conlt[i__ + 1]) / 2.f * (temp[i__ + 1] - 
		    temp[i__]) / dxa;
	    vtb = -(conlt[i__] + conlt[i__ - 1]) / 2.f * (temp[i__] - temp[
		    i__ - 1]) / dxb;
	    v[i__] += (vta * dxb + vtb * dxa) / (dxa + dxb);
	}
	vv[i__] = 0.f;
	if (*lvapor) {
	    vva = -(convh[i__] + convh[i__ + 1]) / 2.f * (hnew[i__ + 1] - 
		    hnew[i__]) / dxa;
	    vvb = -(convh[i__] + convh[i__ - 1]) / 2.f * (hnew[i__] - hnew[
		    i__ - 1]) / dxb;
	    vva -= (convt[i__] + convt[i__ + 1]) / 2.f * (temp[i__ + 1] - 
		    temp[i__]) / dxa;
	    vvb -= (convt[i__] + convt[i__ - 1]) / 2.f * (temp[i__] - temp[
		    i__ - 1]) / dxb;
	    vv[i__] = (vva * dxb + vvb * dxa) / (dxa + dxb);
	}
/* L11: */
    }
    dx1 = x[2] - x[1];
    if (*ldensity) {
	fre = (fro_(&c__1, &conc[(conc_dim1 << 1) + 1]) + fro_(&c__1, &conc[
		conc_dim1 + 1])) / 2.f;
    }
    if (*lcentrif) {
	grav = *cosalf * (*radius + (r__1 = (x[2] + x[1]) / 2.f, dabs(r__1)));
    }
    v[1] = -(con[1] + con[2]) / 2.f * ((hnew[2] - hnew[1]) / dx1 + fre * grav)
	     + dx1 / 2.f * (fre * (thnew[1] - thold[1]) / *dt + sink[1]);
    if (*lwtdep) {
	v[1] -= (conlt[1] + conlt[2]) / 2.f * (temp[2] - temp[1]) / dx1;
    }
    vv[1] = 0.f;
    if (*lvapor) {
	vv[1] = -(convh[1] + convh[2]) / 2.f * (hnew[2] - hnew[1]) / dx1 - (
		convt[1] + convt[2]) / 2.f * (temp[2] - temp[1]) / dx1 + dx1 /
		 2.f * (thvnew[1] - thvold[1]) / *dt;
    }
    return 0;
} /* veloc_ */

/* *********************************************************************** */
/* Subroutine */ int hyster_(integer *numnp, integer *nmat, real *hold, 
	integer *matnum, real *pard, real *parw, real *thnew, real *thold, 
	integer *kappa, real *aths, real *thrr, real *cono, real *conr, real *
	aks, integer *kappao, real *ah, real *ak, integer *ihyst, integer *
	imodel, real *tolth)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__, m;
    extern doublereal fk_(integer *, real *, real *), fs_(integer *, real *, 
	    real *);
    static real ks, kw, rr, ksd, thr, ths, sew, ksw, thsd, thsw;

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
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       Check for reversal */
	kappao[i__] = kappa[i__];
	if ((thnew[i__] - thold[i__]) * kappa[i__] >= -(*tolth) / 1.f) {
	    goto L11;
	}
	kappa[i__] = -kappa[i__];
	m = matnum[i__];
	thr = pard[m * 11 + 1];
	thsd = pard[m * 11 + 2];
	thsw = parw[m * 11 + 2];
	ksd = pard[m * 11 + 5];
	ksw = parw[m * 11 + 5];
/*       Update Ths and Ks for wetting scanning curve */
	if (kappa[i__] == 1) {
	    if (thsw >= thsd * .999f) {
		ths = thsd;
	    } else {
		rr = 1.f / (thsd - thsw) - 1.f / (thsd - thr);
/* Eq. 8 */
		ths = thsd - (thsd - thold[i__]) / (rr * (thsd - thold[i__]) 
			+ 1.f);
/* Eq. 8 */
	    }
	    if (ksw >= ksd * .999f) {
		ks = ksd;
	    } else {
		rr = 1.f / (ksd - ksw) - 1.f / ksd;
/* Eq. 13 */
		ks = ksd - (ksd - cono[i__]) / (rr * (ksd - cono[i__]) + 1.f);
/* Eq. 13 */
	    }
	}
/*       Update parameters for scanning curve */
	if (kappa[i__] == 1) {
/* Wetting */
	    aths[i__] = 1.f;
	    r__1 = hold[i__] / ah[i__];
	    sew = fs_(imodel, &r__1, &parw[m * 11 + 1]);
	    if (sew < .999f) {
		aths[i__] = (thold[i__] - ths) / (1.f - sew) / (thr - thsw);
	    }

	    thrr[i__] = ths - aths[i__] * (thsw - thr);
/* Eq. 7a */
	    aks[i__] = 1.f;
	    conr[i__] = 0.f;
	    if (*ihyst == 2) {
		r__1 = hold[i__] / ah[i__];
		kw = ak[i__] * fk_(imodel, &r__1, &parw[m * 11 + 1]);
		if (kw < ksw * .999f) {
		    aks[i__] = (cono[i__] - ks) / (kw - ksw);
		}
/* Eq. 12 */
		conr[i__] = ks - aks[i__] * ksw;
/* Eq. 12a */
	    }
	} else {
/* Drying */
	    r__1 = hold[i__] / ah[i__];
	    aths[i__] = (thold[i__] - thr) / fs_(imodel, &r__1, &pard[m * 11 
		    + 1]) / (thsd - thr);
/* Eq. 5 */
	    thrr[i__] = thr;
	    aks[i__] = 1.f;
	    conr[i__] = 0.f;
	    if (*ihyst == 2) {
		r__1 = hold[i__] / ah[i__];
		aks[i__] = cono[i__] / fk_(imodel, &r__1, &pard[m * 11 + 1]) /
			 ak[i__];
	    }
/* Eq. 10 */
	}
L11:
	;
    }
    return 0;
} /* hyster_ */

/* *********************************************************************** */
/*     To calculate isothermal vapor hydraulic conductivity, and */
/*     thermal vapor and liquid hydraulic conductivities */
/* Subroutine */ int convapor_(integer *n, integer *nmat, integer *matnum, 
	real *hnew, real *temp, real *con, real *theta, real *ths, real *
	conlt, real *convt, real *convh, real *xconv, real *tconv, logical *
	lvapor, real *hcrita, integer *ienhanc)
{
    /* Initialized data */

    static real gwt = 7.f;
    static real diff0 = 2.12e-5f;
    static real g = 9.81f;
    static real xmol = .018015f;
    static real r__ = 8.314f;
    static real fc = .02f;
    static real gamma0 = 71.89f;

    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static real h__;
    static integer i__, m;
    static real t, hr, eta, tau, row, diff, rovs, rovs1, gamma, difft, conlh, 
	    tkelv, drovs, tkelv1, dgamma, thetaa, thetas;
    static logical llimit;

/*     ConVh - Conductivity for vapor phase due to gradient of h [m/s] */
/*     ConVT - Conductivity for vapor phase due to gradient of T [m2/s/K] */
/*     ConLT - Conductivity for liquid phase due to gradient of T [m2/s/K] */
/*     Gwt   - Gain factor [-] */
/*     Gamma - surface tension [N/m,J/m2], [g/s2] */
/*     dGamma - derivative of surface tension versus temperature [g/s2/K] */
/*     Diff0 - diffusivity of water vapor in air [m2/s] (2.12e-5) */
/*     Tau   - tortuosity [-] */
/*     row   - density of soil water [kg/m3] */
/*     rovs  - saturated vapor density [kg/m3] */
/*     drovs - derivative of saturated vapor density versus temp [kg/m3/K] */
/*     g     - gravitational acceleration [m/s2] (9.81) */
/*     xMol  - molecular weight of water [kg/mol] (0.018015) */
/*     R     - universal gas constant [J/mol/K] (8.314) */
/*     Hr    - relative humidity [-] */
/*     eta   - enhancement factor [-] */
/*     fc    - mass fraction of clay in soil (0.02) */
/*     T     - temperature [C] */
    /* Parameter adjustments */
    --convh;
    --convt;
    --conlt;
    --theta;
    --con;
    --temp;
    --hnew;
    --matnum;
    --ths;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__ = hnew[i__] / *xconv;
/* Conversion to m */
	conlh = con[i__] / *xconv * *tconv;
/* Conversion to m/s */
	llimit = FALSE_;
	if (hnew[*n] < *hcrita * .99f) {
	    llimit = TRUE_;
	}
	t = temp[i__];
	m = matnum[i__];
	thetas = ths[m];
	gamma = 75.6f - t * .1425f - t * 2.38e-4f * t;
	dgamma = -.1425f - t * 4.79e-4f;
	conlt[i__] = conlh * h__ * gwt * dgamma / gamma0;
	if (*lvapor) {
	    tkelv = t + 273.15f;
/* Computing 2nd power */
	    r__1 = tkelv / 273.15f;
	    difft = diff0 * (r__1 * r__1);
	    thetaa = thetas - theta[i__];
	    d__1 = (doublereal) thetaa;
/* Computing 2nd power */
	    r__1 = thetas;
	    tau = pow_dd(&d__1, &c_b65) / (r__1 * r__1);
/* Millington & Quirk */
	    diff = tau * thetaa * difft;
/* Computing 2nd power */
	    r__1 = t - 4.f;
/* Computing 3rd power */
	    r__2 = t - 4.f;
	    row = (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 
		    3.79e-8f) * 1e3f;
	    rovs = exp(31.3716f - 6014.79f / tkelv - tkelv * .00792495f) * 
		    .001f / tkelv;
	    hr = exp(h__ * xmol * g / r__ / tkelv);
	    if (llimit) {
		hr = 1e-6f;
	    }
/* ##runs fast */
	    convh[i__] = diff / row * rovs * xmol * g / r__ / tkelv * hr;
	    tkelv1 = tkelv + 1.f;
	    rovs1 = exp(31.3716f - 6014.79f / tkelv1 - tkelv1 * .00792495f) * 
		    .001f / tkelv1;
	    drovs = rovs1 - rovs;
	    eta = 1.f;
	    if (*ienhanc == 1) {
/* Computing 4th power */
		r__1 = (2.6f / sqrt(fc) + 1.f) * theta[i__] / thetas, r__1 *= 
			r__1;
		eta = theta[i__] * 3.f / thetas + 9.5f - exp(-(r__1 * r__1)) *
			 8.5f;
	    }
	    convt[i__] = diff / row * eta * hr * drovs;
/*         Conversions to HYDRUS units */
	    convh[i__] = convh[i__] * *xconv / *tconv;
	    convt[i__] = convt[i__] * *xconv * *xconv / *tconv;
	}
	conlt[i__] = conlt[i__] * *xconv * *xconv / *tconv;
/* L11: */
    }
    return 0;
} /* convapor_ */

/* ************************************************************************ */
/* Subroutine */ int vaporcontent_(integer *numnp, integer *nmat, real *theta,
	 real *thetav, real *temp, real *hnew, integer *matnum, real *ths, 
	real *xconv)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static real g, h__;
    static integer i__, m;
    static real r__, hr, rov, row, xmol, rovs, tkelv;

/*     g     - gravitational acceleration [m/s2] (9.81) */
/*     xMol  - molecular weight of water [kg/mol] (0.018015) */
/*     R     - universal gas constant [J/mol/K] (8.314) */
/*     ThetaV- volumetric vapor content expressed as an equivalent water content */
/*     row   - density of soil water [kg/m3] */
/*     rovs  - saturated vapor density [kg/m3] */
/*     rov   - vapor density [kg/m3] */
/*     Hr    - relative humidity [-] */
    /* Parameter adjustments */
    --matnum;
    --hnew;
    --temp;
    --thetav;
    --theta;
    --ths;

    /* Function Body */
    g = 9.81f;
    xmol = .018015f;
    r__ = 8.314f;
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__ = hnew[i__] / *xconv;
/* Conversion to m */
	m = matnum[i__];
	tkelv = temp[i__] + 273.15f;
	rovs = exp(31.3716f - 6014.79f / tkelv - tkelv * .00792495f) * .001f /
		 tkelv;
	hr = exp(h__ * xmol * g / r__ / tkelv);
	rov = rovs * hr;
/* Computing 2nd power */
	r__1 = temp[i__] - 4.f;
/* Computing 3rd power */
	r__2 = temp[i__] - 4.f;
	row = (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 3.79e-8f)
		 * 1e3f;
	thetav[i__] = rov * (ths[m] - theta[i__]) / row;
/* L11: */
    }
    return 0;
} /* vaporcontent_ */

/* *********************************************************************** */
/*     To calculate isothermal vapor hydraulic conductivity */
doublereal convh_(real *hnew, real *theta, real *ths, real *xconv, real *
	tconv)
{
    /* Initialized data */

    static real diff0 = 2.12e-5f;
    static real g = 9.81f;
    static real xmol = .018015f;
    static real r__ = 8.314f;

    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static real h__, hr, tau, row, diff, temp, rovs, difft, tkelv, thetaa;

/*     ConVh - Conductivity for vapor phase due to gradient of h [m/s] */
/*     Diff0 - diffusivity of water vapor in air [m2/s] (2.12e-5) */
/*     Tau   - tortuosity [-] */
/*     row   - density of soil water [kg/m3] */
/*     rovs  - saturated vapor density [kg/m3] */
/*     g     - gravitational acceleration [m/s2] (9.81) */
/*     xMol  - molecular weight of water [kg/mol] (0.018015) */
/*     R     - universal gas constant [J/mol/K] (8.314) */
/*     Hr    - relative humidity [-] */
/*     Temp  - temperature [C] */
    temp = 20.f;
    h__ = *hnew / *xconv;
/* Conversion to m */
    tkelv = temp + 273.15f;
/* Computing 2nd power */
    r__1 = tkelv / 273.15f;
    difft = diff0 * (r__1 * r__1);
    thetaa = *ths - *theta;
    if (thetaa > 0.f) {
	d__1 = (doublereal) thetaa;
/* Computing 2nd power */
	r__1 = *ths;
	tau = pow_dd(&d__1, &c_b65) / (r__1 * r__1);
    }
/* Millington & */
    diff = tau * thetaa * difft;
/* Computing 2nd power */
    r__1 = temp - 4.f;
/* Computing 3rd power */
    r__2 = temp - 4.f;
    row = (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 3.79e-8f) * 
	    1e3f;
    rovs = exp(31.3716f - 6014.79f / tkelv - tkelv * .00792495f) * .001f / 
	    tkelv;
    hr = exp(h__ * xmol * g / r__ / tkelv);
    ret_val = diff / row * rovs * xmol * g / r__ / tkelv * hr;
    ret_val = ret_val * *xconv / *tconv;
/* Conversions to HYDRUS units */
    return ret_val;
} /* convh_ */

/* ************************************************************************ */
doublereal fro_(integer *ikod, real *conc)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real a1, a2, a3, a4;

/*     Ratio of bulk densities (dynamic viscosities) at given and zero */
/*     concentrations */
    ret_val = 1.f;
    if (*ikod == 1) {
/* bulk density */
	a1 = 1.f;
	a2 = .75f;
	a3 = 0.f;
	a4 = 0.f;
    } else if (*ikod == 2) {
/* dynamic viscosity */
	a1 = 1.f;
	a2 = 0.f;
	a3 = 0.f;
	a4 = 0.f;
    }
/* Computing 2nd power */
    r__1 = *conc;
/* Computing 3rd power */
    r__2 = *conc;
    ret_val = a1 + a2 * *conc + a3 * (r__1 * r__1) + a4 * (r__2 * (r__2 * 
	    r__2));
    return ret_val;
} /* fro_ */

/* ************************************************************************ */
/* Subroutine */ int dualpor_(integer *n, integer *nmat, integer *matnum, 
	integer *idualpor, real *thold, real *ths, real *thr, real *thnewim, 
	real *tholdim, real *pard, real *sinkim, real *dt, integer *imodel, 
	real *hnew, real *hcrita, real *x, real *wtransf)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static real h__;
    static integer i__, m;
    extern doublereal fh_(integer *, real *, real *), fk_(integer *, real *, 
	    real *);
    static real se, par[10], seim, condf, condm, deltath, trmaxim;

    /* Parameter adjustments */
    --x;
    --hnew;
    --sinkim;
    --tholdim;
    --thnewim;
    --thold;
    --matnum;
    pard -= 12;
    --thr;
    --ths;

    /* Function Body */
    *wtransf = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
/* Computing MIN */
	r__1 = 1.f, r__2 = (tholdim[i__] - pard[m * 11 + 7]) / (pard[m * 11 + 
		8] - pard[m * 11 + 7]);
	seim = dmin(r__1,r__2);
	if (*idualpor == 1) {
/* Water Content driven */
/* Computing MIN */
	    r__1 = 1.f, r__2 = (thold[i__] - thr[m]) / (ths[m] - thr[m]);
	    se = dmin(r__1,r__2);
	    sinkim[i__] = pard[m * 11 + 9] * (se - seim);
	    deltath = (se - seim) / (ths[m] - thr[m] + pard[m * 11 + 8] - 
		    pard[m * 11 + 7]) * (ths[m] - thr[m]) * (pard[m * 11 + 8] 
		    - pard[m * 11 + 7]);
	} else if (*idualpor == 2) {
/* Pressure head driven */
	    h__ = fh_(imodel, &seim, &pard[m * 11 + 7]);
	    par[0] = pard[m * 11 + 7];
/* Variable coeff. (functi */
	    par[1] = pard[m * 11 + 8];
	    par[2] = pard[m * 11 + 9];
	    par[3] = pard[m * 11 + 10];
	    par[4] = pard[m * 11 + 11];
	    par[5] = pard[m * 11 + 6];
	    condm = fk_(imodel, &h__, par);
	    condf = fk_(imodel, &hnew[i__], par);
	    sinkim[i__] = (condm + condf) * .5f * (hnew[i__] - h__);
	    if (i__ == *n && (r__1 = *hcrita - hnew[*n], dabs(r__1)) < *
		    hcrita * -.001f && (r__2 = *hcrita - h__, dabs(r__2)) < *
		    hcrita * -.01f) {
		sinkim[i__] = 0.f;
	    }
	}
	if (sinkim[i__] > 0.f) {
	    trmaxim = (pard[m * 11 + 8] - tholdim[i__]) / *dt;
	    if (*idualpor == 1) {
/* Computing MIN */
		r__1 = deltath / *dt;
		trmaxim = dmin(r__1,trmaxim);
	    }
	    if (sinkim[i__] > trmaxim) {
		sinkim[i__] = trmaxim;
	    }
	}
	if (sinkim[i__] < 0.f) {
	    trmaxim = -(ths[m] - thold[i__]) / *dt;
	    if (*idualpor == 1) {
/* Computing MAX */
		r__1 = deltath / *dt;
		trmaxim = dmax(r__1,trmaxim);
	    }
	    if (sinkim[i__] < trmaxim) {
		sinkim[i__] = trmaxim;
	    }
	}
/* Computing MAX */
/* Computing MIN */
	r__3 = tholdim[i__] + sinkim[i__] * *dt, r__4 = pard[m * 11 + 8];
	r__1 = dmin(r__3,r__4), r__2 = pard[m * 11 + 7];
	thnewim[i__] = dmax(r__1,r__2);
	if (i__ >= 2) {
	    *wtransf += (sinkim[i__ - 1] + sinkim[i__]) / 2.f * (x[i__] - x[
		    i__ - 1]);
	}
/* L11: */
    }
    return 0;
} /* dualpor_ */

/* *********************************************************************** */
/* Subroutine */ int update_(integer *numnp, logical *lwat, logical *lchem, 
	logical *ltemp, logical *lvapor, integer *idualpor, logical *lextrap, 
	real *dt, real *dtold, real *htemp, real *hnew, real *hold, real *
	thold, real *thnew, real *vold, real *vnew, real *thvold, real *
	thvnew, real *vvold, real *vvnew, real *tholdim, real *thnewim, real *
	tempo, real *tempn, real *rtop, real *xconv, real *consmax, integer *
	kodtop, integer *kodbot)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ibot;
    static logical lsat;
    static integer itop;

    /* Parameter adjustments */
    --tempn;
    --tempo;
    --thnewim;
    --tholdim;
    --vvnew;
    --vvold;
    --thvnew;
    --thvold;
    --vnew;
    --vold;
    --thnew;
    --thold;
    --hold;
    --hnew;
    --htemp;

    /* Function Body */
    lsat = TRUE_;
    ibot = 1;
    if (*kodbot > 0) {
	ibot = 2;
    }
    itop = *numnp;
    if (*kodtop > 0) {
	itop = *numnp - 1;
    }
    i__1 = itop;
    for (i__ = ibot; i__ <= i__1; ++i__) {
	if (*lwat) {
	    if (*lextrap && hnew[i__] < 0.f && hold[i__] < 0.f) {
		htemp[i__] = hnew[i__] + (hnew[i__] - hold[i__]) * *dt / *
			dtold;
	    } else {
		htemp[i__] = hnew[i__];
	    }
	    hold[i__] = hnew[i__];
	    hnew[i__] = htemp[i__];
	}
/* L11: */
    }
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*lwat) {
	    thold[i__] = thnew[i__];
	    if (*ltemp || *lchem) {
		vold[i__] = vnew[i__];
	    }
	    if (*lvapor) {
		thvold[i__] = thvnew[i__];
	    }
	    if (*lvapor) {
		vvold[i__] = vvnew[i__];
	    }
	    if (*idualpor > 0) {
		tholdim[i__] = thnewim[i__];
	    }
	    if (hnew[i__] < 0.f) {
		lsat = FALSE_;
	    }
	}
	if (*ltemp) {
	    tempo[i__] = tempn[i__];
	}
/* L12: */
    }
    if (*lwat && lsat && (*kodtop == -1 || *kodtop == 4) && *rtop >= -(*
	    consmax)) {
	if (*kodtop == 4) {
	    *kodtop = -4;
	}
	hnew[*numnp] = *xconv * -.005f;
/* 0.5 cm */
	htemp[*numnp] = *xconv * -.005f;
	hold[*numnp] = *xconv * -.005f;
    }
    return 0;
} /* update_ */

