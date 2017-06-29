/* SOLUTE.f -- translated by f2c (version 12.02.01).
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
static integer c__10 = 10;
static integer c__2 = 2;
static integer c__11 = 11;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__12 = 12;
static integer c__5 = 5;
static integer c__13 = 13;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__9 = 9;
static doublereal c_b89 = 2.3333333333333335;
static doublereal c_b91 = 2.6666666666666665;
static doublereal c_b92 = 1.5;
static doublereal c_b116 = 2.;
static doublereal c_b117 = 3.;
static doublereal c_b129 = .33333333333333331;
static doublereal c_b131 = -.66666666666666663;
static doublereal c_b132 = .125;
static doublereal c_b133 = 1.875;
static doublereal c_b134 = 1.2;
static doublereal c_b135 = -.4;
static doublereal c_b173 = 10.;

/* Source file SOLUTE.FOR ||||||||||||||||||||||||||||||||||||||||||||||| */
/*     To assemble and solve the solute transport equation */
/*     Mass-lumping finite elements */
/* Subroutine */ int solute_(integer *n, integer *nmat, integer *ns, integer *
	nsd, real *x, real *dt, doublereal *t, real *tpulse, real *chpar, 
	integer *matnum, real *tho, real *thn, real *vo, real *vn, real *disp,
	 real *epsi, integer *ktopch, real *ctop, integer *kbotch, real *cbot,
	 real *conc, doublereal *b, doublereal *d__, doublereal *e, 
	doublereal *f, real *g0, real *g1, real *retard, real *cvtop, real *
	cvbot, real *cvch0, real *cvch1, logical *lupw, real *wc, real *
	peclet, real *courant, real *dtmaxc, real *tempo, real *tempn, real *
	cnew, real *cprevo, real *ctemp, real *tdep, real *thsat, real *ctola,
	 real *ctolr, integer *iterc, integer *maxitc, real *vcorr, real *
	sorb, real *sorbn, logical *llinear, logical *lequil, logical *lartd, 
	real *pecr, real *q0, real *q1, real *dsurf, real *catm, logical *
	ltort, real *sink, real *crootmax, real *ssink, real *cvchr, logical *
	lmobim, real *cvchim, integer *tlevel, logical *lbact, real *sorb2, 
	real *sorbn2, real *dtmin, real *dtopt, logical *lwat, logical *
	lfiltr, integer *idualpor, real *thoim, real *thnim, real *sinkim, 
	real *strans, integer *itort, real *xconv, real *tconv, logical *
	lvapor, real *rbot, integer *ierr, integer *imoistdep, integer *nmatd,
	 real *dmoist, real *wdep, integer *iconctype, real *beta, logical *
	ldualneq, logical *atmbc, logical *sinkf, logical *lactrsu, real *
	omegas, real *omegaw, real *spot, real *rkm, real *cmin, logical *
	ldensity)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset, conc_dim1, conc_offset, sorb_dim1, 
	    sorb_offset, sorb2_dim1, sorb2_offset, dmoist_dim1, dmoist_dim2, 
	    dmoist_offset, wdep_dim1, wdep_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    extern /* Subroutine */ int sorbconc_(integer *, integer *, integer *, 
	    integer *, real *, logical *, real *, real *, real *, real *, 
	    real *, integer *, logical *, real *, real *, logical *, real *, 
	    integer *, real *, real *, real *, real *, integer *, integer *, 
	    real *, real *, real *, real *, logical *), fluxconc_(integer *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    integer *, real *, real *, real *, real *, logical *, logical *, 
	    integer *, real *, integer *), masstran_(integer *, integer *, 
	    integer *, integer *, integer *, real *, logical *, logical *, 
	    real *, real *, real *, real *, integer *, real *, real *, real *,
	     real *, real *, real *, real *, real *, real *, logical *, real *
	    , real *, logical *, real *, integer *, real *, real *, real *, 
	    logical *, real *, logical *);
    static integer i__, m;
    static real r__, d1, e1, f1, dg, bn, dn, fn;
    static integer js;
    static real tr, tt, alf;
    static integer jjj;
    static real pecl, rmin;
    static integer iter;
    static real cour;
    extern /* Subroutine */ int coeff_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, integer *, real *, real *, real *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, logical *, logical *, logical *, logical *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *, logical *, real *, logical *, logical *, real *, real *, 
	    logical *, integer *, real *, real *, real *, integer *, real *, 
	    real *, integer *, integer *, real *, real *, logical *, logical *
	    );
    static integer level;
    static real dtold;
    static logical lconv;
    static real dtmxc, henry;
    extern /* Subroutine */ int bansol_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static integer nlevel;
    extern /* Subroutine */ int matset_(integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, integer *, integer *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, doublereal *, doublereal *
	    , doublereal *, doublereal *, real *, real *, real *, real *, 
	    real *, real *, integer *, real *, real *, real *, real *, real *,
	     real *, integer *, logical *, integer *, logical *, real *, 
	    logical *);
    static real dsurft;
    static logical lnequil;
    extern /* Subroutine */ int setssnk_(integer *, integer *, integer *, 
	    doublereal *, real *, real *, real *, real *, integer *, real *, 
	    real *, real *, logical *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --beta;
    --strans;
    --sinkim;
    --thnim;
    --thoim;
    --sorbn2;
    --ssink;
    --sink;
    --q1;
    --q0;
    --sorbn;
    --vcorr;
    --ctemp;
    --cprevo;
    --cnew;
    --tempn;
    --tempo;
    --wc;
    --retard;
    --g1;
    --g0;
    --f;
    --e;
    --d__;
    --b;
    --disp;
    --vn;
    --vo;
    --thn;
    --tho;
    --matnum;
    --x;
    --lmobim;
    --thsat;
    --cvchim;
    --cvchr;
    --crootmax;
    --cvch1;
    --cvch0;
    --cvbot;
    --cvtop;
    --cbot;
    --ctop;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    --llinear;
    sorb_dim1 = *nsd;
    sorb_offset = 1 + sorb_dim1;
    sorb -= sorb_offset;
    --tdep;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;
    wdep_dim1 = 2 + *nmatd;
    wdep_offset = 1 + wdep_dim1;
    wdep -= wdep_offset;
    dmoist_dim1 = *nmatd;
    dmoist_dim2 = *nsd;
    dmoist_offset = 1 + dmoist_dim1 * (1 + dmoist_dim2 * 14);
    dmoist -= dmoist_offset;

    /* Function Body */
    alf = 1.f - *epsi;
    *iterc = 1.f;
    nlevel = 2;
    *peclet = 0.f;
    *courant = 0.f;
    *dtmaxc = 1e30f;
    rmin = 1e-30f;
    tr = 293.15f;
    r__ = 8.314f;
/*     Sequential first order decay goes into equilibrium phase (lNEquil=.false.) or nonequilbrium phase (lNEq
uil=.true.) */
    lnequil = FALSE_;
L10:
/*     Loop on species in the chain */
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
	iter = 0;
	jjj = js - 1 << 4;
	cvtop[js] = 0.f;
	cvbot[js] = 0.f;
	cvch0[js] = 0.f;
	cvch1[js] = 0.f;
	cvchr[js] = 0.f;
	cvchim[js] = 0.f;
	if (*t - *tpulse > *dtmin && ! (*atmbc)) {
	    ctop[js] = 0.f;
	    cbot[js] = 0.f;
	}
	if (*kbotch < 0) {
	    if (vo[1] >= 0.f) {
		cvbot[js] = alf * cbot[js] * vo[1];
	    }
	    if (vo[1] < 0.f) {
		cvbot[js] = alf * conc[js + conc_dim1] * vo[1];
	    }
	    if (*lvapor && *rbot == 0.f) {
		cvbot[js] = 0.f;
	    }
	} else if (*kbotch == 0) {
	    cvbot[js] = alf * conc[js + conc_dim1] * vo[1];
	}
	if ((real) (*ktopch) < 0.f && *tlevel != 1) {
	    if (vo[*n] < 0.f) {
		cvtop[js] = alf * ctop[js] * vo[*n];
	    }
	}
	if (*ktopch == -2) {
	    m = matnum[*n];
	    tr = 293.15f;
	    r__ = 8.314f;
	    tt = (tempo[*n] + 273.15f - tr) / r__ / (tempo[*n] + 273.15f) / 
		    tr;
	    dg = chpar[jjj + 6 + m * chpar_dim1] * exp(tdep[jjj + 6] * tt);
	    henry = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 10] * 
		    tt);
	    dsurft = *dsurf * exp(tdep[jjj + 9] * tt);
	    cvtop[js] = cvtop[js] + alf * dg / dsurft * henry * conc[js + *n *
		     conc_dim1] - dg / dsurft * *catm;
	}
	if (! llinear[js]) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		cnew[i__] = conc[js + i__ * conc_dim1];
		if (! (*lequil)) {
		    sorbn[i__] = sorb[js + i__ * sorb_dim1];
		}
		if (*lbact || *ldualneq) {
		    sorbn2[i__] = sorb2[js + i__ * sorb2_dim1];
		}
/* L11: */
	    }
	}
/*       Root Solute Uptake */
	if (*sinkf) {
	    setssnk_(&js, ns, n, t, &x[1], &beta[1], &sink[1], &ssink[1], nsd,
		     &conc[conc_offset], omegaw, &crootmax[js], lactrsu, 
		    omegas, spot, rkm, cmin);
	}
/*       Iterative loop for a nonlinear adsorption isotherm */
L12:
	++iter;
	if (! llinear[js]) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ctemp[i__] = cnew[i__];
/* L13: */
	    }
	}
/*       To construct the matrix equation */
	i__2 = nlevel;
	for (level = 1; level <= i__2; ++level) {
/*         Calculate the dispersion coefficients, retardation factors, source/ */
/*         decay coefficients, Peclet and Courant numbers, upstream weighting */
/*         factors */
	    coeff_(&js, &level, &nlevel, n, nmat, nsd, &x[1], &disp[1], &vo[1]
		    , &vn[1], &tho[1], &thn[1], &thsat[1], &chpar[
		    chpar_offset], &matnum[1], &tempn[1], &tempo[1], &tdep[1],
		     &g0[1], &g1[1], &retard[1], &conc[conc_offset], &cnew[1],
		     &cprevo[1], dt, &pecl, &cour, &dtmxc, &llinear[1], 
		    lequil, lupw, lartd, &iter, &wc[1], &vcorr[1], &sorb[
		    sorb_offset], &sorbn[1], epsi, pecr, &q0[1], &q1[1], 
		    ltort, &ssink[1], &lmobim[1], lbact, &sorb2[sorb2_offset],
		     &sorbn2[1], lfiltr, idualpor, &thoim[1], &thnim[1], &
		    sinkim[1], itort, xconv, tconv, imoistdep, nmatd, &dmoist[
		    dmoist_offset], &wdep[wdep_offset], &lnequil, ldualneq);
	    *peclet = dmax(*peclet,pecl);
	    *courant = dmax(*courant,cour);
	    *dtmaxc = dmin(*dtmaxc,dtmxc);
/*         Set up the matrix equation */
	    matset_(&js, n, ns, nsd, &level, epsi, &alf, dt, kbotch, ktopch, &
		    cbot[1], &ctop[1], &x[1], &tho[1], &thn[1], &vo[1], &vn[1]
		    , &conc[conc_offset], &disp[1], &retard[1], &wc[1], &g0[1]
		    , &g1[1], &b[1], &d__[1], &e[1], &f[1], &e1, &d1, &f1, &
		    bn, &dn, &fn, nmat, &chpar[chpar_offset], &tempo[1], &
		    tempn[1], &tdep[1], &dsurft, catm, &matnum[1], &lmobim[1],
		     idualpor, lvapor, rbot, lbact);
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		if (level == 1) {
		    vo[i__] += vcorr[i__];
		}
		if (level == 2) {
		    vn[i__] += vcorr[i__];
		}
/* L14: */
	    }
/*         Calculate mass-transfer fluxes at the beginning of the time interval */
	    if (level == 1 && iter == 1) {
		masstran_(&js, ns, nsd, n, &matnum[1], &tempo[1], &lmobim[1], 
			lequil, &chpar[chpar_offset], &tdep[1], &sorb[
			sorb_offset], &conc[conc_offset], nmat, &x[1], &cvch0[
			1], &cvch1[1], &cvchr[1], &cvchim[1], &alf, &q0[1], &
			q1[1], &ssink[1], lbact, &tho[1], &sorb2[sorb2_offset]
			, lfiltr, &vo[1], idualpor, &sinkim[1], xconv, tconv, 
			ldualneq, &strans[1], &llinear[1]);
	    }
/* L15: */
	}
/*       Solve matrix equation */
	bansol_(n, &b[1], &d__[1], &e[1], &f[1]);
/*       Test for convergence for nonlinear problem */
	lconv = TRUE_;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (*ns > 1 && iter == 1 || *ldensity) {
		cprevo[i__] = conc[js + i__ * conc_dim1];
	    }
	    if (llinear[js]) {
/* Computing MAX */
		r__1 = (real) f[i__];
		conc[js + i__ * conc_dim1] = dmax(r__1,0.f);
		if (conc[js + i__ * conc_dim1] < 1e-30f && conc[js + i__ * 
			conc_dim1] > 0.f) {
		    conc[js + i__ * conc_dim1] = 0.f;
		}
	    } else {
		cnew[i__] = (real) f[i__];
		if (cnew[i__] < 1e-30f) {
		    cnew[i__] = 0.f;
		}
		if ((r__1 = cnew[i__] - ctemp[i__], dabs(r__1)) > *ctola + *
			ctolr * conc[js + i__ * conc_dim1]) {
		    lconv = FALSE_;
		}
	    }
/* L16: */
	}
	if (! llinear[js]) {
	    if (! lconv) {
		if (iter < *maxitc) {
		    goto L12;
		} else if (*dt > *dtmin && ! (*lwat)) {
/*              ierr=1 */
		    dtold = *dt;
/* Computing MAX */
		    r__1 = *dt / 3.f;
		    *dt = dmax(r__1,*dtmin);
		    *dtopt = *dt;
		    *t = *t - dtold + *dt;
		    goto L10;
		} else {
		    *ierr = 1;
		}
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		conc[js + i__ * conc_dim1] = cnew[i__];
		if (! (*lequil)) {
		    sorb[js + i__ * sorb_dim1] = sorbn[i__];
		}
		if (*lbact || *ldualneq) {
		    sorb2[js + i__ * sorb2_dim1] = sorbn2[i__];
		}
/* L17: */
	    }
	}
/*       Calculate sorbed concentration for linear noneq. adsorption or */
/*       concentration in the imobile water. */
	if (! (*lequil) && llinear[js]) {
	    sorbconc_(&js, nsd, n, &matnum[1], &tempn[1], &lmobim[1], &chpar[
		    chpar_offset], &tdep[1], &sorb[sorb_offset], &conc[
		    conc_offset], dt, nmat, lbact, &thn[1], &sorb2[
		    sorb2_offset], lfiltr, &vn[1], idualpor, &thnim[1], &
		    thoim[1], &sinkim[1], &strans[1], imoistdep, nmatd, &
		    dmoist[dmoist_offset], &wdep[wdep_offset], xconv, tconv, 
		    ldualneq);
	}
/*       Calculate mass-transfer fluxes at the end of the time interval */
	masstran_(&js, ns, nsd, n, &matnum[1], &tempn[1], &lmobim[1], lequil, 
		&chpar[chpar_offset], &tdep[1], &sorb[sorb_offset], &conc[
		conc_offset], nmat, &x[1], &cvch0[1], &cvch1[1], &cvchr[1], &
		cvchim[1], epsi, &q0[1], &q1[1], &ssink[1], lbact, &thn[1], &
		sorb2[sorb2_offset], lfiltr, &vn[1], idualpor, &sinkim[1], 
		xconv, tconv, ldualneq, &strans[1], &llinear[1]);
/*       Set up mass fluxes */
	if (*ktopch < 0) {
	    if (*tlevel != 1) {
		if (vn[*n] < 0.f) {
		    cvtop[js] += *epsi * vn[*n] * ctop[js];
		}
	    } else {
		if (vn[*n] < 0.f) {
		    cvtop[js] += vn[*n] * ctop[js];
		}
	    }
	} else {
	    cvtop[js] = fn - bn * conc[js + (*n - 1) * conc_dim1] - dn * conc[
		    js + *n * conc_dim1];
	}
	if (*ktopch == -2) {
	    m = matnum[*n];
	    tt = (tempn[*n] + 273.15f - tr) / r__ / (tempn[*n] + 273.15f) / 
		    tr;
	    dg = chpar[jjj + 6 + m * chpar_dim1] * exp(tdep[jjj + 6] * tt);
	    henry = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 10] * 
		    tt);
	    cvtop[js] = cvtop[js] + *epsi * dg / dsurft * henry * conc[js + *
		    n * conc_dim1] - dg / dsurft * *catm;
	}
	if (*kbotch < 0) {
	    if (vn[1] >= 0.f) {
		cvbot[js] += *epsi * cbot[js] * vn[1];
	    }
	    if (vn[1] < 0.f) {
		cvbot[js] += *epsi * conc[js + conc_dim1] * vn[1];
	    }
	    if (*lvapor && *rbot == 0.f) {
		cvbot[js] = 0.f;
	    }
	} else if (*kbotch == 0) {
	    cvbot[js] += *epsi * conc[js + conc_dim1] * vn[1];
	} else {
	    cvbot[js] = d1 * conc[js + conc_dim1] + e1 * conc[js + (conc_dim1 
		    << 1)] - f1;
	}
	*iterc = max(*iterc,iter);
	if ((r__1 = cvtop[js], dabs(r__1)) < rmin) {
	    cvtop[js] = 0.f;
	}
	if ((r__1 = cvbot[js], dabs(r__1)) < rmin) {
	    cvbot[js] = 0.f;
	}
/* L18: */
    }
/*     Calculate flux concentrations */
    if (*iconctype == 2) {
	fluxconc_(n, nmat, nsd, &x[1], &vn[1], &thn[1], &thsat[1], &chpar[
		chpar_offset], &matnum[1], &tempn[1], &tdep[1], &conc[
		conc_offset], &cnew[1], ltort, &lmobim[1], idualpor, &thnim[1]
		, &c__1);
    }
    return 0;
} /* solute_ */

/* *********************************************************************** */
/*     Calculate the dispersion coefficients, retardation factors, source/ */
/*     decay coefficients, Peclet and Courant numbers, upstream weighting */
/*     factors */
/* Subroutine */ int coeff_(integer *js, integer *level, integer *nlevel, 
	integer *numnp, integer *nmat, integer *nsd, real *x, real *disp, 
	real *vo, real *vn, real *tho, real *thn, real *thsat, real *chpar, 
	integer *matnum, real *tempn, real *tempo, real *tdep, real *g0, real 
	*g1, real *retard, real *conc, real *cnew, real *cprevo, real *dt, 
	real *peclet, real *courant, real *dtmaxc, logical *llinear, logical *
	lequil, logical *lupw, logical *lartd, integer *iter, real *wc, real *
	vcorr, real *sorb, real *sorbn, real *epsi, real *pecr, real *q0, 
	real *q1, logical *ltort, real *ssink, logical *lmobim, logical *
	lbact, real *sorb2, real *sorbn2, logical *lfiltr, integer *idualpor, 
	real *thoim, real *thnim, real *sinkim, integer *itort, real *xconv, 
	real *tconv, integer *imoistdep, integer *nmatd, real *dmoist, real *
	wdep, logical *lnequil, logical *ldualneq)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset, conc_dim1, conc_offset, sorb_dim1, 
	    sorb_offset, sorb2_dim1, sorb2_offset, dmoist_dim1, dmoist_dim2, 
	    dmoist_offset, wdep_dim1, wdep_offset, i__1;
    real r__1, r__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12;

    /* Local variables */
    extern /* Subroutine */ int blocking_(integer *, real *, real *, real *, 
	    real *, real *, real *);
    static real dretards;
    static integer i__, j, k, m;
    static real r__, v, f1, aa, cc, dc, cg, dg, dp, dw, dx, vj, ro, tr, ss, 
	    tt, cg1;
    static integer jj1;
    static real ss1, ss2;
    static integer jjj;
    static real dks, thg;
    extern doublereal rmd_(integer *, integer *, integer *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *);
    static real thj, dnu, tti, ttj, thw, xks, tto, ttn, xnu, rka1, rka2, rkd1,
	     rkd2, psi1, psi2, f_em__, frac, gamg, cmid, gaml, derk, gams, 
	    smid, fexp, taug, xmug, thwo, xksn, alfa1, xmul, alfa2, xksp, 
	    xkso, xnup, xnuo, gamg1, xmus, xnun, gaml1, gams1, rka1o, rka2o, 
	    rkd1o, rkd2o;
    static integer ipsi1, ipsi2;
    static real smax1, smax2, psi1o, psi2o, gamgi, dconc, omega, gamli, dmobi,
	     gamlo, gamsi, ddexp, sconc, gamso, cprev, fexpp, henry, fexpo, 
	    fexpn, xmugo, ssorb, gamg1i, xmulo, gaml1i, xmuso, gamg1p, gaml1o,
	     gams1i, gaml1p, gams1o, gams1p, smax1o, smax2o, ssorb2, omegao, 
	    gamloi, dsconc, dconcs, omegas, thimob, gamsoi, sconco, sconcp;
    extern /* Subroutine */ int disper_(integer *, integer *, integer *, 
	    integer *, integer *, real *, logical *, logical *, logical *, 
	    logical *, integer *, integer *, integer *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, integer *, logical *);
    static real sconcs, dhenry;
    extern /* Subroutine */ int nequil_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, logical *, logical *, logical *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, integer *,
	     real *, real *, logical *, real *, real *, real *, real *, 
	    logical *, real *, real *);
    static real henryi, henryj, henryn, henryo, henryp;
    extern /* Subroutine */ int pecour_(integer *, integer *, integer *, 
	    integer *, integer *, logical *, logical *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, integer *, real *);
    static real gaml1oi, gaml1pi, gams1oi, gams1pi, flmacro, dretard, dsconcs,
	     thimobo;
    extern /* Subroutine */ int deposit_(real *, real *, real *, real *, real 
	    *, real *, real *, real *, real *, real *, real *);
    static real sconcos, sconcps, courmax;

    /* Fortran I/O blocks */
    static cilist io___161 = { 0, 6, 0, 0, 0 };
    static cilist io___162 = { 0, 6, 0, 0, 0 };
    static cilist io___163 = { 0, 5, 0, 0, 0 };


/*     Inicialization */
    /* Parameter adjustments */
    --sinkim;
    --thnim;
    --thoim;
    --sorbn2;
    --ssink;
    --q1;
    --q0;
    --sorbn;
    --vcorr;
    --wc;
    --cprevo;
    --cnew;
    --retard;
    --g1;
    --g0;
    --tempo;
    --tempn;
    --matnum;
    --thn;
    --tho;
    --vn;
    --vo;
    --disp;
    --x;
    --lmobim;
    --thsat;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    sorb_dim1 = *nsd;
    sorb_offset = 1 + sorb_dim1;
    sorb -= sorb_offset;
    --llinear;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    --tdep;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;
    wdep_dim1 = 2 + *nmatd;
    wdep_offset = 1 + wdep_dim1;
    wdep -= wdep_offset;
    dmoist_dim1 = *nmatd;
    dmoist_dim2 = *nsd;
    dmoist_offset = 1 + dmoist_dim1 * (1 + dmoist_dim2 * 14);
    dmoist -= dmoist_offset;

    /* Function Body */
    jjj = *js - 1 << 4;
    if (*js > 1) {
	jj1 = jjj - 16;
    }
    *peclet = 0.f;
    *courant = 0.f;
    courmax = 1.f;
    *dtmaxc = 1e30f;
    tr = 293.15f;
    r__ = 8.314f;
    for (i__ = *numnp; i__ >= 1; --i__) {
	j = i__ + 1;
	k = i__ - 1;
	m = matnum[i__];
	if (*level == *nlevel) {
	    thw = thn[i__];
	    thwo = tho[i__];
/* Computing MAX */
	    r__1 = 0.f, r__2 = thsat[m] - thw;
	    thg = dmax(r__1,r__2);
	    if (lmobim[m] && *idualpor == 0 || *lbact) {
		thimob = chpar[m * chpar_dim1 + 4];
		thimobo = thimob;
/* Computing MAX */
		r__1 = thw - thimob;
		thw = dmax(r__1,.001f);
	    }
	    if (*idualpor > 0) {
		thimob = thnim[i__];
		thimobo = thoim[i__];
	    }
	    v = vn[i__];
	    if (i__ != *numnp) {
		vj = vn[j];
		thj = thn[j];
		if (lmobim[m] && *idualpor == 0 || *lbact) {
/* Computing MAX */
		    r__1 = thj - thimob;
		    thj = dmax(r__1,.001f);
		}
	    }
	    tt = (tempn[i__] + 273.15f - tr) / r__ / (tempn[i__] + 273.15f) / 
		    tr;
	    if (*js > 1) {
		cprev = conc[*js - 1 + i__ * conc_dim1];
	    }
	} else {
	    thw = tho[i__];
/* Computing MAX */
	    r__1 = 0.f, r__2 = thsat[m] - thw;
	    thg = dmax(r__1,r__2);
	    if (lmobim[m] && *idualpor == 0 || *lbact) {
		thimob = chpar[m * chpar_dim1 + 4];
/* Computing MAX */
		r__1 = thw - thimob;
		thw = dmax(r__1,.001f);
	    }
	    if (*idualpor > 0) {
		thimob = thoim[i__];
	    }
	    v = vo[i__];
	    if (i__ != *numnp) {
		vj = vo[j];
		thj = tho[j];
		if (lmobim[m] && *idualpor == 0 || *lbact) {
/* Computing MAX */
		    r__1 = thj - thimob;
		    thj = dmax(r__1,.001f);
		}
	    }
	    tt = (tempo[i__] + 273.15f - tr) / r__ / (tempo[i__] + 273.15f) / 
		    tr;
	    if (*js > 1) {
		cprev = cprevo[i__];
	    }
	}
/*       Temperature dependence */
	f1 = 1.f;
	ro = chpar[m * chpar_dim1 + 1] * exp(tdep[1] * tt);
	frac = chpar[m * chpar_dim1 + 3] * exp(tdep[3] * tt);
	dw = chpar[jjj + 5 + m * chpar_dim1] * exp(tdep[jjj + 5] * tt);
	dg = chpar[jjj + 6 + m * chpar_dim1] * exp(tdep[jjj + 6] * tt);
	xks = chpar[jjj + 7 + m * chpar_dim1] * exp(tdep[jjj + 7] * tt);
	xnu = chpar[jjj + 8 + m * chpar_dim1] * exp(tdep[jjj + 8] * tt);
	fexp = chpar[jjj + 9 + m * chpar_dim1];
	henry = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 10] * tt);
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__1, &dmoist[dmoist_offset], &
		    c__1, &wdep[wdep_offset], &thw, imoistdep);
	}
	gaml = chpar[jjj + 11 + m * chpar_dim1] * exp(tdep[jjj + 11] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__10, &dmoist[dmoist_offset], &
		    c__1, &wdep[wdep_offset], &thimob, imoistdep);
	}
	gamli = chpar[jjj + 11 + m * chpar_dim1] * exp(tdep[jjj + 11] * tt) * 
		f1;
/* reaction in */
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__2, &dmoist[dmoist_offset], &
		    c__2, &wdep[wdep_offset], &thw, imoistdep);
	}
	gams = chpar[jjj + 12 + m * chpar_dim1] * exp(tdep[jjj + 12] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__11, &dmoist[dmoist_offset], &
		    c__2, &wdep[wdep_offset], &thimob, imoistdep);
	}
	gamsi = chpar[jjj + 12 + m * chpar_dim1] * exp(tdep[jjj + 12] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__3, &dmoist[dmoist_offset], &
		    c__3, &wdep[wdep_offset], &thw, imoistdep);
	}
	gamg = chpar[jjj + 13 + m * chpar_dim1] * exp(tdep[jjj + 13] * tt) * 
		f1;
	gamgi = gamg;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__4, &dmoist[dmoist_offset], &
		    c__4, &wdep[wdep_offset], &thw, imoistdep);
	}
	gaml1 = chpar[jjj + 14 + m * chpar_dim1] * exp(tdep[jjj + 14] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__12, &dmoist[dmoist_offset], &
		    c__4, &wdep[wdep_offset], &thimob, imoistdep);
	}
	gaml1i = chpar[jjj + 14 + m * chpar_dim1] * exp(tdep[jjj + 14] * tt) *
		 f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__5, &dmoist[dmoist_offset], &
		    c__5, &wdep[wdep_offset], &thw, imoistdep);
	}
	gams1 = chpar[jjj + 15 + m * chpar_dim1] * exp(tdep[jjj + 15] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__13, &dmoist[dmoist_offset], &
		    c__5, &wdep[wdep_offset], &thimob, imoistdep);
	}
	gams1i = chpar[jjj + 15 + m * chpar_dim1] * exp(tdep[jjj + 15] * tt) *
		 f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__6, &dmoist[dmoist_offset], &
		    c__6, &wdep[wdep_offset], &thw, imoistdep);
	}
	gamg1 = chpar[jjj + 16 + m * chpar_dim1] * exp(tdep[jjj + 16] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__7, &dmoist[dmoist_offset], &
		    c__7, &wdep[wdep_offset], &thw, imoistdep);
	}
	xmul = chpar[jjj + 17 + m * chpar_dim1] * exp(tdep[jjj + 17] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__8, &dmoist[dmoist_offset], &
		    c__8, &wdep[wdep_offset], &thw, imoistdep);
	}
	xmus = chpar[jjj + 18 + m * chpar_dim1] * exp(tdep[jjj + 18] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__9, &dmoist[dmoist_offset], &
		    c__9, &wdep[wdep_offset], &thw, imoistdep);
	}
	xmug = chpar[jjj + 19 + m * chpar_dim1] * exp(tdep[jjj + 19] * tt) * 
		f1;
	omega = chpar[jjj + 20 + m * chpar_dim1] * exp(tdep[jjj + 20] * tt);
	f_em__ = 1.f;
	if (*ldualneq) {
	    f_em__ = chpar[jjj + 13 + m * chpar_dim1] * exp(tdep[jjj + 13] * 
		    tt);
	}
	if (*ldualneq) {
	    omegas = chpar[jjj + 16 + m * chpar_dim1] * exp(tdep[jjj + 16] * 
		    tt);
	}
	if (*lbact) {
	    dg = 0.f;
	    gamg = 0.f;
	    gamgi = 0.f;
	    gaml1 = 0.f;
	    gaml1i = 0.f;
	    gams1 = 0.f;
	    gams1i = 0.f;
	    gamg1 = 0.f;
	    gamg1i = 0.f;
	    xmul = 0.f;
	    xmus = 0.f;
	    xmug = 0.f;
	    omega = 0.f;
	    smax2 = chpar[jjj + 15 + m * chpar_dim1] * exp(tdep[jjj + 15] * 
		    tt);
	    rka2 = chpar[jjj + 16 + m * chpar_dim1] * exp(tdep[jjj + 16] * tt)
		    ;
	    rkd2 = chpar[jjj + 17 + m * chpar_dim1] * exp(tdep[jjj + 17] * tt)
		    ;
	    smax1 = chpar[jjj + 18 + m * chpar_dim1] * exp(tdep[jjj + 18] * 
		    tt);
	    rka1 = chpar[jjj + 19 + m * chpar_dim1] * exp(tdep[jjj + 19] * tt)
		    ;
	    rkd1 = chpar[jjj + 20 + m * chpar_dim1] * exp(tdep[jjj + 20] * tt)
		    ;
	    ipsi1 = 0;
	    ipsi2 = 0;
	    if (! (*lfiltr)) {
		ipsi2 = (integer) chpar[jjj + 13 + m * chpar_dim1];
	    }
	    if (! (*lfiltr)) {
		ipsi1 = (integer) chpar[jjj + 14 + m * chpar_dim1];
	    }
	    if (ipsi1 == 0 && smax1 > 0.f) {
		ipsi1 = 1;
	    }
	    if (ipsi2 == 0 && smax2 > 0.f) {
		ipsi2 = 1;
	    }
	    if (ipsi1 >= 3 || ipsi2 >= 3) {
		dc = chpar[jjj + 6 + m * chpar_dim1] * exp(tdep[jjj + 6] * tt)
			;
	    }
	    if (ipsi1 == 5 || ipsi2 == 5) {
		aa = chpar[jjj + 15 + m * chpar_dim1];
	    }
	    if (*level == *nlevel) {
		ss1 = sorbn[i__];
		ss2 = sorbn2[i__];
	    } else {
		ss1 = sorb[*js + i__ * sorb_dim1];
		ss2 = sorb2[*js + i__ * sorb2_dim1];
	    }
	    psi1 = 1.f;
	    psi2 = 1.f;
	    if (ipsi1 > 0) {
		blocking_(&ipsi1, &smax1, &psi1, &x[i__], &ss1, &dc, &aa);
	    }
	    if (ipsi2 > 0) {
		blocking_(&ipsi2, &smax2, &psi2, &x[i__], &ss2, &dc, &aa);
	    }
/*         recalculate ka1 and ka2 based on filtration theory */
	    if (*lfiltr) {
		gamg = 0.f;
		gaml1 = 0.f;
		dc = chpar[jjj + 13 + m * chpar_dim1] * exp(tdep[jjj + 13] * 
			tt);
		dp = chpar[jjj + 14 + m * chpar_dim1] * exp(tdep[jjj + 14] * 
			tt);
		alfa1 = rka1;
		alfa2 = rka2;
		deposit_(&rka1, &rka2, &dc, &dp, &alfa1, &alfa2, &thw, &v, &
			tempn[i__], xconv, tconv);
	    }
	}
	if (*js > 1) {
	    xksp = chpar[jj1 + 7 + m * chpar_dim1] * exp(tdep[jj1 + 7] * tt);
	    xnup = chpar[jj1 + 8 + m * chpar_dim1] * exp(tdep[jj1 + 8] * tt);
	    fexpp = chpar[jj1 + 9 + m * chpar_dim1];
/* *exp(TDep(jj1+ 9)*TT) */
	    henryp = chpar[jj1 + 10 + m * chpar_dim1] * exp(tdep[jj1 + 10] * 
		    tt);
	    f1 = 1.f;
	    if (*imoistdep > 0) {
		i__1 = *js - 1;
		f1 = rmd_(nmatd, nsd, &m, &i__1, &c__4, &dmoist[dmoist_offset]
			, &c__4, &wdep[wdep_offset], &thw, imoistdep);
	    }
	    gaml1p = chpar[jj1 + 14 + m * chpar_dim1] * exp(tdep[jj1 + 14] * 
		    tt) * f1;
	    if (*imoistdep > 0) {
		i__1 = *js - 1;
		f1 = rmd_(nmatd, nsd, &m, &i__1, &c__12, &dmoist[
			dmoist_offset], &c__4, &wdep[wdep_offset], &thimob, 
			imoistdep);
	    }
	    gaml1pi = chpar[jj1 + 14 + m * chpar_dim1] * exp(tdep[jj1 + 14] * 
		    tt) * f1;
	    if (*imoistdep > 0) {
		i__1 = *js - 1;
		f1 = rmd_(nmatd, nsd, &m, &i__1, &c__5, &dmoist[dmoist_offset]
			, &c__5, &wdep[wdep_offset], &thw, imoistdep);
	    }
	    gams1p = chpar[jj1 + 15 + m * chpar_dim1] * exp(tdep[jj1 + 15] * 
		    tt) * f1;
	    if (*imoistdep > 0) {
		i__1 = *js - 1;
		f1 = rmd_(nmatd, nsd, &m, &i__1, &c__13, &dmoist[
			dmoist_offset], &c__5, &wdep[wdep_offset], &thimob, 
			imoistdep);
	    }
	    gams1pi = chpar[jj1 + 15 + m * chpar_dim1] * exp(tdep[jj1 + 15] * 
		    tt) * f1;
	    if (*imoistdep > 0) {
		i__1 = *js - 1;
		f1 = rmd_(nmatd, nsd, &m, &i__1, &c__6, &dmoist[dmoist_offset]
			, &c__6, &wdep[wdep_offset], &thw, imoistdep);
	    }
	    gamg1p = chpar[jj1 + 16 + m * chpar_dim1] * exp(tdep[jj1 + 16] * 
		    tt) * f1;
	    if (*lbact) {
		gaml1p = 0.f;
		gaml1pi = 0.f;
		gams1p = 0.f;
		gams1pi = 0.f;
		gamg1p = 0.f;
	    }
	}
	if (*level == *nlevel) {
	    tto = (tempo[i__] + 273.15f - tr) / r__ / (tempo[i__] + 273.15f) /
		     tr;
	    xkso = chpar[jjj + 7 + m * chpar_dim1] * exp(tdep[jjj + 7] * tto);
	    xnuo = chpar[jjj + 8 + m * chpar_dim1] * exp(tdep[jjj + 8] * tto);
	    fexpo = chpar[jjj + 9 + m * chpar_dim1];
	    henryo = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 10] * 
		    tto);
	    f1 = 1.f;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__1, &dmoist[dmoist_offset], &
			c__1, &wdep[wdep_offset], &tho[i__], imoistdep);
	    }
	    gamlo = chpar[jjj + 11 + m * chpar_dim1] * exp(tdep[jjj + 11] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__10, &dmoist[dmoist_offset], 
			&c__1, &wdep[wdep_offset], &thimobo, imoistdep);
	    }
	    gamloi = chpar[jjj + 11 + m * chpar_dim1] * exp(tdep[jjj + 11] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__2, &dmoist[dmoist_offset], &
			c__2, &wdep[wdep_offset], &tho[i__], imoistdep);
	    }
	    gamso = chpar[jjj + 12 + m * chpar_dim1] * exp(tdep[jjj + 12] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__11, &dmoist[dmoist_offset], 
			&c__2, &wdep[wdep_offset], &thimobo, imoistdep);
	    }
	    gamsoi = chpar[jjj + 12 + m * chpar_dim1] * exp(tdep[jjj + 12] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__4, &dmoist[dmoist_offset], &
			c__4, &wdep[wdep_offset], &tho[i__], imoistdep);
	    }
	    gaml1o = chpar[jjj + 14 + m * chpar_dim1] * exp(tdep[jjj + 14] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__12, &dmoist[dmoist_offset], 
			&c__4, &wdep[wdep_offset], &thimobo, imoistdep);
	    }
	    gaml1oi = chpar[jjj + 14 + m * chpar_dim1] * exp(tdep[jjj + 14] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__5, &dmoist[dmoist_offset], &
			c__5, &wdep[wdep_offset], &tho[i__], imoistdep);
	    }
	    gams1o = chpar[jjj + 15 + m * chpar_dim1] * exp(tdep[jjj + 15] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__13, &dmoist[dmoist_offset], 
			&c__5, &wdep[wdep_offset], &thimobo, imoistdep);
	    }
	    gams1oi = chpar[jjj + 15 + m * chpar_dim1] * exp(tdep[jjj + 15] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__7, &dmoist[dmoist_offset], &
			c__7, &wdep[wdep_offset], &tho[i__], imoistdep);
	    }
	    xmulo = chpar[jjj + 17 + m * chpar_dim1] * exp(tdep[jjj + 17] * 
		    tto) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__8, &dmoist[dmoist_offset], &
			c__8, &wdep[wdep_offset], &tho[i__], imoistdep);
	    }
	    xmuso = chpar[jjj + 18 + m * chpar_dim1] * exp(tdep[jjj + 18] * 
		    tto) * f1;
	    omegao = chpar[jjj + 20 + m * chpar_dim1] * exp(tdep[jjj + 20] * 
		    tto);
	    dks = (xks - xkso) / *dt;
	    dnu = (xnu - xnuo) / *dt;
	    ddexp = (fexp - fexpo) / *dt;
	    dhenry = (henry - henryo) / *dt;
	    if (i__ != 1) {
		tti = (tempn[k] + 273.15f - tr) / r__ / (tempn[k] + 273.15f) /
			 tr;
	    }
	    if (i__ != *numnp) {
		ttj = (tempn[j] + 273.15f - tr) / r__ / (tempn[j] + 273.15f) /
			 tr;
	    }
	    if (*lbact) {
		gams1o = 0.f;
		gams1oi = 0.f;
		xmulo = 0.f;
		xmuso = 0.f;
		xmugo = 0.f;
		omegao = 0.f;
		smax2o = chpar[jjj + 15 + m * chpar_dim1] * exp(tdep[jjj + 15]
			 * tto);
		rka2o = chpar[jjj + 16 + m * chpar_dim1] * exp(tdep[jjj + 16] 
			* tto);
		rkd2o = chpar[jjj + 17 + m * chpar_dim1] * exp(tdep[jjj + 17] 
			* tto);
		smax1o = chpar[jjj + 18 + m * chpar_dim1] * exp(tdep[jjj + 18]
			 * tto);
		rka1o = chpar[jjj + 19 + m * chpar_dim1] * exp(tdep[jjj + 19] 
			* tto);
		rkd1o = chpar[jjj + 20 + m * chpar_dim1] * exp(tdep[jjj + 20] 
			* tto);
		ipsi1 = 0;
		ipsi2 = 0;
		if (! (*lfiltr)) {
		    ipsi2 = (integer) chpar[jjj + 13 + m * chpar_dim1];
		}
		if (! (*lfiltr)) {
		    ipsi1 = (integer) chpar[jjj + 14 + m * chpar_dim1];
		}
		if (ipsi1 == 0 && smax1o > 0.f) {
		    ipsi1 = 1;
		}
		if (ipsi2 == 0 && smax2o > 0.f) {
		    ipsi2 = 1;
		}
		if (ipsi1 >= 3 || ipsi2 >= 3) {
		    dc = chpar[jjj + 6 + m * chpar_dim1] * exp(tdep[jjj + 6] *
			     tt);
		}
		if (ipsi1 == 5 || ipsi2 == 5) {
		    aa = chpar[jjj + 15 + m * chpar_dim1];
		}
		psi1o = 1.f;
		psi2o = 1.f;
		if (ipsi1 > 0) {
		    blocking_(&ipsi1, &smax1o, &psi1o, &x[i__], &sorb[*js + 
			    i__ * sorb_dim1], &dc, &aa);
		}
		if (ipsi2 > 0) {
		    blocking_(&ipsi2, &smax2o, &psi2o, &x[i__], &sorb2[*js + 
			    i__ * sorb2_dim1], &dc, &aa);
		}
		if (*lfiltr) {
		    gaml1o = 0.f;
		    dc = chpar[jjj + 13 + m * chpar_dim1] * exp(tdep[jjj + 13]
			     * tto);
		    dp = chpar[jjj + 14 + m * chpar_dim1] * exp(tdep[jjj + 14]
			     * tto);
		    alfa1 = rka1o;
		    alfa2 = rka2o;
		    deposit_(&rka1o, &rka2o, &dc, &dp, &alfa1, &alfa2, &thwo, 
			    &vo[i__], &tempo[i__], xconv, tconv);
		}
	    }
	} else {
	    ttn = (tempn[i__] + 273.15f - tr) / r__ / (tempn[i__] + 273.15f) /
		     tr;
	    xksn = chpar[jjj + 7 + m * chpar_dim1] * exp(tdep[jjj + 7] * ttn);
	    xnun = chpar[jjj + 8 + m * chpar_dim1] * exp(tdep[jjj + 8] * ttn);
	    fexpn = chpar[jjj + 9 + m * chpar_dim1];
/* *exp(TDep(jjj+ 9)*TTN) */
	    henryn = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 10] * 
		    ttn);
	    dks = (xksn - xks) / *dt;
	    dnu = (xnun - xnu) / *dt;
	    ddexp = (fexpn - fexp) / *dt;
	    dhenry = (henryn - henry) / *dt;
	    if (i__ != 1) {
		tti = (tempo[k] + 273.15f - tr) / r__ / (tempo[k] + 273.15f) /
			 tr;
	    }
	    if (i__ != *numnp) {
		ttj = (tempo[j] + 273.15f - tr) / r__ / (tempo[j] + 273.15f) /
			 tr;
	    }
	}
	if (i__ != 1) {
	    henryi = chpar[jjj + 10 + matnum[k] * chpar_dim1] * exp(tdep[jjj 
		    + 10] * tti);
	}
	if (i__ != *numnp) {
	    henryj = chpar[jjj + 10 + matnum[j] * chpar_dim1] * exp(tdep[jjj 
		    + 10] * ttj);
	}
	dsconc = 1.f;
	dconc = 1.f;
	sconcp = 1.f;
	sconc = 1.f;
	sconco = 1.f;
	dretard = 0.f;
	sconcs = 1.f;
	sconcos = 1.f;
	dsconcs = 1.f;
	dconcs = 1.f;
	sconcps = 1.f;
	dretards = 0.f;
/*       Effects of nonlinear adsorption */
	if (! llinear[*js]) {
	    cc = conc[*js + i__ * conc_dim1];
	    cmid = (conc[*js + i__ * conc_dim1] + cnew[i__]) / 2.f;
	    if (*level == *nlevel) {
		cc = cnew[i__];
	    }
	    if (cc > 0.f) {
		d__1 = (doublereal) cc;
		d__2 = (doublereal) (fexp - 1.f);
		d__3 = (doublereal) cc;
		d__4 = (doublereal) fexp;
/* Computing 2nd power */
		r__1 = xnu * pow_dd(&d__3, &d__4) + 1.f;
		dsconc = fexp * pow_dd(&d__1, &d__2) / (r__1 * r__1);
		d__1 = (doublereal) cc;
		d__2 = (doublereal) (fexp - 1.f);
		d__3 = (doublereal) cc;
		d__4 = (doublereal) fexp;
		sconc = pow_dd(&d__1, &d__2) / (xnu * pow_dd(&d__3, &d__4) + 
			1.f);
	    }
	    if (cmid > 0.f) {
		d__1 = (doublereal) cmid;
		d__2 = (doublereal) (fexp - 1.f);
		d__3 = (doublereal) cmid;
		d__4 = (doublereal) fexp;
/* Computing 2nd power */
		r__1 = xnu * pow_dd(&d__3, &d__4) + 1.f;
		dconc = fexp * pow_dd(&d__1, &d__2) / (r__1 * r__1);
		d__1 = (doublereal) cmid;
		d__2 = (doublereal) fexp;
		d__3 = (doublereal) cmid;
		d__4 = (doublereal) fexp;
		d__5 = (doublereal) cmid;
		d__6 = (doublereal) (fexp * 2.f);
		d__7 = (doublereal) cmid;
		d__8 = (doublereal) fexp;
/* Computing 2nd power */
		r__1 = xnu * pow_dd(&d__7, &d__8) + 1.f;
		d__9 = (doublereal) cmid;
		d__10 = (doublereal) fexp;
		d__11 = (doublereal) cmid;
		d__12 = (doublereal) fexp;
/* Computing 2nd power */
		r__2 = xnu * pow_dd(&d__11, &d__12) + 1.f;
		dretard = pow_dd(&d__1, &d__2) / (xnu * pow_dd(&d__3, &d__4) 
			+ 1.f) * dks - xks * pow_dd(&d__5, &d__6) / (r__1 * 
			r__1) * dnu + xks * log(cmid) * pow_dd(&d__9, &d__10) 
			/ (r__2 * r__2) * ddexp;
	    }
	    if (*level == *nlevel && ! (*lequil) && conc[*js + i__ * 
		    conc_dim1] > 0.f) {
		d__1 = (doublereal) conc[*js + i__ * conc_dim1];
		d__2 = (doublereal) (fexpo - 1.f);
		d__3 = (doublereal) conc[*js + i__ * conc_dim1];
		d__4 = (doublereal) fexpo;
		sconco = pow_dd(&d__1, &d__2) / (xnuo * pow_dd(&d__3, &d__4) 
			+ 1.f);
	    }
	    if (lmobim[m] || *idualpor > 0) {
/* mobile-immobile mode */
		ss = sorb[*js + i__ * sorb_dim1];
		smid = (sorb[*js + i__ * sorb_dim1] + sorbn[i__]) / 2.f;
		if (*level == *nlevel) {
		    ss = sorbn[i__];
		}
		if (ss > 0.f) {
		    d__1 = (doublereal) ss;
		    d__2 = (doublereal) (fexp - 1.f);
		    d__3 = (doublereal) ss;
		    d__4 = (doublereal) fexp;
/* Computing 2nd power */
		    r__1 = xnu * pow_dd(&d__3, &d__4) + 1.f;
		    dsconcs = fexp * pow_dd(&d__1, &d__2) / (r__1 * r__1);
		    d__1 = (doublereal) ss;
		    d__2 = (doublereal) (fexp - 1.f);
		    d__3 = (doublereal) ss;
		    d__4 = (doublereal) fexp;
		    sconcs = pow_dd(&d__1, &d__2) / (xnu * pow_dd(&d__3, &
			    d__4) + 1.f);
		}
		if (smid > 0.f) {
		    d__1 = (doublereal) smid;
		    d__2 = (doublereal) (fexp - 1.f);
		    d__3 = (doublereal) smid;
		    d__4 = (doublereal) fexp;
/* Computing 2nd power */
		    r__1 = xnu * pow_dd(&d__3, &d__4) + 1.f;
		    dconcs = fexp * pow_dd(&d__1, &d__2) / (r__1 * r__1);
		    d__1 = (doublereal) smid;
		    d__2 = (doublereal) fexp;
		    d__3 = (doublereal) smid;
		    d__4 = (doublereal) fexp;
		    d__5 = (doublereal) smid;
		    d__6 = (doublereal) (fexp * 2.f);
		    d__7 = (doublereal) smid;
		    d__8 = (doublereal) fexp;
/* Computing 2nd power */
		    r__1 = xnu * pow_dd(&d__7, &d__8) + 1.f;
		    d__9 = (doublereal) smid;
		    d__10 = (doublereal) fexp;
		    d__11 = (doublereal) smid;
		    d__12 = (doublereal) fexp;
/* Computing 2nd power */
		    r__2 = xnu * pow_dd(&d__11, &d__12) + 1.f;
		    dretards = pow_dd(&d__1, &d__2) / (xnu * pow_dd(&d__3, &
			    d__4) + 1.f) * dks - xks * pow_dd(&d__5, &d__6) / 
			    (r__1 * r__1) * dnu + xks * log(smid) * pow_dd(&
			    d__9, &d__10) / (r__2 * r__2) * ddexp;
		}
		if (*level == *nlevel && ! (*lequil) && sorb[*js + i__ * 
			sorb_dim1] > 0.f) {
		    d__1 = (doublereal) sorb[*js + i__ * sorb_dim1];
		    d__2 = (doublereal) (fexpo - 1.f);
		    d__3 = (doublereal) sorb[*js + i__ * sorb_dim1];
		    d__4 = (doublereal) fexpo;
		    sconcos = pow_dd(&d__1, &d__2) / (xnuo * pow_dd(&d__3, &
			    d__4) + 1.f);
		}
	    }
	} else {
	    if (conc[*js + i__ * conc_dim1] > 0.f) {
		dretard = conc[*js + i__ * conc_dim1] * dks;
	    }
	    if (lmobim[m] || *idualpor > 0) {
/* mobile-immobile mode */
		if (sorb[*js + i__ * sorb_dim1] > 0.f) {
		    dretards = sorb[*js + i__ * sorb_dim1] * dks;
		}
	    }
	}
	if (*js > 1) {
	    if (! llinear[*js - 1]) {
		if (cprev > 0.f) {
		    d__1 = (doublereal) cprev;
		    d__2 = (doublereal) (fexpp - 1.f);
		    d__3 = (doublereal) cprev;
		    d__4 = (doublereal) fexpp;
		    sconcp = pow_dd(&d__1, &d__2) / (xnup * pow_dd(&d__3, &
			    d__4) + 1.f);
		}
		if (sorb[*js - 1 + i__ * sorb_dim1] > 0.f) {
		    d__1 = (doublereal) sorb[*js - 1 + i__ * sorb_dim1];
		    d__2 = (doublereal) (fexpp - 1.f);
		    d__3 = (doublereal) sorb[*js - 1 + i__ * sorb_dim1];
		    d__4 = (doublereal) fexpp;
		    sconcps = pow_dd(&d__1, &d__2) / (xnup * pow_dd(&d__3, &
			    d__4) + 1.f);
		}
	    }
	}
/*       Calculate the retardation factors */
	retard[i__] = (ro * frac * f_em__ * xks * dconc + thg * henry) / thw 
		+ 1.f;
/*       Calculate the dispersion coefficients */
	disper_(&i__, &m, nmat, nsd, numnp, dt, ltort, lartd, lupw, &lmobim[1]
		, idualpor, level, nlevel, &chpar[chpar_offset], &retard[1], &
		disp[1], &thsat[1], &thimob, &thw, &thg, &v, &dw, &dg, &henry,
		 &ro, &frac, &xks, &fexp, &xnu, &cmid, &dsconc, pecr, &taug, 
		itort, lbact);
/*       Calculate the adsorbed concentration on kinetic sites or */
/*       the concentration in an imobile zone, before solving matrix equation */
	if (! (*lequil)) {
	    nequil_(&i__, js, nsd, numnp, nmat, &m, &conc[conc_offset], &sorb[
		    sorb_offset], &sorb2[sorb2_offset], &cnew[1], &sorbn[1], &
		    sorbn2[1], &ssorb, &ssorb2, &lmobim[1], &llinear[1], 
		    lbact, level, nlevel, dt, epsi, &ro, &xks, &xkso, &cc, &
		    sconc, &sconco, &sconcs, &sconcos, &dsconcs, &xmul, &
		    xmulo, &xmus, &xmuso, &dretards, &gamli, &gaml1i, &gamloi,
		     &gaml1oi, &gamsi, &gams1i, &gamsoi, &gams1oi, &omega, &
		    omegao, &rka1, &rka1o, &rka2, &rka2o, &rkd1, &rkd1o, &
		    rkd2, &rkd2o, &thw, &thwo, &psi1, &psi1o, &psi2, &psi2o, &
		    dmobi, &frac, &thimob, &thimobo, idualpor, &sinkim[1], &
		    flmacro, lnequil, &gaml1pi, &gams1pi, &xksp, &sconcps, 
		    ldualneq, &f_em__, &omegas);
	}
/*       Calculate zero-order coefficient g0 */
	g0[i__] = xmul * thw + frac * f_em__ * ro * xmus + thg * xmug - ssink[
		i__];
	q0[i__] = xmul * thw + ro * xmus + thg * xmug;
	if (! (*lequil)) {
	    if ((lmobim[m] || *idualpor > 0) && ! (*lbact)) {
		g0[i__] += omega * ssorb;
		if (*idualpor > 0 && sinkim[i__] <= 0.f) {
		    g0[i__] -= flmacro;
		}
		if (*ldualneq) {
		    g0[i__] += omegas * ro * ssorb2;
		}
	    } else if (! (*lbact)) {
		g0[i__] += omega * ro * ssorb;
	    } else if (*lbact) {
		g0[i__] = g0[i__] + rkd1 * ro * ssorb + rkd2 * ro * ssorb2;
	    }
	}
	if (*js > 1) {
	    cg = cprev * (gaml1p * thw + ro * frac * f_em__ * xksp * gams1p * 
		    sconcp + thg * henryp * gamg1p);
	    cg1 = cg;
	    if (! (*lequil)) {
		if ((lmobim[m] || *idualpor > 0) && ! (*lbact)) {
		    aa = sorb[*js - 1 + i__ * sorb_dim1] * (thimob * gaml1pi 
			    + (1.f - frac) * ro * gams1pi * xksp * sconcps);
		    if (! (*lnequil)) {
			cg += aa;
		    }
		    cg1 += aa;
		    if (*ldualneq) {
			aa = gams1pi * ro * sorb2[*js - 1 + i__ * sorb2_dim1];
			if (! (*lnequil)) {
			    cg += aa;
			}
			cg1 += aa;
		    }
		} else if (! (*lbact)) {
		    aa = gams1pi * ro * sorb[*js - 1 + i__ * sorb_dim1];
		    if (! (*lnequil)) {
			cg += aa;
		    }
		    cg1 += aa;
		} else if (*lbact) {
		    s_wsle(&io___161);
		    do_lio(&c__9, &c__1, "Attachment/dettachment model is im"
			    "plemented \r only for one solute", (ftnlen)65);
		    e_wsle();
		    s_wsle(&io___162);
		    do_lio(&c__9, &c__1, "Press Enter to continue", (ftnlen)
			    23);
		    e_wsle();
		    s_rsle(&io___163);
		    e_rsle();
		    s_stop("", (ftnlen)0);
		}
	    }
	    g0[i__] += cg;
	    q0[i__] += cg1;
	}
	if (cmid > 0.f) {
	    g0[i__] -= ro * frac * f_em__ * dretard;
	}
/*       Calculate first-order coefficient g1 */
	g1[i__] = -(gaml + gaml1) * thw - (gams + gams1) * ro * frac * f_em__ 
		* xks * sconc - (gamg + gamg1) * thg * henry;
/*        if(Level.eq.NLevel) g1(i)=g1(i)-ThG*dHenry-Henry*(ThWO-ThW)/dt */
	if (! (*lequil)) {
	    if ((lmobim[m] || *idualpor > 0) && ! (*lbact)) {
/* mobile- */
		g1[i__] -= omega;
		if (*idualpor > 0 && sinkim[i__] > 0.f) {
		    g1[i__] -= sinkim[i__];
		}
		if (*level == *nlevel && llinear[*js]) {
		    g1[i__] += omega * *dt * omega / dmobi;
		}
		if (*ldualneq) {
		    g1[i__] -= omegas * ro * frac * (1.f - f_em__) * sconc * 
			    xks;
		    if (*level == *nlevel && llinear[*js]) {
			g1[i__] += omegas * ro * (*dt * omegas * frac * (1.f 
				- f_em__) * xks / (*dt * (omegas + gamsi + 
				gams1i) + 2.f));
		    }
		}
	    } else if (! (*lbact)) {
/* two-sit */
		g1[i__] -= omega * ro * (1.f - frac) * sconc * xks;
		if (*level == *nlevel && llinear[*js]) {
		    g1[i__] += omega * ro * (*dt * omega * (1.f - frac) * xks 
			    / (*dt * (omega + gamsi + gams1i) + 2.f));
		}
	    } else if (*lbact) {
/* filtrat */
		g1[i__] -= thw * (rka1 * psi1 + rka2 * psi2);
		if (*level == *nlevel && llinear[*js]) {
		    g1[i__] += *dt * thw * (rkd1 * rka1 / (*dt * (rkd1 + 
			    gamsi) + 2.f) + rkd2 * rka2 / (*dt * (rkd2 + 
			    gamsi) + 2.f));
		}
	    }
	}
	q1[i__] = (-(gaml + gaml1) * thw - (gams + gams1) * ro * frac * 
		f_em__ * xks * sconc - (gamg + gamg1) * thg * henry) * conc[*
		js + i__ * conc_dim1];
	if (! (*lequil)) {
	    if ((lmobim[m] || *idualpor > 0) && ! (*lbact)) {
		q1[i__] -= sorb[*js + i__ * sorb_dim1] * (thimob * (gamli + 
			gaml1i) + (1.f - frac) * ro * xks * sconcs * (gamsi + 
			gams1i));
		if (*ldualneq) {
		    q1[i__] -= (gamsi + gams1i) * ro * sorb2[*js + i__ * 
			    sorb2_dim1];
		}
	    } else if (! (*lbact)) {
		q1[i__] -= (gamsi + gams1i) * ro * sorb[*js + i__ * sorb_dim1]
			;
	    } else if (*lbact) {
		q1[i__] -= ro * (gamsi + gams1i) * (sorb[*js + i__ * 
			sorb_dim1] + sorb2[*js + i__ * sorb2_dim1]);
	    }
	}
/*       Velocity corrections */
	if (i__ == 1) {
	    dx = x[2] - x[1];
	    derk = (henryj - henry) / dx;
	} else if (i__ == *numnp) {
	    dx = x[*numnp] - x[*numnp - 1];
	    derk = (henry - henryi) / dx;
	} else {
	    dx = (x[j] - x[k]) / 2.f;
	    derk = (henryj - henryi) / dx;
	}
	vcorr[i__] = thg * dg * taug * derk;
	if (*level == 1) {
	    vo[i__] -= vcorr[i__];
	}
	if (*level == *nlevel) {
	    vn[i__] -= vcorr[i__];
	}
/*       Calculate the maximum local Peclet and Courant numbers */
	pecour_(&i__, &j, numnp, level, nlevel, lupw, lartd, dt, &x[1], &v, &
		wc[1], &thw, &vj, &thj, &disp[1], &retard[1], peclet, courant,
		 &courmax, pecr, dtmaxc, iter, epsi);
/* L11: */
    }
    return 0;
} /* coeff_ */

/* *********************************************************************** */
/* Subroutine */ int matset_(integer *js, integer *n, integer *ns, integer *
	nsd, integer *level, real *epsi, real *alf, real *dt, integer *kbotch,
	 integer *ktopch, real *cbot, real *ctop, real *x, real *tho, real *
	thn, real *vo, real *vn, real *conc, real *disp, real *retard, real *
	wc, real *g0, real *g1, doublereal *b, doublereal *d__, doublereal *e,
	 doublereal *f, real *e1, real *d1, real *f1, real *bn, real *dn, 
	real *fn, integer *nmat, real *chpar, real *tempo, real *tempn, real *
	tdep, real *dsurf, real *catm, integer *matnum, logical *lmobim, 
	integer *idualpor, logical *lvapor, real *rbot, logical *lbact)
{
    /* System generated locals */
    integer conc_dim1, conc_offset, chpar_dim1, chpar_offset, i__1;
    real r__1;

    /* Local variables */
    static integer i__, m;
    static real r__, a1, b1, f2, dg, fe, dx, tr, tt;
    static integer jjj;
    static real henry, thimob;

    /* Fortran I/O blocks */
    static cilist io___169 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    --matnum;
    --tempn;
    --tempo;
    --f;
    --e;
    --d__;
    --b;
    --g1;
    --g0;
    --wc;
    --retard;
    --disp;
    --vn;
    --vo;
    --thn;
    --tho;
    --x;
    --ctop;
    --cbot;
    --tdep;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    --lmobim;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	if (lmobim[m] && *idualpor == 0 || *lbact) {
	    thimob = chpar[m * chpar_dim1 + 4];
	    if (thimob > tho[i__]) {
		s_wsle(&io___169);
		do_lio(&c__9, &c__1, "Warning !!! ThImob > Theta", (ftnlen)26)
			;
		e_wsle();
	    }
/* Computing MAX */
	    r__1 = thn[i__] - thimob;
	    thn[i__] = dmax(r__1,.001f);
/* Computing MAX */
	    r__1 = tho[i__] - thimob;
	    tho[i__] = dmax(r__1,.001f);
	}
/* L10: */
    }
/*     Lower boundary condition */
    b1 = x[2] - x[1];
    if (*level == 1) {
	*f1 = conc[*js + conc_dim1] * (b1 / 2.f / *dt * tho[1] * retard[1] + *
		alf * (-(tho[1] * disp[1] + tho[2] * disp[2]) / b1 / 2.f - ((
		wc[1] * 3.f + 2.f) * vo[1] + vo[2]) / 6.f + b1 / 12.f * (g1[1]
		 * 3.f + g1[2]))) + conc[*js + (conc_dim1 << 1)] * *alf * ((
		tho[1] * disp[1] + tho[2] * disp[2]) / b1 / 2.f - (vo[1] + (
		2.f - wc[1] * 3.f) * vo[2]) / 6.f + b1 / 12.f * (g1[1] + g1[2]
		)) + *alf * b1 / 6.f * (g0[1] * 2.f + g0[2]);
/*       3. type  BC */
	if (*kbotch == -1) {
	    f[1] = *f1 + *alf * cbot[*js] * vo[1];
	}
    } else {
	*e1 = *epsi * (-(thn[1] * disp[1] + thn[2] * disp[2]) / b1 / 2.f + (
		vn[1] + (2.f - wc[1] * 3.f) * vn[2]) / 6.f - b1 / 12.f * (g1[
		1] + g1[2]));
	*d1 = b1 / 2.f / *dt * thn[1] * retard[1] + *epsi * ((thn[1] * disp[1]
		 + thn[2] * disp[2]) / b1 / 2.f + ((wc[1] * 3.f + 2.f) * vn[1]
		 + vn[2]) / 6.f - b1 / 12.f * (g1[1] * 3.f + g1[2]));
	f2 = *epsi * b1 / 6.f * (g0[1] * 2.f + g0[2]);
	*f1 += f2;
/*       1.type BC */
	if (*kbotch == 1) {
	    d__[1] = 1.f;
	    e[1] = 0.f;
	    f[1] = cbot[*js];
	}
/*       3. type  BC */
	if (*kbotch == -1) {
	    if (vn[1] > 0.f || *lvapor && *rbot == 0.f) {
		e[1] = *e1;
		d__[1] = *d1;
		f[1] = f[1] + f2 + *epsi * cbot[*js] * vn[1];
	    } else {
		d__[1] = -1.f;
		e[1] = 1.f;
		f[1] = 0.f;
	    }
	}
/*       Free drainage */
	if (*kbotch == 0) {
	    d__[1] = -1.f;
	    e[1] = 1.f;
	    f[1] = 0.f;
	}
    }
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	a1 = b1;
	b1 = x[i__ + 1] - x[i__];
	dx = (x[i__ + 1] - x[i__ - 1]) / 2.f;
	if (*level == 1) {
	    f[i__] = conc[*js + (i__ - 1) * conc_dim1] * *alf * ((tho[i__ - 1]
		     * disp[i__ - 1] + tho[i__] * disp[i__]) / a1 / 2.f + ((
		    wc[i__ - 1] * 3.f + 2.f) * vo[i__ - 1] + vo[i__]) / 6.f + 
		    a1 / 12.f * (g1[i__ - 1] + g1[i__])) + conc[*js + i__ * 
		    conc_dim1] * (dx / *dt * tho[i__] * retard[i__] + *alf * (
		    -(tho[i__ - 1] * disp[i__ - 1] + tho[i__] * disp[i__]) / 
		    a1 / 2.f - (tho[i__ + 1] * disp[i__ + 1] + tho[i__] * 
		    disp[i__]) / b1 / 2.f - (vo[i__ + 1] + (wc[i__ - 1] + wc[
		    i__]) * 3.f * vo[i__] - vo[i__ - 1]) / 6.f + (a1 * (g1[
		    i__ - 1] + g1[i__] * 3.f) + b1 * (g1[i__] * 3.f + g1[i__ 
		    + 1])) / 12.f)) + conc[*js + (i__ + 1) * conc_dim1] * *
		    alf * ((tho[i__ + 1] * disp[i__ + 1] + tho[i__] * disp[
		    i__]) / b1 / 2.f - (vo[i__] + (2.f - wc[i__] * 3.f) * vo[
		    i__ + 1]) / 6.f + b1 / 12.f * (g1[i__] + g1[i__ + 1])) + *
		    alf * (a1 * (g0[i__ - 1] + g0[i__] * 2.f) + b1 * (g0[i__] 
		    * 2.f + g0[i__ + 1])) / 6.f;
	} else {
	    b[i__] = *epsi * (-(thn[i__ - 1] * disp[i__ - 1] + thn[i__] * 
		    disp[i__]) / a1 / 2.f - ((wc[i__ - 1] * 3.f + 2.f) * vn[
		    i__ - 1] + vn[i__]) / 6.f - a1 / 12.f * (g1[i__ - 1] + g1[
		    i__]));
	    d__[i__] = dx / *dt * thn[i__] * retard[i__] + *epsi * ((thn[i__ 
		    - 1] * disp[i__ - 1] + thn[i__] * disp[i__]) / a1 / 2.f + 
		    (thn[i__ + 1] * disp[i__ + 1] + thn[i__] * disp[i__]) / 
		    b1 / 2.f + (vn[i__ + 1] + (wc[i__ - 1] + wc[i__]) * 3.f * 
		    vn[i__] - vn[i__ - 1]) / 6.f - (a1 * (g1[i__ - 1] + g1[
		    i__] * 3.f) + b1 * (g1[i__] * 3.f + g1[i__ + 1])) / 12.f);
	    e[i__] = *epsi * (-(thn[i__ + 1] * disp[i__ + 1] + thn[i__] * 
		    disp[i__]) / b1 / 2.f + (vn[i__] + (2.f - wc[i__] * 3.f) *
		     vn[i__ + 1]) / 6.f - b1 / 12.f * (g1[i__] + g1[i__ + 1]))
		    ;
	    f[i__] += *epsi * (a1 * (g0[i__ - 1] + g0[i__] * 2.f) + b1 * (g0[
		    i__] * 2.f + g0[i__ + 1])) / 6.f;
	}
/* L11: */
    }
/*     Upper boundary condition */
    if (*level == 1) {
	*fn = conc[*js + (*n - 1) * conc_dim1] * *alf * ((tho[*n - 1] * disp[*
		n - 1] + tho[*n] * disp[*n]) / b1 / 2.f + ((wc[*n - 1] * 3.f 
		+ 2.f) * vo[*n - 1] + vo[*n]) / 6.f + b1 / 12.f * (g1[*n - 1] 
		+ g1[*n])) + conc[*js + *n * conc_dim1] * (b1 / 2.f / *dt * 
		tho[*n] * retard[*n] + *alf * (-(tho[*n - 1] * disp[*n - 1] + 
		tho[*n] * disp[*n]) / b1 / 2.f + (vo[*n - 1] + (2.f - wc[*n - 
		1] * 3.f) * vo[*n]) / 6.f + b1 / 12.f * (g1[*n - 1] + g1[*n] *
		 3))) + *alf * b1 / 6.f * (g0[*n - 1] + g0[*n] * 2.f);
/*       3. type BC */
	if (*ktopch <= 0) {
	    f[*n] = *fn;
	    if (vo[*n] < 0.f) {
		f[*n] -= *alf * vo[*n] * ctop[*js];
	    }
	    if (*ktopch == -2) {
		m = matnum[*n];
		tr = 293.15f;
		r__ = 8.314f;
		jjj = *js - 1 << 4;
		tt = (tempo[*n] + 273.15f - tr) / r__ / (tempo[*n] + 273.15f) 
			/ tr;
		dg = chpar[jjj + 6 + m * chpar_dim1] * exp(tdep[jjj + 6] * tt)
			;
		henry = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 10] 
			* tt);
		f[*n] = f[*n] - *alf * dg / *dsurf * henry * conc[*js + *n * 
			conc_dim1] + dg / *dsurf * *catm;
	    }
	}
    } else {
	*bn = *epsi * (-(thn[*n - 1] * disp[*n - 1] + thn[*n] * disp[*n]) / 
		b1 / 2.f - ((wc[*n - 1] * 3.f + 2.f) * vn[*n - 1] + vn[*n]) / 
		6.f - b1 / 12.f * (g1[*n - 1] + g1[*n]));
	*dn = b1 / 2.f / *dt * thn[*n] * retard[*n] + *epsi * ((thn[*n - 1] * 
		disp[*n - 1] + thn[*n] * disp[*n]) / b1 / 2.f - (vn[*n - 1] + 
		(2.f - wc[*n - 1] * 3.f) * vn[*n]) / 6.f - b1 / 12.f * (g1[*n 
		- 1] + g1[*n] * 3.f));
	fe = *epsi * b1 / 6.f * (g0[*n - 1] + g0[*n] * 2.f);
	*fn += fe;
/*       1. type BC */
	if (*ktopch > 0) {
	    b[*n] = 0.f;
	    d__[*n] = 1.f;
	    f[*n] = ctop[*js];
/*       3. type BC */
	} else {
	    b[*n] = *bn;
	    d__[*n] = *dn;
	    f[*n] += fe;
	    if (vn[*n] < 0.f) {
		f[*n] -= *epsi * vn[*n] * ctop[*js];
	    }
	    if (*ktopch == -2) {
		m = matnum[*n];
		tr = 293.15f;
		r__ = 8.314f;
		jjj = *js - 1 << 4;
		tt = (tempn[*n] + 273.15f - tr) / r__ / (tempn[*n] + 273.15f) 
			/ tr;
		dg = chpar[jjj + 6 + m * chpar_dim1] * exp(tdep[jjj + 6] * tt)
			;
		henry = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 10] 
			* tt);
		d__[*n] += *epsi * dg / *dsurf * henry;
	    }
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	if (lmobim[m] && *idualpor == 0 || *lbact) {
	    thimob = chpar[m * chpar_dim1 + 4];
	    thn[i__] += thimob;
	    tho[i__] += thimob;
	}
/* L12: */
    }
    return 0;
} /* matset_ */

/* ************************************************************************ */
/*     Solve matrix equation */
/* Subroutine */ int bansol_(integer *n, doublereal *a, doublereal *b, 
	doublereal *c__, doublereal *f)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;

    /* Parameter adjustments */
    --f;
    --c__;
    --b;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	b[i__] -= a[i__] * c__[i__ - 1] / b[i__ - 1];
	f[i__] -= a[i__] * f[i__ - 1] / b[i__ - 1];
/* L11: */
    }
    f[*n] /= b[*n];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = *n - i__ + 1;
	f[j] = (f[j] - c__[j] * f[j + 1]) / b[j];
/* L12: */
    }
    return 0;
} /* bansol_ */

/* ************************************************************************ */
/*     Calculate the maximum local Peclet and Courant numbers */
/* Subroutine */ int pecour_(integer *i__, integer *j, integer *numnp, 
	integer *level, integer *nlevel, logical *lupw, logical *lartd, real *
	dt, real *x, real *v, real *wc, real *thw, real *vj, real *thj, real *
	disp, real *retard, real *peclet, real *courant, real *courmax, real *
	pecr, real *dtmaxc, integer *iter, real *epsi)
{
    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static real dd, dx, vv, pe2, vv1, pec, rmin, cour, vmax, cour1, dtmax;

    /* Parameter adjustments */
    --retard;
    --disp;
    --wc;
    --x;

    /* Function Body */
    if (*i__ != *numnp) {
	dx = x[*j] - x[*i__];
	vv = 0.f;
	if (*thw > 1e-6f && *thj > 1e-6f) {
	    vv = (dabs(*v) / *thw + dabs(*vj) / *thj) / 2.f;
	}
	vv1 = 0.f;
	if (*thw > 1e-6f && *thj > 1e-6f) {
	    vv1 = (*v / *thw + *vj / *thj) / 2.f;
	}
	dd = (disp[*i__] + disp[*j]) / 2.f;
	if (*level == *nlevel) {
	    pec = 99999.f;
	    dtmax = 1e30f;
/*          vMax=amax1(abs(v)/ThW,abs(vj)/Thj) */
	    vmax = (dabs(*v) + dabs(*vj)) / (*thw + *thj);
/* Computing MIN */
	    r__1 = retard[*i__], r__2 = retard[*j];
	    rmin = dmin(r__1,r__2);
	    if (dd > 0.f) {
		pec = dabs(vv) * dx / dd;
	    }
	    cour = vmax * *dt / dx / rmin;
	    *peclet = dmax(*peclet,pec);
	    *courant = dmax(*courant,cour);
	    cour1 = *courmax;
	    if (! (*lupw) && ! (*lartd)) {
		if (pec != 99999.f) {
/* Computing MIN */
		    r__1 = 1.f, r__2 = *pecr / dmax(.5f,pec);
		    cour1 = dmin(r__1,r__2);
		}
	    }
	    if (*epsi < 1.f && vmax > 1e-20f) {
		dtmax = cour1 * dx * rmin / vmax;
	    }
/*         the von Neumann time step limit */
/*          RThE=(ThW+thj)/2.*RMin */
/*          if(abs(DD).gt.1.e-20)dtMax=amin1(dtMax,10.*RThE*dx*dx/2./DD) */
	    *dtmaxc = dmin(*dtmaxc,dtmax);
/*       Calculate upstream weighting factors */
	} else if (*lupw && *iter == 1) {
	    pe2 = 11.f;
	    if (dd > 0.f) {
		pe2 = dx * vv1 / dd / 2.f;
	    }
	    if (dabs(vv) < 1e-30f) {
		wc[*i__] = 0.f;
	    } else if (dabs(pe2) > 10.f) {
		if (vv1 > 0.f) {
		    wc[*i__] = 1.f;
		}
		if (vv1 < 0.f) {
		    wc[*i__] = -1.f;
		}
	    } else {
		wc[*i__] = 1.f / ((exp(pe2) - exp(-pe2)) / (exp(pe2) + exp(
			-pe2))) - 1.f / pe2;
/* Computing MIN */
		r__1 = 1.f, r__2 = wc[*i__];
		wc[*i__] = dmin(r__1,r__2);
/* Computing MAX */
		r__1 = -1.f, r__2 = wc[*i__];
		wc[*i__] = dmax(r__1,r__2);
	    }
	}
    }
    return 0;
} /* pecour_ */

/* ************************************************************************ */
/*     Calculate the dispersion coefficients */
/* Subroutine */ int disper_(integer *i__, integer *m, integer *nmat, integer 
	*nsd, integer *numnp, real *dt, logical *ltort, logical *lartd, 
	logical *lupw, logical *lmobim, integer *idualpor, integer *level, 
	integer *nlevel, real *chpar, real *retard, real *disp, real *thsat, 
	real *thimob, real *thw, real *thg, real *v, real *dw, real *dg, real 
	*henry, real *ro, real *frac, real *xks, real *fexp, real *xnu, real *
	cmid, real *dsconc, real *pecr, real *taug, integer *itort, logical *
	lbact)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset;
    real r__1, r__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static real dd, fi, ths, dpom, tauw;

    /* Parameter adjustments */
    --thsat;
    --lmobim;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;
    --disp;
    --retard;

    /* Function Body */
    if (*ltort) {
	ths = thsat[*m];
	if (lmobim[*m] && *idualpor == 0 || *lbact) {
/* Computing MAX */
	    r__1 = thsat[*m] - *thimob;
	    ths = dmax(r__1,.001f);
	}
	if (*idualpor > 0) {
	    ths = thsat[*m] + *thimob;
	}
	if (*itort == 0) {
	    d__1 = (doublereal) (*thw);
/* Computing 2nd power */
	    r__1 = ths;
	    tauw = pow_dd(&d__1, &c_b89) / (r__1 * r__1);
	    d__1 = (doublereal) (*thg);
/* Computing 2nd power */
	    r__1 = ths;
	    *taug = pow_dd(&d__1, &c_b89) / (r__1 * r__1);
	} else {
	    d__1 = (doublereal) (*thw / ths);
	    tauw = pow_dd(&d__1, &c_b91) * .66f;
	    d__1 = (doublereal) (*thg);
	    *taug = pow_dd(&d__1, &c_b92) / ths;
	}
    } else {
	tauw = 1.f;
	*taug = 1.f;
    }
    disp[*i__] = chpar[*m * chpar_dim1 + 2] * dabs(*v) / *thw + *dw * tauw + *
	    thg / *thw * *dg * *henry * *taug;
    if (! (*lartd) && ! (*lupw)) {
	fi = 0.f;
	if (*cmid > 0.f) {
	    d__1 = (doublereal) (*cmid);
	    d__2 = (doublereal) (*fexp - 1.f);
	    d__3 = (doublereal) (*cmid);
	    d__4 = (doublereal) (*fexp);
/* Computing 2nd power */
	    r__1 = *xnu * pow_dd(&d__3, &d__4) + 1.f;
	    d__5 = (doublereal) (*cmid);
	    d__6 = (doublereal) (*fexp);
	    fi = *thw * 6.f * *ro * *xks * pow_dd(&d__1, &d__2) * (*fexp / (
		    r__1 * r__1) - 1.f / (*xnu * pow_dd(&d__5, &d__6) + 1.f));
	}
/* Computing MAX */
	r__1 = *dt / (*thw * 6.f * (*thw + *ro * *frac * *xks * *dsconc + *
		thg * *henry) + fi);
	dpom = dmax(r__1,0.f);
	if (*level != *nlevel) {
	    disp[*i__] += *v * *v * dpom;
	} else {
/* Computing MAX */
	    r__1 = disp[*i__] - *v * *v * dpom, r__2 = disp[*i__] / 2.f;
	    disp[*i__] = dmax(r__1,r__2);
	}
    }
    if (*lartd) {
	dd = 0.f;
	if (*pecr != 0.f && dabs(*v) > 1e-15f) {
	    dd = *v * *v * *dt / *thw / *thw / retard[*i__] / *pecr;
	}
	if (dd > disp[*i__]) {
	    disp[*i__] = dd;
	}
    }
    return 0;
} /* disper_ */

/* ************************************************************************ */
/*     Calculate the adsorbed concentration on kinetic sites or */
/*     the concentration in an imobile zone, before solving matrix equation */
/* Subroutine */ int nequil_(integer *i__, integer *js, integer *nsd, integer 
	*numnp, integer *nmat, integer *m, real *conc, real *sorb, real *
	sorb2, real *cnew, real *sorbn, real *sorbn2, real *ssorb, real *
	ssorb2, logical *lmobim, logical *llinear, logical *lbact, integer *
	level, integer *nlevel, real *dt, real *epsi, real *ro, real *xks, 
	real *xkso, real *cc, real *sconc, real *sconco, real *sconcs, real *
	sconcos, real *dsconcs, real *xmul, real *xmulo, real *xmus, real *
	xmuso, real *dretards, real *gaml, real *gaml1, real *gamlo, real *
	gaml1o, real *gams, real *gams1, real *gamso, real *gams1o, real *
	omega, real *omegao, real *rka1, real *rka1o, real *rka2, real *rka2o,
	 real *rkd1, real *rkd1o, real *rkd2, real *rkd2o, real *thw, real *
	thwo, real *psi1, real *psi1o, real *psi2, real *psi2o, real *dmobi, 
	real *frac, real *thimob, real *thimobo, integer *idualpor, real *
	sinkim, real *flmacro, logical *lnequil, real *gaml1pi, real *gams1pi,
	 real *xksp, real *sconcps, logical *ldualneq, real *f_em__, real *
	omegas)
{
    /* System generated locals */
    integer conc_dim1, conc_offset, sorb_dim1, sorb_offset, sorb2_dim1, 
	    sorb2_offset;

    /* Local variables */
    static real cs, amobi, bmobi, emobi, gmobi0, bmobio, dtheta, emobio;

    /* Parameter adjustments */
    --llinear;
    --sinkim;
    --sorbn2;
    --sorbn;
    --cnew;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    sorb_dim1 = *nsd;
    sorb_offset = 1 + sorb_dim1;
    sorb -= sorb_offset;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    --lmobim;

    /* Function Body */
    *flmacro = 0.f;
    if (*idualpor > 0) {
/* mobile-immobile model */
	if (sinkim[*i__] > 0.f) {
	    *flmacro = sinkim[*i__] * conc[*js + *i__ * conc_dim1];
	} else {
	    *flmacro = sinkim[*i__] * sorb[*js + *i__ * sorb_dim1];
	}
    }
    *ssorb = sorb[*js + *i__ * sorb_dim1];
    *ssorb2 = sorb2[*js + *i__ * sorb2_dim1];
    if (*level == *nlevel) {
/*       mobile-immobile model */
	if ((lmobim[*m] || *idualpor > 0) && ! (*lbact)) {
	    amobi = (*thimob + *thimobo) / 2.f + (1.f - *frac) * *ro * *xks * 
		    *dsconcs;
	    dtheta = *thimob - *thimobo;
	    emobi = *thimob * *xmul + (1.f - *frac) * *ro * *xmus - (1.f - *
		    frac) * *ro * *dretards;
	    emobio = *thimobo * *xmulo + (1.f - *frac) * *ro * *xmuso - (1.f 
		    - *frac) * *ro * *dretards;
	    emobi += *flmacro;
	    emobio += *flmacro;
	    if (*lnequil) {
		cs = sorb[*js - 1 + *i__ * sorb_dim1] * (*thimob * *gaml1pi + 
			(1.f - *frac) * *ro * *gams1pi * *xksp * *sconcps);
		emobi += cs;
		emobio += cs;
	    }
	    bmobi = *thimob * (*gaml + *gaml1) + (1.f - *frac) * *ro * (*gams 
		    + *gams1) * *xks * *sconcs;
	    bmobio = *thimobo * (*gamlo + *gaml1o) + (1.f - *frac) * *ro * (*
		    gamso + *gams1o) * *xkso * *sconcos;
	    if (llinear[*js]) {
		*dmobi = amobi * 2.f + *dt * (*omega + bmobi) + dtheta;
		gmobi0 = (amobi * 2.f - *dt * (*omegao + bmobio) - dtheta) / *
			dmobi;
		sorb[*js + *i__ * sorb_dim1] = sorb[*js + *i__ * sorb_dim1] * 
			gmobi0 + *dt * (*omegao * conc[*js + *i__ * conc_dim1]
			 + emobi + emobio) / *dmobi;
		*ssorb = sorb[*js + *i__ * sorb_dim1];
	    } else {
		sorbn[*i__] = sorb[*js + *i__ * sorb_dim1] + *dt / amobi * (*
			epsi * (*omega * (cnew[*i__] - sorbn[*i__]) - bmobi * 
			sorbn[*i__] + emobi) + (1.f - *epsi) * (*omegao * (
			conc[*js + *i__ * conc_dim1] - sorb[*js + *i__ * 
			sorb_dim1]) - bmobio * sorb[*js + *i__ * sorb_dim1] + 
			emobio));
		*ssorb = sorbn[*i__];
	    }
	    if (*ldualneq) {
		cs = 0.f;
		if (*lnequil) {
		    cs = *gams1pi * sorb2[*js - 1 + *i__ * sorb2_dim1] * *dt *
			     2.f;
		}
		if (llinear[*js]) {
		    sorb2[*js + *i__ * sorb2_dim1] = ((2.f - (*omegas + *
			    gamso + *gams1o) * *dt) * sorb2[*js + *i__ * 
			    sorb2_dim1] + *dt * *frac * (1.f - *f_em__) * *
			    omegas * *xkso * conc[*js + *i__ * conc_dim1] + *
			    dt * (1.f - *f_em__) * (*xmuso + *xmus) + cs) / (*
			    dt * (*omegas + *gams + *gams1) + 2.f);
		    *ssorb2 = sorb2[*js + *i__ * sorb2_dim1];
		} else {
		    sorbn2[*i__] = sorb2[*js + *i__ * sorb2_dim1] + *dt * (*
			    epsi * (*omegas * (*frac * (1.f - *f_em__) * *
			    sconc * *xks * *cc - sorbn2[*i__]) - (*gams + *
			    gams1) * sorbn2[*i__] + (1.f - *f_em__) * *xmus) 
			    + (1.f - *epsi) * (*omegas * (*frac * (1.f - *
			    f_em__) * *sconco * *xkso * conc[*js + *i__ * 
			    conc_dim1] - *ssorb2) - (*gamso + *gams1o) * *
			    ssorb2 + (1.f - *f_em__) * *xmuso));
		    *ssorb2 = sorbn2[*i__];
		}
	    }
/*       two-site sorption model */
	} else if (! (*lbact)) {
	    cs = 0.f;
	    if (*lnequil) {
		cs = *gams1pi * sorb[*js - 1 + *i__ * sorb_dim1] * *dt * 2.f;
	    }
	    if (llinear[*js]) {
		sorb[*js + *i__ * sorb_dim1] = ((2.f - (*omegao + *gamso + *
			gams1o) * *dt) * sorb[*js + *i__ * sorb_dim1] + *dt * 
			(1.f - *frac) * *omegao * *xkso * conc[*js + *i__ * 
			conc_dim1] + *dt * (1.f - *frac) * (*xmuso + *xmus) + 
			cs) / (*dt * (*omega + *gams + *gams1) + 2.f);
		*ssorb = sorb[*js + *i__ * sorb_dim1];
	    } else {
		sorbn[*i__] = sorb[*js + *i__ * sorb_dim1] + *dt * (*epsi * (*
			omega * ((1.f - *frac) * *sconc * *xks * *cc - sorbn[*
			i__]) - (*gams + *gams1) * sorbn[*i__] + (1.f - *frac)
			 * *xmus) + (1.f - *epsi) * (*omegao * ((1.f - *frac) 
			* *sconco * *xkso * conc[*js + *i__ * conc_dim1] - *
			ssorb) - (*gamso + *gams1o) * *ssorb + (1.f - *frac) *
			 *xmuso));
		*ssorb = sorbn[*i__];
	    }
/*       filtration model */
	} else if (*lbact) {
	    if (llinear[*js]) {
		sorb[*js + *i__ * sorb_dim1] = ((2.f - *dt * (*rkd1o + *gamso 
			+ *gams1o)) * sorb[*js + *i__ * sorb_dim1] + *dt * *
			rka1o * *thw * conc[*js + *i__ * conc_dim1] / *ro) / (
			*dt * (*rkd1 + *gams + *gams1) + 2.f);
		*ssorb = sorb[*js + *i__ * sorb_dim1];
		sorb2[*js + *i__ * sorb2_dim1] = ((2.f - *dt * (*rkd2o + *
			gamso + *gams1o)) * sorb2[*js + *i__ * sorb2_dim1] + *
			dt * *rka2o * *thw * conc[*js + *i__ * conc_dim1] / *
			ro) / (*dt * (*rkd2 + *gams + *gams1) + 2.f);
		*ssorb2 = sorb2[*js + *i__ * sorb2_dim1];
	    } else {
		sorbn[*i__] = sorb[*js + *i__ * sorb_dim1] + *dt * (*epsi * (*
			rka1 * *thw / *ro * *psi1 * *cc - (*rkd1 + *gams + *
			gams1) * sorbn[*i__]) + (1.f - *epsi) * (*rka1o * *
			thwo / *ro * *psi1o * conc[*js + *i__ * conc_dim1] - (
			*rkd1o + *gamso + *gams1o) * sorb[*js + *i__ * 
			sorb_dim1]));
		*ssorb = sorbn[*i__];
		sorbn2[*i__] = sorb2[*js + *i__ * sorb2_dim1] + *dt * (*epsi *
			 (*rka2 * *thw / *ro * *psi2 * *cc - (*rkd2 + *gams + 
			*gams1) * sorbn2[*i__]) + (1.f - *epsi) * (*rka2o * *
			thwo / *ro * *psi2o * conc[*js + *i__ * conc_dim1] - (
			*rkd2o + *gamso + *gams1o) * sorb2[*js + *i__ * 
			sorb2_dim1]));
		*ssorb2 = sorbn2[*i__];
	    }
	}
    }
    return 0;
} /* nequil_ */

/* ************************************************************************ */
/*     Calculate sorbed concentration for linear noneq. adsorption or */
/*     concentration in the imobile water. */
/*     At the end of the time step. After solving matrix equation */
/* Subroutine */ int sorbconc_(integer *js, integer *nsd, integer *n, integer 
	*matnum, real *temp, logical *lmobim, real *chpar, real *tdep, real *
	sorb, real *conc, real *dt, integer *nmat, logical *lbact, real *thw, 
	real *sorb2, logical *lfiltr, real *veloc, integer *idualpor, real *
	thim, real *thoim, real *sinkim, real *strans, integer *imoistdep, 
	integer *nmatd, real *dmoist, real *wdep, real *xconv, real *tconv, 
	logical *ldualneq)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset, conc_dim1, conc_offset, sorb_dim1, 
	    sorb_offset, sorb2_dim1, sorb2_offset, dmoist_dim1, dmoist_dim2, 
	    dmoist_offset, wdep_dim1, wdep_offset, i__1;

    /* Local variables */
    static integer i__, m;
    static real r__, f1, dc, dp, ro, tr, tt;
    static integer jjj;
    extern doublereal rmd_(integer *, integer *, integer *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *);
    static real xks, rka1, rka2, rkd1, rkd2, f_em__, frac, gaml, gams, alfa1, 
	    alfa2, gaml1, gams1, amobi, omega, bmobi, dmobi, theta, dtheta, 
	    omegas, thimob, flmacro, thimobo;
    extern /* Subroutine */ int deposit_(real *, real *, real *, real *, real 
	    *, real *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --tdep;
    --strans;
    --sinkim;
    --thoim;
    --thim;
    --veloc;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    --thw;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    sorb_dim1 = *nsd;
    sorb_offset = 1 + sorb_dim1;
    sorb -= sorb_offset;
    --temp;
    --matnum;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;
    --lmobim;
    wdep_dim1 = 2 + *nmatd;
    wdep_offset = 1 + wdep_dim1;
    wdep -= wdep_offset;
    dmoist_dim1 = *nmatd;
    dmoist_dim2 = *nsd;
    dmoist_offset = 1 + dmoist_dim1 * (1 + dmoist_dim2 * 14);
    dmoist -= dmoist_offset;

    /* Function Body */
    tr = 293.15f;
    r__ = 8.314f;
    jjj = *js - 1 << 4;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	tt = (temp[i__] + 273.15f - tr) / r__ / (temp[i__] + 273.15f) / tr;
	frac = chpar[m * chpar_dim1 + 3] * exp(tdep[3] * tt);
	xks = chpar[jjj + 7 + m * chpar_dim1] * exp(tdep[jjj + 7] * tt);
	f1 = 1.f;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__11, &dmoist[dmoist_offset], &
		    c__2, &wdep[wdep_offset], &thw[i__], imoistdep);
	}
	gams = chpar[jjj + 12 + m * chpar_dim1] * exp(tdep[jjj + 12] * tt) * 
		f1;
	if (*imoistdep > 0) {
	    f1 = rmd_(nmatd, nsd, &m, js, &c__13, &dmoist[dmoist_offset], &
		    c__5, &wdep[wdep_offset], &thw[i__], imoistdep);
	}
	gams1 = chpar[jjj + 15 + m * chpar_dim1] * exp(tdep[jjj + 15] * tt) * 
		f1;
	omega = chpar[jjj + 20 + m * chpar_dim1] * exp(tdep[jjj + 20] * tt);
	if ((lmobim[m] || *idualpor > 0) && ! (*lbact)) {
/* mobile-im */
	    ro = chpar[m * chpar_dim1 + 1] * exp(tdep[1] * tt);
	    if (*idualpor == 0) {
		thimob = chpar[m * chpar_dim1 + 4] * exp(tdep[4] * tt);
		thimobo = thimob;
	    } else if (*idualpor > 0) {
		thimob = thim[i__];
		thimobo = thoim[i__];
	    }
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__10, &dmoist[dmoist_offset], 
			&c__1, &wdep[wdep_offset], &thimob, imoistdep);
	    }
	    gaml = chpar[jjj + 11 + m * chpar_dim1] * exp(tdep[jjj + 11] * tt)
		     * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__12, &dmoist[dmoist_offset], 
			&c__4, &wdep[wdep_offset], &thimob, imoistdep);
	    }
	    gaml1 = chpar[jjj + 14 + m * chpar_dim1] * exp(tdep[jjj + 14] * 
		    tt) * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__11, &dmoist[dmoist_offset], 
			&c__2, &wdep[wdep_offset], &thimob, imoistdep);
	    }
	    gams = chpar[jjj + 12 + m * chpar_dim1] * exp(tdep[jjj + 12] * tt)
		     * f1;
	    if (*imoistdep > 0) {
		f1 = rmd_(nmatd, nsd, &m, js, &c__13, &dmoist[dmoist_offset], 
			&c__5, &wdep[wdep_offset], &thimob, imoistdep);
	    }
	    gams1 = chpar[jjj + 15 + m * chpar_dim1] * exp(tdep[jjj + 15] * 
		    tt) * f1;
	    dtheta = thimob - thimobo;
	    amobi = (thimob + thimobo) / 2.f + (1.f - frac) * ro * xks;
	    bmobi = thimob * (gaml + gaml1) + (gams + gams1) * ro * (1.f - 
		    frac) * xks;
	    dmobi = amobi * 2.f + *dt * (omega + bmobi) + dtheta;
	    sorb[*js + i__ * sorb_dim1] += *dt * omega * conc[*js + i__ * 
		    conc_dim1] / dmobi;
	    flmacro = 0.f;
	    if (*idualpor > 0) {
		if (sinkim[i__] > 0.f) {
		    flmacro = sinkim[i__] * conc[*js + i__ * conc_dim1];
		} else {
		    flmacro = sinkim[i__] * sorb[*js + i__ * sorb_dim1];
		}
	    }
	    if (*js == 1) {
		strans[i__] = omega * (conc[*js + i__ * conc_dim1] - sorb[*js 
			+ i__ * sorb_dim1]) + flmacro;
	    }
	    if (*ldualneq) {
		f_em__ = chpar[jjj + 13 + m * chpar_dim1] * exp(tdep[jjj + 13]
			 * tt);
		omegas = chpar[jjj + 16 + m * chpar_dim1] * exp(tdep[jjj + 16]
			 * tt);
		sorb2[*js + i__ * sorb2_dim1] += *dt * omegas * frac * (1.f - 
			f_em__) * xks * conc[*js + i__ * conc_dim1] / (*dt * (
			omegas + gams + gams1) + 2.f);
	    }
	} else if (! (*lbact)) {
/* two-site sorption mod */
	    sorb[*js + i__ * sorb_dim1] += *dt * omega * (1.f - frac) * xks * 
		    conc[*js + i__ * conc_dim1] / (*dt * (omega + gams + 
		    gams1) + 2.f);
	} else if (*lbact) {
/* filtration model */
	    ro = chpar[m * chpar_dim1 + 1] * exp(tdep[1] * tt);
	    thimob = chpar[m * chpar_dim1 + 4];
	    theta = thw[i__] - thimob;
	    gams = chpar[jjj + 12 + m * chpar_dim1] * exp(tdep[jjj + 12] * tt)
		    ;
	    gams1 = 0.f;
	    rka1 = chpar[jjj + 19 + m * chpar_dim1] * exp(tdep[jjj + 19] * tt)
		    ;
	    rkd1 = chpar[jjj + 20 + m * chpar_dim1] * exp(tdep[jjj + 20] * tt)
		    ;
	    rka2 = chpar[jjj + 16 + m * chpar_dim1] * exp(tdep[jjj + 16] * tt)
		    ;
	    rkd2 = chpar[jjj + 17 + m * chpar_dim1] * exp(tdep[jjj + 17] * tt)
		    ;
	    if (*lfiltr) {
		dc = chpar[jjj + 13 + m * chpar_dim1] * exp(tdep[jjj + 13] * 
			tt);
		dp = chpar[jjj + 14 + m * chpar_dim1] * exp(tdep[jjj + 14] * 
			tt);
		alfa1 = rka1;
		alfa2 = rka2;
		deposit_(&rka1, &rka2, &dc, &dp, &alfa1, &alfa2, &theta, &
			veloc[i__], &temp[i__], xconv, tconv);
	    }
	    sorb[*js + i__ * sorb_dim1] += *dt * rka1 * theta * conc[*js + 
		    i__ * conc_dim1] / ro / (*dt * (rkd1 + gams + gams1) + 
		    2.f);
	    sorb2[*js + i__ * sorb2_dim1] += *dt * rka2 * theta * conc[*js + 
		    i__ * conc_dim1] / ro / (*dt * (rkd2 + gams + gams1) + 
		    2.f);
	}
/* L11: */
    }
    return 0;
} /* sorbconc_ */

/* ************************************************************************ */
/*     Calculate mass-transfer fluxes at the end of the time interval */
/* Subroutine */ int masstran_(integer *js, integer *ns, integer *nsd, 
	integer *n, integer *matnum, real *temp, logical *lmobim, logical *
	lequil, real *chpar, real *tdep, real *sorb, real *conc, integer *
	nmat, real *x, real *cvch0, real *cvch1, real *cvchr, real *cvchim, 
	real *epsi, real *q0, real *q1, real *ssink, logical *lbact, real *
	theta, real *sorb2, logical *lfiltr, real *veloc, integer *idualpor, 
	real *sinkim, real *xconv, real *tconv, logical *ldualneq, real *
	strans, logical *llinear)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset, conc_dim1, conc_offset, sorb_dim1, 
	    sorb_offset, sorb2_dim1, sorb2_offset, i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    extern /* Subroutine */ int blocking_(integer *, real *, real *, real *, 
	    real *, real *, real *);
    static real flmacroi, flmacroj;
    static integer i__, j;
    static real r__, aa, dc, dp;
    static integer mi, mj;
    static real dx, tr, cci, ccj;
    static integer jjj;
    static real roi, roj, tti, ttj, rka1, rka2, xksi, xksj, xnui, xnuj, alfa1,
	     alfa2, rka1i, rka1j, rka2i, rkd1i, rkd1j, rka2j, rkd2i, rkd2j;
    static integer ipsi1, ipsi2;
    static real psi1i, psi1j, psi2i, psi2j, f_emi__, f_emj__, fraci, fracj, 
	    fexpi, fexpj, smax1i, smax1j, smax2i, smax2j, omegai, omegaj, 
	    thetai, thetaj, sorbei, sorbej, omegasi, omegasj, thimobi, 
	    thimobj;
    extern /* Subroutine */ int deposit_(real *, real *, real *, real *, real 
	    *, real *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --llinear;
    --cvchim;
    --cvchr;
    --cvch1;
    --cvch0;
    --tdep;
    --strans;
    --sinkim;
    --veloc;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    --theta;
    --ssink;
    --q1;
    --q0;
    --x;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    sorb_dim1 = *nsd;
    sorb_offset = 1 + sorb_dim1;
    sorb -= sorb_offset;
    --temp;
    --matnum;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;
    --lmobim;

    /* Function Body */
    tr = 293.15f;
    r__ = 8.314f;
    jjj = *js - 1 << 4;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mi = matnum[i__];
	tti = (temp[i__] + 273.15f - tr) / r__ / (temp[i__] + 273.15f) / tr;
	if (i__ == *n) {
	    goto L10;
	}
	j = i__ + 1;
	dx = x[j] - x[i__];
	cvch0[*js] += *epsi * dx * (q0[i__] + q0[j]) / 2.f;
	cvch1[*js] += *epsi * dx * (q1[i__] + q1[j]) / 2.f;
	cvchr[*js] += *epsi * dx * (ssink[i__] + ssink[j]) / 2.f;
	if (! (*lequil)) {
	    mj = matnum[j];
	    ttj = (temp[j] + 273.15f - tr) / r__ / (temp[j] + 273.15f) / tr;
	    omegai = chpar[jjj + 20 + mi * chpar_dim1] * exp(tdep[jjj + 20] * 
		    tti);
	    omegaj = chpar[jjj + 20 + mj * chpar_dim1] * exp(tdep[jjj + 20] * 
		    ttj);
/*         mobile-immobile model */
	    if ((lmobim[mi] || *idualpor > 0) && ! (*lbact)) {
		flmacroi = 0.f;
		flmacroj = 0.f;
		if (*idualpor > 0) {
		    if (sinkim[i__] > 0.f) {
			flmacroi = sinkim[i__] * conc[*js + i__ * conc_dim1];
		    } else {
			flmacroi = sinkim[i__] * sorb[*js + i__ * sorb_dim1];
		    }
		    if (sinkim[j] > 0.f) {
			flmacroj = sinkim[j] * conc[*js + j * conc_dim1];
		    } else {
			flmacroj = sinkim[j] * sorb[*js + j * sorb_dim1];
		    }
		}
		cvchim[*js] += *epsi * dx / 2.f * (omegai * (conc[*js + i__ * 
			conc_dim1] - sorb[*js + i__ * sorb_dim1]) + omegaj * (
			conc[*js + j * conc_dim1] - sorb[*js + j * sorb_dim1])
			 + flmacroi + flmacroj);
		if (*js == 1 && ! llinear[*js]) {
		    strans[i__] = *epsi * (omegai * (conc[*js + i__ * 
			    conc_dim1] - sorb[*js + i__ * sorb_dim1]) + 
			    flmacroi);
		}
		if (*ldualneq) {
		    roi = chpar[mi * chpar_dim1 + 1] * exp(tdep[1] * tti);
		    fraci = chpar[mi * chpar_dim1 + 3] * exp(tdep[3] * tti);
		    xksi = chpar[jjj + 7 + mi * chpar_dim1] * exp(tdep[jjj + 
			    7] * tti);
		    xnui = chpar[jjj + 8 + mi * chpar_dim1] * exp(tdep[jjj + 
			    8] * tti);
		    fexpi = chpar[jjj + 9 + mi * chpar_dim1];
/* *exp(TDep(jjj+ 9)*TTi) */
		    roj = chpar[mj * chpar_dim1 + 1] * exp(tdep[1] * ttj);
		    fracj = chpar[mj * chpar_dim1 + 3] * exp(tdep[3] * ttj);
		    xksj = chpar[jjj + 7 + mj * chpar_dim1] * exp(tdep[jjj + 
			    7] * ttj);
		    xnuj = chpar[jjj + 8 + mj * chpar_dim1] * exp(tdep[jjj + 
			    8] * ttj);
		    fexpj = chpar[jjj + 9 + mj * chpar_dim1];
/* *exp(TDep(jjj+ 9)*TTj) */
		    f_emi__ = chpar[jjj + 13 + mi * chpar_dim1] * exp(tdep[
			    jjj + 13] * tti);
		    f_emj__ = chpar[jjj + 13 + mj * chpar_dim1] * exp(tdep[
			    jjj + 13] * ttj);
		    omegasi = chpar[jjj + 16 + mi * chpar_dim1] * exp(tdep[
			    jjj + 16] * tti);
		    omegasj = chpar[jjj + 16 + mj * chpar_dim1] * exp(tdep[
			    jjj + 16] * ttj);
		    cci = conc[*js + i__ * conc_dim1];
		    ccj = conc[*js + j * conc_dim1];
		    sorbei = 0.f;
		    sorbej = 0.f;
		    if (cci > 0.f) {
			d__1 = (doublereal) cci;
			d__2 = (doublereal) fexpi;
			d__3 = (doublereal) cci;
			d__4 = (doublereal) fexpi;
			sorbei = fraci * (1.f - f_emi__) * xksi * pow_dd(&
				d__1, &d__2) / (xnui * pow_dd(&d__3, &d__4) + 
				1.f);
		    }
		    if (ccj > 0.f) {
			d__1 = (doublereal) ccj;
			d__2 = (doublereal) fexpj;
			d__3 = (doublereal) ccj;
			d__4 = (doublereal) fexpj;
			sorbej = fracj * (1.f - f_emj__) * xksj * pow_dd(&
				d__1, &d__2) / (xnuj * pow_dd(&d__3, &d__4) + 
				1.f);
		    }
		    cvchim[*js] += *epsi * dx / 2.f * (roi * omegasi * (
			    sorbei - sorb2[*js + i__ * sorb2_dim1]) + roj * 
			    omegasj * (sorbej - sorb2[*js + j * sorb2_dim1]));
		    if (*js == 1 && ! llinear[*js]) {
			strans[i__] += *epsi * roi * omegasi * (sorbei - 
				sorb2[*js + i__ * sorb2_dim1]);
		    }
		}
/*         two-site sorption model */
	    } else if (! (*lbact)) {
		roi = chpar[mi * chpar_dim1 + 1] * exp(tdep[1] * tti);
		fraci = chpar[mi * chpar_dim1 + 3] * exp(tdep[3] * tti);
		xksi = chpar[jjj + 7 + mi * chpar_dim1] * exp(tdep[jjj + 7] * 
			tti);
		xnui = chpar[jjj + 8 + mi * chpar_dim1] * exp(tdep[jjj + 8] * 
			tti);
		fexpi = chpar[jjj + 9 + mi * chpar_dim1];
/* *exp(TDep(jjj+ 9)*TTi) */
		roj = chpar[mj * chpar_dim1 + 1] * exp(tdep[1] * ttj);
		fracj = chpar[mj * chpar_dim1 + 3] * exp(tdep[3] * ttj);
		xksj = chpar[jjj + 7 + mj * chpar_dim1] * exp(tdep[jjj + 7] * 
			ttj);
		xnuj = chpar[jjj + 8 + mj * chpar_dim1] * exp(tdep[jjj + 8] * 
			ttj);
		fexpj = chpar[jjj + 9 + mj * chpar_dim1];
/* *exp(TDep(jjj+ 9)*TTj) */
		cci = conc[*js + i__ * conc_dim1];
		ccj = conc[*js + j * conc_dim1];
		sorbei = 0.f;
		sorbej = 0.f;
		if (cci > 0.f) {
		    d__1 = (doublereal) cci;
		    d__2 = (doublereal) fexpi;
		    d__3 = (doublereal) cci;
		    d__4 = (doublereal) fexpi;
		    sorbei = (1.f - fraci) * xksi * pow_dd(&d__1, &d__2) / (
			    xnui * pow_dd(&d__3, &d__4) + 1.f);
		}
		if (ccj > 0.f) {
		    d__1 = (doublereal) ccj;
		    d__2 = (doublereal) fexpj;
		    d__3 = (doublereal) ccj;
		    d__4 = (doublereal) fexpj;
		    sorbej = (1.f - fracj) * xksj * pow_dd(&d__1, &d__2) / (
			    xnuj * pow_dd(&d__3, &d__4) + 1.f);
		}
		cvchim[*js] += *epsi * dx / 2.f * (roi * omegai * (sorbei - 
			sorb[*js + i__ * sorb_dim1]) + roj * omegaj * (sorbej 
			- sorb[*js + j * sorb_dim1]));
		if (*js == 1) {
		    strans[i__] = *epsi * roi * omegai * (sorbei - sorb[*js + 
			    i__ * sorb_dim1]);
		}
/*         filtration model */
	    } else if (*lbact) {
		roi = chpar[mi * chpar_dim1 + 1] * exp(tdep[1] * tti);
		roj = chpar[mj * chpar_dim1 + 1] * exp(tdep[1] * ttj);
		smax1i = chpar[jjj + 18 + mi * chpar_dim1] * exp(tdep[jjj + 
			18] * tti);
		smax1j = chpar[jjj + 18 + mj * chpar_dim1] * exp(tdep[jjj + 
			18] * ttj);
		rka1i = chpar[jjj + 19 + mi * chpar_dim1] * exp(tdep[jjj + 19]
			 * tti);
		rka1j = chpar[jjj + 19 + mj * chpar_dim1] * exp(tdep[jjj + 19]
			 * ttj);
		rkd1i = chpar[jjj + 20 + mi * chpar_dim1] * exp(tdep[jjj + 20]
			 * tti);
		rkd1j = chpar[jjj + 20 + mj * chpar_dim1] * exp(tdep[jjj + 20]
			 * ttj);
		smax2i = chpar[jjj + 15 + mi * chpar_dim1] * exp(tdep[jjj + 
			15] * tti);
		smax2j = chpar[jjj + 15 + mj * chpar_dim1] * exp(tdep[jjj + 
			15] * ttj);
		rka2i = chpar[jjj + 16 + mi * chpar_dim1] * exp(tdep[jjj + 16]
			 * tti);
		rka2j = chpar[jjj + 16 + mj * chpar_dim1] * exp(tdep[jjj + 16]
			 * ttj);
		rkd2i = chpar[jjj + 17 + mi * chpar_dim1] * exp(tdep[jjj + 17]
			 * tti);
		rkd2j = chpar[jjj + 17 + mj * chpar_dim1] * exp(tdep[jjj + 17]
			 * ttj);
		thimobi = chpar[mi * chpar_dim1 + 4];
		thimobj = chpar[mj * chpar_dim1 + 4];
		thetai = theta[i__] - thimobi;
		thetaj = theta[j] - thimobj;
		ipsi1 = 0;
		ipsi2 = 0;
		if (! (*lfiltr)) {
		    ipsi2 = (integer) chpar[jjj + 13 + mi * chpar_dim1];
		}
		if (! (*lfiltr)) {
		    ipsi1 = (integer) chpar[jjj + 14 + mj * chpar_dim1];
		}
		if (ipsi1 == 0 && smax1i > 0.f) {
		    ipsi1 = 1;
		}
		if (ipsi2 == 0 && smax2i > 0.f) {
		    ipsi2 = 1;
		}
		if (ipsi1 >= 3 || ipsi2 >= 3) {
		    dc = chpar[jjj + 6 + mi * chpar_dim1] * exp(tdep[jjj + 6] 
			    * tti);
		}
		if (ipsi1 == 5 || ipsi2 == 5) {
		    aa = chpar[jjj + 15 + mi * chpar_dim1];
		}
		psi1i = 1.f;
		psi1j = 1.f;
		psi2i = 1.f;
		psi2j = 1.f;
		if (ipsi1 > 0) {
		    blocking_(&ipsi1, &smax1i, &psi1i, &x[i__], &sorb[*js + 
			    i__ * sorb_dim1], &dc, &aa);
		    blocking_(&ipsi1, &smax1j, &psi1j, &x[j], &sorb[*js + j * 
			    sorb_dim1], &dc, &aa);
		}
		if (ipsi2 > 0) {
		    blocking_(&ipsi2, &smax2i, &psi2i, &x[i__], &sorb2[*js + 
			    i__ * sorb2_dim1], &dc, &aa);
		    blocking_(&ipsi2, &smax2j, &psi2j, &x[j], &sorb2[*js + j *
			     sorb2_dim1], &dc, &aa);
		}
		if (*lfiltr) {
		    dc = chpar[jjj + 13 + mi * chpar_dim1] * exp(tdep[jjj + 
			    13] * tti);
		    dp = chpar[jjj + 14 + mi * chpar_dim1] * exp(tdep[jjj + 
			    14] * tti);
		    alfa1 = rka1i;
		    alfa2 = rka2i;
		    deposit_(&rka1, &rka2, &dc, &dp, &alfa1, &alfa2, &thetai, 
			    &veloc[i__], &temp[i__], xconv, tconv);
		    rka1i = rka1;
		    rka1j = rka1;
		    rka2i = rka2;
		    rka2j = rka2;
		}
		cvchim[*js] += *epsi * dx / 2.f * (conc[*js + i__ * conc_dim1]
			 * thetai * (psi1i * rka1i + psi2i * rka2i) + conc[*
			js + j * conc_dim1] * thetaj * (psi1j * rka1j + psi2j 
			* rka2j) - roi * (sorb[*js + i__ * sorb_dim1] * rkd1i 
			+ sorb2[*js + i__ * sorb2_dim1] * rkd2i) - roj * (
			sorb[*js + j * sorb_dim1] * rkd1j + sorb2[*js + j * 
			sorb2_dim1] * rkd2j));
		if (*js == 1) {
		    strans[i__] = *epsi * (conc[*js + i__ * conc_dim1] * 
			    thetai * (psi1i * rka1i + psi2i * rka2i) - roi * (
			    sorb[*js + i__ * sorb_dim1] * rkd1i + sorb2[*js + 
			    i__ * sorb2_dim1] * rkd2i));
		}
	    }
	}
L10:
/* L11: */
	;
    }
    return 0;
} /* masstran_ */

/* ************************************************************************ */
/*     Calculate flux concentration for the first solute */
/* Subroutine */ int fluxconc_(integer *numnp, integer *nmat, integer *nsd, 
	real *x, real *v, real *theta, real *thsat, real *chpar, integer *
	matnum, real *temp, real *tdep, real *conc, real *concf, logical *
	ltort, logical *lmobim, integer *idualpor, real *thim, integer *js)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset, conc_dim1, conc_offset, i__1;
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, m;
    static real r__, dg, dw, tr, tt, qw;
    static integer jjj;
    static real thg, ths, thw, taug, disp, tauw, cgrad, henry, thimob;

    /* Parameter adjustments */
    --thim;
    --concf;
    --temp;
    --matnum;
    --theta;
    --v;
    --x;
    --lmobim;
    --thsat;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    --tdep;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;

    /* Function Body */
    jjj = *js - 1 << 4;
    tr = 293.15f;
    r__ = 8.314f;
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	thw = theta[i__];
/* Computing MAX */
	r__1 = 0.f, r__2 = thsat[m] - thw;
	thg = dmax(r__1,r__2);
	if (lmobim[m]) {
	    if (*idualpor == 0) {
		thimob = chpar[m * chpar_dim1 + 4];
/* Computing MAX */
		r__1 = thw - thimob;
		thw = dmax(r__1,.001f);
	    } else if (*idualpor > 0) {
		thimob = thim[i__];
	    }
	}
	tt = (temp[i__] + 273.15f - tr) / r__ / (temp[i__] + 273.15f) / tr;
	dw = chpar[jjj + 5 + m * chpar_dim1] * exp(tdep[jjj + 5] * tt);
	dg = chpar[jjj + 6 + m * chpar_dim1] * exp(tdep[jjj + 6] * tt);
	henry = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 10] * tt);
	if (*ltort) {
	    ths = thsat[m];
	    if (lmobim[m] && *idualpor == 0) {
/* Computing MAX */
		r__1 = thsat[m] - thimob;
		ths = dmax(r__1,.001f);
	    }
	    if (*idualpor > 0) {
		ths = thsat[m] + thimob;
	    }
	    d__1 = (doublereal) thw;
/* Computing 2nd power */
	    r__1 = ths;
	    tauw = pow_dd(&d__1, &c_b89) / (r__1 * r__1);
	    d__1 = (doublereal) thg;
/* Computing 2nd power */
	    r__1 = ths;
	    taug = pow_dd(&d__1, &c_b89) / (r__1 * r__1);
	} else {
	    tauw = 1.f;
	    taug = 1.f;
	}
	qw = v[i__];
	disp = chpar[m * chpar_dim1 + 2] * dabs(qw) / thw + dw * tauw + thg / 
		thw * dg * henry * taug;
	cgrad = 0.f;
	if (i__ == 1) {
	    cgrad = (conc[*js + (i__ + 1) * conc_dim1] - conc[*js + i__ * 
		    conc_dim1]) / (x[i__ + 1] - x[i__]);
	} else if (i__ == *numnp) {
	    cgrad = (conc[*js + i__ * conc_dim1] - conc[*js + (i__ - 1) * 
		    conc_dim1]) / (x[i__] - x[i__ - 1]);
	} else {
	    cgrad = (conc[*js + (i__ + 1) * conc_dim1] - conc[*js + (i__ - 1) 
		    * conc_dim1]) / (x[i__ + 1] - x[i__ - 1]);
	}
	concf[i__] = conc[*js + i__ * conc_dim1];
	if (qw != 0.f) {
	    concf[i__] = conc[*js + i__ * conc_dim1] - disp * thw / qw * 
		    cgrad;
	}
/* L11: */
    }
    return 0;
} /* fluxconc_ */

/* ************************************************************************ */
/*     Calculate blocking coefficient for the attachment process */
/* Subroutine */ int blocking_(integer *ipsi, real *smax, real *psi, real *x, 
	real *ss, real *dc, real *smax2)
{
    /* System generated locals */
    real r__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static real binf, minf, sinf, const__;

    *psi = 1.f;
    if (*ipsi == 1) {
	if (*smax > 0.f) {
	    *psi = 1.f - *ss / *smax;
	}
    } else if (*ipsi == 2) {
	if (*smax > 0.f) {
/* Computing MAX */
	    d__1 = (doublereal) (*ss);
	    d__2 = (doublereal) (*smax);
	    r__1 = pow_dd(&d__1, &d__2);
	    *psi = dmax(r__1,*psi);
	}
    } else if (*ipsi == 3) {
	binf = 1.f / *smax;
	sinf = .546f;
	minf = *dc;
	const__ = sinf * binf * *ss;
	if (*ss <= *smax * .8f) {
	    d__1 = (doublereal) const__;
	    d__2 = (doublereal) const__;
	    *psi = 1.f - const__ * 4.f + pow_dd(&d__1, &c_b116) * 3.08f + 
		    pow_dd(&d__2, &c_b117) * 1.4069f;
	}
	if (*ss > *smax * .8f) {
	    d__1 = (doublereal) (1.f - binf * *ss);
	    d__2 = (doublereal) minf;
	    d__3 = (doublereal) binf;
	    *psi = pow_dd(&d__1, &c_b117) / (pow_dd(&d__2, &c_b116) * 2.f * 
		    pow_dd(&d__3, &c_b117));
	}
    } else if (*ipsi == 4) {
	if (*smax > 0.f && *dc > 0.f) {
	    d__1 = (doublereal) ((dabs(*x) + *dc) / *dc);
	    d__2 = (doublereal) (-(*smax));
	    *psi = pow_dd(&d__1, &d__2);
	}
    } else if (*ipsi == 5) {
	if (*smax > 0.f && *dc > 0.f) {
	    d__1 = (doublereal) ((dabs(*x) + *dc) / *dc);
	    d__2 = (doublereal) (-(*smax));
	    *psi = pow_dd(&d__1, &d__2);
	}
	if (*smax2 > 0.f) {
	    *psi *= 1.f - *ss / *smax2;
	}
    }
    return 0;
} /* blocking_ */

/* ************************************************************************ */
/*     Calculate the deposition coefficient for the bacteria transport, */
/*     All calculations within this subroutines are in meters and seconds */
/*     Conversions are needed */
/* Subroutine */ int deposit_(real *ka1, real *ka2, real *dc1, real *dp1, 
	real *alfa1, real *alfa2, real *theta, real *q, real *temp, real *
	xconv, real *tconv)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;
    doublereal d__1, d__2;

    /* Local variables */
    static real g, h__, dc, bk, dp, as, pi, mu, n_g__, eta, n_r__, rof, rop, 
	    n_pe__, n_lo__, gamma, veloc, e_diff__, e_grav__, pveloc, 
	    e_inter__;

    /* Fortran I/O blocks */
    static cilist io___325 = { 0, 6, 0, 0, 0 };
    static cilist io___326 = { 0, 6, 0, 0, 0 };
    static cilist io___327 = { 0, 5, 0, 0, 0 };


/*     Ka       - deposition coefficient (output) [1/T] */
/*     Dc       - diameter of the sand grains (m) */
/*     Dp       - diameter of the bacteria (0.95 microm) (m) */
/*     Alfa     - sticking efficiency (-) */
/*     Theta    - porosity (-) */
/*     q        - Darcy flux [L/T] */
/*     Temp     - Temperature in Celcius */
    dc = *dc1 / *xconv;
    dp = *dp1 / *xconv;
    if (dp <= 0.f && dc <= 0.f) {
	s_wsle(&io___325);
	do_lio(&c__9, &c__1, "Both Dp and Dc are equal to zero !!!", (ftnlen)
		36);
	e_wsle();
	s_wsle(&io___326);
	do_lio(&c__9, &c__1, "Press Enter to continue", (ftnlen)23);
	e_wsle();
	s_rsle(&io___327);
	e_rsle();
	s_stop("", (ftnlen)0);
    }
    pi = 3.1415f;
/* Ludolf's number */
    mu = 9.3e-4f;
/* fluid viscosity (Pa s) */
    bk = 1.38048e-23f;
/* Boltzman constatnt (J/K) */
    h__ = 1e-20f;
/* Hamaker constant (J) */
    g = 9.81f;
/* gravitational acceleration (m/s2) */
    rop = 1080.f;
/* bacterial density (kg/m3) */
    rof = 998.f;
/* fluid density (kg/m3) */
    veloc = (r__1 = *q / *xconv * *tconv, dabs(r__1));
/* absolute value of Darcy flux (converted */
    pveloc = veloc / *theta;
/* pore velocity (converted to m/s) */
    dc = *dc1 / *xconv;
/* conversion to m */
    dp = *dp1 / *xconv;
/* conversion to m */
    if (veloc > 0.f) {
	d__1 = (doublereal) (1.f - *theta);
	gamma = pow_dd(&d__1, &c_b129);
/* Computing 5th power */
	r__1 = gamma, r__2 = r__1, r__1 *= r__1;
/* Computing 5th power */
	r__3 = gamma, r__4 = r__3, r__3 *= r__3;
/* Computing 6th power */
	r__5 = gamma, r__5 *= r__5;
	as = (1.f - r__2 * (r__1 * r__1)) * 2.f / (2.f - gamma * 3.f + r__4 * 
		(r__3 * r__3) * 3.f - r__5 * (r__5 * r__5) * 2.f);
/* Corr */
	n_pe__ = pi * 3.f * mu * dp * dc * veloc / (bk * (*temp + 273.15f));
/* Peclet number */
	d__1 = (doublereal) as;
	d__2 = (doublereal) n_pe__;
	e_diff__ = pow_dd(&d__1, &c_b129) * 4.f * pow_dd(&d__2, &c_b131);
/* removal by diffus */
/* Computing 2nd power */
	r__1 = dp;
	n_lo__ = h__ * 4.f / (pi * 9.f * mu * (r__1 * r__1) * veloc);
/* London number */
	n_r__ = dp / dc;
/* Interception numb */
	d__1 = (doublereal) n_lo__;
	d__2 = (doublereal) n_r__;
	e_inter__ = as * pow_dd(&d__1, &c_b132) * pow_dd(&d__2, &c_b133);
/* removal intercept */
/* Computing 2nd power */
	r__1 = dp;
	n_g__ = g * (rop - rof) * (r__1 * r__1) / (mu * 18.f * veloc);
/* Gravitation numbe */
	d__1 = (doublereal) n_g__;
	d__2 = (doublereal) n_r__;
	e_grav__ = as * .00338f * pow_dd(&d__1, &c_b134) * pow_dd(&d__2, &
		c_b135);
/*                                                      sedimentation */
/* removal by gravit */
    } else {
	e_diff__ = 0.f;
	e_inter__ = 0.f;
	e_grav__ = 0.f;
    }
    eta = e_diff__ + e_inter__ + e_grav__;
/*     Original Filtration Theory */
/* single-collector */
    *ka1 = (1.f - *theta) * 3.f / 2.f / dc * eta * *alfa1 * pveloc;
    *ka1 /= *tconv;
    *ka2 = (1.f - *theta) * 3.f / 2.f / dc * eta * *alfa2 * pveloc;
    *ka2 /= *tconv;
    return 0;
} /* deposit_ */

/* ************************************************************************ */
/*     Nonequilibrium phase is initially in equilibrium with liquid phase */
/* Subroutine */ int noneqinit_(integer *numnp, integer *nsd, integer *ns, 
	integer *nmat, integer *matnum, real *tdep, real *temp, real *chpar, 
	real *conc, real *sorb, logical *llinear, logical *lmobim, integer *
	idualpor, logical *lbact, real *sorb2, real *theta)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset, conc_dim1, conc_offset, sorb_dim1, 
	    sorb_offset, sorb2_dim1, sorb2_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, m;
    static real r__, cc;
    static integer js;
    static real ro, tr, tt;
    static integer jjj;
    static real xks, xnu, rka1, rka2, rkd1, rkd2, xks1, xks2, frac, fexp, 
	    sconc;

    /* Parameter adjustments */
    --theta;
    --temp;
    --matnum;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    --llinear;
    sorb_dim1 = *nsd;
    sorb_offset = 1 + sorb_dim1;
    sorb -= sorb_offset;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    --tdep;
    --lmobim;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;

    /* Function Body */
    tr = 293.15f;
    r__ = 8.314f;
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
	jjj = js - 1 << 4;
	i__2 = *numnp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    m = matnum[i__];
	    tt = (temp[i__] + 273.15f - tr) / r__ / (temp[i__] + 273.15f) / 
		    tr;
	    if (lmobim[m] || *idualpor > 0) {
		sorb[js + i__ * sorb_dim1] = conc[js + i__ * conc_dim1];
	    } else {
		ro = chpar[m * chpar_dim1 + 1] * exp(tdep[1] * tt);
		frac = chpar[m * chpar_dim1 + 3] * exp(tdep[3] * tt);
		xks = chpar[jjj + 7 + m * chpar_dim1] * exp(tdep[jjj + 7] * 
			tt);
		xnu = chpar[jjj + 8 + m * chpar_dim1] * exp(tdep[jjj + 8] * 
			tt);
		fexp = chpar[jjj + 9 + m * chpar_dim1];
/* *exp(TDep(jjj+ 9)*TT) */
		sconc = 1.f;
		cc = conc[js + i__ * conc_dim1];
		if (! llinear[js] && cc > 0.f) {
		    d__1 = (doublereal) cc;
		    d__2 = (doublereal) (fexp - 1.f);
		    d__3 = (doublereal) cc;
		    d__4 = (doublereal) fexp;
		    sconc = pow_dd(&d__1, &d__2) / (xnu * pow_dd(&d__3, &d__4)
			     + 1.f);
		}
		sorb[js + i__ * sorb_dim1] = (1.f - frac) * sconc * xks * cc;
		if (*lbact) {
		    frac = 0.f;
		    rka1 = chpar[jjj + 19 + m * chpar_dim1] * exp(tdep[jjj + 
			    19] * tt);
		    rkd1 = chpar[jjj + 20 + m * chpar_dim1] * exp(tdep[jjj + 
			    20] * tt);
		    xks1 = theta[i__] * rka1 / ro / rkd1;
		    sorb[js + i__ * sorb_dim1] = xks1 * conc[js + i__ * 
			    conc_dim1];
		    rka2 = chpar[jjj + 16 + m * chpar_dim1] * exp(tdep[jjj + 
			    16] * tt);
		    rkd2 = chpar[jjj + 17 + m * chpar_dim1] * exp(tdep[jjj + 
			    17] * tt);
		    xks2 = theta[i__] * rka1 / ro / rkd1;
		    sorb2[js + i__ * sorb2_dim1] = xks2 * conc[js + i__ * 
			    conc_dim1];
		}
	    }
/* L11: */
	}
/* L12: */
    }
    return 0;
} /* noneqinit_ */

/* *********************************************************************** */
/*     Subroutine calculating accesible water content for colloids and colloid velocity */
/*     Only for steady-state water flow and homogeneous soil profile */
/* Subroutine */ int exclusion_(integer *numnp, integer *nmat, integer *nsd, 
	real *par, real *chpar, real *thnew, real *vnew, real *thold, real *
	vold)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset, i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal sw_c_eff__;
    static integer i__, m;
    static real xl;
    static doublereal sw, r_c__;
    static real thc;
    static doublereal thr, ths, krw, swr, vgm1, vgn1, vgn2, vgm2, pc_c__;
    static real th_c__;
    static doublereal kr_c__;
    static real velc;
    static doublereal sw_c__, alpha;
    static real porvelsolute, porvelcolloid;
    static doublereal sw_eff__;

    /* Fortran I/O blocks */
    static cilist io___379 = { 0, 6, 0, 0, 0 };
    static cilist io___380 = { 0, 6, 0, 0, 0 };
    static cilist io___381 = { 0, 5, 0, 0, 0 };


    /* Parameter adjustments */
    --vold;
    --thold;
    --vnew;
    --thnew;
    par -= 12;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;

    /* Function Body */
    m = 1;
    th_c__ = chpar[m * chpar_dim1 + 4];
/* Water Content from which colloids are excluded */
    if (dabs(th_c__) < 1e-20f) {
	return 0;
    }
    chpar[m * chpar_dim1 + 4] = 0.f;
    thr = par[m * 11 + 1];
    ths = par[m * 11 + 2];
    swr = thr / ths;
    sw_c__ = th_c__ / ths;
    alpha = par[m * 11 + 3];
    vgn1 = par[m * 11 + 4];
    xl = par[m * 11 + 6];
    vgm1 = 1.f - 1.f / vgn1;
    vgn2 = vgn1 + 1;
    vgm2 = 1.f - 2.f / vgn2;
    if (sw_c__ < 0. || sw_c__ > 1.f) {
	s_wsle(&io___379);
	do_lio(&c__9, &c__1, "Problem with size exclusion!", (ftnlen)28);
	e_wsle();
	s_wsle(&io___380);
	do_lio(&c__9, &c__1, "Press Enter to continue", (ftnlen)23);
	e_wsle();
	s_rsle(&io___381);
	e_rsle();
    }
/*     Accessible Water Content to Colloid */
    thc = thnew[1] - ths * sw_c__;
    sw = thnew[1] / ths;
    sw_eff__ = (sw - swr) / (1.f - swr);
    sw_c_eff__ = (sw_c__ - swr) / (1.f - swr);
/*     Colloid Permeability according to Burdine Model */
    if (sw_eff__ > sw_c_eff__) {
	d__2 = 1.f / vgm2;
	d__1 = 1.f - pow_dd(&sw_c_eff__, &d__2);
	d__4 = 1.f / vgm2;
	d__3 = 1.f - pow_dd(&sw_eff__, &d__4);
	kr_c__ = pow_dd(&sw_eff__, &c_b116) * (pow_dd(&d__1, &vgm2) - pow_dd(&
		d__3, &vgm2));
    } else {
	kr_c__ = 0.f;
    }
/*     Mualem Water Relative Permeability */
    if (sw_eff__ > 0.f) {
	d__1 = (doublereal) xl;
	d__4 = 1 / vgm1;
	d__3 = 1 - pow_dd(&sw_eff__, &d__4);
	d__2 = 1 - pow_dd(&d__3, &vgm1);
	krw = pow_dd(&sw_eff__, &d__1) * pow_dd(&d__2, &c_b116);
	velc = vnew[1] * kr_c__ / krw;
    } else {
	krw = 0.f;
	velc = vnew[1];
    }
/*     Convert sw_c to r_c */
    d__2 = -1.f / vgm1;
    d__1 = pow_dd(&sw_c_eff__, &d__2) - 1;
    d__3 = 1.f / vgn1;
    pc_c__ = 1.f / alpha * pow_dd(&d__1, &d__3);
    r_c__ = 144.f / (pc_c__ * 981.f);
    porvelsolute = vnew[1] / thnew[1];
    porvelcolloid = velc / thc;
/*     write(*,*) "r_c microns", 10000.*r_c */
    i__1 = *numnp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	thnew[i__] = thc;
	thold[i__] = thc;
	vnew[i__] = velc;
	vold[i__] = velc;
/* L11: */
    }
    return 0;
} /* exclusion_ */

/* *********************************************************************** */
/*     Reads parameters for the function expressing reaction rate dependence */
/*     on the water content */
/* Subroutine */ int moistdepin_(char *cdatapath, char *cfilename, integer *
	nmat, integer *nmatd, integer *ns, integer *nsd, real *dmoist, 
	integer *imoistdep, ftnlen cdatapath_len, ftnlen cfilename_len)
{
    /* System generated locals */
    address a__1[2];
    integer dmoist_dim1, dmoist_dim2, dmoist_offset, i__1[2], i__2, i__3, 
	    i__4;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    extern integer len_trim__(char *, ftnlen);
    static integer i__, m, js, ilengthpath, jreact;

    /* Fortran I/O blocks */
    static cilist io___395 = { 1, 15, 0, 0, 0 };
    static cilist io___397 = { 1, 15, 0, 0, 0 };
    static cilist io___399 = { 1, 15, 0, 0, 0 };
    static cilist io___401 = { 1, 15, 0, 0, 0 };
    static cilist io___403 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    dmoist_dim1 = *nmatd;
    dmoist_dim2 = *nsd;
    dmoist_offset = 1 + dmoist_dim1 * (1 + dmoist_dim2 * 14);
    dmoist -= dmoist_offset;

    /* Function Body */
    ilengthpath = len_trim__(cdatapath, (ftnlen)260);
/* Writing concatenation */
    i__1[0] = ilengthpath, a__1[0] = cdatapath;
    i__1[1] = 11, a__1[1] = "MoistDep.in";
    s_cat(cfilename, a__1, i__1, &c__2, (ftnlen)260);
    o__1.oerr = 1;
    o__1.ounit = 15;
    o__1.ofnmlen = 260;
    o__1.ofnm = cfilename;
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__2 = f_open(&o__1);
    if (i__2 != 0) {
	goto L901;
    }
    i__2 = s_rsle(&io___395);
    if (i__2 != 0) {
	goto L902;
    }
    i__2 = e_rsle();
    if (i__2 != 0) {
	goto L902;
    }
    i__2 = *nmat;
    for (m = 1; m <= i__2; ++m) {
	i__3 = s_rsle(&io___397);
	if (i__3 != 0) {
	    goto L902;
	}
	i__3 = e_rsle();
	if (i__3 != 0) {
	    goto L902;
	}
	i__3 = *ns;
	for (js = 1; js <= i__3; ++js) {
	    i__4 = s_rsle(&io___399);
	    if (i__4 != 0) {
		goto L902;
	    }
	    i__4 = e_rsle();
	    if (i__4 != 0) {
		goto L902;
	    }
	    for (jreact = 1; jreact <= 13; ++jreact) {
		i__4 = s_rsle(&io___401);
		if (i__4 != 0) {
		    goto L902;
		}
		for (i__ = 1; i__ <= 6; ++i__) {
		    i__4 = do_lio(&c__4, &c__1, (char *)&dmoist[m + (js + (
			    jreact + i__ * 13) * dmoist_dim2) * dmoist_dim1], 
			    (ftnlen)sizeof(real));
		    if (i__4 != 0) {
			goto L902;
		    }
		}
		i__4 = e_rsle();
		if (i__4 != 0) {
		    goto L902;
		}
/* L11: */
	    }
/* L12: */
	}
/* L13: */
    }
    cl__1.cerr = 0;
    cl__1.cunit = 15;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
/*     Error opening an input file */
L901:
    if (*imoistdep == 2) {
	*imoistdep = 0;
    }
    return 0;
/*     Error reading from an input file */
L902:
    if (*imoistdep == 2) {
	*imoistdep = 0;
    }
    s_wsle(&io___403);
    do_lio(&c__9, &c__1, "Error reading from an input file MoistDep.in !!!!", 
	    (ftnlen)49);
    e_wsle();
    cl__1.cerr = 0;
    cl__1.cunit = 15;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* moistdepin_ */

/* *********************************************************************** */
doublereal rmd_(integer *nmatd, integer *nsd, integer *m, integer *js, 
	integer *jreact, real *dmoist, integer *ireact, real *wdep, real *
	theta, integer *imoistdep)
{
    /* System generated locals */
    integer dmoist_dim1, dmoist_dim2, dmoist_offset, wdep_dim1, wdep_offset;
    real ret_val;
    doublereal d__1, d__2;

    /* Local variables */
    static integer jjj;
    static real theta0, theta1, theta2, theta3, reacmin0, reacmin1;

/*     Function expressing reaction rate dependence on the water content */
/*     ReacMin0  - relative minimum rate of reaction at low water contents */
/*     Theta0    - water content at which reaction rate start increasing */
/*     Theta1    - water content at which reaction rate stops increasing */
/*     Theta2    - water content at which reaction rate start decreasing */
/*     Theta3    - water content at which reaction rate stops decreasing */
/*     ReacMin1  - relative minimum rate of reaction at high water contents */
/*     If theta2=theta3=thetaS -> Anaerobic process */
/*     If theta0=theta1=0      -> Aerobic process */
/*     If theta2=0 -> no reduction */
    /* Parameter adjustments */
    wdep_dim1 = 2 + *nmatd;
    wdep_offset = 1 + wdep_dim1;
    wdep -= wdep_offset;
    dmoist_dim1 = *nmatd;
    dmoist_dim2 = *nsd;
    dmoist_offset = 1 + dmoist_dim1 * (1 + dmoist_dim2 * 14);
    dmoist -= dmoist_offset;

    /* Function Body */
    ret_val = 1.f;
    if (*imoistdep == 2) {
	if (*jreact == 0) {
	    return ret_val;
	}
	reacmin0 = dmoist[*m + (*js + (*jreact + 13) * dmoist_dim2) * 
		dmoist_dim1];
	theta0 = dmoist[*m + (*js + (*jreact + 26) * dmoist_dim2) * 
		dmoist_dim1];
	theta1 = dmoist[*m + (*js + (*jreact + 39) * dmoist_dim2) * 
		dmoist_dim1];
	theta2 = dmoist[*m + (*js + (*jreact + 52) * dmoist_dim2) * 
		dmoist_dim1];
	theta3 = dmoist[*m + (*js + (*jreact + 65) * dmoist_dim2) * 
		dmoist_dim1];
	reacmin1 = dmoist[*m + (*js + (*jreact + 78) * dmoist_dim2) * 
		dmoist_dim1];
	if (dabs(theta2) < .001f) {
	    return ret_val;
	}
	if (*theta <= theta0) {
	    ret_val = reacmin0;
	} else if (*theta <= theta1) {
	    ret_val = reacmin0 + (*theta - theta0) / (theta1 - theta0) * (1.f 
		    - reacmin0);
	} else if (*theta <= theta2) {
	    ret_val = 1.f;
	} else if (*theta <= theta3) {
	    ret_val = reacmin1 + (*theta - theta3) / (theta2 - theta3) * (1.f 
		    - reacmin1);
	} else {
	    ret_val = reacmin1;
	}
    } else if (*imoistdep == 1) {
/* Walker's formula */
	jjj = (*js - 1) * 9;
	if (wdep[*m + 2 + (jjj + *ireact) * wdep_dim1] > *theta && wdep[*m + 
		2 + (jjj + *ireact) * wdep_dim1] > 0.f) {
	    d__1 = (doublereal) (*theta / wdep[*m + 2 + (jjj + *ireact) * 
		    wdep_dim1]);
	    d__2 = (doublereal) wdep[(jjj + *ireact) * wdep_dim1 + 1];
	    ret_val = pow_dd(&d__1, &d__2);
	}
    }
    return ret_val;
} /* rmd_ */

/* *********************************************************************** */
doublereal rmd1_(integer *nmatd, integer *nsd, integer *m, integer *js, 
	integer *jreact, real *dmoist, real *theta)
{
    /* System generated locals */
    integer dmoist_dim1, dmoist_dim2, dmoist_offset;
    real ret_val, r__1;

    /* Local variables */
    static real rtype, theta0, theta1, reacmin;

/*     Function expressing reaction rate dependence on the water content */
/*     Type  =1 or -1: Rate increases or decreases with water content, respectively) */
/*     Theta0 - water content at which reaction rate start increasing or decreasing */
/*     Theta1 - water content at which reaction rate stops increasing or decreasing */
/*     ReacMin	- relative minimum rate of reaction */
    /* Parameter adjustments */
    dmoist_dim1 = *nmatd;
    dmoist_dim2 = *nsd;
    dmoist_offset = 1 + dmoist_dim1 * (1 + dmoist_dim2 * 10);
    dmoist -= dmoist_offset;

    /* Function Body */
    rtype = dmoist[*m + (*js + (*jreact + 9) * dmoist_dim2) * dmoist_dim1];
    theta0 = dmoist[*m + (*js + (*jreact + 18) * dmoist_dim2) * dmoist_dim1];
    theta1 = dmoist[*m + (*js + (*jreact + 27) * dmoist_dim2) * dmoist_dim1];
    reacmin = dmoist[*m + (*js + (*jreact + 36) * dmoist_dim2) * dmoist_dim1];
    ret_val = 1.f;
    if (dabs(rtype) < .1f) {
	return ret_val;
    }
    if (rtype > 0.f) {
/* increasing rate */
	if (*theta >= theta1) {
	    ret_val = 1.f;
	} else if (*theta <= theta0) {
	    ret_val = reacmin;
	} else {
	    if ((r__1 = theta1 - theta0, dabs(r__1)) > 0.f) {
		ret_val = reacmin + (*theta - theta0) / (theta1 - theta0) * (
			1.f - reacmin);
	    }
	}
    } else {
/* decreasing rate */
	if (*theta >= theta1) {
	    ret_val = reacmin;
	} else if (*theta <= theta0) {
	    ret_val = 1.f;
	} else {
	    if ((r__1 = theta1 - theta0, dabs(r__1)) > 0.f) {
		ret_val = reacmin + (*theta - theta1) / (theta0 - theta1) * (
			1.f - reacmin);
	    }
	}
    }
    return ret_val;
} /* rmd1_ */

/* ************************************************************************ */
/*     Distribute mass into different phases */
/* Subroutine */ int massinit_(integer *numnp, integer *nsd, integer *ns, 
	integer *nmat, integer *matnum, real *tdep, real *temp, real *chpar, 
	real *conc, real *theta, real *thetaim, real *thsat, logical *llinear,
	 logical *lbact)
{
    /* System generated locals */
    integer chpar_dim1, chpar_offset, conc_dim1, conc_offset, i__1, i__2;

    /* Local variables */
    static integer i__, m;
    static real r__;
    static integer js;
    static real tr, tt;
    static integer jjj;
    static real par[10], rka1, rka2, rkd1, rkd2, xks1, xks2;
    static integer npar;
    extern doublereal cinit_(real *, real *, integer *);

    /* Parameter adjustments */
    --thetaim;
    --theta;
    --temp;
    --matnum;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    --tdep;
    --llinear;
    --thsat;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;

    /* Function Body */
    tr = 293.15f;
    r__ = 8.314f;
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
	jjj = js - 1 << 4;
	i__2 = *numnp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (conc[js + i__ * conc_dim1] > 0.f) {
		m = matnum[i__];
		tt = (temp[i__] + 273.15f - tr) / r__ / (temp[i__] + 273.15f) 
			/ tr;
		par[0] = chpar[m * chpar_dim1 + 1] * exp(tdep[1] * tt);
/* ro */
		par[1] = chpar[m * chpar_dim1 + 3] * exp(tdep[3] * tt);
/* frac */
		par[2] = chpar[jjj + 7 + m * chpar_dim1] * exp(tdep[jjj + 7] *
			 tt);
/* xKs */
		par[3] = chpar[jjj + 8 + m * chpar_dim1] * exp(tdep[jjj + 8] *
			 tt);
/* xNu */
		par[4] = chpar[jjj + 9 + m * chpar_dim1];
/* *exp(TDep(jjj+ 9)*TT) !fExp */
		par[5] = chpar[jjj + 10 + m * chpar_dim1] * exp(tdep[jjj + 7] 
			* tt);
/* xKH */
		par[6] = theta[i__] + thetaim[i__];
		par[7] = thsat[m] - par[6];
		if (*lbact) {
		    rka1 = chpar[jjj + 19 + m * chpar_dim1] * exp(tdep[jjj + 
			    19] * tt);
		    rkd1 = chpar[jjj + 20 + m * chpar_dim1] * exp(tdep[jjj + 
			    20] * tt);
		    xks1 = theta[i__] * rka1 / par[0] / rkd1;
		    rka2 = chpar[jjj + 16 + m * chpar_dim1] * exp(tdep[jjj + 
			    16] * tt);
		    rkd2 = chpar[jjj + 17 + m * chpar_dim1] * exp(tdep[jjj + 
			    17] * tt);
		    xks2 = theta[i__] * rka1 / par[0] / rkd1;
		    par[2] = xks1 + xks2;
		}
		if (llinear[js]) {
		    conc[js + i__ * conc_dim1] /= par[6] + par[0] * par[2] + 
			    par[7] * par[5];
		} else {
		    conc[js + i__ * conc_dim1] = cinit_(&conc[js + i__ * 
			    conc_dim1], par, &npar);
		}
	    }
/* L11: */
	}
/* L12: */
    }
    return 0;
} /* massinit_ */

/* *********************************************************************** */
/*     Evaluate Liquid concentration from the total solute mass */
doublereal cinit_(real *xmass, real *par, integer *npar)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real x1, x2, xb1, xb2;
    extern /* Subroutine */ int zbrak1_(real *, real *, real *, real *, real *
	    , real *, integer *);
    extern doublereal zbrent1_(real *, real *, real *, real *, integer *);

    /* Parameter adjustments */
    --par;

    /* Function Body */
    x1 = .001f;
    x2 = 1e3f;
    zbrak1_(&x1, &x2, &xb1, &xb2, xmass, &par[1], npar);
    ret_val = zbrent1_(&xb1, &xb2, xmass, &par[1], npar);
    return ret_val;
} /* cinit_ */

/* *********************************************************************** */
doublereal solmass_(real *conc, real *xmass, real *par, integer *npar)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static real ro, xkh, xks, xnu, frac, fexp, theta, ymass, thetaa;

/*     Calculate total solute mass for concentration Conc */
    /* Parameter adjustments */
    --par;

    /* Function Body */
    ro = par[1];
    frac = par[2];
    xks = par[3];
    xnu = par[4];
    fexp = par[5];
    xkh = par[6];
    theta = par[7];
    thetaa = par[8];
    d__1 = (doublereal) (*conc);
    d__2 = (doublereal) fexp;
    d__3 = (doublereal) (*conc);
    d__4 = (doublereal) fexp;
    ymass = theta * *conc + ro * xks * pow_dd(&d__1, &d__2) / (xnu * pow_dd(&
	    d__3, &d__4) + 1.f) + thetaa * xkh * *conc;
    ret_val = *xmass - ymass;
    return ret_val;
} /* solmass_ */

/* *********************************************************************** */
/*     Bracketing of the root, Numerical recepies (345) */
/* Subroutine */ int zbrak1_(real *x1, real *x2, real *xb1, real *xb2, real *
	xmass, real *par, integer *npar)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static real fc;
    static integer nb;
    static real fp, dx2;
    static integer nbb;
    static real dlh;
    extern doublereal solmass_(real *, real *, real *, integer *);

    /* Parameter adjustments */
    --par;

    /* Function Body */
    nbb = 1;
    nb = 1000;
    dlh = (r_lg10(x2) - r_lg10(x1)) / (nb - 1);
    fp = solmass_(x1, xmass, &par[1], npar);
    i__1 = nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dx2 = r_lg10(x1) + i__ * dlh;
	d__1 = (doublereal) dx2;
	*x2 = pow_dd(&c_b173, &d__1);
	fc = solmass_(x2, xmass, &par[1], npar);
	if (fc * fp < 0.f) {
	    ++nbb;
	    *xb1 = *x1;
	    *xb2 = *x2;
	    return 0;
	}
	fp = fc;
	*x1 = *x2;
	if (nbb == nb) {
	    return 0;
	}
/* L11: */
    }
    return 0;
} /* zbrak1_ */

/* *********************************************************************** */
/*     Brent method of finding root that lies between x1 and x2, */
/*     Numerical recepies (354) */
doublereal zbrent1_(real *x1, real *x2, real *xmass, real *par, integer *npar)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;

    /* Local variables */
    static real a, b, c__, d__, e, p, q, r__, s, fa, fb, fc, xm, tol1;
    static integer iter;
    extern doublereal solmass_(real *, real *, real *, integer *);

    /* Parameter adjustments */
    --par;

    /* Function Body */
    a = *x1;
    b = *x2;
    fa = solmass_(&a, xmass, &par[1], npar);
    fb = solmass_(&b, xmass, &par[1], npar);
    if (fb * fa > 0.f) {
	s_paus("Root must be bracketed for ZBRENT1.", (ftnlen)35);
    }
    fc = fb;
    for (iter = 1; iter <= 100; ++iter) {
	if (fb * fc > 0.f) {
	    c__ = a;
	    fc = fa;
	    d__ = b - a;
	    e = d__;
	}
	if (dabs(fc) < dabs(fb)) {
	    a = b;
	    b = c__;
	    c__ = a;
	    fa = fb;
	    fb = fc;
	    fc = fa;
	}
	tol1 = dabs(b) * 5.9999999999999995e-8f + 4.9999999999999998e-7f;
	xm = (c__ - b) * .5f;
	if (dabs(xm) <= tol1 || fb == 0.f) {
	    ret_val = b;
	    return ret_val;
	}
	if (dabs(e) >= tol1 && dabs(fa) > dabs(fb)) {
	    s = fb / fa;
	    if (a == c__) {
		p = xm * 2.f * s;
		q = 1.f - s;
	    } else {
		q = fa / fc;
		r__ = fb / fc;
		p = s * (xm * 2.f * q * (q - r__) - (b - a) * (r__ - 1.f));
		q = (q - 1.f) * (r__ - 1.f) * (s - 1.f);
	    }
	    if (p > 0.f) {
		q = -q;
	    }
	    p = dabs(p);
/* Computing MIN */
	    r__3 = xm * 3.f * q - (r__1 = tol1 * q, dabs(r__1)), r__4 = (r__2 
		    = e * q, dabs(r__2));
	    if (p * 2.f < dmin(r__3,r__4)) {
		e = d__;
		d__ = p / q;
	    } else {
		d__ = xm;
		e = d__;
	    }
	} else {
	    d__ = xm;
	    e = d__;
	}
	a = b;
	fa = fb;
	if (dabs(d__) > tol1) {
	    b += d__;
	} else {
	    b += r_sign(&tol1, &xm);
	}
	fb = solmass_(&b, xmass, &par[1], npar);
/* L11: */
    }
    s_paus("ZBRENT1 exceeding maximum iterations.", (ftnlen)37);
    ret_val = b;
    return ret_val;
} /* zbrent1_ */

