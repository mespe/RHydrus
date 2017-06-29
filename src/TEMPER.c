/* TEMPER.f -- translated by f2c (version 12.02.01).
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

static integer c__9 = 9;
static integer c__1 = 1;

/* Source file TEMPER.FOR ||||||||||||||||||||||||||||||||||||||||||||||| */
/*     Calculation of heat transport */
/* Subroutine */ int temper_(integer *n, integer *nmat, real *x, real *dt, 
	doublereal *t, integer *matnum, real *tempo, real *tempn, real *tpar, 
	real *ampl, doublereal *b, doublereal *d__, doublereal *e, doublereal 
	*f, real *vold, real *vnew, real *thold, real *thnew, real *cap, real 
	*cond, real *sink, real *tperiod, integer *ktopt, real *ttop, integer 
	*kbott, real *tbot, logical *lvapor, real *thvold, real *thvnew, real 
	*vvold, real *vvnew, real *g0, logical *lenbal, real *heatfl, real *
	xconv, real *tconv, real *dtmaxt, integer *icampbell, integer *itemp)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;

    /* Local variables */
    static real dtempmax;
    static integer i__;
    extern /* Subroutine */ int setupheat_(integer *, integer *, logical *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    real *, real *, real *, real *, real *, real *);
    static real pi, cw, cv, xlat;
    static integer level;
    static real dtemp, ttopa;
    extern /* Subroutine */ int bansol_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *), tempadj_(integer *, real *, real *, 
	    real *, real *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, integer *, integer *, real *
	    , real *, real *, real *, real *, real *), tempcap_(integer *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    integer *, real *, real *);
    extern doublereal xlatent_(real *);

    /* Parameter adjustments */
    --g0;
    --vvnew;
    --vvold;
    --thvnew;
    --thvold;
    --sink;
    --cond;
    --cap;
    --thnew;
    --thold;
    --vnew;
    --vold;
    --f;
    --e;
    --d__;
    --b;
    --tempn;
    --tempo;
    --matnum;
    --x;
    tpar -= 11;

    /* Function Body */
    cw = tpar[19];
    pi = 3.141592654f;
    ttopa = *ttop;
    if (*tperiod > 0.f) {
	ttopa = *ttop + *ampl * sin(pi * 2.f * (real) (*t) / *tperiod - pi * 
		7.f / 12.f);
    }
/*     Upper Flux BC */
    if (*ktopt < 0 && ! (*lenbal)) {
	*heatfl = cw * ttopa * (vnew[*n] + vold[*n]) / 2.f;
	if (*lvapor) {
	    cv = 1.8e6f / *xconv / *tconv / *tconv;
	    xlat = xlatent_(&tempn[*n]) / *xconv / *tconv / *tconv;
	    *heatfl = *heatfl + cv * ttopa * (vvnew[*n] + vvold[*n]) / 2.f + 
		    xlat * (vvnew[*n] + vvold[*n]) / 2.f;
	}
    }
    for (level = 1; level <= 2; ++level) {
	if (level == 1) {
	    tempcap_(n, nmat, &matnum[1], &thold[1], &vold[1], &cap[1], &cond[
		    1], &tpar[11], icampbell, xconv, tconv);
	} else {
	    tempcap_(n, nmat, &matnum[1], &thnew[1], &vnew[1], &cap[1], &cond[
		    1], &tpar[11], icampbell, xconv, tconv);
	}
	setupheat_(n, &level, lenbal, kbott, ktopt, tbot, &ttopa, heatfl, dt, 
		&x[1], &cw, &b[1], &d__[1], &e[1], &f[1], &tempo[1], &cap[1], 
		&cond[1], &vnew[1], &vold[1], &sink[1]);
/* L11: */
    }
/*     Adjust matrix for vapor flow effects */
    if (*lvapor) {
	tempadj_(n, &x[1], dt, &tempn[1], &tempo[1], &b[1], &d__[1], &e[1], &
		f[1], &vvold[1], &vvnew[1], &g0[1], ktopt, kbott, &ttopa, 
		tbot, &thvnew[1], &thvold[1], xconv, tconv);
    }
/*     Solve matrix equation */
    bansol_(n, &b[1], &d__[1], &e[1], &f[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tempn[i__] = (real) f[i__];
/* L12: */
    }
/*     Max time step */
    *dtmaxt = 1e30f;
    cv = 1.8e6f / *xconv / *tconv / *tconv;
    if (*lenbal) {
	if (*itemp == 0) {
	    *dtmaxt = 0.f;
/* Computing MAX */
	    r__2 = .1f, r__3 = tempn[*n];
	    *dtmaxt = (r__1 = cap[*n] * (x[*n] - x[*n - 1]) / (dmax(r__2,r__3)
		     * (cw * vnew[*n] + cv * vvnew[*n])), dabs(r__1));
	}
	dtempmax = 2.f;
	dtemp = (r__1 = tempn[*n] - tempo[*n], dabs(r__1));
/* Computing MIN */
	r__1 = *dtmaxt, r__2 = dtempmax / dmax(.1f,dtemp) * *dt;
	*dtmaxt = dmin(r__1,r__2);
    }
    return 0;
} /* temper_ */

/* *********************************************************************** */
/* Subroutine */ int setupheat_(integer *n, integer *level, logical *lenbal, 
	integer *kbott, integer *ktopt, real *tbot, real *ttopa, real *heatfl,
	 real *dt, real *x, real *cw, doublereal *b, doublereal *d__, 
	doublereal *e, doublereal *f, real *tempo, real *cap, real *cond, 
	real *vnew, real *vold, real *sink)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static real dx, dxa, dxb;

    /* Parameter adjustments */
    --sink;
    --vold;
    --vnew;
    --cond;
    --cap;
    --tempo;
    --f;
    --e;
    --d__;
    --b;
    --x;

    /* Function Body */
    dx = x[2] - x[1];
    if (*kbott > 0) {
	d__[1] = 1.f;
	e[1] = 0.f;
	f[1] = *tbot;
    } else if (*kbott < 0) {
	if (*level == 2) {
	    d__[1] = dx / 2.f / *dt * cap[1] + (cond[1] + cond[2]) / dx / 4.f 
		    + *cw * (vnew[1] * 2.f + vnew[2]) / 12.f + dx / 24.f * *
		    cw * (sink[1] * 3.f + sink[2]);
	    e[1] = -(cond[1] + cond[2]) / 4.f / dx + *cw * (vnew[2] * 2.f + 
		    vnew[1]) / 12.f + dx / 24.f * *cw * (sink[1] + sink[2]);
	} else {
	    f[1] = tempo[1] * (dx / 2.f / *dt * cap[1] - (cond[1] + cond[2]) /
		     dx / 4.f - *cw * (vold[1] * 2.f + vold[2]) / 12.f - dx / 
		    24.f * *cw * (sink[1] * 3.f + sink[2])) + tempo[2] * ((
		    cond[1] + cond[2]) / 4.f / dx - *cw * (vold[2] * 2.f + 
		    vold[1]) / 12.f - dx / 24.f * *cw * (sink[1] + sink[2])) 
		    + *tbot * *cw * (vnew[1] + vold[1]) / 2.f;
	}
    } else {
	d__[1] = -1.f;
	e[1] = 1.f;
	f[1] = 0.f;
    }
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dxa = x[i__] - x[i__ - 1];
	dxb = x[i__ + 1] - x[i__];
	dx = (x[i__ + 1] - x[i__ - 1]) / 2.f;
	if (*level == 2) {
	    b[i__] = -(cond[i__] + cond[i__ - 1]) / 4.f / dxa - *cw * (vnew[
		    i__] + vnew[i__ - 1] * 2.f) / 12.f + dxa / 24.f * *cw * (
		    sink[i__ - 1] + sink[i__]);
	    d__[i__] = (cond[i__ - 1] + cond[i__]) / 4.f / dxa + (cond[i__] + 
		    cond[i__ + 1]) / 4.f / dxb + dx / *dt * cap[i__] + *cw * (
		    vnew[i__ + 1] - vnew[i__ - 1]) / 12.f + dxa / 24.f * *cw *
		     (sink[i__ - 1] + sink[i__] * 3.f) + dxb / 24.f * *cw * (
		    sink[i__] * 3.f + sink[i__ + 1]);
	    e[i__] = -(cond[i__] + cond[i__ + 1]) / 4.f / dxb + *cw * (vnew[
		    i__ + 1] * 2.f + vnew[i__]) / 12.f + dxb / 24 * *cw * (
		    sink[i__ + 1] + sink[i__]);
	} else {
	    f[i__] = tempo[i__ - 1] * ((cond[i__] + cond[i__ - 1]) / 4.f / 
		    dxa + *cw * (vold[i__] + vold[i__ - 1] * 2.f) / 12.f - 
		    dxa / 24.f * *cw * (sink[i__ - 1] + sink[i__])) + tempo[
		    i__] * (-(*cw) * (vold[i__ + 1] - vold[i__ - 1]) / 12.f + 
		    dx / *dt * cap[i__] - (cond[i__ + 1] + cond[i__]) / 4.f / 
		    dxb - (cond[i__] + cond[i__ - 1]) / 4.f / dxa - dxa / 
		    24.f * *cw * (sink[i__ - 1] + sink[i__] * 3.f) - dxb / 
		    24.f * *cw * (sink[i__] * 3.f + sink[i__ + 1])) + tempo[
		    i__ + 1] * ((cond[i__ + 1] + cond[i__]) / 4.f / dxb - *cw 
		    * (vold[i__ + 1] * 2.f + vold[i__]) / 12.f - dxb / 24.f * 
		    *cw * (sink[i__ + 1] + sink[i__]));
	}
/* L12: */
    }
    if (*ktopt > 0) {
	b[*n] = 0.f;
	d__[*n] = 1.f;
	f[*n] = *ttopa;
    } else if (*ktopt < 0) {
	dx = x[*n] - x[*n - 1];
	if (*level == 2) {
	    b[*n] = -(cond[*n] + cond[*n - 1]) / 4.f / dx - *cw * (vnew[*n] + 
		    vnew[*n - 1] * 2.f) / 12.f + dx / 24.f * *cw * (sink[*n - 
		    1] + sink[*n]);
	    d__[*n] = dx / 2.f / *dt * cap[*n] + (cond[*n - 1] + cond[*n]) / 
		    4.f / dx - *cw * (vnew[*n] * 2.f + vnew[*n - 1]) / 12.f + 
		    dx / 24.f * *cw * (sink[*n - 1] + sink[*n] * 3.f);
	} else {
	    f[*n] = tempo[*n - 1] * ((cond[*n] + cond[*n - 1]) / 4.f / dx + *
		    cw * (vold[*n] + vold[*n - 1] * 2.f) / 12.f - dx / 24.f * 
		    *cw * (sink[*n - 1] + sink[*n])) + tempo[*n] * (dx / 2.f /
		     *dt * cap[*n] - (cond[*n - 1] + cond[*n]) / 4.f / dx + *
		    cw * (vold[*n] * 2.f + vold[*n - 1]) / 12.f - dx / 24.f * 
		    *cw * (sink[*n - 1] + sink[*n] * 3.f));
	    if (! (*lenbal)) {
		f[*n] -= *heatfl;
		f[*n] -= dmin(0.f,*heatfl);
	    } else {
		f[*n] -= *heatfl;
	    }
	}
    }
    return 0;
} /* setupheat_ */

/* *********************************************************************** */
/* Subroutine */ int tempadj_(integer *n, real *x, real *dt, real *tempn, 
	real *tempo, doublereal *b, doublereal *d__, doublereal *e, 
	doublereal *f, real *vvold, real *vvnew, real *g0, integer *ktopt, 
	integer *kbott, real *ttopa, real *tbot, real *thvnew, real *thvold, 
	real *xconv, real *tconv)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static real cv, dx, dxa, dxb, lat;
    static integer ilevel;
    static real vvgrad, thvgrad;
    extern doublereal xlatent_(real *);

/*     Cv  - volumetric specific heat of vapor [J/m3/K,kg/m/s2/K] */
/*     Lat - volumetric latent heat of vaporization of water [J/m3,kg/m/s2] */
    /* Parameter adjustments */
    --thvold;
    --thvnew;
    --g0;
    --vvnew;
    --vvold;
    --f;
    --e;
    --d__;
    --b;
    --tempo;
    --tempn;
    --x;

    /* Function Body */
    cv = 1.8e6f / *xconv / *tconv / *tconv;
    for (ilevel = 1; ilevel <= 2; ++ilevel) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (ilevel == 1) {
		if (i__ == 1) {
		    vvgrad = (vvold[i__ + 1] - vvold[i__]) / (x[i__ + 1] - x[
			    i__]);
		} else if (i__ == *n) {
		    vvgrad = (vvold[i__] - vvold[i__ - 1]) / (x[i__] - x[i__ 
			    - 1]);
		} else {
		    vvgrad = (vvold[i__ + 1] - vvold[i__ - 1]) / (x[i__ + 1] 
			    - x[i__ - 1]) * 2.f;
		}
		lat = xlatent_(&tempo[i__]) / *xconv / *tconv / *tconv;
	    } else {
		if (i__ == 1) {
		    vvgrad = (vvnew[i__ + 1] - vvnew[i__]) / (x[i__ + 1] - x[
			    i__]);
		} else if (i__ == *n) {
		    vvgrad = (vvnew[i__] - vvnew[i__ - 1]) / (x[i__] - x[i__ 
			    - 1]);
		} else {
		    vvgrad = (vvnew[i__ + 1] - vvnew[i__ - 1]) / (x[i__ + 1] 
			    - x[i__ - 1]) * 2.f;
		}
		lat = xlatent_(&tempn[i__]) / *xconv / *tconv / *tconv;
	    }
	    thvgrad = (thvnew[i__] - thvold[i__]) / *dt;
	    g0[i__] = -lat * (vvgrad + thvgrad);
/* L11: */
	}
	dx = x[2] - x[1];
	if (*kbott < 0) {
	    if (ilevel == 1) {
		f[1] = f[1] + tempo[1] * (-cv * (vvold[1] * 2.f + vvold[2]) / 
			12.f) + tempo[2] * (-cv * (vvold[2] * 2.f + vvold[1]) 
			/ 12.f) + dx / 12.f * (g0[1] * 2.f + g0[2]) + *tbot * 
			cv * (vvnew[1] + vvold[1]) / 2.f;
	    } else {
		d__[1] += cv * (vvnew[1] * 2.f + vvnew[2]) / 12.f;
		e[1] += cv * (vvnew[2] * 2.f + vvnew[1]) / 12.f;
		f[1] += dx / 12.f * (g0[1] * 2.f + g0[2]);
	    }
	}
	i__1 = *n - 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    dxa = x[i__] - x[i__ - 1];
	    dxb = x[i__ + 1] - x[i__];
	    dx = (x[i__ + 1] - x[i__ - 1]) / 2.f;
	    if (ilevel == 1) {
		f[i__] = f[i__] + tempo[i__ - 1] * (cv * (vvold[i__] + vvold[
			i__ - 1] * 2.f) / 12.f) + tempo[i__] * (-cv * (vvold[
			i__ + 1] - vvold[i__ - 1]) / 12.f) + tempo[i__ + 1] * 
			(-cv * (vvold[i__ + 1] * 2.f + vvold[i__]) / 12.f) + 
			dxa * (g0[i__ - 1] + g0[i__] * 2.f) / 12.f + dxb * (
			g0[i__] * 2.f + g0[i__ + 1]) / 12.f;
	    } else {
		b[i__] -= cv * (vvnew[i__] + vvnew[i__ - 1] * 2.f) / 12.f;
		d__[i__] += cv * (vvnew[i__ + 1] - vvnew[i__ - 1]) / 12.f;
		e[i__] += cv * (vvnew[i__ + 1] * 2.f + vvnew[i__]) / 12.f;
		f[i__] = f[i__] + dxa * (g0[i__ - 1] + g0[i__] * 2.f) / 12.f 
			+ dxb * (g0[i__] * 2.f + g0[i__ + 1]) / 12.f;
	    }
/* L12: */
	}
	if (*ktopt < 0) {
	    dx = x[*n] - x[*n - 1];
	    if (ilevel == 1) {
		f[*n] = f[*n] + tempo[*n - 1] * (cv * (vvold[*n] + vvold[*n - 
			1] * 2.f) / 12.f) + tempo[*n] * (cv * (vvold[*n] * 
			2.f + vvold[*n - 1]) / 12.f) + dx / 12.f * (g0[*n - 1]
			 + g0[*n] * 2.f);
		if (*ktopt == -1) {
		    f[*n] -= *ttopa * cv * (vvnew[*n] + vvold[*n]) / 2.f;
		}
	    } else {
		b[*n] -= cv * (vvnew[*n] + vvnew[*n - 1] * 2.f) / 12.f;
		d__[*n] -= cv * (vvnew[*n] * 2.f + vvnew[*n - 1]) / 12.f;
		f[*n] += dx / 12.f * (g0[*n - 1] + g0[*n] * 2.f);
	    }
	}
/* L13: */
    }
    return 0;
} /* tempadj_ */

/* *********************************************************************** */
/* Subroutine */ int tempcap_(integer *n, integer *nmat, integer *matnum, 
	real *theta, real *veloc, real *cap, real *cond, real *tpar, integer *
	icampbell, real *xconv, real *tconv)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, m;
    static real v, aa, bb, cc, dd, ee, xc, th, xlamb;

    /* Fortran I/O blocks */
    static cilist io___27 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 5, 0, 0, 0 };


    /* Parameter adjustments */
    --cond;
    --cap;
    --veloc;
    --theta;
    --matnum;
    tpar -= 11;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = matnum[i__];
	th = theta[i__];
	v = veloc[i__];
	cap[i__] = tpar[m * 10 + 7] * tpar[m * 10 + 1] + tpar[m * 10 + 8] * 
		tpar[m * 10 + 2] + tpar[m * 10 + 9] * th;
	if (cap[i__] == 0.f) {
	    s_wsle(&io___27);
	    do_lio(&c__9, &c__1, "Heat capacity is equal to zero", (ftnlen)30)
		    ;
	    e_wsle();
	    s_rsle(&io___28);
	    e_rsle();
	    s_stop("", (ftnlen)0);
	}
/* Computing MAX */
	r__1 = 0.f, r__2 = tpar[m * 10 + 4] + tpar[m * 10 + 5] * th + tpar[m *
		 10 + 6] * sqrt(th);
	cond[i__] = dmax(r__1,r__2);
	if (*icampbell == 1) {
/*         TPar(1,M) - the volume fraction of solids */
/*         TPar(4,M) - the volume fraction of quartz */
/*         TPar(5,M) - the volume fraction of other minerals */
/*         TPar(6,M) - the volume fraction of clay */
	    aa = (tpar[m * 10 + 4] * 1.73f + .57f + tpar[m * 10 + 5] * .93f) /
		     (1.f - tpar[m * 10 + 4] * .74f - tpar[m * 10 + 5] * .49f)
		     - tpar[m * 10 + 1] * 2.8f * (1.f - tpar[m * 10 + 1]);
	    bb = tpar[m * 10 + 1] * 2.8f;
/* Computing MAX */
	    r__1 = .005f, r__2 = tpar[m * 10 + 6];
	    xc = dmax(r__1,r__2);
	    cc = 2.6f / sqrt(xc) + 1.f;
/* Computing 2nd power */
	    r__1 = tpar[m * 10 + 1];
	    dd = r__1 * r__1 * .7f + .03f;
	    ee = 4.f;
	    d__1 = (doublereal) (cc * th);
	    d__2 = (doublereal) ee;
	    xlamb = aa + bb * th - (aa - dd) * exp(-pow_dd(&d__1, &d__2));
/* Computing MAX */
	    r__1 = 0.f, r__2 = xlamb * *xconv / *tconv / *tconv / *tconv;
	    cond[i__] = dmax(r__1,r__2);
	}
	cond[i__] += tpar[m * 10 + 9] * tpar[m * 10 + 3] * dabs(v);
/* L11: */
    }
    return 0;
} /* tempcap_ */

/* *********************************************************************** */
doublereal xlatent_(real *temp)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real row, xlw;

/*     Function calculating the volumetric latent heat of vaporization of */
/*     water [J/m3,kg/m/s2] */
/*     row   - density of soil water [kg/m3] */
/*     xLw   - latent heat of vaporization of water [J/kg,ML2/T2] */
/* Computing 2nd power */
    r__1 = *temp - 4.f;
/* Computing 3rd power */
    r__2 = *temp - 4.f;
    row = (1.f - r__1 * r__1 * 7.37e-6f + r__2 * (r__2 * r__2) * 3.79e-8f) * 
	    1e3f;
    xlw = 2.501e6f - *temp * 2369.2f;
    ret_val = row * xlw;
    return ret_val;
} /* xlatent_ */

