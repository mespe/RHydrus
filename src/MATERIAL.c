/* MATERIAL.f -- translated by f2c (version 12.02.01).
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

static doublereal c_b2 = 1e300;
static integer c__1010 = 1010;
static integer c__3 = 3;
static integer c__1 = 1;
static doublereal c_b9 = 10.;
static doublereal c_b13 = 2.;
static doublereal c_b15 = 6.2831853080000002;
static doublereal c_b16 = .5;
static integer c__10 = 10;
static integer c__5 = 5;
static real c_b53 = 0.f;
static integer c__2 = 2;

/* Source file MATERIAL.FOR ||||||||||||||||||||||||||||||||||||||||||||| */
/*     iModel = 0: van Genuchten */
/*              1: modified van Genuchten (Vogel and Cislerova) */
/*              2: Brooks and Corey */
/*              3: van Genuchte with air entry value of 2 cm */
/*              4: log-normal (Kosugi) */
/*              5: dual-porosity (Durner) */
/*             10: fractal model (Shlomo Orr) */
/*     lAltern = VG model with alfa and n different for retention curve */
/*               and hydraulic conductivity function. (iModel=1) */
/*               One needs to comment out line 250 in input2.for */
/*     lUnBound = Unbound n and m in VG-M function. */
/*                one needs to uncheck m=Par(6) in other routines as well. */
doublereal fk_(integer *imodel, real *h__, real *par)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static logical lunbound;
    static doublereal d__, m, n, t, m2, n2, w1, w2, aa, bb, ha, hh, qa, hk, 
	    kk, qe, ax, mn, qk, kr, ks, qm, hs, mt, ex, qr, qs, dw, sk1, sk2, 
	    sv1, sw1, sw2, sv2, qee, ffq, qek, wcl, dlg1, dlg2, alfa, dlgc, 
	    beta;
    extern doublereal binc_(doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal bpar, qeek, ffqk, hmin, qees, dlgw;
    static integer ppar;
    static doublereal term, alfa2;
    extern doublereal gamma_(doublereal *);
    static doublereal dlgks;
    extern doublereal qnorm_(doublereal *);
    static doublereal lambda, rdenom, rnumer;
    static logical laltern;

    /* Fortran I/O blocks */
    static cilist io___31 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    --par;

    /* Function Body */
    qr = par[1];
    qs = par[2];
    alfa = par[3];
    n = par[4];
    ks = dmax(par[5],1e-37f);
    bpar = par[6];
    if (*imodel == 0 || *imodel == 1 || *imodel == 3) {
/*        BPar=.5d0 */
/* VG and mod */
	ppar = 2;
	if (*imodel == 0 || *imodel == 3) {
	    qm = qs;
	    qa = qr;
	    qk = qs;
	    kk = ks;
	} else if (*imodel == 1) {
	    laltern = FALSE_;
	    if (! laltern) {
		qm = par[7];
		qa = par[8];
		qk = par[9];
		kk = par[10];
	    } else {
		qm = qs;
		qa = qr;
		qk = qs;
		kk = ks;
		alfa = par[9];
		n = par[10];
	    }
	}
	if (*imodel == 3) {
	    qm = par[7];
	}
	m = 1. - 1. / n;
	lunbound = FALSE_;
	if (lunbound) {
	    m = par[6];
	    bpar = .5;
	}
	d__1 = 1. / n;
	hmin = -pow_dd(&c_b2, &d__1) / max(alfa,1.);
/* Computing MAX */
	d__1 = (doublereal) (*h__);
	hh = max(d__1,hmin);
/* Computing MIN */
	d__1 = (qs - qa) / (qm - qa);
	qees = min(d__1,.999999999999999);
/* Computing MIN */
	d__1 = (qk - qa) / (qm - qa);
	qeek = min(d__1,qees);
	d__2 = -1. / m;
	d__1 = pow_dd(&qees, &d__2) - 1.;
	d__3 = 1. / n;
	hs = -1. / alfa * pow_dd(&d__1, &d__3);
	d__2 = -1. / m;
	d__1 = pow_dd(&qeek, &d__2) - 1.;
	d__3 = 1. / n;
	hk = -1. / alfa * pow_dd(&d__1, &d__3);
	if ((doublereal) (*h__) < hk) {
	    if (! lunbound) {
/* m=1-1/n */
		d__2 = -alfa * hh;
		d__1 = pow_dd(&d__2, &n) + 1.;
		d__3 = -m;
		qee = pow_dd(&d__1, &d__3);
		qe = (qm - qa) / (qs - qa) * qee;
		qek = (qm - qa) / (qs - qa) * qeek;
		d__2 = 1. / m;
		d__1 = 1. - pow_dd(&qee, &d__2);
		ffq = 1. - pow_dd(&d__1, &m);
		d__2 = 1. / m;
		d__1 = 1. - pow_dd(&qeek, &d__2);
		ffqk = 1. - pow_dd(&d__1, &m);
		if (ffq <= 0.) {
		    d__1 = 1. / m;
		    ffq = m * pow_dd(&qee, &d__1);
		}
		d__1 = qe / qek;
		d__2 = ffq / ffqk;
		kr = pow_dd(&d__1, &bpar) * pow_di(&d__2, &ppar) * kk / ks;
		if (*imodel == 0) {
		    kr = pow_dd(&qe, &bpar) * pow_di(&ffq, &ppar);
		}
/*           Gardner's model {K=Ks*[exp(-a*h)]} */
/*            Kr=dexp(-BPar*dble(h)) */
/* Computing MAX */
		d__1 = ks * kr;
		ret_val = (real) max(d__1,1e-37);
	    } else {
/* unbounded n and m */
		mn = m * n;
		mt = 1.;
/*           Calculate complete Beta function */
		aa = m + mt / n;
		bb = 1. - mt / n;
		if (bb <= .004) {
		    s_wsle(&io___31);
		    do_lio(&c__3, &c__1, (char *)&c__1010, (ftnlen)sizeof(
			    integer));
		    e_wsle();
/* L1010: */
		    s_stop("", (ftnlen)0);
		}
		d__1 = m + 1.f;
		beta = gamma_(&aa) * gamma_(&bb) / gamma_(&d__1);
/* Computing MAX */
		d__1 = 2. / (m + 2.);
		wcl = max(d__1,.2);
		dlgks = d_lg10(&ks);
/*           Water content */
		ax = alfa * (-(*h__));
		if (ax < 1e-20) {
		    qe = 1.;
		} else {
		    ex = n * d_lg10(&ax);
		    if (ex < -10.) {
			qe = 1.;
		    } else if (ex < 10.) {
			d__1 = pow_dd(&ax, &n) + 1.f;
			d__2 = -m;
			qe = pow_dd(&d__1, &d__2);
		    } else {
			ex = m * ex;
			if (ex < 30.) {
			    d__1 = -m * n;
			    qe = pow_dd(&ax, &d__1);
			} else {
			    qe = 0.;
			}
		    }
		}
/*           Conductivity */
		if (qe <= 1e-10) {
		    ret_val = 1e-37f;
		} else if (qe > .999999) {
		    ret_val = ks;
		} else {
		    dlgw = d_lg10(&qe);
		    dlg2 = 3. - mt + bpar + 2. / mn;
		    dlgc = dlg2 * dlgw + dlgks;
		    if (dlgc > -37. && dlgw > m * -15.) {
			d__1 = 1. / m;
			dw = pow_dd(&qe, &d__1);
			if (dw < 1e-6) {
			    d__1 = n / (beta * (mn + mt));
			    dlg1 = (3. - mt) * d_lg10(&d__1);
			    dlgc += dlg1;
			    ret_val = pow_dd(&c_b9, &dlgc);
			    return ret_val;
			}
			if (qe - wcl <= 0.) {
			    term = binc_(&dw, &aa, &bb, &beta);
			} else {
			    d__1 = 1. - dw;
			    term = 1. - binc_(&d__1, &bb, &aa, &beta);
			}
			kr = pow_dd(&qe, &bpar) * term;
			if (mt < 1.5) {
			    kr *= term;
			}
			dlgc = d_lg10(&kr) + dlgks;
		    }
		    dlgc = max(-37.,dlgc);
		    ret_val = pow_dd(&c_b9, &dlgc);
		}
	    }
	}
	if ((doublereal) (*h__) >= hk && (doublereal) (*h__) < hs) {
	    kr = (1. - kk / ks) / (hs - hk) * ((doublereal) (*h__) - hs) + 1.;
	    ret_val = (real) (ks * kr);
	}
	if ((doublereal) (*h__) >= hs) {
	    ret_val = (real) ks;
	}
    } else if (*imodel == 2) {
/*        BPar=1.d0 */
/* Brooks and Cores */
	lambda = 2.;
/*  !=2 for Mualem Model, =1.5 for Burdine model */
	hs = -1. / alfa;
	if (*h__ < hs) {
	    d__1 = -alfa * *h__;
	    d__2 = n * (bpar + lambda) + 2.;
	    kr = 1. / pow_dd(&d__1, &d__2);
/* Computing MAX */
	    d__1 = ks * kr;
	    ret_val = (real) max(d__1,1e-37);
	} else {
	    ret_val = (real) ks;
	}
    } else if (*imodel == 4) {
/* Log-normal model */
	hs = 0.;
	if (*h__ < hs) {
	    d__1 = log(-(*h__) / alfa) / n;
	    qee = qnorm_(&d__1);
	    d__1 = log(-(*h__) / alfa) / n + n;
	    t = qnorm_(&d__1);
	    kr = pow_dd(&qee, &bpar) * t * t;
/* Computing MAX */
	    d__1 = ks * kr;
	    ret_val = (real) max(d__1,1e-37);
	} else {
	    ret_val = (real) ks;
	}
    } else if (*imodel == 5) {
/* Dual-porosity model */
	w2 = par[7];
	alfa2 = par[8];
	n2 = par[9];
	m = 1. - 1. / n;
	m2 = 1. - 1. / n2;
	w1 = 1. - w2;
	d__2 = -alfa * *h__;
	d__1 = pow_dd(&d__2, &n) + 1.;
	d__3 = -m;
	sw1 = w1 * pow_dd(&d__1, &d__3);
	d__2 = -alfa2 * *h__;
	d__1 = pow_dd(&d__2, &n2) + 1.;
	d__3 = -m2;
	sw2 = w2 * pow_dd(&d__1, &d__3);
	qe = sw1 + sw2;
	d__1 = -alfa * *h__;
	d__2 = n - 1;
	sv1 = pow_dd(&d__1, &d__2);
	d__1 = -alfa2 * *h__;
	d__2 = n2 - 1;
	sv2 = pow_dd(&d__1, &d__2);
	d__2 = -alfa * *h__;
	d__1 = pow_dd(&d__2, &n) + 1.;
	d__3 = -m;
	sk1 = w1 * alfa * (1. - sv1 * pow_dd(&d__1, &d__3));
	d__2 = -alfa2 * *h__;
	d__1 = pow_dd(&d__2, &n2) + 1.;
	d__3 = -m2;
	sk2 = w2 * alfa2 * (1. - sv2 * pow_dd(&d__1, &d__3));
	rnumer = sk1 + sk2;
	rdenom = w1 * alfa + w2 * alfa2;
	if (rdenom != 0.f) {
/* Computing 2nd power */
	    d__1 = rnumer / rdenom;
	    kr = pow_dd(&qe, &bpar) * (d__1 * d__1);
	}
/* Computing MAX */
	d__1 = ks * kr;
	ret_val = (real) max(d__1,1e-37);
    } else if (*imodel == -1) {
/* Shlomo Orr */
	ha = alfa;
	d__ = n;
	kr = 1.f;
	if (-(*h__) > ha) {
	    d__2 = -ha / *h__;
	    d__3 = 3.f - d__;
	    d__1 = 1.f - (1.f - pow_dd(&d__2, &d__3)) / (1.f - qr);
	    d__4 = d__ / (3.f - d__);
	    kr = pow_dd(&d__1, &d__4);
	}
/* Computing MAX */
	d__1 = ks * kr;
	ret_val = (real) max(d__1,1e-37);
    }
    return ret_val;
} /* fk_ */

/* *********************************************************************** */
doublereal fc_(integer *imodel, real *h__, real *par)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal d__, m, n, t, c1, c2, m2, n2, w1, w2, ha, hh, qa, hs, 
	    qm, qr, qs, c1a, c1b, c2a, c2b, alfa, hmin, qees, alfa2;

    /* Parameter adjustments */
    --par;

    /* Function Body */
    qr = par[1];
    qs = par[2];
    alfa = par[3];
    n = par[4];
    if (*imodel == 0 || *imodel == 1 || *imodel == 3) {
	if (*imodel == 0 || *imodel == 3) {
	    qm = qs;
	    qa = qr;
	} else if (*imodel == 1) {
	    qm = par[7];
	    qa = par[8];
	}
	if (*imodel == 3) {
	    qm = par[7];
	}
	m = 1. - 1. / n;
/*        m=Par(6) */
	d__1 = 1. / n;
	hmin = -pow_dd(&c_b2, &d__1) / max(alfa,1.);
/* Computing MAX */
	d__1 = (doublereal) (*h__);
	hh = max(d__1,hmin);
/* Computing MIN */
	d__1 = (qs - qa) / (qm - qa);
	qees = min(d__1,.999999999999999);
	d__2 = -1. / m;
	d__1 = pow_dd(&qees, &d__2) - 1.;
	d__3 = 1. / n;
	hs = -1. / alfa * pow_dd(&d__1, &d__3);
	if ((doublereal) (*h__) < hs) {
	    d__2 = -alfa * hh;
	    d__1 = pow_dd(&d__2, &n) + 1.;
	    d__3 = -m - 1.;
	    c1 = pow_dd(&d__1, &d__3);
	    d__1 = -hh;
	    d__2 = n - 1.;
	    c2 = (qm - qa) * m * n * pow_dd(&alfa, &n) * pow_dd(&d__1, &d__2) 
		    * c1;
	    ret_val = (real) max(c2,1e-37);
	    return ret_val;
	} else {
	    ret_val = 0.f;
	}
    } else if (*imodel == 2) {
	hs = -1. / alfa;
	if (*h__ < hs) {
	    d__1 = -n;
	    d__2 = (doublereal) (-(*h__));
	    d__3 = -n - 1.;
	    c2 = (qs - qr) * n * pow_dd(&alfa, &d__1) * pow_dd(&d__2, &d__3);
	    ret_val = (real) max(c2,1e-37);
	} else {
	    ret_val = 0.f;
	}
    } else if (*imodel == 4) {
	hs = 0.;
	if (*h__ < hs) {
	    d__1 = log(-(*h__) / alfa);
	    t = exp(pow_dd(&d__1, &c_b13) * -1. / (pow_dd(&n, &c_b13) * 2.));
	    c2 = (qs - qr) / pow_dd(&c_b15, &c_b16) / n / (-(*h__)) * t;
	    ret_val = (real) max(c2,1e-37);
	} else {
	    ret_val = 0.f;
	}
    } else if (*imodel == 5) {
	w2 = par[7];
	alfa2 = par[8];
	n2 = par[9];
	m = 1. - 1. / n;
	m2 = 1. - 1. / n2;
	w1 = 1. - w2;
	d__2 = -alfa * *h__;
	d__1 = pow_dd(&d__2, &n) + 1.;
	d__3 = -m - 1.;
	c1a = pow_dd(&d__1, &d__3);
	d__2 = -alfa2 * *h__;
	d__1 = pow_dd(&d__2, &n2) + 1.;
	d__3 = -m2 - 1.;
	c1b = pow_dd(&d__1, &d__3);
	d__1 = (doublereal) (-(*h__));
	d__2 = n - 1.;
	c2a = (qs - qr) * m * n * pow_dd(&alfa, &n) * pow_dd(&d__1, &d__2) * 
		c1a * w1;
	d__1 = (doublereal) (-(*h__));
	d__2 = n2 - 1.;
	c2b = (qs - qr) * m2 * n2 * pow_dd(&alfa2, &n2) * pow_dd(&d__1, &d__2)
		 * c1b * w2;
	ret_val = c2a + c2b;
    } else if (*imodel == -1) {
/* Shlomo Orr */
	ha = alfa;
	d__ = n;
	if (-(*h__) < ha) {
	    c1 = 0.f;
	} else {
	    d__1 = 3.f - d__;
	    d__2 = (doublereal) (-(*h__));
	    d__3 = d__ - 4.f;
	    c1 = pow_dd(&ha, &d__1) * -1.f * (d__ - 3.f) * pow_dd(&d__2, &
		    d__3);
	}
	ret_val = (real) max(c1,1e-37);
    }
    return ret_val;
} /* fc_ */

/* *********************************************************************** */
doublereal fq_(integer *imodel, real *h__, real *par)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static doublereal d__, m, n, m2, n2, w1, w2, ha, hh, qa, qe, hs, qm, qr, 
	    qs, sw1, sw2, qee, alfa, hmin, qees, alfa2;
    extern doublereal qnorm_(doublereal *);

    /* Parameter adjustments */
    --par;

    /* Function Body */
    qr = par[1];
    qs = par[2];
    alfa = par[3];
    n = par[4];
    if (*imodel == 0 || *imodel == 1 || *imodel == 3) {
	if (*imodel == 0 || *imodel == 3) {
	    qm = qs;
	    qa = qr;
	} else if (*imodel == 1) {
	    qm = par[7];
	    qa = par[8];
	}
	if (*imodel == 3) {
	    qm = par[7];
	}
	m = 1. - 1. / n;
/*        m=Par(6) */
	d__1 = 1. / n;
	hmin = -pow_dd(&c_b2, &d__1) / max(alfa,1.);
/* Computing MAX */
	d__1 = (doublereal) (*h__);
	hh = max(d__1,hmin);
/* Computing MIN */
	d__1 = (qs - qa) / (qm - qa);
	qees = min(d__1,.999999999999999);
	d__2 = -1. / m;
	d__1 = pow_dd(&qees, &d__2) - 1.;
	d__3 = 1. / n;
	hs = -1. / alfa * pow_dd(&d__1, &d__3);
	if ((doublereal) (*h__) < hs) {
	    d__2 = -alfa * hh;
	    d__1 = pow_dd(&d__2, &n) + 1.;
	    d__3 = -m;
	    qee = pow_dd(&d__1, &d__3);
/* Computing MAX */
	    d__1 = qa + (qm - qa) * qee;
	    ret_val = (real) max(d__1,1e-37);
	    return ret_val;
	} else {
	    ret_val = (real) qs;
	}
    } else if (*imodel == 2) {
	hs = -1. / alfa;
	if (*h__ < hs) {
	    d__1 = -alfa * *h__;
	    d__2 = -n;
	    qee = pow_dd(&d__1, &d__2);
/* Computing MAX */
	    d__1 = qr + (qs - qr) * qee;
	    ret_val = (real) max(d__1,1e-37);
	} else {
	    ret_val = (real) qs;
	}
    } else if (*imodel == 4) {
	hs = 0.;
	if (*h__ < hs) {
	    d__1 = log(-(*h__) / alfa) / n;
	    qee = qnorm_(&d__1);
/* Computing MAX */
	    d__1 = qr + (qs - qr) * qee;
	    ret_val = (real) max(d__1,1e-37);
	} else {
	    ret_val = (real) qs;
	}
    } else if (*imodel == 5) {
	w2 = par[7];
	alfa2 = par[8];
	n2 = par[9];
	m = 1. - 1. / n;
	m2 = 1. - 1. / n2;
	w1 = 1. - w2;
	d__2 = -alfa * *h__;
	d__1 = pow_dd(&d__2, &n) + 1.;
	d__3 = -m;
	sw1 = w1 * pow_dd(&d__1, &d__3);
	d__2 = -alfa2 * *h__;
	d__1 = pow_dd(&d__2, &n2) + 1.;
	d__3 = -m2;
	sw2 = w2 * pow_dd(&d__1, &d__3);
	qe = sw1 + sw2;
/* Computing MAX */
	d__1 = qr + (qs - qr) * qe;
	ret_val = (real) max(d__1,1e-37);
    } else if (*imodel == -1) {
/* Shlomo Orr */
	ha = alfa;
	d__ = n;
	ret_val = qs;
	if (ha > 0.f) {
/* Computing MAX */
/* Computing MIN */
	    d__4 = -(*h__) / ha;
	    d__5 = d__ - 3.f;
	    d__2 = qs, d__3 = qs + pow_dd(&d__4, &d__5) - 1.f;
	    d__1 = min(d__2,d__3);
	    ret_val = max(d__1,qr);
	}
    }
    return ret_val;
} /* fq_ */

/* *********************************************************************** */
doublereal fh_(integer *imodel, real *qe, real *par)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal d__, h__, m, n, p, x, y, m2, n2, w1, w2, ha, qa, th, qm,
	     qr, qs, qee, alfa, qeem, hmin, alfa2;
    extern doublereal xmualem_(doublereal *, real *, integer *);

    /* Parameter adjustments */
    --par;

    /* Function Body */
    qr = par[1];
    qs = par[2];
    alfa = par[3];
    n = par[4];
    if (*imodel == 0 || *imodel == 1 || *imodel == 3) {
	if (*imodel == 0 || *imodel == 3) {
	    qm = qs;
	    qa = qr;
	} else if (*imodel == 1) {
	    qm = par[7];
	    qa = par[8];
	}
	if (*imodel == 3) {
	    qm = par[7];
	}
	m = 1. - 1. / n;
/*        m=Par(6) */
	d__1 = 1. / n;
	hmin = -pow_dd(&c_b2, &d__1) / max(alfa,1.);
	d__2 = -alfa * hmin;
	d__1 = pow_dd(&d__2, &n) + 1.;
	d__3 = -m;
	qeem = pow_dd(&d__1, &d__3);
/* Computing MIN */
/* Computing MAX */
	d__2 = *qe * (qs - qa) / (qm - qa);
	d__1 = max(d__2,qeem);
	qee = min(d__1,.999999999999999);
/* Computing MAX */
	d__3 = -1. / m;
	d__2 = pow_dd(&qee, &d__3) - 1.;
	d__4 = 1. / n;
	d__1 = -1. / alfa * pow_dd(&d__2, &d__4);
	ret_val = (real) max(d__1,-1e37);
    } else if (*imodel == 2) {
/* Computing MAX */
	d__2 = (doublereal) dmax(*qe,1e-10f);
	d__3 = -1. / n;
	d__1 = -1. / alfa * pow_dd(&d__2, &d__3);
	ret_val = (real) max(d__1,-1e37);
    } else if (*imodel == 4) {
	if (*qe > .9999f) {
	    ret_val = 0.f;
	} else if (*qe < 1e-5f) {
	    ret_val = -1e8f;
	} else {
	    y = *qe * 2.;
	    if (y < 1.f) {
		p = sqrt(-log(y / 2.));
	    }
	    if (y >= 1.f) {
		p = sqrt(-log(1 - y / 2.));
	    }
/* Computing 3rd power */
	    d__1 = p;
/* Computing 2nd power */
	    d__2 = p;
/* Computing 3rd power */
	    d__3 = p;
/* Computing 4th power */
	    d__4 = p, d__4 *= d__4;
	    x = p - (p * .9425908f + 1.881796f + d__1 * (d__1 * d__1) * 
		    .0546028f) / (p * 2.356868f + 1.f + d__2 * d__2 * 
		    .3087091f + d__3 * (d__3 * d__3) * .0937563f + d__4 * 
		    d__4 * .021914f);
	    if (y >= 1.f) {
		x = -x;
	    }
	    ret_val = (real) (-alfa * exp(sqrt(2.f) * n * x));
	}
    } else if (*imodel == 5) {
	w2 = par[7];
	alfa2 = par[8];
	n2 = par[9];
	m = 1. - 1. / n;
	m2 = 1. - 1. / n2;
	w1 = 1. - w2;
	qee = *qe;
	if (qee > .9999) {
	    ret_val = 0.f;
	} else if (qee < 1e-5) {
	    ret_val = -1e8f;
	} else {
	    h__ = xmualem_(&qee, &par[1], &c__10);
	    ret_val = (real) max(h__,-1e37);
	}
    } else if (*imodel == -1) {
/* Shlomo Orr */
	ha = alfa;
	d__ = n;
	th = qr + (qs - qr) * *qe;
	h__ = 0.f;
	if (d__ != 3.f) {
	    d__1 = 1.f - qs + th;
	    d__2 = 1.f / (d__ - 3.f);
	    h__ = -ha * pow_dd(&d__1, &d__2);
	}
	ret_val = (real) max(h__,-1e37);
    }
    return ret_val;
} /* fh_ */

/* *********************************************************************** */
doublereal fs_(integer *imodel, real *h__, real *par)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static doublereal d__, m, n, m2, n2, w1, w2, ha, hh, qa, qe, hs, qm, qr, 
	    qs, sw1, sw2, qee, alfa, hmin, qees, alfa2;
    extern doublereal qnorm_(doublereal *);

    /* Parameter adjustments */
    --par;

    /* Function Body */
    qr = par[1];
    qs = par[2];
    alfa = par[3];
    n = par[4];
    if (*imodel == 0 || *imodel == 1 || *imodel == 3) {
	if (*imodel == 0 || *imodel == 3) {
	    qm = qs;
	    qa = qr;
	} else if (*imodel == 1) {
	    qm = par[7];
	    qa = par[8];
	}
	if (*imodel == 3) {
	    qm = par[7];
	}
	m = 1. - 1. / n;
/*        m=Par(6) */
/* Computing MIN */
	d__1 = (qs - qa) / (qm - qa);
	qees = min(d__1,.999999999999999);
	d__2 = -1. / m;
	d__1 = pow_dd(&qees, &d__2) - 1.;
	d__3 = 1. / n;
	hs = -1. / alfa * pow_dd(&d__1, &d__3);
	if (*h__ < hs) {
	    d__1 = 1.f / n;
	    hmin = -pow_dd(&c_b2, &d__1) / max(alfa,1.);
/* Computing MAX */
	    d__1 = (doublereal) (*h__);
	    hh = max(d__1,hmin);
	    d__2 = -alfa * hh;
	    d__1 = pow_dd(&d__2, &n) + 1.;
	    d__3 = -m;
	    qee = pow_dd(&d__1, &d__3);
	    qe = qee * (qm - qa) / (qs - qa);
	    ret_val = (real) max(qe,1e-37);
	} else {
	    ret_val = 1.f;
	}
    } else if (*imodel == 2) {
	hs = -1. / alfa;
	if (*h__ < hs) {
	    d__1 = -alfa * *h__;
	    d__2 = -n;
	    qe = pow_dd(&d__1, &d__2);
	    ret_val = (real) max(qe,1e-37);
	} else {
	    ret_val = 1.f;
	}
    } else if (*imodel == 4) {
	hs = 0.;
	if (*h__ < hs) {
	    d__1 = log(-(*h__) / alfa) / n;
	    qee = qnorm_(&d__1);
	    ret_val = (real) max(qee,1e-37);
	} else {
	    ret_val = 1.f;
	}
    } else if (*imodel == 5) {
	w2 = par[7];
	alfa2 = par[8];
	n2 = par[9];
	m = 1. - 1. / n;
	m2 = 1. - 1. / n2;
	w1 = 1. - w2;
	d__2 = -alfa * *h__;
	d__1 = pow_dd(&d__2, &n) + 1.;
	d__3 = -m;
	sw1 = w1 * pow_dd(&d__1, &d__3);
	d__2 = -alfa2 * *h__;
	d__1 = pow_dd(&d__2, &n2) + 1.;
	d__3 = -m2;
	sw2 = w2 * pow_dd(&d__1, &d__3);
	qe = sw1 + sw2;
	ret_val = (real) max(qe,1e-37);
    } else if (*imodel == -1) {
/* Shlomo Orr */
	ha = alfa;
	d__ = n;
	ret_val = 1.f;
	if (ha > 0.f) {
/* Computing MAX */
/* Computing MIN */
	    d__4 = -(*h__) / ha;
	    d__5 = d__ - 3.f;
	    d__2 = 1.f, d__3 = (qs + pow_dd(&d__4, &d__5) - 1.f - qr) / (qs - 
		    qr);
	    d__1 = min(d__2,d__3);
	    ret_val = max(d__1,0.);
	}
    }
    return ret_val;
} /* fs_ */

/* *********************************************************************** */
doublereal fkq_(integer *imodel, real *th, real *par)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__, m, n, s, qa, kk, qe, kr, ks, qm, qk, qr, qs, qx, 
	    sx, qee, ffq, qek, alfa, bpar, qeek, ffqk, qees;
    static integer ppar;

    /* Parameter adjustments */
    --par;

    /* Function Body */
    qr = par[1];
    qs = par[2];
    alfa = par[3];
    n = par[4];
    ks = par[5];
    bpar = par[6];
    if (*imodel == 0 || *imodel == 1 || *imodel == 3) {
/* VG and mo */
	ppar = 2;
	if (*imodel == 0 || *imodel == 3) {
	    qm = qs;
	    qa = qr;
	    qk = qs;
	    kk = ks;
	} else if (*imodel == 1) {
	    qm = par[7];
	    qa = par[8];
	    qk = par[9];
	    kk = par[10];
	}
	if (*imodel == 3) {
	    qm = par[7];
	}
	m = 1. - 1. / n;
/* Computing MIN */
	d__1 = (qs - qa) / (qm - qa);
	qees = min(d__1,.999999999999999);
/* Computing MIN */
	d__1 = (qk - qa) / (qm - qa);
	qeek = min(d__1,qees);
	if ((doublereal) (*th) < qk) {
	    qee = ((doublereal) (*th) - qa) / (qm - qa);
	    qe = (qm - qa) / (qs - qa) * qee;
	    qek = (qm - qa) / (qs - qa) * qeek;
	    d__2 = 1. / m;
	    d__1 = 1. - pow_dd(&qee, &d__2);
	    ffq = 1. - pow_dd(&d__1, &m);
	    d__2 = 1. / m;
	    d__1 = 1. - pow_dd(&qeek, &d__2);
	    ffqk = 1. - pow_dd(&d__1, &m);
	    if (ffq <= 0.) {
		d__1 = 1. / m;
		ffq = m * pow_dd(&qee, &d__1);
	    }
	    d__1 = qe / qek;
	    d__2 = ffq / ffqk;
	    kr = pow_dd(&d__1, &bpar) * pow_di(&d__2, &ppar) * kk / ks;
/* Computing MAX */
	    d__1 = ks * kr;
	    ret_val = (real) max(d__1,1e-37);
	}
	if ((doublereal) (*th) >= qs) {
	    ret_val = (real) ks;
	}
    } else if (*imodel == -1) {
/* Shlomo Orr */
	d__ = n;
	kr = 1.f;
	qx = 0.f;
	if (d__ != 3.f) {
	    qx = qr + (1.f - qs) * 2.f / (d__ / (3.f - d__) - 2.f);
	}
	if (*th > qx) {
	    s = *th / qs;
	    if (d__ != 3.f) {
		d__1 = 1.f - qs * (1 - s) / (1.f - qr);
		d__2 = d__ / (3.f - d__);
		kr = pow_dd(&d__1, &d__2);
	    }
	} else {
	    sx = qx / qs;
	    if (d__ != 3.f) {
		d__1 = 1.f - qs * (1 - sx) / (1.f - qr);
		d__2 = d__ / (3.f - d__);
		kr = pow_dd(&d__1, &d__2);
		if (qx > qr) {
/* Computing 2nd power */
		    d__1 = *th - qr;
/* Computing 2nd power */
		    d__2 = qx - qr;
		    kr = kr * (d__1 * d__1) / (d__2 * d__2);
		}
	    }
	}
/* Computing MAX */
	d__1 = ks * kr;
	ret_val = (real) max(d__1,1e-37);
    }
    return ret_val;
} /* fkq_ */

/* *********************************************************************** */
doublereal fks_(integer *imodel, real *s, real *par)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__, m, n, qa, kk, qe, se;
    static real th;
    static doublereal kr, ks, qm, qk, qr, qs, qx, sx, qee, ffq, qek, fkq, 
	    alfa, bpar, qeek, ffqk, qees;
    static integer ppar;

    /* Parameter adjustments */
    --par;

    /* Function Body */
    qr = par[1];
    qs = par[2];
    alfa = par[3];
    n = par[4];
    ks = par[5];
    bpar = par[6];
    if (*imodel == 0 || *imodel == 1 || *imodel == 3) {
/* VG and modi */
	ppar = 2;
	if (*imodel == 0 || *imodel == 3) {
	    qm = qs;
	    qa = qr;
	    qk = qs;
	    kk = ks;
	} else if (*imodel == 1) {
	    qm = par[7];
	    qa = par[8];
	    qk = par[9];
	    kk = par[10];
	}
	if (*imodel == 3) {
	    qm = par[7];
	}
	m = 1. - 1. / n;
/* Computing MIN */
	d__1 = (qs - qa) / (qm - qa);
	qees = min(d__1,.999999999999999);
/* Computing MIN */
	d__1 = (qk - qa) / (qm - qa);
	qeek = min(d__1,qees);
	th = *s * (qs - qr) + qr;
	if ((doublereal) th < qk) {
	    qee = ((doublereal) th - qa) / (qm - qa);
	    qe = (qm - qa) / (qs - qa) * qee;
	    qek = (qm - qa) / (qs - qa) * qeek;
	    d__2 = 1. / m;
	    d__1 = 1. - pow_dd(&qee, &d__2);
	    ffq = 1. - pow_dd(&d__1, &m);
	    d__2 = 1. / m;
	    d__1 = 1. - pow_dd(&qeek, &d__2);
	    ffqk = 1. - pow_dd(&d__1, &m);
	    if (ffq <= 0.) {
		d__1 = 1. / m;
		ffq = m * pow_dd(&qee, &d__1);
	    }
	    d__1 = qe / qek;
	    d__2 = ffq / ffqk;
	    kr = pow_dd(&d__1, &bpar) * pow_di(&d__2, &ppar) * kk / ks;
/* Computing MAX */
	    d__1 = ks * kr;
	    ret_val = (real) max(d__1,1e-37);
	}
	if ((doublereal) th >= qs) {
	    ret_val = (real) ks;
	}
    } else if (*imodel == -1) {
/* Shlomo Orr */
	d__ = n;
	kr = 1.f;
	qx = 0.f;
	th = *s * (qs - qr) + qr;
	if (d__ != 3.f) {
	    qx = qr + (1.f - qs) * 2.f / (d__ / (3.f - d__) - 2.f);
	}
	if (th > qx) {
	    se = th / qs;
	    if (d__ != 3.f) {
		d__1 = 1.f - qs * (1 - se) / (1.f - qr);
		d__2 = d__ / (3.f - d__);
		kr = pow_dd(&d__1, &d__2);
	    }
	} else {
	    sx = qx / qs;
	    if (d__ != 3.f) {
		d__1 = 1.f - qs * (1 - sx) / (1.f - qr);
		d__2 = d__ / (3.f - d__);
		kr = pow_dd(&d__1, &d__2);
		if (qx > qr) {
/* Computing 2nd power */
		    d__1 = th - qr;
/* Computing 2nd power */
		    d__2 = qx - qr;
		    kr = kr * (d__1 * d__1) / (d__2 * d__2);
		}
	    }
	}
/* Computing MAX */
	d__1 = ks * kr;
	fkq = (real) max(d__1,1e-37);
    }
    return ret_val;
} /* fks_ */

/* *********************************************************************** */
doublereal qnorm_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal t, z__, erfc;

    z__ = (d__1 = *x / pow_dd(&c_b13, &c_b16), abs(d__1));
    t = 1.f / (z__ * .5f + 1.f);
    erfc = t * exp(-z__ * z__ - 1.26551223f + t * (t * (t * (t * (t * (t * (t 
	    * (t * (t * .17087277f - .82215223f) + 1.48851587f) - 1.13520398f)
	     + .27886807f) - .18628806f) + .09678418f) + .37409196f) + 
	    1.00002368f));
    if (*x < 0.f) {
	erfc = 2.f - erfc;
    }
    ret_val = erfc / 2.f;
    return ret_val;
} /* qnorm_ */

/* *********************************************************************** */
/*     Evaluate h for given theta_e for dual-porosity function */
doublereal xmualem_(doublereal *se, real *par, integer *npar)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal x1, x2, xb1, xb2, hhh;
    extern /* Subroutine */ int zbrak_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, real *, integer *);
    extern doublereal zbrent_(doublereal *, doublereal *, doublereal *, real *
	    , integer *);

    /* Parameter adjustments */
    --par;

    /* Function Body */
    x1 = -1e-6f;
    x2 = -1e6f;
    zbrak_(&x1, &x2, &xb1, &xb2, se, &par[1], npar);
    hhh = zbrent_(&xb1, &xb2, se, &par[1], npar);
    ret_val = hhh;
/* for calculation hh */
    if (hhh != 0.f) {
/*        xMualem=1./hhh ! for integration */
    } else {
	s_paus("xMualem: h is equal to zero!", (ftnlen)28);
    }
    return ret_val;
} /* xmualem_ */

/* *********************************************************************** */
doublereal doublepor_(doublereal *hh, doublereal *se, real *par, integer *
	npar)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static doublereal w1, w2, rm, rn, rm2, rn2, sw1, sw2, wcr, wcs, rwc, 
	    alpha, alpha2;

/*     Double porosity function - for evaluation of h for given theta_e */
    /* Parameter adjustments */
    --par;

    /* Function Body */
    wcr = par[1];
    wcs = par[2];
    alpha = par[3];
    rn = par[4];
    rm = 1.f - 1.f / rn;
    w2 = par[7];
    w1 = 1.f - w2;
    alpha2 = par[8];
    rn2 = par[9];
    rm2 = 1.f - 1.f / rn2;
    d__2 = -alpha * *hh;
    d__1 = pow_dd(&d__2, &rn) + 1.f;
    d__3 = -rm;
    sw1 = w1 * pow_dd(&d__1, &d__3);
    d__2 = -alpha2 * *hh;
    d__1 = pow_dd(&d__2, &rn2) + 1.f;
    d__3 = -rm2;
    sw2 = w2 * pow_dd(&d__1, &d__3);
    rwc = sw1 + sw2;
    ret_val = *se - rwc;
    return ret_val;
} /* doublepor_ */

/* *********************************************************************** */
/*     Bracketing of the root, Numerical recepies (345) */
/* Subroutine */ int zbrak_(doublereal *x1, doublereal *x2, doublereal *xb1, 
	doublereal *xb2, doublereal *se, real *par, integer *npar)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, n;
    extern doublereal doublepor_(doublereal *, doublereal *, real *, integer *
	    );
    static doublereal fc;
    static integer nb;
    static doublereal fp, dx2;
    static integer nbb;
    static doublereal dlh;

    /* Parameter adjustments */
    --par;

    /* Function Body */
    nb = 1;
    nbb = nb;
    nb = 0;
    n = 1000;
    d__1 = -(*x2);
    d__2 = -(*x1);
    dlh = (d_lg10(&d__1) - d_lg10(&d__2)) / (n - 1);
    fp = doublepor_(x1, se, &par[1], npar);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = -(*x1);
	dx2 = d_lg10(&d__1) + i__ * dlh;
	*x2 = -pow_dd(&c_b9, &dx2);
	fc = doublepor_(x2, se, &par[1], npar);
	if (fc * fp < 0.f) {
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
} /* zbrak_ */

/* *********************************************************************** */
/*     Brent method of finding root that lies between x1 and x2, */
/*     Numerical recepies (354) */
doublereal zbrent_(doublereal *x1, doublereal *x2, doublereal *se, real *par, 
	integer *npar)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal a, b, c__, d__, e, p, q, r__, s;
    extern doublereal doublepor_(doublereal *, doublereal *, real *, integer *
	    );
    static doublereal fa, fb, fc, xm, tol1;
    static integer iter;

    /* Parameter adjustments */
    --par;

    /* Function Body */
    a = *x1;
    b = *x2;
    fa = doublepor_(&a, se, &par[1], npar);
    fb = doublepor_(&b, se, &par[1], npar);
    if (fb * fa > 0.f) {
	s_paus("Root must be bracketed for ZBRENT.", (ftnlen)34);
    }
    fc = fb;
    for (iter = 1; iter <= 100; ++iter) {
	if (fb * fc > 0.f) {
	    c__ = a;
	    fc = fa;
	    d__ = b - a;
	    e = d__;
	}
	if (abs(fc) < abs(fb)) {
	    a = b;
	    b = c__;
	    c__ = a;
	    fa = fb;
	    fb = fc;
	    fc = fa;
	}
	tol1 = abs(b) * 5.9999999999999995e-8 + 4.9999999999999998e-7;
	xm = (c__ - b) * .5f;
	if (abs(xm) <= tol1 || fb == 0.f) {
	    ret_val = b;
	    return ret_val;
	}
	if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
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
	    p = abs(p);
/* Computing MIN */
	    d__3 = xm * 3.f * q - (d__1 = tol1 * q, abs(d__1)), d__4 = (d__2 =
		     e * q, abs(d__2));
	    if (p * 2.f < min(d__3,d__4)) {
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
	if (abs(d__) > tol1) {
	    b += d__;
	} else {
	    b += d_sign(&tol1, &xm);
	}
	fb = doublepor_(&b, se, &par[1], npar);
/* L11: */
    }
    s_paus("ZBRENT exceeding maximum iterations.", (ftnlen)36);
    ret_val = b;
    return ret_val;
} /* zbrent_ */

/* *********************************************************************** */
doublereal gamma_(doublereal *z__)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal x, y, fy;

/*     Purpose:  To calculate the Gamma function for positive Z */
    if (*z__ < 33.f) {
	goto L11;
    }
    ret_val = 1e36;
    return ret_val;
L11:
    x = *z__;
    ret_val = 1.f;
    if (x - 2.f <= 0.) {
	goto L14;
    } else {
	goto L13;
    }
L12:
    if (x - 2.f <= 0.) {
	goto L16;
    } else {
	goto L13;
    }
L13:
    x += -1.f;
    ret_val *= x;
    goto L12;
L14:
    if ((d__1 = x - 1.f) < 0.) {
	goto L15;
    } else if (d__1 == 0) {
	goto L17;
    } else {
	goto L16;
    }
L15:
    ret_val /= x;
    x += 1.f;
L16:
    y = x - 1.f;
    fy = 1.f - y * (.5771017f - y * (.985854f - y * (.8764218f - y * (
	    .8328212f - y * (.5684729f - y * (.2548205f - y * .0514993f))))));
    ret_val *= fy;
L17:
    return ret_val;
} /* gamma_ */

/* *********************************************************************** */
doublereal binc_(doublereal *x, doublereal *a, doublereal *b, doublereal *
	beta)
{
    /* Initialized data */

    static integer nt = 10;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__, k;
    static doublereal t[200], y, y2;
    static integer nt1;

/*     Purpose: To calculate the incomplete Beta-function */
    nt1 = nt + 1;
    t[0] = -(*a + *b) * *x / (*a + 1.f);
    i__1 = nt;
    for (i__ = 2; i__ <= i__1; i__ += 2) {
	y = (real) (i__ / 2);
	y2 = (real) i__;
	t[i__ - 1] = y * (*b - y) * *x / ((*a + y2 - 1.f) * (*a + y2));
	t[i__] = -(*a + y) * (*a + *b + y) * *x / ((*a + y2) * (*a + y2 + 1.f)
		);
/* L11: */
    }
    ret_val = 1.f;
    i__1 = nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = nt1 - i__;
	ret_val = t[k - 1] / ret_val + 1.f;
/* L12: */
    }
    d__1 = 1.f - *x;
    ret_val = pow_dd(x, a) * pow_dd(&d__1, b) / (ret_val * *a * *beta);
    return ret_val;
} /* binc_ */

/* *********************************************************************** */
/* Subroutine */ int qromb_(real *a, real *b, real *ss, integer *imodel, real 
	*par)
{
    /* Local variables */
    static real h__[21];
    static integer j;
    static real s[21], dss;
    extern /* Subroutine */ int trapzd_(real *, real *, real *, integer *, 
	    integer *, real *), polint_(real *, real *, integer *, real *, 
	    real *, real *);

    /* Parameter adjustments */
    --par;

    /* Function Body */
    h__[0] = 1.f;
    for (j = 1; j <= 20; ++j) {
	trapzd_(a, b, &s[j - 1], &j, imodel, &par[1]);
	if (j >= 5) {
	    polint_(&h__[j - 5], &s[j - 5], &c__5, &c_b53, ss, &dss);
	    if (dabs(dss) <= dabs(*ss) * 1e-6f) {
		return 0;
	    }
	}
	s[j] = s[j - 1];
	h__[j] = h__[j - 1] * .25f;
/* L11: */
    }
    s_paus("too many steps in qromb", (ftnlen)23);
    return 0;
} /* qromb_ */

/* ********************************************************************** */
/* Subroutine */ int trapzd_(real *a, real *b, real *s, integer *n, integer *
	imodel, real *par)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static real x;
    extern doublereal fh_(integer *, real *, real *);
    static integer it;
    static real del, tnm, sum;

    /* Parameter adjustments */
    --par;

    /* Function Body */
    if (*n == 1) {
	*s = (*b - *a) * .5f * (fh_(imodel, a, &par[1]) + fh_(imodel, b, &par[
		1]));
    } else {
	i__1 = *n - 2;
	it = pow_ii(&c__2, &i__1);
	tnm = (real) it;
	del = (*b - *a) / tnm;
	x = *a + del * .5f;
	sum = 0.f;
	i__1 = it;
	for (j = 1; j <= i__1; ++j) {
	    sum += fh_(imodel, &x, &par[1]);
	    x += del;
/* L11: */
	}
	*s = (*s + (*b - *a) * sum / tnm) * .5f;
    }
    return 0;
} /* trapzd_ */

/* *********************************************************************** */
/* Subroutine */ int polint_(real *xa, real *ya, integer *n, real *x, real *y,
	 real *dy)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static real c__[10], d__[10];
    static integer i__, m;
    static real w, ho, hp;
    static integer ns;
    static real dif, den, dift;

    /* Parameter adjustments */
    --ya;
    --xa;

    /* Function Body */
    ns = 1;
    dif = (r__1 = *x - xa[1], dabs(r__1));
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dift = (r__1 = *x - xa[i__], dabs(r__1));
	if (dift < dif) {
	    ns = i__;
	    dif = dift;
	}
	c__[i__ - 1] = ya[i__];
	d__[i__ - 1] = ya[i__];
/* L11: */
    }
    *y = ya[ns];
    --ns;
    i__1 = *n - 1;
    for (m = 1; m <= i__1; ++m) {
	i__2 = *n - m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ho = xa[i__] - *x;
	    hp = xa[i__ + m] - *x;
	    w = c__[i__] - d__[i__ - 1];
	    den = ho - hp;
	    if (den == 0.f) {
		s_paus("failure in polint", (ftnlen)17);
	    }
	    den = w / den;
	    d__[i__ - 1] = hp * den;
	    c__[i__ - 1] = ho * den;
/* L12: */
	}
	if (ns << 1 < *n - m) {
	    *dy = c__[ns];
	} else {
	    *dy = d__[ns - 1];
	    --ns;
	}
	*y += *dy;
/* L13: */
    }
    return 0;
} /* polint_ */

