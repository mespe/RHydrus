/* OUTPUT.f -- translated by f2c (version 12.02.01).
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
static integer c__9 = 9;
static doublereal c_b604 = .93196644920782856;

/* Source file OUTPUT.FOR ||||||||||||||||||||||||||||||||||||||||||||||| */
/* Subroutine */ int tlinf_(integer *n, real *con, real *x, real *cosalf, 
	doublereal *t, real *dt, integer *iterw, integer *iterc, integer *
	tlevel, real *rtop, real *rroot, real *vroot, real *hnew, real *hroot,
	 real *cumq, integer *itcum, integer *kodtop, integer *kodbot, 
	logical *convgf, logical *lwat, logical *lchem, real *croot, integer *
	ns, integer *nsd, real *conc, real *cvtop, real *cvbot, real *cvch0, 
	real *cvch1, real *peclet, real *courant, real *wcumt, real *wcuma, 
	real *ccumt, real *ccuma, real *cumch, real *thnew, real *thold, real 
	*sink, logical *lscreen, integer *ierr, real *cvchr, real *cvchim, 
	logical *lprint, logical *lvapor, logical *lwtdep, real *conlt, real *
	convh, real *convt, real *temp, real *rsoil, real *prec, integer *
	nprstep, real *thvold, real *thvnew, real *xconv, integer *idualpor, 
	real *sinkim, real *wtransf, logical *ldensity, real *snowlayer, 
	logical *lcentrif, real *radius, real *thnewim, logical *wlayer, real 
	*hcrits, logical *lend, logical *lfluxout, integer *jprint, real *
	vtop, real *vbot, logical *lflux, integer *nobs, integer *node, real *
	vnew, real *cnew, real *ctop)
{
    /* Format strings */
    static char fmt_110[] = "(/\002         Time ItW   ItCum  vTop    SvTop "
	    "   SvRoot   SvBot   \002,\002 hTop hRoot hBot\002/)";
    static char fmt_122[] = "(f14.4,i3,i7,4e9.2,f7.2,2f5.2)";
    static char fmt_120[] = "(f13.4,i3,i7,4e9.2,f8.1,2f6.0)";
    static char fmt_123[] = "(e14.7,i3,i7,4e9.2,f7.2,2f5.2)";
    static char fmt_121[] = "(e14.7,i3,i7,4e9.2,f7.0,2f6.0)";
    static char fmt_130[] = "(/\002       Time          rTop        rRoot   "
	    "     vTop         vRoot        vBot       sum(rTop)   sum(rRoot)"
	    "    sum(vTop)   sum(vRoot)    sum(vBot)      hTop         hRoot "
	    "       hBot        RunOff    sum(RunOff)     Volume     sum(Infi"
	    "l)    sum(Evap) TLevel Cum(WTrans)  SnowLayer\002/\002        [T"
	    "]         [L/T]        [L/T]        [L/T]        [L/T]        [L"
	    "/T]         [L]          [L]          [L]         [L]           "
	    "[L]         [L]           [L]         [L]          [L/T]        "
	    " [L]          [L]          [L]          [L]\002/)";
    static char fmt_111[] = "(/\002          time     PrecipP      EvaporP  "
	    "   FluxTopP     FluxTopA     InfiltrA      EvaporA      TranspP "
	    "     TranspA      FluxBot       RunOff     Sum(PrecP)   Sum(Evap"
	    "P)  Sum(FlTopP)  Sum(FlTopA)   Sum(InfA)   Sum(EvapA)  Sum(Trans"
	    "P)  Sum(TransA)  Sum(FluxBot)  Sum(RunOff)    Storage   SnowLaye"
	    "r      vTopW        vTopV\002/)";
    static char fmt_150[] = "(//\002    TLevel      Time           dt     It"
	    "rW ItrC    ItC\002,\002um  KodT  KodB Converg  Peclet   Couran"
	    "t\002/)";
    static char fmt_160[] = "(\002 All solute fluxes and cumulative solute f"
	    "luxes are positive into the region\002//\002       Time         "
	    "cvTop        cvBot      Sum(cvTop)   Sum(cvBot)     cvCh0       "
	    " cvCh1         cTop        cRoot         cBot        cvRoot    S"
	    "um(cvRoot)  Sum(cvNEql) TLevel      cGWL        cRunOff   Sum(cR"
	    "unOff)    (cv(i),    Sum(cv(i)), i=1,NObs)\002/\002        [T]  "
	    "      [M/L2/T]     [M/L2/T]      [M/L2]       [M/L2]       [M/L2"
	    "]      [M/L2]        [M/L3]      [M/L3]        [M/L3]      [M/L2"
	    "/T]      [M/L2]       [M/L2]              [M/L3]        [M/L2]  "
	    "    [M/L3]      [M/L2/T]      [M/L2]\002)";
    static char fmt_140[] = "(//\002    TLevel      Time          dt      It"
	    "er    ItCum  K\002,\002odT  KodB  Convergency\002/)";
    static char fmt_170[] = "(f13.4,11e13.5,2e13.5,5e13.5,i7,e13.5,f11.3)";
    static char fmt_171[] = "(e14.8,11e13.5,2e13.5,5e13.5,i7,e13.5,f11.3)";
    static char fmt_180[] = "(i9,e15.7,e13.5,2i5,i9,2i6,l6,2f10.3)";
    static char fmt_190[] = "(f14.4,12e13.5,i8,e13.5,8e13.5)";
    static char fmt_191[] = "(e15.8,12e13.5,i8,e13.5,8e13.5)";
    static char fmt_200[] = "(i9,e15.7,e13.5,i5,i9,2i6,l6,2f10.3)";
    static char fmt_222[] = "(f14.4,24e13.5)";

    /* System generated locals */
    integer conc_dim1, conc_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4, r__5;

    /* Local variables */
    static integer i__, j, m;
    static real vb, dx;
    static integer js;
    static real vt, dx1, fre;
    extern doublereal fro_(integer *, real *);
    static real dxn, cgwl[15], dgwl, grav;
    static integer mobs;
    static real revap, vnewi, vtopv, vtopw;
    static integer ibreak;
    static real rinfil, volume, crunoff[15], vnewimi, vrunoff;

    /* Fortran I/O blocks */
    static cilist io___24 = { 0, 6, 0, fmt_110, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_122, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_120, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_123, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_121, 0 };
    static cilist io___29 = { 1, 71, 0, fmt_130, 0 };
    static cilist io___30 = { 0, 44, 0, fmt_111, 0 };
    static cilist io___31 = { 1, 70, 0, fmt_150, 0 };
    static cilist io___32 = { 1, 0, 0, fmt_160, 0 };
    static cilist io___33 = { 1, 70, 0, fmt_140, 0 };
    static cilist io___34 = { 1, 71, 0, fmt_170, 0 };
    static cilist io___35 = { 1, 71, 0, fmt_171, 0 };
    static cilist io___36 = { 1, 70, 0, fmt_180, 0 };
    static cilist io___38 = { 1, 0, 0, fmt_190, 0 };
    static cilist io___39 = { 1, 0, 0, fmt_191, 0 };
    static cilist io___40 = { 1, 70, 0, fmt_200, 0 };
    static cilist io___41 = { 0, 44, 0, fmt_222, 0 };


    /* Parameter adjustments */
    --cnew;
    --vnew;
    --thnewim;
    --sinkim;
    --thvnew;
    --thvold;
    --temp;
    --convt;
    --convh;
    --conlt;
    --sink;
    --thold;
    --thnew;
    --hnew;
    --x;
    --con;
    --cumq;
    --ctop;
    --cvchim;
    --cvchr;
    cumch -= 11;
    --ccuma;
    --ccumt;
    --cvch1;
    --cvch0;
    --cvbot;
    --cvtop;
    --croot;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    --node;

    /* Function Body */
    fre = 1.f;
    grav = *cosalf;
    m = *n - 1;
    if (*ldensity) {
	fre = fro_(&c__1, &conc[*n * conc_dim1 + 1]);
    }
    if (*lcentrif) {
	grav = *cosalf * (*radius + (r__1 = (x[*n] + x[m]) / 2.f, dabs(r__1)))
		;
    }
    dxn = x[*n] - x[m];
    vt = -(con[*n] + con[m]) / 2.f * ((hnew[*n] - hnew[m]) / dxn + grav * fre)
	     - (thnew[*n] - thold[*n]) * fre * dxn / 2.f / *dt - sink[*n] * 
	    dxn / 2.f;
    if (*ldensity) {
	fre = fro_(&c__1, &conc[conc_dim1 + 1]);
    }
    dx1 = x[2] - x[1];
    if (*lcentrif) {
	grav = *cosalf * (*radius + (r__1 = (x[2] + x[1]) / 2.f, dabs(r__1)));
    }
    vb = -(con[1] + con[2]) / 2.f * ((hnew[2] - hnew[1]) / dx1 + grav * fre) 
	    + (thnew[1] - thold[1]) * fre * dx1 / 2.f / *dt + sink[1] * dx1 / 
	    2.f;
    if (*idualpor > 0) {
	vt -= sinkim[*n] * dxn / 2.f;
    }
    if (*lwtdep) {
	vt -= (conlt[*n] + conlt[m]) / 2.f * (temp[*n] - temp[m]) / dxn;
	vb -= (conlt[1] + conlt[2]) / 2.f * (temp[2] - temp[1]) / dx1;
    }
    vtopw = vt;
    if (*lvapor) {
	vt -= (convh[*n] + convh[m]) / 2.f * (hnew[*n] - hnew[m]) / dxn;
	vb -= (convh[1] + convh[2]) / 2.f * (hnew[2] - hnew[1]) / dx1;
	vt = vt - (convt[*n] + convt[m]) / 2.f * (temp[*n] - temp[m]) / dxn - 
		(thvnew[*n] - thvold[*n]) * dxn / 2.f / *dt;
	vb = vb - (convt[1] + convt[2]) / 2.f * (temp[2] - temp[1]) / dx1 + (
		thvnew[1] - thvold[1]) * dx1 / 2.f / *dt;
    }
    vtopv = vt - vtopw;
    *vtop = vt;
    *vbot = vb;
    vrunoff = 0.f;
    if ((! (*wlayer) || *wlayer && hnew[*n] >= *hcrits) && *rtop < 0.f) {
	vrunoff = (r__1 = *rtop - *vtop, dabs(r__1));
    }
    if (vrunoff < 1e-5f) {
	vrunoff = 0.f;
    }
    rinfil = 0.f;
    revap = 0.f;
    if (*vtop < 0.f && (*prec > 0.f || *wlayer && hnew[*n] > 0.f)) {
	rinfil = -(*vtop) + *rsoil;
    }
    if (*vtop >= 0.f && *prec > 0.f) {
	rinfil = *prec;
    }
    if (*vtop > 0.f) {
	revap = *vtop + *prec;
    }
    if (*vtop <= 0.f && *rsoil > 0.f && *prec > 0.f) {
	revap = *rsoil;
    }
    if (*vtop < 0.f && *wlayer && hnew[*n] > 0.f) {
	revap = *rsoil;
    }
    if (*vtop < 0.f && *kodtop > 0) {
	rinfil = -(*vtop);
    }
    cumq[1] += *rtop * *dt;
    cumq[2] += *rroot * *dt;
    cumq[3] += *vtop * *dt;
    cumq[4] += *vroot * *dt;
    cumq[5] += *vbot * *dt;
    cumq[6] += vrunoff * *dt;
    cumq[7] += rinfil * *dt;
    cumq[8] += revap * *dt;
    if (*lfluxout) {
	cumq[9] += *prec * *dt;
	cumq[10] += *rsoil * *dt;
    }
    cumq[11] += *wtransf * *dt;
    *wcumt += (*vbot - *vtop - *vroot) * *dt;
    *wcuma += (dabs(*vbot) + dabs(*vtop) + dabs(*vroot)) * *dt;
    if (*lchem) {
	i__1 = *ns;
	for (js = 1; js <= i__1; ++js) {
	    cumch[js * 10 + 1] -= cvtop[js] * *dt;
	    cumch[js * 10 + 2] += cvbot[js] * *dt;
	    cumch[js * 10 + 3] += cvch0[js] * *dt;
	    cumch[js * 10 + 4] += cvch1[js] * *dt;
	    cumch[js * 10 + 5] += cvchr[js] * *dt;
	    cumch[js * 10 + 6] += cvchim[js] * *dt;
	    ccumt[js] += (cvtop[js] - cvbot[js] - cvch0[js] - cvch1[js] + 
		    cvchr[js]) * *dt;
	    ccuma[js] += ((r__1 = cvbot[js], dabs(r__1)) + (r__2 = cvtop[js], 
		    dabs(r__2)) + (r__3 = cvch0[js], dabs(r__3)) + (r__4 = 
		    cvch1[js], dabs(r__4)) + (r__5 = cvchr[js], dabs(r__5))) *
		     *dt;
	    if (*lflux) {
		if (js == 1) {
/* using flux concentration (available) */
		    if (*nobs >= 1) {
			cumch[js * 10 + 7] += vnew[node[1]] * cnew[node[1]] * 
				*dt;
		    }
		    if (*nobs >= 2) {
			cumch[js * 10 + 8] += vnew[node[2]] * cnew[node[2]] * 
				*dt;
		    }
		    if (*nobs >= 3) {
			cumch[js * 10 + 9] += vnew[node[3]] * cnew[node[3]] * 
				*dt;
		    }
		} else {
/* using resident concentration (flux conc u */
		    if (*nobs >= 1) {
			cumch[js * 10 + 7] += vnew[node[1]] * conc[js + node[
				1] * conc_dim1] * *dt;
		    }
		    if (*nobs >= 2) {
			cumch[js * 10 + 8] += vnew[node[2]] * conc[js + node[
				2] * conc_dim1] * *dt;
		    }
		    if (*nobs >= 3) {
			cumch[js * 10 + 9] += vnew[node[3]] * conc[js + node[
				3] * conc_dim1] * *dt;
		    }
		}
	    }
	    crunoff[js - 1] = vrunoff * ctop[js];
	    cumch[js * 10 + 10] += crunoff[js - 1] * *dt;
/*         Average GWL concentration */
	    cgwl[js - 1] = 0.f;
	    dgwl = 0.f;
	    ibreak = 0;
	    i__2 = *n - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		j = i__ + 1;
		dx = x[j] - x[i__];
		if (hnew[j] > 0.f && ibreak == 0) {
		    cgwl[js - 1] += (conc[js + i__ * conc_dim1] + conc[js + j 
			    * conc_dim1]) / 2.f * dx;
		    dgwl += dx;
		} else {
		    ibreak = 1;
		}
/* L14: */
	    }
	    if (dgwl > 0.f) {
		cgwl[js - 1] /= dgwl;
	    }
/* L11: */
	}
    }
    volume = 0.f;
    for (i__ = *n - 1; i__ >= 1; --i__) {
	j = i__ + 1;
	dx = x[j] - x[i__];
	vnewi = dx * (thnew[i__] + thnew[j]) / 2.f;
	if (*ldensity) {
	    r__1 = conc[i__ * conc_dim1 + 1];
	    r__2 = conc[j * conc_dim1 + 1];
	    vnewi = dx * (thnew[i__] * fro_(&c__1, &r__1) + thnew[j] * fro_(&
		    c__1, &r__2)) / 2.f;
	}
	volume += vnewi;
	if ((real) (*idualpor) > 0.f) {
	    vnewimi = dx * (thnewim[i__] + thnewim[j]) / 2.f;
/*          if(lDensity) VNewImi=dx*(ThNewIm(i)*fRo(1,(Sorb(1,i)))+ */
/*     !                             ThNewIm(j)*fRo(1,(Sorb(1,j))))/2. */
	    volume += vnewimi;
	}
/* L10: */
    }
    if (*lscreen) {
	if (*jprint == 1 && (r__1 = (real) ((*tlevel + *nprstep * 20 - 1) / 
		20 / *nprstep) - (*tlevel + *nprstep * 20 - 1) / (real) (*
		nprstep * 20), dabs(r__1)) < 1e-4f) {
	    s_wsfe(&io___24);
	    e_wsfe();
	}
	if (*jprint == 1) {
	    if (*t < 9999999.f) {
		if (*xconv == 1.f) {
		    s_wsfe(&io___25);
		    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&(*iterw), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*itcum), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*vtop), (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[3], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[4], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[5], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&hnew[*n], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&(*hroot), (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&hnew[1], (ftnlen)sizeof(real));
		    e_wsfe();
		} else {
		    s_wsfe(&io___26);
		    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&(*iterw), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*itcum), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*vtop), (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[3], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[4], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[5], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&hnew[*n], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&(*hroot), (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&hnew[1], (ftnlen)sizeof(real));
		    e_wsfe();
		}
	    } else {
		if (*xconv == 1.f) {
		    s_wsfe(&io___27);
		    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&(*iterw), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*itcum), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*vtop), (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[3], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[4], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[5], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&hnew[*n], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&(*hroot), (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&hnew[1], (ftnlen)sizeof(real));
		    e_wsfe();
		} else {
		    s_wsfe(&io___28);
		    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&(*iterw), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*itcum), (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*vtop), (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[3], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[4], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&cumq[5], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&hnew[*n], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&(*hroot), (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&hnew[1], (ftnlen)sizeof(real));
		    e_wsfe();
		}
	    }
	}
    }
    if (*tlevel == 1 && *lprint) {
	i__1 = s_wsfe(&io___29);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	if (*lfluxout) {
	    s_wsfe(&io___30);
	    e_wsfe();
	}
	if (*lchem) {
	    i__1 = s_wsfe(&io___31);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = *ns;
	    for (js = 1; js <= i__1; ++js) {
		io___32.ciunit = js + 80;
		i__2 = s_wsfe(&io___32);
		if (i__2 != 0) {
		    goto L903;
		}
		i__2 = e_wsfe();
		if (i__2 != 0) {
		    goto L903;
		}
/* L12: */
	    }
	} else {
	    i__1 = s_wsfe(&io___33);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	}
    }
    if (*lprint && *jprint == 1) {
	if (*lwat || *tlevel == 1 || *lend) {
	    if (*t < 9999999.f) {
/*            write(71,170,err=901) t,rTop,vTopW,vTop,vTopV,vBot,(CumQ(i), */
		i__1 = s_wsfe(&io___34);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*rtop), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*rroot), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*vtop), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*vroot), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*vbot), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		for (i__ = 1; i__ <= 5; ++i__) {
		    i__1 = do_fio(&c__1, (char *)&cumq[i__], (ftnlen)sizeof(
			    real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = do_fio(&c__1, (char *)&hnew[*n], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*hroot), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&hnew[1], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&vrunoff, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cumq[6], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&volume, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cumq[7], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cumq[8], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*tlevel), (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cumq[11], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*snowlayer), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    } else {
		i__1 = s_wsfe(&io___35);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*rtop), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*rroot), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*vtop), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*vroot), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*vbot), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		for (i__ = 1; i__ <= 5; ++i__) {
		    i__1 = do_fio(&c__1, (char *)&cumq[i__], (ftnlen)sizeof(
			    real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = do_fio(&c__1, (char *)&hnew[*n], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*hroot), (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&hnew[1], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&vrunoff, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cumq[6], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&volume, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cumq[7], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cumq[8], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*tlevel), (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cumq[11], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*snowlayer), (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	}
	if (*lchem) {
	    i__1 = s_wsfe(&io___36);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*tlevel), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*dt), (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*iterw), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*iterc), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*itcum), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*kodtop), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*kodbot), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*convgf), (ftnlen)sizeof(logical));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*peclet), (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*courant), (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = *ns;
	    for (js = 1; js <= i__1; ++js) {
		mobs = *nobs;
		if (! (*lflux)) {
		    mobs = 0;
		}
		if (*t < 99999999.f) {
		    io___38.ciunit = js + 80;
		    i__2 = s_wsfe(&io___38);
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__2 != 0) {
			goto L903;
		    }
		    r__1 = -cvtop[js];
		    i__2 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cvbot[js], (ftnlen)sizeof(
			    real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 1], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 2], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 3], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 4], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&conc[js + *n * conc_dim1], (
			    ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&croot[js], (ftnlen)sizeof(
			    real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&conc[js + conc_dim1], (
			    ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cvchr[js], (ftnlen)sizeof(
			    real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 5], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 6], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&(*tlevel), (ftnlen)sizeof(
			    integer));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cgwl[js - 1], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&crunoff[js - 1], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 10], (
			    ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__3 = min(mobs,3);
		    for (j = 1; j <= i__3; ++j) {
			r__2 = vnew[node[j]] * conc[js + node[j] * conc_dim1];
			i__2 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__2 != 0) {
			    goto L903;
			}
			i__2 = do_fio(&c__1, (char *)&cumch[j + 6 + js * 10], 
				(ftnlen)sizeof(real));
			if (i__2 != 0) {
			    goto L903;
			}
		    }
		    i__2 = e_wsfe();
		    if (i__2 != 0) {
			goto L903;
		    }
		} else {
		    io___39.ciunit = js + 80;
		    i__2 = s_wsfe(&io___39);
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__2 != 0) {
			goto L903;
		    }
		    r__1 = -cvtop[js];
		    i__2 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cvbot[js], (ftnlen)sizeof(
			    real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 1], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 2], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 3], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 4], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&conc[js + *n * conc_dim1], (
			    ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&croot[js], (ftnlen)sizeof(
			    real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&conc[js + conc_dim1], (
			    ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cvchr[js], (ftnlen)sizeof(
			    real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 5], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 6], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&(*tlevel), (ftnlen)sizeof(
			    integer));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cgwl[js - 1], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&crunoff[js - 1], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__2 = do_fio(&c__1, (char *)&cumch[js * 10 + 10], (
			    ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L903;
		    }
		    i__3 = min(mobs,3);
		    for (j = 1; j <= i__3; ++j) {
			r__2 = vnew[node[j]] * conc[js + node[j] * conc_dim1];
			i__2 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__2 != 0) {
			    goto L903;
			}
			i__2 = do_fio(&c__1, (char *)&cumch[j + 6 + js * 10], 
				(ftnlen)sizeof(real));
			if (i__2 != 0) {
			    goto L903;
			}
		    }
		    i__2 = e_wsfe();
		    if (i__2 != 0) {
			goto L903;
		    }
		}
/* L13: */
	    }
	} else {
	    i__1 = s_wsfe(&io___40);
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*tlevel), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*dt), (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*iterw), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*itcum), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*kodtop), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*kodbot), (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*convgf), (ftnlen)sizeof(logical));
	    if (i__1 != 0) {
		goto L902;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L902;
	    }
	}
    }
    if (*lprint && *lfluxout && *jprint == 1) {
	s_wsfe(&io___41);
	do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*prec), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*rsoil), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*rtop), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*vtop), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rinfil, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&revap, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*rroot), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*vroot), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*vbot), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vrunoff, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[9], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[10], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[3], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[7], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[8], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[2], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[4], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[5], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cumq[6], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&volume, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*snowlayer), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vtopw, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vtopv, (ftnlen)sizeof(real));
	e_wsfe();
    }
    return 0;
/*     Error when writing into an output file */
L901:
    *ierr = 1;
    return 0;
L902:
    *ierr = 2;
    return 0;
L903:
    *ierr = 3;
    return 0;
/* 120   format('+',f12.3,2i3,i6,4e9.2,3f6.0)  ! writing at one line */
} /* tlinf_ */

/* *********************************************************************** */
/* Subroutine */ int alinf_(doublereal *t, real *cumq, real *hnewn, real *
	hroot, real *hnew1, integer *alevel, integer *ierr)
{
    /* Format strings */
    static char fmt_110[] = "(//\002   Time         sum(rTop)     sum(rRoot)"
	    "    sum(vTop)     sum(vRoot)     sum(vBot)    hTop       hRoot  "
	    "    hBot      A-level\002/\002    [T]           [L]           [L"
	    "]           [L]           [L]            [L]        [L]         "
	    "[L]       [L] \002/)";
    static char fmt_120[] = "(f12.5,5e14.6,3f11.3,i8)";
    static char fmt_121[] = "(e14.7,5e14.6,3f11.3,i8)";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___42 = { 1, 72, 0, fmt_110, 0 };
    static cilist io___43 = { 1, 72, 0, fmt_120, 0 };
    static cilist io___45 = { 1, 72, 0, fmt_121, 0 };


    /* Parameter adjustments */
    --cumq;

    /* Function Body */
    if (*alevel == 1) {
	i__1 = s_wsfe(&io___42);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    if (*t < 999999.f) {
	i__1 = s_wsfe(&io___43);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L901;
	}
	for (i__ = 1; i__ <= 5; ++i__) {
	    i__1 = do_fio(&c__1, (char *)&cumq[i__], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = do_fio(&c__1, (char *)&(*hnewn), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*hroot), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*hnew1), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*alevel), (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
    } else {
	i__1 = s_wsfe(&io___45);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L901;
	}
	for (i__ = 1; i__ <= 5; ++i__) {
	    i__1 = do_fio(&c__1, (char *)&cumq[i__], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = do_fio(&c__1, (char *)&(*hnewn), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*hroot), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*hnew1), (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*alevel), (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    return 0;
/*     Error when writing into an output file */
L901:
    *ierr = 1;
    return 0;
} /* alinf_ */

/* *********************************************************************** */
/* Subroutine */ int subreg_(integer *n, integer *nmat, integer *nlay, real *
	hnew, real *thn, real *tho, real *x, integer *matnum, integer *laynum,
	 doublereal *t, real *dt, real *cosalf, real *con, logical *lchem, 
	real *conc, real *chpar, integer *plevel, real *ths, real *wcumt, 
	real *wcuma, real *ccumt, real *ccuma, real *wvoli, real *cvoli, real 
	*watin, real *solin, logical *lwat, logical *ltemp, real *temp, real *
	tpar, real *tdep, integer *ns, integer *nsd, real *sorb, logical *
	llinear, logical *lequil, logical *lmobim, integer *ierr, real *
	subvol, real *area, logical *lprint, logical *lbact, real *sorb2, 
	logical *lvapor, real *thvold, real *thvnew, logical *lwtdep, real *
	conlt, real *convh, real *convt, integer *idualpor, real *thnewim, 
	real *tholdim, logical *ldensity, logical *lcentrif, real *radius, 
	logical *ldualneq, real *cprev)
{
    /* Format strings */
    static char fmt_110[] = "(/\002-----------------------------------------"
	    "-----------------\002/\002 Time       [T]\002,f14.4/\002--------"
	    "--------------------------------------------------\002)";
    static char fmt_111[] = "(/\002-----------------------------------------"
	    "-----------------\002/\002 Time       [T]\002,e15.8/\002--------"
	    "--------------------------------------------------\002)";
    static char fmt_120[] = "(\002 Sub-region num.               \002,9(i7,6"
	    "x))";
    static char fmt_130[] = "(\002------------------------------------------"
	    "----------------\002)";
    static char fmt_140[] = "(\002 Area     [L]      \002,e13.5,9e13.5)";
    static char fmt_150[] = "(\002 W-volume [L]      \002,e13.5,9e13.5)";
    static char fmt_151[] = "(\002 W-volumeI[L]      \002,e13.5,9e13.5)";
    static char fmt_160[] = "(\002 In-flow  [L/T]    \002,e13.5,9e13.5)";
    static char fmt_170[] = "(\002 h Mean   [L]      \002,e13.5,9e13.5)";
    static char fmt_180[] = "(\002 HeatVol  [M/T2]   \002,e13.5,10e13.5)";
    static char fmt_190[] = "(\002 tMean    [K]      \002,f13.3,10f13.3)";
    static char fmt_200[] = "(\002 ConcVol  [M/L2] \002,i1,1x,e13.5,10e13.5)";
    static char fmt_210[] = "(\002 cMean    [M/L3] \002,i1,1x,e13.5,10e13.5)";
    static char fmt_201[] = "(\002 ConcVolIm[M/L2] \002,i1,1x,e13.5,10e13.5)";
    static char fmt_211[] = "(\002 cMeanIm  [M/L3] \002,i1,1x,e13.5,10e13.5)";
    static char fmt_202[] = "(\002 SorbVolIm[M/L2] \002,i1,1x,e13.5,10e13.5)";
    static char fmt_212[] = "(\002 sMeanIm  [-]    \002,i1,1x,e13.5,10e13.5)";
    static char fmt_203[] = "(\002 SorbVolIm2[M/L2]\002,i1,1x,e13.5,10e13.5)";
    static char fmt_220[] = "(\002 Top Flux [L/T]    \002,e13.5/\002 Bot Flu"
	    "x [L/T]    \002,e13.5)";
    static char fmt_230[] = "(\002 WatBalT  [L]      \002,e13.5)";
    static char fmt_240[] = "(\002 WatBalR  [%]      \002,f13.3)";
    static char fmt_250[] = "(\002 CncBalT  [M]    \002,i1,1x,e13.5)";
    static char fmt_260[] = "(\002 CncBalR  [%]    \002,i1,1x,f13.3)";

    /* System generated locals */
    integer conc_dim1, conc_offset, chpar_dim1, chpar_offset, sorb_dim1, 
	    sorb_offset, sorb2_dim1, sorb2_offset, i__1, i__2, i__3;
    real r__1, r__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static real consubim[110]	/* was [11][10] */, convolim[11], volumeim;
    static integer i__, j;
    static real consubim2[110]	/* was [11][10] */, r__, convolim2[11], c1, 
	    c2, s1, s2, v1, cc, ce, he;
    static integer mi, mj;
    static real te, dx;
    static integer js;
    static real vn, tr, tt, ww, dx1, cel, fre;
    static integer jjj, lay;
    extern doublereal fro_(integer *, real *);
    static real dxn, tti, ttj, f_em__, ceim, thgi, thgj, grav, ctot[11], atot,
	     thwi, thwj, subt[10], htot, xksi, xksj, xnui, tvol, xnuj, ttot, 
	    cmean[110]	/* was [11][10] */, cbalr, cbalt, hmean[10], deltc, 
	    tmean[10], cnewi, wbalr, wbalt, deltw, voldi, tnewe, fexpi, vnewi,
	     fexpj, change, subcha[10], consub[110]	/* was [11][10] */, 
	    henryi, convol[11], ctotim[11], henryj, volume, cmeanim[110]	
	    /* was [11][10] */, thimobi, thimobj, cnewiim, voldimi, vnewimi, 
	    cnewiim2;

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 6, 0, 0, 0 };
    static cilist io___121 = { 1, 76, 0, fmt_110, 0 };
    static cilist io___122 = { 1, 76, 0, fmt_111, 0 };
    static cilist io___123 = { 1, 76, 0, fmt_120, 0 };
    static cilist io___124 = { 1, 76, 0, fmt_130, 0 };
    static cilist io___125 = { 1, 76, 0, fmt_140, 0 };
    static cilist io___126 = { 1, 76, 0, fmt_150, 0 };
    static cilist io___127 = { 1, 76, 0, fmt_151, 0 };
    static cilist io___128 = { 1, 76, 0, fmt_160, 0 };
    static cilist io___129 = { 1, 76, 0, fmt_170, 0 };
    static cilist io___130 = { 1, 76, 0, fmt_180, 0 };
    static cilist io___131 = { 1, 76, 0, fmt_190, 0 };
    static cilist io___132 = { 1, 76, 0, fmt_200, 0 };
    static cilist io___133 = { 1, 76, 0, fmt_210, 0 };
    static cilist io___134 = { 1, 76, 0, fmt_201, 0 };
    static cilist io___135 = { 1, 76, 0, fmt_211, 0 };
    static cilist io___136 = { 1, 76, 0, fmt_202, 0 };
    static cilist io___137 = { 1, 76, 0, fmt_212, 0 };
    static cilist io___138 = { 1, 76, 0, fmt_203, 0 };
    static cilist io___139 = { 1, 76, 0, fmt_220, 0 };
    static cilist io___141 = { 1, 76, 0, fmt_230, 0 };
    static cilist io___144 = { 1, 76, 0, fmt_240, 0 };
    static cilist io___146 = { 1, 76, 0, fmt_250, 0 };
    static cilist io___149 = { 1, 76, 0, fmt_260, 0 };
    static cilist io___150 = { 1, 76, 0, fmt_130, 0 };


    /* Parameter adjustments */
    --cprev;
    --tholdim;
    --thnewim;
    --convt;
    --convh;
    --conlt;
    --thvnew;
    --thvold;
    --temp;
    --solin;
    --watin;
    --con;
    --laynum;
    --matnum;
    --x;
    --tho;
    --thn;
    --hnew;
    --lmobim;
    tpar -= 11;
    --ths;
    --llinear;
    --cvoli;
    --ccuma;
    --ccumt;
    sorb2_dim1 = *nsd;
    sorb2_offset = 1 + sorb2_dim1;
    sorb2 -= sorb2_offset;
    sorb_dim1 = *nsd;
    sorb_offset = 1 + sorb_dim1;
    sorb -= sorb_offset;
    --tdep;
    chpar_dim1 = (*nsd << 4) + 4;
    chpar_offset = 1 + chpar_dim1;
    chpar -= chpar_offset;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;
    --subvol;
    --area;

    /* Function Body */
    fre = 1.f;
    grav = *cosalf;
    atot = 0.f;
    tr = 293.15f;
    r__ = 8.314f;
    if (*lwat || *plevel == 0) {
	volume = 0.f;
	volumeim = 0.f;
	change = 0.f;
	htot = 0.f;
	deltw = 0.f;
    }
    if (*ltemp) {
	ttot = 0.f;
	tvol = 0.f;
    }
    if (*lchem) {
	if (*ns > 11) {
	    s_wsle(&io___58);
	    do_lio(&c__9, &c__1, "Dimensions in the Subreg subroutine need t"
		    "o be increased !", (ftnlen)58);
	    e_wsle();
	    s_stop("", (ftnlen)0);
	}
	i__1 = *ns;
	for (js = 1; js <= i__1; ++js) {
	    ctot[js - 1] = 0.f;
	    convol[js - 1] = 0.f;
	    if (! (*lequil)) {
		convolim[js - 1] = 0.f;
	    }
	    if (! (*lequil)) {
		ctotim[js - 1] = 0.f;
	    }
	    if (*lbact || *ldualneq) {
		convolim2[js - 1] = 0.f;
	    }
/* L11: */
	}
	deltc = 0.f;
    }
    i__1 = *nlay;
    for (lay = 1; lay <= i__1; ++lay) {
	area[lay] = 0.f;
	if (*lwat || *plevel == 0) {
	    subvol[lay] = 0.f;
	    subcha[lay - 1] = 0.f;
	    hmean[lay - 1] = 0.f;
	}
	if (*ltemp) {
	    subt[lay - 1] = 0.f;
	    tmean[lay - 1] = 0.f;
	}
	if (*lchem) {
	    i__2 = *ns;
	    for (js = 1; js <= i__2; ++js) {
		consub[js + lay * 11 - 12] = 0.f;
		cmean[js + lay * 11 - 12] = 0.f;
		if (! (*lequil)) {
		    consubim[js + lay * 11 - 12] = 0.f;
		}
		if (*lbact || *ldualneq) {
		    consubim2[js + lay * 11 - 12] = 0.f;
		}
		if (! (*lequil)) {
		    cmeanim[js + lay * 11 - 12] = 0.f;
		}
/* L12: */
	    }
	}
/* L13: */
    }
    for (i__ = *n - 1; i__ >= 1; --i__) {
	j = i__ + 1;
	cel = 0.f;
	mi = matnum[i__];
	mj = matnum[j];
	lay = laynum[i__];
	dx = x[j] - x[i__];
	area[lay] += dx;
	atot += dx;
	tt = (temp[i__] + temp[j]) / 2.f + 273.15f;
	if (*lwat || *plevel == 0) {
	    he = (hnew[i__] + hnew[j]) / 2.f;
	    vnewi = dx * (thn[i__] + thn[j]) / 2.f;
	    voldi = dx * (tho[i__] + tho[j]) / 2.f;
	    if (*ldensity) {
		r__1 = conc[i__ * conc_dim1 + 1];
		r__2 = conc[j * conc_dim1 + 1];
		vnewi = dx * (thn[i__] * fro_(&c__1, &r__1) + thn[j] * fro_(&
			c__1, &r__2)) / 2.f;
		r__1 = cprev[i__];
		r__2 = cprev[j];
		voldi = dx * (tho[i__] * fro_(&c__1, &r__1) + tho[j] * fro_(&
			c__1, &r__2)) / 2.f;
	    }
	    if (*lvapor) {
		vnewi += dx * (thvnew[i__] + thvnew[j]) / 2.f;
		voldi += dx * (thvold[i__] + thvold[j]) / 2.f;
	    }
	    volume += vnewi;
	    change += (vnewi - voldi) / *dt;
	    subcha[lay - 1] += (vnewi - voldi) / *dt;
	    subvol[lay] += vnewi;
	    htot += he * dx;
	    hmean[lay - 1] += he * dx;
	    if (*idualpor > 0) {
		vnewimi = dx * (thnewim[i__] + thnewim[j]) / 2.f;
		voldimi = dx * (tholdim[i__] + tholdim[j]) / 2.f;
		if (*ldensity) {
		    r__1 = sorb[i__ * sorb_dim1 + 1];
		    r__2 = sorb[j * sorb_dim1 + 1];
		    vnewimi = dx * (thnewim[i__] * fro_(&c__1, &r__1) + 
			    thnewim[j] * fro_(&c__1, &r__2)) / 2.f;
		    r__1 = sorb[i__ * sorb_dim1 + 1];
		    r__2 = sorb[j * sorb_dim1 + 1];
		    voldimi = dx * (tholdim[i__] * fro_(&c__1, &r__1) + 
			    tholdim[j] * fro_(&c__1, &r__2)) / 2.f;
		}
		volumeim += dx * (thnewim[i__] + thnewim[j]) / 2.f;
		change += (vnewimi - voldimi) / *dt;
		subcha[lay - 1] += (vnewimi - voldimi) / *dt;
		subvol[lay] += vnewimi;
	    }
	}
	if (*ltemp) {
	    te = (temp[i__] + temp[j]) / 2.f;
	    tnewe = dx * ((temp[i__] + 273.15f) * (tpar[mi * 10 + 1] * tpar[
		    mi * 10 + 7] + tpar[mi * 10 + 2] * tpar[mi * 10 + 8] + 
		    tpar[mi * 10 + 9] * thn[i__]) + (temp[j] + 273.15f) * (
		    tpar[mj * 10 + 1] * tpar[mj * 10 + 7] + tpar[mj * 10 + 2] 
		    * tpar[mj * 10 + 8] + tpar[mj * 10 + 9] * thn[j])) / 2.f;
	    tvol += tnewe;
	    subt[lay - 1] += tnewe;
	    ttot += te * dx;
	    tmean[lay - 1] += te * dx;
	}
	if (*lchem) {
	    i__1 = *ns;
	    for (js = 1; js <= i__1; ++js) {
		jjj = js - 1 << 4;
		ce = (conc[js + i__ * conc_dim1] + conc[js + j * conc_dim1]) /
			 2.f;
		tti = (temp[i__] + 273.15f - tr) / r__ / tt / tr;
		xksi = chpar[jjj + 7 + mi * chpar_dim1] * exp(tdep[jjj + 7] * 
			tti);
		xnui = chpar[jjj + 8 + mi * chpar_dim1] * exp(tdep[jjj + 8] * 
			tti);
		fexpi = chpar[jjj + 9 + mi * chpar_dim1] * exp(tdep[jjj + 9] *
			 tti);
		henryi = chpar[jjj + 10 + mi * chpar_dim1] * exp(tdep[jjj + 
			10] * tti);
		ttj = (temp[j] + 273.15f - tr) / r__ / tt / tr;
		xksj = chpar[jjj + 7 + mj * chpar_dim1] * exp(tdep[jjj + 7] * 
			ttj);
		xnuj = chpar[jjj + 8 + mj * chpar_dim1] * exp(tdep[jjj + 8] * 
			ttj);
		fexpj = chpar[jjj + 9 + mj * chpar_dim1] * exp(tdep[jjj + 9] *
			 ttj);
		henryj = chpar[jjj + 10 + mj * chpar_dim1] * exp(tdep[jjj + 
			10] * ttj);
		c1 = 1.f;
		c2 = 1.f;
		if (! llinear[js]) {
		    if (conc[js + i__ * conc_dim1] > 0.f) {
			d__1 = (doublereal) conc[js + i__ * conc_dim1];
			d__2 = (doublereal) (fexpi - 1.f);
			d__3 = (doublereal) conc[js + i__ * conc_dim1];
			d__4 = (doublereal) fexpi;
			c1 = pow_dd(&d__1, &d__2) / (xnui * pow_dd(&d__3, &
				d__4) + 1.f);
		    }
		    if (conc[js + j * conc_dim1] > 0.f) {
			d__1 = (doublereal) conc[js + j * conc_dim1];
			d__2 = (doublereal) (fexpj - 1.f);
			d__3 = (doublereal) conc[js + j * conc_dim1];
			d__4 = (doublereal) fexpj;
			c2 = pow_dd(&d__1, &d__2) / (xnuj * pow_dd(&d__3, &
				d__4) + 1.f);
		    }
		}
		thwi = thn[i__];
		thwj = thn[j];
		thimobi = chpar[mi * chpar_dim1 + 4];
		thimobj = chpar[mj * chpar_dim1 + 4];
/* Computing MAX */
		r__1 = 0.f, r__2 = ths[mi] - thwi;
		thgi = dmax(r__1,r__2);
/* Computing MAX */
		r__1 = 0.f, r__2 = ths[mj] - thwj;
		thgj = dmax(r__1,r__2);
		if (*idualpor > 0) {
		    thimobi = thnewim[i__];
		    thimobj = thnewim[j];
/*              ThGi=amax1(0.,ths(Mi)-ThWi+thSIm(Mi)-ThImobi) */
/*              ThGj=amax1(0.,ths(Mj)-ThWj+thSIm(Mj)-ThImobj) */
		}
		if (lmobim[mi] && *idualpor == 0 || *lbact) {
/* Computing MAX */
		    r__1 = thwi - thimobi;
		    thwi = dmax(r__1,.001f);
		}
		if (lmobim[mj] && *idualpor == 0 || *lbact) {
/* Computing MAX */
		    r__1 = thwj - thimobj;
		    thwj = dmax(r__1,.001f);
		}
		f_em__ = 1.f;
		if (*ldualneq) {
		    f_em__ = chpar[jjj + 13 + mi * chpar_dim1];
		}
		cnewi = dx / 2.f * (conc[js + i__ * conc_dim1] * (thwi + 
			f_em__ * chpar[mi * chpar_dim1 + 3] * chpar[mi * 
			chpar_dim1 + 1] * xksi * c1 + thgi * henryi) + conc[
			js + j * conc_dim1] * (thwj + f_em__ * chpar[mj * 
			chpar_dim1 + 3] * chpar[mj * chpar_dim1 + 1] * xksj * 
			c2 + thgj * henryj));
		convol[js - 1] += cnewi;
		consub[js + lay * 11 - 12] += cnewi;
		if (! (*lequil)) {
		    if (lmobim[mi] || *idualpor > 0) {
/* mobile-immobile mod */
			s1 = 1.f;
			s2 = 1.f;
			if (! llinear[js]) {
			    if (sorb[js + i__ * sorb_dim1] > 0.f) {
				d__1 = (doublereal) sorb[js + i__ * sorb_dim1]
					;
				d__2 = (doublereal) (fexpi - 1.f);
				d__3 = (doublereal) sorb[js + i__ * sorb_dim1]
					;
				d__4 = (doublereal) fexpi;
				s1 = pow_dd(&d__1, &d__2) / (xnui * pow_dd(&
					d__3, &d__4) + 1.f);
			    }
			    if (sorb[js + j * sorb_dim1] > 0.f) {
				d__1 = (doublereal) sorb[js + j * sorb_dim1];
				d__2 = (doublereal) (fexpj - 1.f);
				d__3 = (doublereal) sorb[js + j * sorb_dim1];
				d__4 = (doublereal) fexpj;
				s2 = pow_dd(&d__1, &d__2) / (xnuj * pow_dd(&
					d__3, &d__4) + 1.f);
			    }
			}
			cnewiim = dx / 2.f * (sorb[js + i__ * sorb_dim1] * (
				thimobi + (1.f - chpar[mi * chpar_dim1 + 3]) *
				 chpar[mi * chpar_dim1 + 1] * xksi * s1) + 
				sorb[js + j * sorb_dim1] * (thimobj + (1.f - 
				chpar[mj * chpar_dim1 + 3]) * chpar[mj * 
				chpar_dim1 + 1] * xksj * s2));
			if (*ldualneq) {
			    cnewiim2 = dx / 2.f * (chpar[mi * chpar_dim1 + 1] 
				    * sorb2[js + i__ * sorb2_dim1] + chpar[mj 
				    * chpar_dim1 + 1] * sorb2[js + j * 
				    sorb2_dim1]);
			    convolim2[js - 1] += cnewiim2;
			    consubim2[js + lay * 11 - 12] += cnewiim2;
			}
		    } else {
/* two-site sorption m */
			cnewiim = dx / 2.f * (chpar[mi * chpar_dim1 + 1] * 
				sorb[js + i__ * sorb_dim1] + chpar[mj * 
				chpar_dim1 + 1] * sorb[js + j * sorb_dim1]);
			if (*lbact) {
			    cnewiim2 = dx / 2.f * (chpar[mi * chpar_dim1 + 1] 
				    * sorb2[js + i__ * sorb2_dim1] + chpar[mj 
				    * chpar_dim1 + 1] * sorb2[js + j * 
				    sorb2_dim1]);
			    convolim2[js - 1] += cnewiim2;
			    consubim2[js + lay * 11 - 12] += cnewiim2;
			}
		    }
		    convolim[js - 1] += cnewiim;
		    consubim[js + lay * 11 - 12] += cnewiim;
		    ceim = (sorb[js + i__ * sorb_dim1] + sorb[js + j * 
			    sorb_dim1]) / 2.f;
		    cmeanim[js + lay * 11 - 12] += ceim * dx;
		    ctotim[js - 1] += ceim * dx;
		}
		ctot[js - 1] += ce * dx;
		cmean[js + lay * 11 - 12] += ce * dx;
		if (js == 1) {
		    cel = cnewi;
		}
/* L14: */
	    }
	}
	if (*plevel == 0) {
	    if (*lwat) {
		watin[i__] = vnewi;
	    }
	    if (*lwat && *idualpor > 0) {
		watin[i__] = vnewi + vnewimi;
	    }
	    if (*lchem) {
		solin[i__] = cel;
	    }
	} else {
	    if (*lwat) {
		if (*idualpor <= 0) {
		    deltw += (r__1 = watin[i__] - vnewi, dabs(r__1));
		}
		if (*idualpor > 0) {
		    deltw += (r__1 = watin[i__] - vnewi - vnewimi, dabs(r__1))
			    ;
		}
	    }
	    if (*lchem) {
		deltc += (r__1 = solin[i__] - cel, dabs(r__1));
	    }
	}
/* L15: */
    }
    i__1 = *nlay;
    for (lay = 1; lay <= i__1; ++lay) {
	if (area[lay] > 0.f) {
	    if (*lwat || *plevel == 0) {
		hmean[lay - 1] /= area[lay];
	    }
	    if (*ltemp) {
		tmean[lay - 1] /= area[lay];
	    }
	    i__2 = *ns;
	    for (js = 1; js <= i__2; ++js) {
		if (*lchem) {
		    cmean[js + lay * 11 - 12] /= area[lay];
		    if (! (*lequil)) {
			cmeanim[js + lay * 11 - 12] /= area[lay];
		    }
		}
/* L16: */
	    }
	}
/* L17: */
    }
    if (atot > 0.f) {
	if (*lwat || *plevel == 0) {
	    htot /= atot;
	}
	if (*ltemp) {
	    ttot /= atot;
	}
    }
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
	if (*lchem && atot > 0.f) {
	    ctot[js - 1] /= atot;
	    if (! (*lequil)) {
		ctotim[js - 1] /= atot;
	    }
	}
/* L18: */
    }
    if (*ldensity) {
	fre = fro_(&c__1, &conc[conc_dim1 + 1]);
    }
    if (*lcentrif) {
	grav = *cosalf * (*radius + (r__1 = (x[2] + x[1]) / 2.f, dabs(r__1)));
    }
    dx1 = x[2] - x[1];
    v1 = -(con[1] + con[2]) / 2.f * ((hnew[2] - hnew[1]) / dx1 + grav * fre);
    if (*ldensity) {
	fre = fro_(&c__1, &conc[*n * conc_dim1 + 1]);
    }
    if (*lcentrif) {
	grav = *cosalf * (*radius + (r__1 = (x[*n] + x[*n - 1]) / 2.f, dabs(
		r__1)));
    }
    dxn = x[*n] - x[*n - 1];
    vn = -(con[*n] + con[*n - 1]) / 2 * ((hnew[*n] - hnew[*n - 1]) / dxn + 
	    grav * fre);
    if (*lwtdep) {
	v1 -= (conlt[1] + conlt[2]) / 2.f * (temp[2] - temp[1]) / dx1;
	vn -= (conlt[*n] + conlt[*n - 1]) / 2.f * (temp[*n] - temp[*n - 1]) / 
		dxn;
    }
    if (*lvapor) {
	v1 -= (convh[1] + convh[2]) / 2.f * (hnew[2] - hnew[1]) / dx1;
	vn -= (convh[*n] + convh[*n - 1]) / 2.f * (hnew[*n] - hnew[*n - 1]) / 
		dxn;
	v1 -= (convt[1] + convt[2]) / 2.f * (temp[2] - temp[1]) / dx1;
	vn -= (convt[*n] + convt[*n - 1]) / 2.f * (temp[*n] - temp[*n - 1]) / 
		dxn;
    }
    if (*lprint) {
	if (*t < 99999999.f) {
	    i__1 = s_wsfe(&io___121);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	} else {
	    i__1 = s_wsfe(&io___122);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = s_wsfe(&io___123);
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *nlay;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_wsfe(&io___124);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = s_wsfe(&io___125);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&atot, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L901;
	}
	i__2 = *nlay;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = do_fio(&c__1, (char *)&area[i__], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	if (*lwat || *plevel == 0) {
	    i__1 = s_wsfe(&io___126);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&volume, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__2 = *nlay;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = do_fio(&c__1, (char *)&subvol[i__], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	    if (*idualpor > 0) {
		i__1 = s_wsfe(&io___127);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&volumeim, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	    i__1 = s_wsfe(&io___128);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&change, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__2 = *nlay;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = do_fio(&c__1, (char *)&subcha[i__ - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_wsfe(&io___129);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&htot, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__2 = *nlay;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = do_fio(&c__1, (char *)&hmean[i__ - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	if (*ltemp) {
	    i__1 = s_wsfe(&io___130);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&tvol, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__2 = *nlay;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = do_fio(&c__1, (char *)&subt[i__ - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = s_wsfe(&io___131);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&ttot, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__2 = *nlay;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = do_fio(&c__1, (char *)&tmean[i__ - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	if (*lchem) {
	    i__1 = *ns;
	    for (js = 1; js <= i__1; ++js) {
		i__2 = s_wsfe(&io___132);
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(integer));
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = do_fio(&c__1, (char *)&convol[js - 1], (ftnlen)sizeof(
			real));
		if (i__2 != 0) {
		    goto L901;
		}
		i__3 = *nlay;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__2 = do_fio(&c__1, (char *)&consub[js + i__ * 11 - 12], 
			    (ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L901;
		    }
		}
		i__2 = e_wsfe();
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = s_wsfe(&io___133);
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(integer));
		if (i__2 != 0) {
		    goto L901;
		}
		i__2 = do_fio(&c__1, (char *)&ctot[js - 1], (ftnlen)sizeof(
			real));
		if (i__2 != 0) {
		    goto L901;
		}
		i__3 = *nlay;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__2 = do_fio(&c__1, (char *)&cmean[js + i__ * 11 - 12], (
			    ftnlen)sizeof(real));
		    if (i__2 != 0) {
			goto L901;
		    }
		}
		i__2 = e_wsfe();
		if (i__2 != 0) {
		    goto L901;
		}
		if (! (*lequil)) {
		    if (lmobim[1] || *idualpor > 0) {
			i__2 = s_wsfe(&io___134);
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(
				integer));
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&convolim[js - 1], (
				ftnlen)sizeof(real));
			if (i__2 != 0) {
			    goto L901;
			}
			i__3 = *nlay;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    i__2 = do_fio(&c__1, (char *)&consubim[js + i__ * 
				    11 - 12], (ftnlen)sizeof(real));
			    if (i__2 != 0) {
				goto L901;
			    }
			}
			i__2 = e_wsfe();
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = s_wsfe(&io___135);
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(
				integer));
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&ctotim[js - 1], (ftnlen)
				sizeof(real));
			if (i__2 != 0) {
			    goto L901;
			}
			i__3 = *nlay;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    i__2 = do_fio(&c__1, (char *)&cmeanim[js + i__ * 
				    11 - 12], (ftnlen)sizeof(real));
			    if (i__2 != 0) {
				goto L901;
			    }
			}
			i__2 = e_wsfe();
			if (i__2 != 0) {
			    goto L901;
			}
		    } else {
			i__2 = s_wsfe(&io___136);
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(
				integer));
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&convolim[js - 1], (
				ftnlen)sizeof(real));
			if (i__2 != 0) {
			    goto L901;
			}
			i__3 = *nlay;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    i__2 = do_fio(&c__1, (char *)&consubim[js + i__ * 
				    11 - 12], (ftnlen)sizeof(real));
			    if (i__2 != 0) {
				goto L901;
			    }
			}
			i__2 = e_wsfe();
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = s_wsfe(&io___137);
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(
				integer));
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&ctotim[js - 1], (ftnlen)
				sizeof(real));
			if (i__2 != 0) {
			    goto L901;
			}
			i__3 = *nlay;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    i__2 = do_fio(&c__1, (char *)&cmeanim[js + i__ * 
				    11 - 12], (ftnlen)sizeof(real));
			    if (i__2 != 0) {
				goto L901;
			    }
			}
			i__2 = e_wsfe();
			if (i__2 != 0) {
			    goto L901;
			}
		    }
		}
		if (*lbact || *ldualneq) {
		    i__2 = s_wsfe(&io___138);
		    if (i__2 != 0) {
			goto L901;
		    }
		    i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(integer))
			    ;
		    if (i__2 != 0) {
			goto L901;
		    }
		    i__2 = do_fio(&c__1, (char *)&convolim2[js - 1], (ftnlen)
			    sizeof(real));
		    if (i__2 != 0) {
			goto L901;
		    }
		    i__3 = *nlay;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__2 = do_fio(&c__1, (char *)&consubim2[js + i__ * 11 
				- 12], (ftnlen)sizeof(real));
			if (i__2 != 0) {
			    goto L901;
			}
		    }
		    i__2 = e_wsfe();
		    if (i__2 != 0) {
			goto L901;
		    }
		}
/* L19: */
	    }
	}
	if (*lwat || *plevel == 0) {
	    i__1 = s_wsfe(&io___139);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&vn, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
    }
/*     Mass balance calculation */
    if (*plevel == 0) {
	*wvoli = volume;
	if (*idualpor > 0) {
	    *wvoli = volume + volumeim;
	}
	if (*lchem) {
	    i__1 = *ns;
	    for (js = 1; js <= i__1; ++js) {
		cvoli[js] = convol[js - 1];
		if (! (*lequil)) {
		    cvoli[js] += convolim[js - 1];
		}
		if (*lbact || *ldualneq) {
		    cvoli[js] += convolim2[js - 1];
		}
/* L20: */
	    }
	}
    } else {
	if (*lwat) {
	    wbalt = volume - *wvoli - *wcumt;
	    if (*idualpor > 0) {
		wbalt = volume - *wvoli - *wcumt + volumeim;
	    }
	    if (*lprint) {
		i__1 = s_wsfe(&io___141);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&wbalt, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	    ww = dmax(deltw,*wcuma);
	    if (ww > 1e-25f) {
		wbalr = dabs(wbalt) / ww * 100.f;
		if (*lprint) {
		    i__1 = s_wsfe(&io___144);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&wbalr, (ftnlen)sizeof(real))
			    ;
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
	    }
	}
	if (*lchem) {
	    i__1 = *ns;
	    for (js = 1; js <= i__1; ++js) {
		cbalt = convol[js - 1] - cvoli[js] + ccumt[js];
		if (! (*lequil)) {
		    cbalt += convolim[js - 1];
		}
		if (*lbact || *ldualneq) {
		    cbalt += convolim2[js - 1];
		}
		if (*lprint) {
		    i__2 = s_wsfe(&io___146);
		    if (i__2 != 0) {
			goto L901;
		    }
		    i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(integer))
			    ;
		    if (i__2 != 0) {
			goto L901;
		    }
		    i__2 = do_fio(&c__1, (char *)&cbalt, (ftnlen)sizeof(real))
			    ;
		    if (i__2 != 0) {
			goto L901;
		    }
		    i__2 = e_wsfe();
		    if (i__2 != 0) {
			goto L901;
		    }
		}
/* Computing MAX */
		r__1 = deltc, r__2 = ccuma[js];
		cc = dmax(r__1,r__2);
		if (cc > 1e-25f) {
		    cbalr = dabs(cbalt) / cc * 100.f;
		    if (*lprint) {
			i__2 = s_wsfe(&io___149);
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&js, (ftnlen)sizeof(
				integer));
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = do_fio(&c__1, (char *)&cbalr, (ftnlen)sizeof(
				real));
			if (i__2 != 0) {
			    goto L901;
			}
			i__2 = e_wsfe();
			if (i__2 != 0) {
			    goto L901;
			}
		    }
		}
/* L21: */
	    }
	}
    }
    if (*lprint) {
	i__1 = s_wsfe(&io___150);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    return 0;
/*     Error when writing into an output file */
L901:
    *ierr = 1;
    return 0;
} /* subreg_ */

/* ********************************************************************** */
/* Subroutine */ int nodout_(integer *n, integer *nmat, real *hnew, real *thn,
	 real *con, real *x, real *xsurf, real *cosalf, doublereal *tprint, 
	integer *matnum, real *cap, real *bxz, real *sink, real *cons, 
	integer *ns, integer *nsd, real *conc, real *temp, real *sorb, 
	integer *kappa, logical *lbact, real *sorb2, logical *lvapor, logical 
	*lwtdep, real *conlt, real *convt, real *convh, real *tholdt, real *
	dt, integer *idualpor, real *thnewim, real *sinkim, real *strans, 
	logical *ldensity, logical *lcentrif, real *radius, logical *
	lvaporout, logical *ldualneq, integer *ierr)
{
    /* Format strings */
    static char fmt_110[] = "(//\002 Time:\002,f14.4//)";
    static char fmt_111[] = "(//\002 Time:\002,e15.8//)";
    static char fmt_112[] = "(\002 Node      Depth      Head Moisture       "
	    "K          C         \002,\002Flux        Sink         Kappa   v"
	    "/KsTop   Temp\002/\002           [L]        [L]    [-]        [L"
	    "/T]      [1/L]      [\002,\002L/T]        [1/T]         [-]     "
	    " [-]      [C]\002/)";
    static char fmt_113[] = "(\002 Node      Depth      Head Moisture    WTr"
	    "ans    Im.Moist.     \002,\002Flux       STrans        Kappa   v"
	    "/KsTop   Temp\002/\002           [L]        [L]    [-]        [1"
	    "/T]       [-]       [\002,\002L/T]     [M/L*3/T]        [-]     "
	    " [-]      [C]\002/)";
    static char fmt_114[] = "(\002 Node      Depth      Head Moisture       "
	    "K          C         \002,\002Flux        Sink         Kappa   v"
	    "/KsTop   Temp   Conc(1..NS) Sorb(1...NS)\002/\002           [L] "
	    "       [L]    [-]        [L/T]      [1/L]      [\002,\002L/T]   "
	    "     [1/T]         [-]      [-]      [C]      [M/L*3]\002/)";
    static char fmt_115[] = "(\002 Node      Depth      Head Moisture    WTr"
	    "ans    Im.Moist.     \002,\002Flux       STrans        Kappa   v"
	    "/KsTop   Temp   Conc(1..NS) Sorb(1...NS)\002/\002           [L] "
	    "       [L]    [-]        [1/T]       [-]       [\002,\002L/T]   "
	    "  [M/L*3/T]        [-]      [-]      [C]      [M/L*3]\002/)";
    static char fmt_116[] = "(\002 Node    Depth        Con        ConLT    "
	    "    ConVh        ConVT       vLiquid       vVapor       vTotal  "
	    "     vVapIso     vVapTerm\002/\002          [L]        [L/T]    "
	    " [L2/K/T]       [L/T]      [L2/K/T]       [L/T]         [L/T]   "
	    "     [L/T]        [L/T]        [L/T]\002/)";
    static char fmt_140[] = "(i4,1x,f10.4,1x,9e13.5)";
    static char fmt_120[] = "(i4,1x,f10.4,1x,f11.3,1x,f6.4,1x,4e12.4,i8,1x,e"
	    "11.3,f8.2,30e12.4)";
    static char fmt_130[] = "(i4,1x,f10.4,1x,e11.4,1x,f6.4,1x,4e12.4,i8,2x,f"
	    "10.3,f8.2,30e12.4)";

    /* System generated locals */
    integer conc_dim1, conc_offset, sorb_dim1, sorb_offset, sorb2_dim1, 
	    sorb2_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__, n1;
    static real va, vb, dx;
    static integer js;
    static real vi, xs, fre, dxa, dxb;
    static integer mat;
    extern doublereal fro_(integer *, real *);
    static real vat, vbt, vva, vvb, vvi, grav, vvii, vvat, vvbt, vvti, consn;

    /* Fortran I/O blocks */
    static cilist io___154 = { 1, 75, 0, fmt_110, 0 };
    static cilist io___155 = { 1, 45, 0, fmt_110, 0 };
    static cilist io___156 = { 1, 75, 0, fmt_111, 0 };
    static cilist io___157 = { 1, 45, 0, fmt_111, 0 };
    static cilist io___158 = { 1, 75, 0, fmt_112, 0 };
    static cilist io___159 = { 1, 75, 0, fmt_113, 0 };
    static cilist io___160 = { 1, 75, 0, fmt_114, 0 };
    static cilist io___161 = { 1, 75, 0, fmt_115, 0 };
    static cilist io___162 = { 1, 45, 0, fmt_116, 0 };
    static cilist io___182 = { 1, 45, 0, fmt_140, 0 };
    static cilist io___183 = { 1, 75, 0, fmt_120, 0 };
    static cilist io___185 = { 1, 75, 0, fmt_120, 0 };
    static cilist io___186 = { 1, 75, 0, fmt_120, 0 };
    static cilist io___187 = { 1, 75, 0, fmt_120, 0 };
    static cilist io___188 = { 1, 75, 0, fmt_120, 0 };
    static cilist io___189 = { 1, 75, 0, fmt_130, 0 };
    static cilist io___190 = { 1, 75, 0, fmt_130, 0 };
    static cilist io___191 = { 1, 75, 0, fmt_130, 0 };
    static cilist io___192 = { 1, 75, 0, fmt_130, 0 };
    static cilist io___193 = { 1, 75, 0, fmt_130, 0 };
    static cilist io___194 = { 1, 75, 0, "('end')", 0 };
    static cilist io___195 = { 1, 45, 0, "('end')", 0 };


    /* Parameter adjustments */
    --strans;
    --sinkim;
    --thnewim;
    --convh;
    --convt;
    --conlt;
    --kappa;
    --temp;
    --sink;
    --bxz;
    --cap;
    --matnum;
    --x;
    --con;
    --thn;
    --hnew;
    --cons;
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
    fre = 1.f;
    grav = *cosalf;
    xs = *xsurf + *radius;
    if (*tprint < 99999999.f) {
	i__1 = s_wsfe(&io___154);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*tprint), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	if (*lvaporout && *lvapor) {
	    i__1 = s_wsfe(&io___155);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*tprint), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
    } else {
	i__1 = s_wsfe(&io___156);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = do_fio(&c__1, (char *)&(*tprint), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
	if (*lvaporout && *lvapor) {
	    i__1 = s_wsfe(&io___157);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&(*tprint), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
    }
    if (*ns == 0) {
	if (*idualpor == 0) {
	    i__1 = s_wsfe(&io___158);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	} else {
	    i__1 = s_wsfe(&io___159);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
    } else {
	if (*idualpor == 0) {
	    i__1 = s_wsfe(&io___160);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	} else {
	    i__1 = s_wsfe(&io___161);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
    }
    if (*lvaporout && *lvapor) {
	i__1 = s_wsfe(&io___162);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    for (i__ = *n; i__ >= 1; --i__) {
	mat = matnum[i__];
	if (i__ == 1) {
	    if (*ldensity) {
		fre = fro_(&c__1, &conc[conc_dim1 + 1]);
	    }
	    if (*lcentrif) {
		grav = *cosalf * (*radius + (r__1 = (x[2] + x[1]) / 2.f, dabs(
			r__1)));
	    }
	    dx = x[2] - x[1];
	    vi = -(con[1] + con[2]) / 2.f * ((hnew[2] - hnew[1]) / dx + grav *
		     fre);
	    if (*lwtdep) {
		vi -= (conlt[1] + conlt[2]) / 2.f * (temp[2] - temp[1]) / dx;
	    }
	    vvi = 0.f;
	    if (*lvapor) {
		vvii = -(convh[1] + convh[2]) / 2.f * (hnew[2] - hnew[1]) / 
			dx;
		vvti = -(convt[1] + convt[2]) / 2.f * (temp[2] - temp[1]) / 
			dx;
		vvi = vvii + vvti;
	    }
	} else if (i__ == *n) {
	    consn = cons[mat] * bxz[*n];
	    n1 = *n - 1;
	    dx = x[*n] - x[*n - 1];
	    if (*ldensity) {
		fre = fro_(&c__1, &conc[*n * conc_dim1 + 1]);
	    }
	    if (*lcentrif) {
		grav = *cosalf * (*radius + (r__1 = (x[*n] + x[n1]) / 2.f, 
			dabs(r__1)));
	    }
	    vi = -(con[*n] + con[n1]) / 2.f * ((hnew[*n] - hnew[n1]) / dx + 
		    grav * fre) - (thn[*n] - *tholdt) * fre * dx / 2.f / *dt 
		    - sink[*n] * dx / 2.f;
	    if (*lwtdep) {
		vi -= (conlt[*n] + conlt[n1]) / 2.f * (temp[*n] - temp[n1]) / 
			dx;
	    }
	    vvi = 0.f;
	    if (*lvapor) {
		vvii = -(convh[*n] + convh[n1]) / 2.f * (hnew[*n] - hnew[n1]) 
			/ dx;
		vvti = -(convt[*n] + convt[n1]) / 2.f * (temp[*n] - temp[n1]) 
			/ dx;
		vvi = vvti + vvii;
	    }
	} else {
	    dxa = x[i__ + 1] - x[i__];
	    dxb = x[i__] - x[i__ - 1];
	    if (*ldensity) {
		fre = (fro_(&c__1, &conc[i__ * conc_dim1 + 1]) + fro_(&c__1, &
			conc[(i__ + 1) * conc_dim1 + 1])) / 2.f;
	    }
	    if (*lcentrif) {
		grav = *cosalf * (*radius + (r__1 = (x[i__ + 1] + x[i__]) / 
			2.f, dabs(r__1)));
	    }
	    va = -(con[i__] + con[i__ + 1]) / 2.f * ((hnew[i__ + 1] - hnew[
		    i__]) / dxa + grav * fre);
	    if (*ldensity) {
		fre = (fro_(&c__1, &conc[i__ * conc_dim1 + 1]) + fro_(&c__1, &
			conc[(i__ - 1) * conc_dim1 + 1])) / 2.f;
	    }
	    if (*lcentrif) {
		grav = *cosalf * (*radius + (r__1 = (x[i__] + x[i__ - 1]) / 
			2.f, dabs(r__1)));
	    }
	    vb = -(con[i__] + con[i__ - 1]) / 2.f * ((hnew[i__] - hnew[i__ - 
		    1]) / dxb + grav * fre);
	    vi = (va * dxa + vb * dxb) / (dxa + dxb);
	    if (*lwtdep) {
		vat = -(conlt[i__] + conlt[i__ + 1]) / 2.f * (temp[i__ + 1] - 
			temp[i__]) / dxa;
		vbt = -(conlt[i__] + conlt[i__ - 1]) / 2.f * (temp[i__] - 
			temp[i__ - 1]) / dxb;
		vi += (vat * dxa + vbt * dxb) / (dxa + dxb);
	    }
	    vvi = 0.f;
	    if (*lvapor) {
		vva = -(convh[i__] + convh[i__ + 1]) / 2.f * (hnew[i__ + 1] - 
			hnew[i__]) / dxa;
		vvb = -(convh[i__] + convh[i__ - 1]) / 2.f * (hnew[i__] - 
			hnew[i__ - 1]) / dxb;
		vvii = (vva * dxa + vvb * dxb) / (dxa + dxb);
		vvat = -(convt[i__] + convt[i__ + 1]) / 2.f * (temp[i__ + 1] 
			- temp[i__]) / dxa;
		vvbt = -(convt[i__] + convt[i__ - 1]) / 2.f * (temp[i__] - 
			temp[i__ - 1]) / dxb;
		vvti = (vvat * dxa + vvbt * dxb) / (dxa + dxb);
		vvi = vvii + vvti;
	    }
	}
	if (*lvaporout && *lvapor) {
	    i__1 = s_wsfe(&io___182);
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__2 = *n - i__ + 1;
	    i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L901;
	    }
	    r__1 = x[i__] - xs;
	    i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&con[i__], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&conlt[i__], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&convh[i__], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&convt[i__], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&vi, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&vvi, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    r__2 = vi + vvi;
	    i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&vvii, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = do_fio(&c__1, (char *)&vvti, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L901;
	    }
	    i__1 = e_wsfe();
	    if (i__1 != 0) {
		goto L901;
	    }
	}
	if (hnew[i__] > -9.9e5f) {
	    if (! (*lbact)) {
		if (*idualpor == 0) {
		    if (! (*ldualneq)) {
			i__1 = s_wsfe(&io___183);
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *n - i__ + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__1 = x[i__] - xs;
			i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&con[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&cap[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			r__2 = vi + vvi;
			i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&sink[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)
				sizeof(integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__3 = vi / consn;
			i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (js = 1; js <= i__3; ++js) {
			    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * 
				    conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__4 = *ns;
			for (js = 1; js <= i__4; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * 
				    sorb_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__1 = e_wsfe();
			if (i__1 != 0) {
			    goto L901;
			}
		    } else {
			i__1 = s_wsfe(&io___185);
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *n - i__ + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__1 = x[i__] - xs;
			i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&con[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&cap[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			r__2 = vi + vvi;
			i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&sink[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)
				sizeof(integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__3 = vi / consn;
			i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (js = 1; js <= i__3; ++js) {
			    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * 
				    conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__4 = *ns;
			for (js = 1; js <= i__4; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * 
				    sorb_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__5 = *ns;
			for (js = 1; js <= i__5; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb2[js + i__ * 
				    sorb2_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__1 = e_wsfe();
			if (i__1 != 0) {
			    goto L901;
			}
		    }
		} else {
		    if (! (*ldualneq)) {
			i__1 = s_wsfe(&io___186);
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *n - i__ + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__1 = x[i__] - xs;
			i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&sinkim[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thnewim[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			r__2 = vi + vvi;
			i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&strans[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)
				sizeof(integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__3 = vi / consn;
			i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (js = 1; js <= i__3; ++js) {
			    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * 
				    conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__4 = *ns;
			for (js = 1; js <= i__4; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * 
				    sorb_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__1 = e_wsfe();
			if (i__1 != 0) {
			    goto L901;
			}
		    } else {
			i__1 = s_wsfe(&io___187);
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *n - i__ + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__1 = x[i__] - xs;
			i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&sinkim[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thnewim[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			r__2 = vi + vvi;
			i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&strans[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)
				sizeof(integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__3 = vi / consn;
			i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (js = 1; js <= i__3; ++js) {
			    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * 
				    conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__4 = *ns;
			for (js = 1; js <= i__4; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * 
				    sorb_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__5 = *ns;
			for (js = 1; js <= i__5; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb2[js + i__ * 
				    sorb2_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__1 = e_wsfe();
			if (i__1 != 0) {
			    goto L901;
			}
		    }
		}
	    } else {
		i__1 = s_wsfe(&io___188);
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *n - i__ + 1;
		i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		if (i__1 != 0) {
		    goto L901;
		}
		r__1 = x[i__] - xs;
		i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&con[i__], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cap[i__], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		r__2 = vi + vvi;
		i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&sink[i__], (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L901;
		}
		r__3 = vi / consn;
		i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__3 = *ns;
		for (js = 1; js <= i__3; ++js) {
		    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * conc_dim1], 
			    (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__4 = *ns;
		for (js = 1; js <= i__4; ++js) {
		    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * sorb_dim1], 
			    (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__5 = *ns;
		for (js = 1; js <= i__5; ++js) {
		    i__1 = do_fio(&c__1, (char *)&sorb2[js + i__ * sorb2_dim1]
			    , (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	} else {
	    if (! (*lbact)) {
		if (*idualpor == 0) {
		    if (! (*ldualneq)) {
			i__1 = s_wsfe(&io___189);
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *n - i__ + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__1 = x[i__] - xs;
			i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&con[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&cap[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			r__2 = vi + vvi;
			i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&sink[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)
				sizeof(integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__3 = vi / consn;
			i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (js = 1; js <= i__3; ++js) {
			    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * 
				    conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__4 = *ns;
			for (js = 1; js <= i__4; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * 
				    sorb_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__1 = e_wsfe();
			if (i__1 != 0) {
			    goto L901;
			}
		    } else {
			i__1 = s_wsfe(&io___190);
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *n - i__ + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__1 = x[i__] - xs;
			i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&con[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&cap[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			r__2 = vi + vvi;
			i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&sink[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)
				sizeof(integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__3 = vi / consn;
			i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (js = 1; js <= i__3; ++js) {
			    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * 
				    conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__4 = *ns;
			for (js = 1; js <= i__4; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * 
				    sorb_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__5 = *ns;
			for (js = 1; js <= i__5; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb2[js + i__ * 
				    sorb2_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__1 = e_wsfe();
			if (i__1 != 0) {
			    goto L901;
			}
		    }
		} else {
		    if (! (*ldualneq)) {
			i__1 = s_wsfe(&io___191);
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *n - i__ + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__1 = x[i__] - xs;
			i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&sinkim[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thnewim[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			r__2 = vi + vvi;
			i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&strans[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)
				sizeof(integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__3 = vi / consn;
			i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (js = 1; js <= i__3; ++js) {
			    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * 
				    conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__4 = *ns;
			for (js = 1; js <= i__4; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * 
				    sorb_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__1 = e_wsfe();
			if (i__1 != 0) {
			    goto L901;
			}
		    } else {
			i__1 = s_wsfe(&io___192);
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *n - i__ + 1;
			i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(
				integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__1 = x[i__] - xs;
			i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&sinkim[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&thnewim[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			r__2 = vi + vvi;
			i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&strans[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)
				sizeof(integer));
			if (i__1 != 0) {
			    goto L901;
			}
			r__3 = vi / consn;
			i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(
				real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (js = 1; js <= i__3; ++js) {
			    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * 
				    conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__4 = *ns;
			for (js = 1; js <= i__4; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * 
				    sorb_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__5 = *ns;
			for (js = 1; js <= i__5; ++js) {
			    i__1 = do_fio(&c__1, (char *)&sorb2[js + i__ * 
				    sorb2_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
			i__1 = e_wsfe();
			if (i__1 != 0) {
			    goto L901;
			}
		    }
		}
	    } else {
		i__1 = s_wsfe(&io___193);
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *n - i__ + 1;
		i__1 = do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		if (i__1 != 0) {
		    goto L901;
		}
		r__1 = x[i__] - xs;
		i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&hnew[i__], (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&thn[i__], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&con[i__], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&cap[i__], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		r__2 = vi + vvi;
		i__1 = do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&sink[i__], (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&kappa[i__], (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L901;
		}
		r__3 = vi / consn;
		i__1 = do_fio(&c__1, (char *)&r__3, (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&temp[i__], (ftnlen)sizeof(real))
			;
		if (i__1 != 0) {
		    goto L901;
		}
		i__3 = *ns;
		for (js = 1; js <= i__3; ++js) {
		    i__1 = do_fio(&c__1, (char *)&conc[js + i__ * conc_dim1], 
			    (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__4 = *ns;
		for (js = 1; js <= i__4; ++js) {
		    i__1 = do_fio(&c__1, (char *)&sorb[js + i__ * sorb_dim1], 
			    (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__5 = *ns;
		for (js = 1; js <= i__5; ++js) {
		    i__1 = do_fio(&c__1, (char *)&sorb2[js + i__ * sorb2_dim1]
			    , (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	}
/* L11: */
    }
    i__1 = s_wsfe(&io___194);
    if (i__1 != 0) {
	goto L901;
    }
    i__1 = e_wsfe();
    if (i__1 != 0) {
	goto L901;
    }
    if (*lvaporout && *lvapor) {
	i__1 = s_wsfe(&io___195);
	if (i__1 != 0) {
	    goto L901;
	}
	i__1 = e_wsfe();
	if (i__1 != 0) {
	    goto L901;
	}
    }
    return 0;
/*     Error when writing into an output file */
L901:
    *ierr = 1;
    return 0;
} /* nodout_ */

/* ********************************************************************** */
/* Subroutine */ int obsnod_(doublereal *t, integer *n, integer *nobs, 
	integer *ns, integer *nsd, integer *node, real *conc, real *hnew, 
	real *thnew, real *tempn, logical *lchem, real *thnewim, real *vnew, 
	real *vvnew, logical *lflux, integer *ierr)
{
    /* Format strings */
    static char fmt_102[] = "(2x,f14.4,100(f12.2,f8.4,e11.3,2x))";
    static char fmt_100[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,2x))";
    static char fmt_103[] = "(x,e15.8,100(f12.2,f8.4,e11.3,2x))";
    static char fmt_101[] = "(x,e15.8,100(f12.2,f8.4,f9.3,2x))";
    static char fmt_112[] = "(2x,f14.4,100(f12.2,f8.4,2e12.4,2x))";
    static char fmt_122[] = "(2x,f14.4,100(f12.2,f8.4,3e12.4,2x))";
    static char fmt_132[] = "(2x,f14.4,100(f12.2,f8.4,4e12.4,2x))";
    static char fmt_140[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,4e12.4,2x))";
    static char fmt_150[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,5e12.4,2x))";
    static char fmt_160[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,6e12.4,2x))";
    static char fmt_170[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,7e12.4,2x))";
    static char fmt_180[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,8e12.4,2x))";
    static char fmt_190[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,9e12.4,2x))";
    static char fmt_200[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,10e12.4,2x))";
    static char fmt_210[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,11e12.4,2x))";
    static char fmt_220[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,12e12.4,2x))";
    static char fmt_110[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,e12.4,2x))";
    static char fmt_120[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,2e12.4,2x))";
    static char fmt_130[] = "(2x,f14.4,100(f12.2,f8.4,f9.3,3e12.4,2x))";
    static char fmt_111[] = "(x,e15.8,100(f12.2,f8.4,f9.3,e12.4,2x))";
    static char fmt_121[] = "(x,e15.8,100(f12.2,f8.4,f9.3,2e12.4,2x))";
    static char fmt_131[] = "(x,e15.8,100(f12.2,f8.4,f9.3,3e12.4,2x))";
    static char fmt_141[] = "(x,e15.8,100(f12.2,f8.4,f9.3,4e12.4,2x))";
    static char fmt_151[] = "(x,e15.8,100(f12.2,f8.4,f9.3,5e12.4,2x))";
    static char fmt_161[] = "(x,e15.8,100(f12.2,f8.4,f9.3,6e12.4,2x))";
    static char fmt_171[] = "(x,e15.8,100(f12.2,f8.4,f9.3,7e12.4,2x))";
    static char fmt_181[] = "(x,e15.8,100(f12.2,f8.4,f9.3,8e12.4,2x))";
    static char fmt_191[] = "(x,e15.8,100(f12.2,f8.4,f9.3,9e12.4,2x))";
    static char fmt_201[] = "(x,e15.8,100(f12.2,f8.4,f9.3,10e12.4,2x))";
    static char fmt_211[] = "(x,e15.8,100(f12.2,f8.4,f9.3,11e12.4,2x))";
    static char fmt_221[] = "(x,e15.8,100(f12.2,f8.4,f9.3,12e12.4,2x))";

    /* System generated locals */
    integer conc_dim1, conc_offset, i__1, i__2, i__3;
    real r__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static real th[1000], eca[10];
    static logical lec;
    static real ecw, thw;

    /* Fortran I/O blocks */
    static cilist io___198 = { 1, 77, 0, fmt_102, 0 };
    static cilist io___199 = { 1, 77, 0, fmt_100, 0 };
    static cilist io___200 = { 1, 77, 0, fmt_103, 0 };
    static cilist io___201 = { 1, 77, 0, fmt_101, 0 };
    static cilist io___206 = { 1, 77, 0, fmt_112, 0 };
    static cilist io___208 = { 1, 77, 0, fmt_122, 0 };
    static cilist io___209 = { 1, 77, 0, fmt_132, 0 };
    static cilist io___210 = { 1, 77, 0, fmt_140, 0 };
    static cilist io___211 = { 1, 77, 0, fmt_150, 0 };
    static cilist io___212 = { 1, 77, 0, fmt_160, 0 };
    static cilist io___213 = { 1, 77, 0, fmt_170, 0 };
    static cilist io___214 = { 1, 77, 0, fmt_180, 0 };
    static cilist io___215 = { 1, 77, 0, fmt_190, 0 };
    static cilist io___216 = { 1, 77, 0, fmt_200, 0 };
    static cilist io___217 = { 1, 77, 0, fmt_210, 0 };
    static cilist io___218 = { 1, 77, 0, fmt_220, 0 };
    static cilist io___219 = { 1, 77, 0, fmt_110, 0 };
    static cilist io___220 = { 1, 77, 0, fmt_120, 0 };
    static cilist io___221 = { 1, 77, 0, fmt_130, 0 };
    static cilist io___222 = { 1, 77, 0, fmt_140, 0 };
    static cilist io___223 = { 1, 77, 0, fmt_150, 0 };
    static cilist io___224 = { 1, 77, 0, fmt_160, 0 };
    static cilist io___225 = { 1, 77, 0, fmt_170, 0 };
    static cilist io___226 = { 1, 77, 0, fmt_180, 0 };
    static cilist io___227 = { 1, 77, 0, fmt_190, 0 };
    static cilist io___228 = { 1, 77, 0, fmt_200, 0 };
    static cilist io___229 = { 1, 77, 0, fmt_210, 0 };
    static cilist io___230 = { 1, 77, 0, fmt_220, 0 };
    static cilist io___231 = { 1, 77, 0, fmt_111, 0 };
    static cilist io___232 = { 1, 77, 0, fmt_121, 0 };
    static cilist io___233 = { 1, 77, 0, fmt_131, 0 };
    static cilist io___234 = { 1, 77, 0, fmt_141, 0 };
    static cilist io___235 = { 1, 77, 0, fmt_151, 0 };
    static cilist io___236 = { 1, 77, 0, fmt_161, 0 };
    static cilist io___237 = { 1, 77, 0, fmt_171, 0 };
    static cilist io___238 = { 1, 77, 0, fmt_181, 0 };
    static cilist io___239 = { 1, 77, 0, fmt_191, 0 };
    static cilist io___240 = { 1, 77, 0, fmt_201, 0 };
    static cilist io___241 = { 1, 77, 0, fmt_211, 0 };
    static cilist io___242 = { 1, 77, 0, fmt_221, 0 };
    static cilist io___243 = { 1, 77, 0, fmt_111, 0 };
    static cilist io___244 = { 1, 77, 0, fmt_121, 0 };
    static cilist io___245 = { 1, 77, 0, fmt_131, 0 };
    static cilist io___246 = { 1, 77, 0, fmt_141, 0 };
    static cilist io___247 = { 1, 77, 0, fmt_151, 0 };
    static cilist io___248 = { 1, 77, 0, fmt_161, 0 };
    static cilist io___249 = { 1, 77, 0, fmt_171, 0 };
    static cilist io___250 = { 1, 77, 0, fmt_181, 0 };
    static cilist io___251 = { 1, 77, 0, fmt_191, 0 };
    static cilist io___252 = { 1, 77, 0, fmt_201, 0 };
    static cilist io___253 = { 1, 77, 0, fmt_211, 0 };
    static cilist io___254 = { 1, 77, 0, fmt_221, 0 };


    /* Parameter adjustments */
    --vvnew;
    --vnew;
    --thnewim;
    --tempn;
    --thnew;
    --hnew;
    --node;
    conc_dim1 = *nsd;
    conc_offset = 1 + conc_dim1;
    conc -= conc_offset;

    /* Function Body */
    i__1 = *nobs;
    for (i__ = 1; i__ <= i__1; ++i__) {
	th[i__ - 1] = thnew[node[i__]] + thnewim[node[i__]];
/* L10: */
    }
    if (! (*lchem)) {
	if (*t < 99999999.f) {
	    if (*lflux) {
		i__1 = s_wsfe(&io___198);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *nobs;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)sizeof(
			    real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    r__1 = vnew[node[i__]] + vvnew[node[i__]];
		    i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    } else {
		i__1 = s_wsfe(&io___199);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *nobs;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)sizeof(
			    real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	} else {
	    if (*lflux) {
		i__1 = s_wsfe(&io___200);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *nobs;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)sizeof(
			    real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    r__1 = vnew[node[i__]] + vvnew[node[i__]];
		    i__1 = do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    } else {
		i__1 = s_wsfe(&io___201);
		if (i__1 != 0) {
		    goto L901;
		}
		i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal)
			);
		if (i__1 != 0) {
		    goto L901;
		}
		i__2 = *nobs;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)sizeof(
			    real));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (ftnlen)
			    sizeof(real));
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		i__1 = e_wsfe();
		if (i__1 != 0) {
		    goto L901;
		}
	    }
	}
    } else {
	lec = FALSE_;
	if (lec) {
	    i__1 = *nobs;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		d__1 = (doublereal) (conc[node[i__] * conc_dim1 + 1] / 
			.008465f);
		ecw = pow_dd(&d__1, &c_b604);
		thw = thnew[node[i__]];
		eca[i__ - 1] = ecw * 1.45f * thw * thw + thw * .102f;
/* L11: */
	    }
	}
	if (*t < 99999999.f) {
	    if (*lflux) {
		if (*ns == 1) {
		    i__1 = s_wsfe(&io___206);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 2) {
		    i__1 = s_wsfe(&io___208);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 3) {
		    i__1 = s_wsfe(&io___209);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 4) {
		    i__1 = s_wsfe(&io___210);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 5) {
		    i__1 = s_wsfe(&io___211);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 6) {
		    i__1 = s_wsfe(&io___212);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 7) {
		    i__1 = s_wsfe(&io___213);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 8) {
		    i__1 = s_wsfe(&io___214);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 9) {
		    i__1 = s_wsfe(&io___215);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 10) {
		    i__1 = s_wsfe(&io___216);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 11) {
		    i__1 = s_wsfe(&io___217);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 12) {
		    i__1 = s_wsfe(&io___218);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
	    } else {
		if (*ns == 1) {
		    i__1 = s_wsfe(&io___219);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 2) {
		    i__1 = s_wsfe(&io___220);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 3) {
		    i__1 = s_wsfe(&io___221);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 4) {
		    i__1 = s_wsfe(&io___222);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 5) {
		    i__1 = s_wsfe(&io___223);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 6) {
		    i__1 = s_wsfe(&io___224);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 7) {
		    i__1 = s_wsfe(&io___225);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 8) {
		    i__1 = s_wsfe(&io___226);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 9) {
		    i__1 = s_wsfe(&io___227);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 10) {
		    i__1 = s_wsfe(&io___228);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 11) {
		    i__1 = s_wsfe(&io___229);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 12) {
		    i__1 = s_wsfe(&io___230);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
	    }
	} else {
	    if (*lflux) {
		if (*ns == 1) {
		    i__1 = s_wsfe(&io___231);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 2) {
		    i__1 = s_wsfe(&io___232);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 3) {
		    i__1 = s_wsfe(&io___233);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 4) {
		    i__1 = s_wsfe(&io___234);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 5) {
		    i__1 = s_wsfe(&io___235);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 6) {
		    i__1 = s_wsfe(&io___236);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 7) {
		    i__1 = s_wsfe(&io___237);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 8) {
		    i__1 = s_wsfe(&io___238);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 9) {
		    i__1 = s_wsfe(&io___239);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 10) {
		    i__1 = s_wsfe(&io___240);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 11) {
		    i__1 = s_wsfe(&io___241);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 12) {
		    i__1 = s_wsfe(&io___242);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&vnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
	    } else {
		if (*ns == 1) {
		    i__1 = s_wsfe(&io___243);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 2) {
		    i__1 = s_wsfe(&io___244);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 3) {
		    i__1 = s_wsfe(&io___245);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 4) {
		    i__1 = s_wsfe(&io___246);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 5) {
		    i__1 = s_wsfe(&io___247);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 6) {
		    i__1 = s_wsfe(&io___248);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 7) {
		    i__1 = s_wsfe(&io___249);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 8) {
		    i__1 = s_wsfe(&io___250);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 9) {
		    i__1 = s_wsfe(&io___251);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 10) {
		    i__1 = s_wsfe(&io___252);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 11) {
		    i__1 = s_wsfe(&io___253);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__2 = *nobs;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__3 = *ns;
			for (j = 1; j <= i__3; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
		if (*ns == 12) {
		    i__1 = s_wsfe(&io___254);
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__1 = do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(
			    doublereal));
		    if (i__1 != 0) {
			goto L901;
		    }
		    i__3 = *nobs;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__1 = do_fio(&c__1, (char *)&hnew[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&th[i__ - 1], (ftnlen)
				sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__1 = do_fio(&c__1, (char *)&tempn[node[i__]], (
				ftnlen)sizeof(real));
			if (i__1 != 0) {
			    goto L901;
			}
			i__2 = *ns;
			for (j = 1; j <= i__2; ++j) {
			    i__1 = do_fio(&c__1, (char *)&conc[j + node[i__] *
				     conc_dim1], (ftnlen)sizeof(real));
			    if (i__1 != 0) {
				goto L901;
			    }
			}
		    }
		    i__1 = e_wsfe();
		    if (i__1 != 0) {
			goto L901;
		    }
		}
	    }
	}
    }
    return 0;
/*     Error when writing into an output file */
L901:
    *ierr = 1;
    return 0;
} /* obsnod_ */

