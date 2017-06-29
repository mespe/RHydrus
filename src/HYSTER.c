/* HYSTER.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    doublereal alphad, alphai, xn, xm, xxm, sm, sarwi;
    integer ipath, ihyst;
    doublereal xnw, xmw, xxmw;
} properties_;

#define properties_1 properties_

struct {
    doublereal rhsw[7007]	/* was [1001][7] */, rasw[7007]	/* was [1001][
	    7] */, sarw[1001], rswaw[1001], phsw[1001];
    integer mpsw[1001], jjh[1001], ipsw[1001];
} glob_;

#define glob_1 glob_

struct {
    doublereal sw, sraw, sat, rsw, raw, permw;
    integer mpsw1, jjjh, ipsw1, ilsw;
    doublereal esw, esat, asw, dww;
} local_;

#define local_1 local_

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b89 = .5;
static doublereal c_b90 = 2.;

/* Source file HYSTER.FOR |||||||||||||||||||||||||||||||||||||||||||||| */
/*     Written/revised by RJ Lenhard, Nov 2004 */
/*     Program to evaluate hysteretic water-wet K-S-P routines */
/* Subroutine */ int hyst_(integer *numnp, integer *nmatd, real *pard, real *
	parw, integer *matnum, integer *kappa, real *hnew, real *hold, real *
	theta, real *con, real *cap, integer *ikappa, integer *ikod)
{
    /* Format strings */
    static char fmt_1000[] = "(//)";
    static char fmt_1002[] = "(3x,\002CONVERGED RESULTS\002)";
    static char fmt_1003[] = "(/5x,\002RASW-1 = \002,f5.3,5x,\002RASW-2 ="
	    " \002,f5.3,5x,\002RASW-3 = \002,f5.3,5x,\002RASW-4 = \002,f5.3,/"
	    "5x,\002RASW-5 = \002,f5.3,5x,\002RASW-6 = \002,f5.3,5x,\002RASW-"
	    "7 = \002,f5.3)";
    static char fmt_1004[] = "(/5x,\002RHSW-1 = \002,f5.1,5x,\002RHSW-2 ="
	    " \002,f5.1,5x,\002RHSW-3 = \002,f5.1,5x,\002RHSW-4 = \002,f5.1,/"
	    "5x,\002RHSW-5 = \002,f5.1,5x,\002RHSW-6 = \002,f5.1,5x,\002RHSW-"
	    "7 = \002,f5.1)";
    static char fmt_1005[] = "(//5x,\002PHSW = \002,f7.1,5x,\002SARW = \002,"
	    "f5.3,5x,\002SRAW = \002,f5.3,5x,\002SARWI = \002,f5.3,/5x,\002RS"
	    "WAW = \002,f5.3,4x,\002RSW = \002,f5.3,6x,\002JJJH = \002,i1,9x"
	    ",\002JJH = \002,i1,/5x,\002IPSW = \002,i1,9x,\002IPSW1 = \002,i1"
	    ",8x,\002MPSW = \002,i1,9x,\002MPSW1 = \002,i1,/5x,\002RAW = \002"
	    ",f5.3)";
    static char fmt_1006[] = "(//5x,\002ASW = \002,f5.3,6x,\002ESW = \002,f5"
	    ".3,6x,\002SW = \002,f5.3,5x,\002ILSW = \002,i1,/5x\002ESAT = "
	    "\002,f5.3,5x,\002SAT = \002,f5.3,6x,\002PERMW = \002,e12.3)";
    static char fmt_1007[] = "(3x,\002**************************************"
	    "********\002)";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer m, n;
    extern /* Subroutine */ int hysterini_(integer *);
    static doublereal ha, hw;
    extern /* Subroutine */ int updatehyst_(integer *, doublereal *, 
	    doublereal *), path_(integer *, doublereal *, doublereal *);
    static logical lprint;

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___7 = { 0, 6, 0, fmt_1002, 0 };
    static cilist io___8 = { 0, 6, 0, fmt_1003, 0 };
    static cilist io___9 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___10 = { 0, 6, 0, fmt_1005, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_1006, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_1007, 0 };


/*     SW  - actual water content */
/*     ESW - effective water content */
/*     ASW - apparent water content */
/*     PERMW - relative conductivity */
/*     ALPHAD = 'alpha for main drainage' */
/*     ALPHAI = 'alpha for main imbibition' */
/*     XN = 'VG n parameter' for drying */
/*     XNW= 'VG n parameter' for wetting */
/*     XM = 1.d+0 - 1.d+0/XN */
/*     XXM= 1.d+0/XM */
/*     SM = 'Sr' */
    /* Parameter adjustments */
    --cap;
    --con;
    --theta;
    --hold;
    --hnew;
    --kappa;
    --matnum;
    parw -= 12;
    pard -= 12;

    /* Function Body */
    lprint = FALSE_;
/*     iPath - number of reversal points the code remembers (<=7) */
/*     =1: Fluid entrapment only, non hysteretic, i.e., alphaw=alphad */
    properties_1.ipath = 7;
/*     iHyst - Location of the initial condition */
/*     =1: Main drainage curve */
/*     =2: Main inhibition curve */
/*     =3: Primary drainage curve (with maximum entrapped air) */
    properties_1.ihyst = 1;
    if (*ikod == 1) {
	if (*ikappa == -1) {
	    properties_1.ihyst = 1;
	}
	if (*ikappa == 1) {
	    properties_1.ihyst = 2;
	}
	hysterini_(numnp);
    }
    i__1 = *numnp;
    for (n = 1; n <= i__1; ++n) {
	m = matnum[n];
	properties_1.alphad = pard[m * 11 + 3];
	properties_1.alphai = parw[m * 11 + 3];
	properties_1.xn = pard[m * 11 + 4];
	properties_1.xm = 1. - 1. / properties_1.xn;
	properties_1.xxm = 1. / properties_1.xm;
/*       different n values for the drying curve */
	properties_1.xnw = parw[m * 11 + 4];
	properties_1.xmw = 1. - 1. / properties_1.xnw;
	properties_1.xxmw = 1. / properties_1.xmw;
/*       Residual water saturation */
	properties_1.sm = pard[m * 11 + 1] / pard[m * 11 + 2];
/*        sm=ParD(1,m)/(ParD(2,m)-ParD(1,m)) */
/*       Max amount of air that get traped on main inhibition branch (saturation) */
/*        SARWI=(ParD(2,m)-ParW(2,m))/ParD(2,m) */
	properties_1.sarwi = (pard[m * 11 + 2] - parw[m * 11 + 2]) / (pard[m *
		 11 + 2] - pard[m * 11 + 1]);
/*       Initialize arrays for hysteretic saturations */
	hw = hnew[n];
	ha = 0.f;
	glob_1.phsw[n - 1] = -hold[n];
/*       subprogram for hysteretic routines */
	path_(&n, &hw, &ha);
/*       after convergence */
	if (*ikod == 1 || *ikod == 3) {
	    updatehyst_(&n, &hw, &ha);
	}
	theta[n] = local_1.sw * pard[m * 11 + 2];
/*        Theta(n)=SW*(ParD(2,m)-ParD(1,m))+ParD(1,m) */
	con[n] = pard[m * 11 + 5] * (real) local_1.permw;
	cap[n] = (real) local_1.dww * pard[m * 11 + 2];
	if (glob_1.jjh[n - 1] == 0) {
	    kappa[n] = -1;
	}
	if (glob_1.jjh[n - 1] == 1) {
	    kappa[n] = 1;
	}
	if (lprint && *ikod == 3) {
	    s_wsfe(&io___6);
	    e_wsfe();
	    s_wsfe(&io___7);
	    e_wsfe();
	    s_wsfe(&io___8);
	    do_fio(&c__1, (char *)&glob_1.rasw[0], (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&glob_1.rasw[1001], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rasw[2002], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rasw[3003], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rasw[4004], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rasw[5005], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rasw[6006], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    s_wsfe(&io___9);
	    do_fio(&c__1, (char *)&glob_1.rhsw[0], (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&glob_1.rhsw[1001], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rhsw[2002], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rhsw[3003], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rhsw[4004], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rhsw[5005], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rhsw[6006], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    s_wsfe(&io___10);
	    do_fio(&c__1, (char *)&glob_1.phsw[0], (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&glob_1.sarw[0], (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&local_1.sraw, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&properties_1.sarwi, (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&glob_1.rswaw[0], (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&local_1.rsw, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&local_1.jjjh, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&glob_1.jjh[0], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&glob_1.ipsw[0], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&local_1.ipsw1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&glob_1.mpsw[0], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&local_1.mpsw1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&local_1.raw, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    s_wsfe(&io___11);
	    do_fio(&c__1, (char *)&local_1.asw, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&local_1.esw, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&local_1.sw, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&local_1.ilsw, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&local_1.esat, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&local_1.sat, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&local_1.permw, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    s_wsfe(&io___12);
	    e_wsfe();
	}
/* L11: */
    }
    return 0;
/* L1001: */
} /* hyst_ */

/* *********************************************************************** */
/*     Initialize arrays for hysteretic saturations. */
/*     Written/revised by RJ Lenhard, Nov. 2004 */
/* Subroutine */ int hysterini_(integer *n)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, ij, ji;
    static doublereal one, zero;

/*     Initializing reversal points, depending on number of saturation */
/*     paths being considered (odd versus even number) */
/*      IRES = restart option */
/*      IF( IRES .GE. 1) GO TO 210 */
    zero = 0.;
    one = 1.;
    i__1 = properties_1.ipath;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (properties_1.ipath - (i__ << 1) == 0) {
	    goto L110;
	}
/* L100: */
    }
    goto L160;
L110:
    i__1 = properties_1.ipath;
    for (ij = 4; ij <= i__1; ij += 2) {
	i__2 = *n;
	for (ji = 1; ji <= i__2; ++ji) {
	    glob_1.rasw[ji + ij * 1001 - 1002] = one;
/* L120: */
	}
/* L130: */
    }
    i__1 = properties_1.ipath - 1;
    for (ij = 3; ij <= i__1; ij += 2) {
	i__2 = *n;
	for (ji = 1; ji <= i__2; ++ji) {
	    glob_1.rasw[ji + ij * 1001 - 1002] = zero;
/* L140: */
	}
/* L150: */
    }
    goto L210;
L160:
    i__1 = properties_1.ipath - 1;
    for (ij = 4; ij <= i__1; ij += 2) {
	i__2 = *n;
	for (ji = 1; ji <= i__2; ++ji) {
	    glob_1.rasw[ji + ij * 1001 - 1002] = one;
/* L170: */
	}
/* L180: */
    }
    i__1 = properties_1.ipath;
    for (ij = 3; ij <= i__1; ij += 2) {
	i__2 = *n;
	for (ji = 1; ji <= i__2; ++ji) {
	    glob_1.rasw[ji + ij * 1001 - 1002] = zero;
/* L190: */
	}
/* L200: */
    }
L210:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	glob_1.rasw[i__ - 1] = one;
	glob_1.rhsw[i__ - 1] = zero;
	glob_1.mpsw[i__ - 1] = 3;
	glob_1.sarw[i__ - 1] = zero;
	glob_1.rswaw[i__ - 1] = one;
/* L220: */
    }
/*     Setting historical values for given initial conditions */
/*     For nonhysteresis and fluid entrapment options */
    if (properties_1.ipath == 1) {
	if (properties_1.sarwi == zero) {
	    goto L260;
	}
	if (properties_1.ihyst == 2) {
	    glob_1.rswaw[i__ - 1] = zero;
	}
	goto L260;
    }
/*     For hysteresis options */
    if (properties_1.ihyst == 1) {
/* ---    Starting from the main drainage branch */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    glob_1.rhsw[i__ + 1000] = zero;
	    glob_1.rasw[i__ + 1000] = one;
	    glob_1.jjh[i__ - 1] = 0;
	    glob_1.phsw[i__ - 1] = zero;
	    glob_1.rswaw[i__ - 1] = one;
	    glob_1.ipsw[i__ - 1] = 1;
/* L230: */
	}
    } else if (properties_1.ihyst == 2) {
/*       Starting from the main imbibition branch */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    glob_1.rhsw[i__ + 1000] = 1e3;
	    glob_1.rasw[i__ + 1000] = zero;
	    glob_1.phsw[i__ - 1] = 1e3;
	    glob_1.jjh[i__ - 1] = 1;
	    glob_1.rswaw[i__ - 1] = zero;
	    glob_1.ipsw[i__ - 1] = 2;
/* L240: */
	}
    } else if (properties_1.ihyst == 3) {
/*       Starting from drainage scanning curve when all depths where */
/*       previously apparent saturated */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    glob_1.rhsw[i__ + 1000] = 1e5;
	    glob_1.rhsw[i__ + 2001] = zero;
	    glob_1.rasw[i__ + 1000] = zero;
	    glob_1.rasw[i__ + 2001] = one;
	    glob_1.phsw[i__ - 1] = zero;
	    glob_1.jjh[i__ - 1] = 0;
	    glob_1.rswaw[i__ - 1] = zero;
	    glob_1.ipsw[i__ - 1] = 3;
/* L250: */
	}
    }
L260:
/*    End of HYINI group */
    return 0;
} /* hysterini_ */

/* *********************************************************************** */
/*     Beginning subprogram for hysteretic routines */
/*     Written/revised by RJ Lenhard, Nov 2004 */
/* Subroutine */ int path_(integer *n, doublereal *hw, doublereal *ha)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal haw, one, zero;
    extern /* Subroutine */ int hawpath_(integer *, doublereal *);

    zero = 0.;
    one = 1.;
/*     Assigning property values */
    if (properties_1.sarwi == zero) {
	local_1.raw = zero;
	local_1.sraw = zero;
    } else {
	local_1.raw = one / properties_1.sarwi - one;
    }
/* Computing MAX */
    d__1 = zero, d__2 = *ha - *hw;
    haw = max(d__1,d__2);
/*     Initializing common block variables */
    local_1.esat = zero;
    local_1.asw = zero;
/*     Setting global variables to local variables */
    local_1.jjjh = glob_1.jjh[*n - 1];
    local_1.ipsw1 = glob_1.ipsw[*n - 1];
    local_1.mpsw1 = glob_1.mpsw[*n - 1];
    local_1.rsw = glob_1.rswaw[*n - 1];
    if (local_1.rsw == 0.f) {
	glob_1.sarw[*n - 1] = properties_1.sarwi;
    }
    local_1.sraw = glob_1.sarw[*n - 1];
/*     Two phase K-S-P relations */
    hawpath_(n, &haw);
    return 0;
/*     End of PATH group */
} /* path_ */

/* *********************************************************************** */
/*     Determining saturation path history & calling k-S-P routines */
/*     Written/revised by RJ Lenhard, Nov 2004 */
/* Subroutine */ int hawpath_(integer *n, doublereal *haw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer l;
    static doublereal one;
    extern /* Subroutine */ int wet_(integer *, doublereal *), dry_(integer *,
	     doublereal *);
    static doublereal zero;
    static integer ipsw2;
    extern /* Subroutine */ int drain_(integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___25 = { 0, 6, 0, "(A)", 0 };
    static cilist io___26 = { 0, 6, 0, "(A)", 0 };
    static cilist io___27 = { 0, 6, 0, "(A)", 0 };


    zero = 0.;
    one = 1.;
/*     For fluid entrapment only option */
    if (properties_1.ipath == 1) {
	drain_(n, haw);
	local_1.ilsw = 0;
	return 0;
    }
/*     For hysteresis option */
/*     Main air-water drainage branch */
L100:
    if (*haw >= glob_1.rhsw[*n + 1000]) {
	drain_(n, haw);
	local_1.ilsw = 0;
	local_1.jjjh = 0;
	local_1.ipsw1 = 1;
	local_1.mpsw1 = 3;
	return 0;
    }
/*     Wetting scanning paths */
L110:
    if ((*haw < glob_1.phsw[*n - 1] && local_1.jjjh == 0 || *haw <= 
	    glob_1.phsw[*n - 1] && local_1.jjjh == 1) && local_1.mpsw1 == 3 ||
	     local_1.mpsw1 == 1 || local_1.jjjh == 1 && *haw == zero) {
/*       On max. sat. path, but passed reveral point and became a drying path */
	if (local_1.ipsw1 == properties_1.ipath && *haw > glob_1.rhsw[*n + 
		properties_1.ipath * 1001 - 1002]) {
	    local_1.mpsw1 = 0;
	    local_1.ipsw1 = properties_1.ipath - 1;
	    goto L100;
	}
/*       Close of path 2 at HAW=0. */
	if (*haw == 0.) {
	    if (local_1.ipsw1 < properties_1.ipath) {
		local_1.mpsw1 = 3;
	    }
	    local_1.ipsw1 = 2;
	    local_1.jjjh = 1;
	    wet_(n, haw);
	    return 0;
	}
/*       Switching from drying to wetting paths */
	if (local_1.jjjh == 0 && local_1.mpsw1 == 3) {
	    ++local_1.ipsw1;
/*         Switched to max. sat. path, which is a wetting path */
	    if (local_1.ipsw1 >= properties_1.ipath) {
		local_1.mpsw1 = 1;
		local_1.ipsw1 = properties_1.ipath;
	    }
	}
/*       Evaluating the current wetting path number and whether */
/*       saturation paths have closed */
	ipsw2 = local_1.ipsw1;
	i__1 = ipsw2 - 1;
	for (l = 0; l <= i__1; l += 2) {
	    if (*haw < glob_1.rhsw[*n + (local_1.ipsw1 - l) * 1001 - 1002] && 
		    *haw > glob_1.rhsw[*n + (local_1.ipsw1 - 1 - l) * 1001 - 
		    1002]) {
		local_1.ipsw1 -= l;
		if (local_1.ipsw1 < properties_1.ipath) {
		    local_1.mpsw1 = 3;
		}
		wet_(n, haw);
		local_1.ilsw = l;
		local_1.jjjh = 1;
		return 0;
	    }
/* L120: */
	}
	s_wsfe(&io___25);
	do_fio(&c__1, "ERROR: Passed Through Wetting Scanning Path Loop With"
		"out Entering subroutine WET2P", (ftnlen)82);
	e_wsfe();
	goto L200;
/*     Drying scanning paths */
    } else if ((*haw > glob_1.phsw[*n - 1] && local_1.jjjh == 1 || *haw >= 
	    glob_1.phsw[*n - 1] && local_1.jjjh == 0) && local_1.mpsw1 == 3 ||
	     local_1.mpsw1 == 0) {
/*       K-S-P Relations corresponding to the main drainage branch */
	if (local_1.ipsw1 == 1) {
	    drain_(n, haw);
	    local_1.ilsw = 0;
	    local_1.jjjh = 0;
	    local_1.ipsw1 = 1;
	    local_1.mpsw1 = 3;
	    return 0;
	}
/*       On max. sat. path, but passed reveral point and became */
/*       a wetting path */
	if (local_1.ipsw1 == properties_1.ipath && *haw < glob_1.rhsw[*n + 
		properties_1.ipath * 1001 - 1002]) {
	    local_1.mpsw1 = 1;
	    local_1.ipsw1 = properties_1.ipath - 1;
	    goto L110;
	}
/*       Switching from wetting to drying paths */
	if (local_1.jjjh == 1 && local_1.mpsw1 == 3) {
	    ++local_1.ipsw1;
	    if (local_1.ipsw1 >= properties_1.ipath) {
		local_1.mpsw1 = 0;
		local_1.ipsw1 = properties_1.ipath;
	    }
	}
/*       Evaluating the current drying path number and whether */
/*       saturation paths have closed */
	ipsw2 = local_1.ipsw1;
	i__1 = ipsw2 - 2;
	for (l = 0; l <= i__1; l += 2) {
	    if (*haw > glob_1.rhsw[*n + (local_1.ipsw1 - l) * 1001 - 1002] && 
		    *haw < glob_1.rhsw[*n + (local_1.ipsw1 - 1 - l) * 1001 - 
		    1002]) {
		local_1.ipsw1 -= l;
		if (local_1.ipsw1 < properties_1.ipath) {
		    local_1.mpsw1 = 3;
		}
		dry_(n, haw);
		local_1.ilsw = l;
		local_1.jjjh = 0;
		return 0;
	    }
/* L140: */
	}
	s_wsfe(&io___26);
	do_fio(&c__1, "ERROR: Passed Through Three-Phase Liquid Wetting Loop"
		" Without Entering subroutine DRY2P", (ftnlen)87);
	e_wsfe();
	goto L200;
    } else {
	s_wsfe(&io___27);
	do_fio(&c__1, "ERROR: Passed Through subroutine HAWPATH Without Ente"
		"ring a Two-Phase Liquid Saturation Path", (ftnlen)92);
	e_wsfe();
    }
L200:
    return 0;
} /* hawpath_ */

/* *********************************************************************** */
/*     Main drainage k-S-P relations */
/*     Options: nonhysteretic, fluid entrapment only, hysteretic */
/*     Written/revised by RJ Lenhard, Nov 2004 */
/* Subroutine */ int drain_(integer *n, doublereal *haw)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal p0, p3, x1, x7, one, dwwd, zero, csatw;

    /* Fortran I/O blocks */
    static cilist io___36 = { 0, 6, 0, "(A,I6)", 0 };


/*     Apparent water saturation for all options */
    zero = 0.;
    one = 1.;
/* Computing MIN */
    d__3 = properties_1.alphad * *haw;
    d__2 = pow_dd(&d__3, &properties_1.xn) + 1.;
    d__4 = -properties_1.xm;
    d__1 = pow_dd(&d__2, &d__4);
    local_1.asw = min(d__1,one);
    if (local_1.asw < local_1.rsw) {
	local_1.rsw = local_1.asw;
	if (properties_1.sarwi != 0.f) {
	    local_1.sraw = (one - local_1.asw) / (one + local_1.raw * (one - 
		    local_1.asw));
	} else {
	    local_1.sraw = 0.;
	}
    }
    if (properties_1.ipath == 1) {
	local_1.esat = local_1.sraw * (local_1.asw - local_1.rsw) / (one - 
		local_1.rsw);
    } else {
	local_1.esat = 0.;
    }
/*     Effective and actual saturations for all options */
    local_1.esw = local_1.asw - local_1.esat;
    local_1.sw = local_1.esw * (one - properties_1.sm) + properties_1.sm;
    local_1.sat = local_1.esat * (one - properties_1.sm);
    if (properties_1.ipath != 1) {
	goto L100;
    }
/*     Derivative terms for nonhysteretic/fluid entrapment options */
    x1 = pow_dd(&local_1.asw, &properties_1.xxm);
    d__1 = one - x1;
    dwwd = (one - properties_1.sm) * properties_1.alphad * (properties_1.xn - 
	    one) * x1 * pow_dd(&d__1, &properties_1.xm);
    if (local_1.esat == 0.f) {
	local_1.dww = dwwd;
    } else {
	csatw = local_1.sraw / (one - local_1.rsw);
/* Computing MAX */
	d__1 = zero, d__2 = one - csatw;
	x7 = max(d__1,d__2);
	local_1.dww = x7 * dwwd;
    }
/*     Permeability terms for nonhysteretic/fluid entrapment options */
    if (local_1.esat == 0.f) {
	csatw = 0.;
	p3 = 0.;
    } else {
	csatw = local_1.sraw / (one - local_1.rsw);
	d__1 = one - pow_dd(&local_1.rsw, &properties_1.xxm);
	p3 = csatw * pow_dd(&d__1, &properties_1.xm);
    }
    d__1 = one - pow_dd(&local_1.asw, &properties_1.xxm);
    p0 = one - (one - csatw) * pow_dd(&d__1, &properties_1.xm);
    d__1 = p0 - p3;
    local_1.permw = pow_dd(&local_1.esw, &c_b89) * pow_dd(&d__1, &c_b90);
    goto L110;
/*     Derivative terms for hysteresis option */
L100:
    x1 = pow_dd(&local_1.asw, &properties_1.xxm);
    d__1 = one - x1;
    local_1.dww = (one - properties_1.sm) * properties_1.alphad * (
	    properties_1.xn - one) * x1 * pow_dd(&d__1, &properties_1.xm);
/*     Permeability terms for hysteresis option */
    d__1 = one - pow_dd(&local_1.esw, &properties_1.xxm);
    p0 = one - pow_dd(&d__1, &properties_1.xm);
    local_1.permw = pow_dd(&local_1.esw, &c_b89) * pow_dd(&p0, &c_b90);
L110:
/*     Relative permeability error message */
    if (local_1.permw < 0. || local_1.permw > one) {
	s_wsfe(&io___36);
	do_fio(&c__1, "ERROR: Water Relative Permeability in Subro-\rutine D"
		"RAIN @ Node ", (ftnlen)64);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfe();
	s_stop("", (ftnlen)0);
    }
/*     End of DRAIN group */
    return 0;
} /* drain_ */

/* *********************************************************************** */
/*     Hysteretic k-S-P relations of scanning drying path */
/*     Written/revised by RJ Lenhard, Nov 2004 */
/* Subroutine */ int dry_(integer *n, doublereal *haw)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal p0, p3, x1, x7, x8, one, swd, zero, csatw, ridswd, 
	    rdiswd;

    /* Fortran I/O blocks */
    static cilist io___48 = { 0, 6, 0, "(A,I6)", 0 };


/*     Saturations of drying scanning curves for hysteresis option */
    zero = 0.;
    one = 1.;
    d__2 = properties_1.alphad * glob_1.rhsw[*n + local_1.ipsw1 * 1001 - 1002]
	    ;
    d__1 = one + pow_dd(&d__2, &properties_1.xn);
    d__3 = -properties_1.xm;
    ridswd = pow_dd(&d__1, &d__3);
    d__2 = properties_1.alphad * glob_1.rhsw[*n + (local_1.ipsw1 - 1) * 1001 
	    - 1002];
    d__1 = one + pow_dd(&d__2, &properties_1.xn);
    d__3 = -properties_1.xm;
    rdiswd = pow_dd(&d__1, &d__3);
    d__2 = properties_1.alphad * *haw;
    d__1 = one + pow_dd(&d__2, &properties_1.xn);
    d__3 = -properties_1.xm;
    swd = pow_dd(&d__1, &d__3);
    local_1.asw = (swd - rdiswd) * (glob_1.rasw[*n + local_1.ipsw1 * 1001 - 
	    1002] - glob_1.rasw[*n + (local_1.ipsw1 - 1) * 1001 - 1002]) / (
	    ridswd - rdiswd) + glob_1.rasw[*n + (local_1.ipsw1 - 1) * 1001 - 
	    1002];
    local_1.asw = min(local_1.asw,one);
    local_1.esat = local_1.sraw * ((local_1.asw - local_1.rsw) / (one - 
	    local_1.rsw));
    local_1.esw = local_1.asw - local_1.esat;
    local_1.sw = local_1.esw * (one - properties_1.sm) + properties_1.sm;
    local_1.sat = local_1.esat * (one - properties_1.sm);
/*     Permeability terms for hysteresis option */
    d__1 = one - pow_dd(&local_1.asw, &properties_1.xxm);
    p0 = one - pow_dd(&d__1, &properties_1.xm);
    csatw = local_1.sraw / (one - local_1.rsw);
    d__1 = one - pow_dd(&local_1.rsw, &properties_1.xxm);
    d__2 = one - pow_dd(&local_1.asw, &properties_1.xxm);
    p3 = pow_dd(&d__1, &properties_1.xm) - pow_dd(&d__2, &properties_1.xm);
    p3 = csatw * p3;
    d__1 = p0 - p3;
    local_1.permw = pow_dd(&local_1.esw, &c_b89) * pow_dd(&d__1, &c_b90);
/*     Derivative terms for hysteresis option */
    x1 = pow_dd(&swd, &properties_1.xxm);
/* Computing MAX */
    d__1 = zero, d__2 = one - csatw;
    x7 = max(d__1,d__2);
    x8 = (glob_1.rasw[*n + local_1.ipsw1 * 1001 - 1002] - glob_1.rasw[*n + (
	    local_1.ipsw1 - 1) * 1001 - 1002]) / (ridswd - rdiswd);
    d__1 = one - x1;
    local_1.dww = (one - properties_1.sm) * properties_1.alphad * (
	    properties_1.xn - one) * x1 * x7 * x8 * pow_dd(&d__1, &
	    properties_1.xm);
/*     Relative permeability error message */
    if (local_1.permw < 0. || local_1.permw > one) {
	s_wsfe(&io___48);
	do_fio(&c__1, "ERROR: Water Relative Permeability in Subroutine DRY "
		"@ Node ", (ftnlen)60);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfe();
	s_stop("", (ftnlen)0);
    }
/*     End of DRY group */
    return 0;
} /* dry_ */

/* *********************************************************************** */
/*     Hysteretic k-S-P relations of scanning wetting path */
/*     Written/revised by RJ Lenhard, Nov 2004 */
/* Subroutine */ int wet_(integer *n, doublereal *haw)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal p0, p3, x1, x7, x8, one, swi, zero, csatw, ridswi, 
	    rdiswi;

    /* Fortran I/O blocks */
    static cilist io___60 = { 0, 6, 0, "(A,I6)", 0 };


/*     Saturations of wetting scanning paths for hysteresis option */
    zero = 0.;
    one = 1.;
    d__2 = properties_1.alphai * glob_1.rhsw[*n + (local_1.ipsw1 - 1) * 1001 
	    - 1002];
    d__1 = one + pow_dd(&d__2, &properties_1.xnw);
    d__3 = -properties_1.xmw;
    ridswi = pow_dd(&d__1, &d__3);
    d__2 = properties_1.alphai * glob_1.rhsw[*n + local_1.ipsw1 * 1001 - 1002]
	    ;
    d__1 = one + pow_dd(&d__2, &properties_1.xnw);
    d__3 = -properties_1.xmw;
    rdiswi = pow_dd(&d__1, &d__3);
    d__2 = properties_1.alphai * *haw;
    d__1 = one + pow_dd(&d__2, &properties_1.xnw);
    d__3 = -properties_1.xmw;
    swi = pow_dd(&d__1, &d__3);
    local_1.asw = (swi - ridswi) * (glob_1.rasw[*n + local_1.ipsw1 * 1001 - 
	    1002] - glob_1.rasw[*n + (local_1.ipsw1 - 1) * 1001 - 1002]) / (
	    rdiswi - ridswi) + glob_1.rasw[*n + (local_1.ipsw1 - 1) * 1001 - 
	    1002];
    local_1.asw = min(local_1.asw,one);
    local_1.esat = local_1.sraw * ((local_1.asw - local_1.rsw) / (one - 
	    local_1.rsw));
    local_1.esw = local_1.asw - local_1.esat;
    local_1.sw = local_1.esw * (one - properties_1.sm) + properties_1.sm;
    local_1.sat = local_1.esat * (one - properties_1.sm);
/*     Permeability terms for hysteresis option */
    d__1 = one - pow_dd(&local_1.asw, &properties_1.xxm);
    p0 = one - pow_dd(&d__1, &properties_1.xm);
    csatw = local_1.sraw / (one - local_1.rsw);
    d__1 = one - pow_dd(&local_1.rsw, &properties_1.xxm);
    d__2 = one - pow_dd(&local_1.asw, &properties_1.xxm);
    p3 = pow_dd(&d__1, &properties_1.xm) - pow_dd(&d__2, &properties_1.xm);
    p3 = csatw * p3;
    d__1 = p0 - p3;
    local_1.permw = pow_dd(&local_1.esw, &c_b89) * pow_dd(&d__1, &c_b90);
/*     Derivative terms for hysteresis option */
    x1 = pow_dd(&swi, &properties_1.xxmw);
/* Computing MAX */
    d__1 = zero, d__2 = one - csatw;
    x7 = max(d__1,d__2);
    x8 = (glob_1.rasw[*n + local_1.ipsw1 * 1001 - 1002] - glob_1.rasw[*n + (
	    local_1.ipsw1 - 1) * 1001 - 1002]) / (rdiswi - ridswi);
    d__1 = one - x1;
    local_1.dww = (one - properties_1.sm) * properties_1.alphai * (
	    properties_1.xnw - one) * x1 * x8 * x7 * pow_dd(&d__1, &
	    properties_1.xmw);
/*     Relative permeability error message */
    if (local_1.permw < zero || local_1.permw > one) {
	s_wsfe(&io___60);
	do_fio(&c__1, "ERROR: Water Relative Permeability in Subroutine WET "
		"@ Node ", (ftnlen)60);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfe();
	s_stop("", (ftnlen)0);
    }
/*     End of WET group */
    return 0;
} /* wet_ */

/* *********************************************************************** */
/*     Sets hysteresis variables following convergence */
/*     Written/revised by RJ Lenhard, Nov 2004 */
/* Subroutine */ int updatehyst_(integer *n, doublereal *hw, doublereal *ha)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer m;
    static doublereal haw, one, zero;

    zero = 0.;
    one = 1.;
/*       HW = 'converged water pressure for the node' */
/*       HA = 'converged air pressure for the node' */
/* Computing MAX */
    d__1 = zero, d__2 = *ha - *hw;
    haw = max(d__1,d__2);
    if (local_1.asw < glob_1.rswaw[*n - 1]) {
	glob_1.rswaw[*n - 1] = local_1.asw;
	if (properties_1.sarwi != zero) {
	    glob_1.sarw[*n - 1] = (one - local_1.esw) / (one + local_1.raw * (
		    one - local_1.esw));
	}
    }
    if (properties_1.ipath == 1) {
	goto L200;
    }
/*       Setting reversal points when HAW = 0 on path 2 */
    if (haw == zero && local_1.ipsw1 != 1) {
	local_1.ipsw1 = 2;
	local_1.jjjh = 1;
	i__1 = properties_1.ipath;
	for (m = local_1.ipsw1 + 3; m <= i__1; m += 2) {
	    glob_1.rasw[*n + m * 1001 - 1002] = zero;
/* L102: */
	}
	i__1 = properties_1.ipath;
	for (m = local_1.ipsw1 + 2; m <= i__1; m += 2) {
	    glob_1.rasw[*n + m * 1001 - 1002] = one;
/* L103: */
	}
	glob_1.rhsw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = haw;
	glob_1.rasw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = local_1.asw;
	glob_1.phsw[*n - 1] = haw;
    }
/*       A drying scanning path has closed - resetting reversal points */
    if (local_1.ilsw > 0 && local_1.jjjh == 0) {
	i__1 = properties_1.ipath;
	for (m = local_1.ipsw1 + 3; m <= i__1; m += 2) {
	    glob_1.rasw[*n + m * 1001 - 1002] = one;
/* L120: */
	}
	i__1 = properties_1.ipath;
	for (m = local_1.ipsw1 + 2; m <= i__1; m += 2) {
	    glob_1.rasw[*n + m * 1001 - 1002] = zero;
/* L121: */
	}
	glob_1.rhsw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = haw;
	glob_1.rasw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = local_1.asw;
	glob_1.phsw[*n - 1] = haw;
/*       A wetting scanning path has closed - resetting reversal points */
    } else if (local_1.ilsw > 0 && local_1.jjjh == 1) {
	i__1 = properties_1.ipath;
	for (m = local_1.ipsw1 + 3; m <= i__1; m += 2) {
	    glob_1.rasw[*n + m * 1001 - 1002] = zero;
/* L122: */
	}
	i__1 = properties_1.ipath;
	for (m = local_1.ipsw1 + 2; m <= i__1; m += 2) {
	    glob_1.rasw[*n + m * 1001 - 1002] = one;
/* L123: */
	}
	glob_1.rhsw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = haw;
	glob_1.rasw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = local_1.asw;
	glob_1.phsw[*n - 1] = haw;
    }
/*       Drying scanning path closed with main drainage - resetting reversals */
    if (haw >= glob_1.rhsw[*n + 1000] && glob_1.rasw[*n + 2001] > zero) {
	local_1.ipsw1 = 1;
	local_1.jjjh = 0;
	i__1 = properties_1.ipath;
	for (m = local_1.ipsw1 + 2; m <= i__1; m += 2) {
	    glob_1.rasw[*n + m * 1001 - 1002] = zero;
/* L124: */
	}
	i__1 = properties_1.ipath;
	for (m = local_1.ipsw1 + 3; m <= i__1; m += 2) {
	    glob_1.rasw[*n + m * 1001 - 1002] = one;
/* L125: */
	}
	glob_1.rhsw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = haw;
	glob_1.rasw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = local_1.asw;
	glob_1.phsw[*n - 1] = haw;
    }
/*       Setting gobal variables from converged local variables */
    glob_1.jjh[*n - 1] = local_1.jjjh;
    glob_1.ipsw[*n - 1] = local_1.ipsw1;
    glob_1.mpsw[*n - 1] = local_1.mpsw1;
    if (local_1.ipsw1 == properties_1.ipath) {
	goto L200;
    }
/*       No scanning paths closed, only continued */
    if (local_1.jjjh == 1 && local_1.asw > glob_1.rasw[*n + (local_1.ipsw1 + 
	    1) * 1001 - 1002]) {
	glob_1.rhsw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = haw;
	glob_1.rasw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = local_1.asw;
	glob_1.phsw[*n - 1] = haw;
    }
    if (local_1.jjjh == 0 && local_1.asw < glob_1.rasw[*n + (local_1.ipsw1 + 
	    1) * 1001 - 1002]) {
	glob_1.rhsw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = haw;
	glob_1.rasw[*n + (local_1.ipsw1 + 1) * 1001 - 1002] = local_1.asw;
	glob_1.phsw[*n - 1] = haw;
    }
    if (local_1.ipsw1 == 1) {
	glob_1.phsw[*n - 1] = glob_1.rhsw[*n + 1000];
    }
L200:
    return 0;
} /* updatehyst_ */

