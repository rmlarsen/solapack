/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#ifndef SOLAPACK_MAXSIZE_H
#define SOLAPACK_MAXSIZE_H

#define MAXRAD    (103)    /* Max radial grid points */
#define MAXTHETA  (53)     /* Max co-latitude grid points */
#define MAXPTS    (MAXRAD*MAXTHETA)  /* Total grid points */
#define MAXNLM    (45001)  /* Max number of splittings */
#define MAXNL     (5001)   /* Max number of multiplets (n,l) */
#define MAXL      (301)    /* Max angular degree l */
#define MAXLM     ((MAXL+1)*(MAXL+2)/2)  /* Max (l,m) pairs */
#define MAXEXP    (36)     /* Max expansion coefficients */
#define MAXITER   (2000)   /* Max Lanczos iterations */

#endif
