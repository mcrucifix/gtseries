#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Declare the C functions from fmft.c */
extern int fmft(int *localnfreq, double *localminfreq, double *localmaxfreq,
                int *localflag, int *localndata, double *localxdata, double *localydata,
                struct component *signal1, struct component *signal2, struct component *signal3);

extern int fastgsec(double *localminfreq, double *localmaxfreq,
                    int *localndata, double *localxdata, double *localydata, double *outfreq);

/* Register routines */

static const R_CMethodDef CEntries[] = {
    {"fmft",    (DL_FUNC) &fmft,    10},
    {"fastgsec",(DL_FUNC) &fastgsec, 6},
    {NULL, NULL, 0}
};

void R_init_gtseries(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
