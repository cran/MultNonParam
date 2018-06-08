#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(signtestperm)(void *, void *, void *, void *, void *);
extern void F77_NAME(aovp)(void *,void *,void *, void *, void *, void *, void *);
extern void F77_NAME(betatestf)(void *,void *,void *, void *);
extern void F77_NAME(probestf)( void *,void *,void *, void *, void *,void *,void *, void *, void *,void *,void *, void *, void *,void *,void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"signtestperm", (DL_FUNC) &F77_NAME(signtestperm), 5},
    {"aovp", (DL_FUNC) &F77_NAME(aovp), 7},
    {"betatestf", (DL_FUNC) &F77_NAME(betatestf), 4},
    {"probestf", (DL_FUNC) &F77_NAME(probestf), 16},
    {NULL, NULL, 0}
};

void R_init_MultNonParam(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
