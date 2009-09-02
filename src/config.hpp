#ifndef CONFIG_H
#define CONFIG_H 1


/* ATLAS */
#ifdef BLAS_ATLAS
#define USE_BLAS 1
#define BLAS_HEADER <cblas.h>
#define BLAS(FUNC) cblas_ ## FUNC
#endif

/* ACML */
#ifdef BLAS_ACML
#define USE_BLAS 1
#define BLAS_HEADER <acml.h>
#define BLAS(FUNC) FUNC
#endif



#endif
















