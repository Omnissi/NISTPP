#pragma once

/*	Real Fast Fourier Transform

void __ogg_fdrffti(int n, double *wsave, int *ifac)
void __ogg_fdrfftf(int n, double *r, double *wsave, int *ifac)
void __ogg_fdrfftb(int n, double *r, double *wsave, int *ifac)

__ogg_fdrffti == initialization
__ogg_fdrfftf == forward transform
__ogg_fdrfftb == backward transform

Parameters are
n == length of sequence
r == sequence to be transformed (input)
== transformed sequence (output)
wsave == work array of length 2n (allocated by caller)
ifac == work array of length 15 (allocated by caller)

	Cosine quarter-wave FFT

void __ogg_fdcosqi(int n, double *wsave, int *ifac)
void __ogg_fdcosqf(int n, double *x, double *wsave, int *ifac)
void __ogg_fdcosqb(int n, double *x, double *wsave, int *ifac)
*/

// ========= Real FFT ====================
void __ogg_fdrffti(int n, double *wsave, int *ifac);
void __ogg_fdrfftf(int n, double *X, double *wsave, int *ifac);
void __ogg_fdrfftb(int n, double *r, double *wsave, int *ifac);

// ========= Cosine quarter-wave FFT=======
void __ogg_fdcosqi(int n, double *wsave, int *ifac);
void __ogg_fdcosqf(int n, double *x, double *wsave, int *ifac);
void __ogg_fdcosqb(int n, double *x, double *wsave, int *ifac);
