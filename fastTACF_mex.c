/* fastTACF_mex.c
   MEX: compute 3x3 ROI AutoCorr, fit f(t)= a/(1+t/b)+c by fast 1D search,
   return Lifetime (b), D, and n maps.

   Compile with:
     mex -O fastTACF_mex.c
*/

#include "mex.h"
#include <math.h>
#include <string.h>

/* Helper macro for indexing */
#define IDX3(r,c,t,rows,cols) ((r) + (c)*(rows) + (t)*(rows)*(cols))

/* Entry point */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* --- Input checks --- */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("fastTACF:usage", "Usage: [Lifetime, D, n] = fastTACF_mex(data, Tau, C)");
    }
    if (nlhs < 1) {
        mexErrMsgIdAndTxt("fastTACF:usage", "Need at least one output.");
    }

    /* Retrieve inputs */
    const mxArray *mxData = prhs[0];
    const mxArray *mxTau  = prhs[1];
    const mxArray *mxC    = prhs[2];

    if (!mxIsDouble(mxData) || mxIsComplex(mxData) || mxGetNumberOfDimensions(mxData) != 3) {
        mexErrMsgIdAndTxt("fastTACF:input", "data must be a real double 3-D array.");
    }
    if (!mxIsDouble(mxTau) || mxIsComplex(mxTau)) {
        mexErrMsgIdAndTxt("fastTACF:input", "Tau must be real double vector.");
    }
    if (!mxIsDouble(mxC) || mxIsComplex(mxC) || mxGetNumberOfElements(mxC) != 1) {
        mexErrMsgIdAndTxt("fastTACF:input", "C must be a real scalar double.");
    }

    mwSize rows  = mxGetDimensions(mxData)[0];
    mwSize cols  = mxGetDimensions(mxData)[1];
    mwSize frames= mxGetDimensions(mxData)[2];

    double *data = mxGetPr(mxData);
    double *Tau  = mxGetPr(mxTau);
    double C     = mxGetScalar(mxC);

    /* ensure Tau length matches frames */
    mwSize TauLen = mxGetNumberOfElements(mxTau);
    if (TauLen != frames) {
        mexErrMsgIdAndTxt("fastTACF:input", "Length of Tau must equal the number of frames in data.");
    }

    /* Prepare outputs: Lifetime, D, n (all rows x cols) */
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double *outLifetime = mxGetPr(plhs[0]);

    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    } else {
        plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL); /* placeholder */
    }
    double *outD = mxGetPr(plhs[1]);

    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    } else {
        plhs[2] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
    double *outn = mxGetPr(plhs[2]);

    /* Temporary buffers (allocated once) */
    double *Y = (double*) mxCalloc(frames, sizeof(double));   /* AutoCorr vector */
    double *R = (double*) mxCalloc(frames, sizeof(double));   /* regressor for given b */

    mwSize planeSize = rows * cols;

    /* Constants for n calculation */
    const double kB = 1.380649e-23;
    const double Tkelvin = 296.15;
    const double viscosity_term = 6.0 * M_PI * 20e-9; /* 20 nm radius? matches your formula */

    /* Precompute some Tau related safe values */
    double Tau_first = (frames >= 2) ? Tau[1] : ( (frames>=1) ? Tau[0] : 1.0 );
    double Tau_last  = Tau[frames-1];

    /* Golden-section constants */
    const double phi = (sqrt(5.0) - 1.0) / 2.0;

    /* Loop all pixels */
    mwSize r, c;
    for (c = 0; c < cols; ++c) {
        for (r = 0; r < rows; ++r) {

            /* compute 3x3 ROI bounds with clamping */
            mwSize r0 = (r >= 1) ? (r-1) : 0;
            mwSize r1 = (r+1 < rows) ? (r+1) : (rows-1);
            mwSize c0 = (c >= 1) ? (c-1) : 0;
            mwSize c1 = (c+1 < cols) ? (c+1) : (cols-1);

            mwSize num_pixels = (r1 - r0 + 1) * (c1 - c0 + 1);

            /* Zero Y buffer */
            {
                mwSize t;
                for (t = 0; t < frames; ++t) Y[t] = 0.0;
            }

            /* Sum temporal traces of pixels in ROI into Y */
            mwSize rr, cc;
            for (cc = c0; cc <= c1; ++cc) {
                for (rr = r0; rr <= r1; ++rr) {
                    mwSize base = rr + cc*rows; /* index of pixel at frame 0 */
                    /* add that pixel's continuous time series */
                    mwSize t;
                    const double *p = data + base;
                    for (t = 0; t < frames; ++t) {
                        Y[t] += p[t * planeSize];
                    }
                }
            }

            /* divide by Numel -> get mean time trace */
            {
                mwSize t;
                double invNum = 1.0 / (double)num_pixels;
                for (t = 0; t < frames; ++t) Y[t] *= invNum;
            }

            /* Quick check: if Y is constant or too small variance, return zeros */
            double meanY = 0.0;
            double varY = 0.0;
            {
                mwSize t;
                for (t=0; t < frames; ++t) meanY += Y[t];
                meanY /= (double)frames;
                for (t=0; t < frames; ++t) {
                    double d = Y[t] - meanY;
                    varY += d*d;
                }
            }
            if (varY < 1e-20) {
                outLifetime[r + c*rows] = 0.0;
                if (outD) outD[r + c*rows] = 0.0;
                if (outn) outn[r + c*rows] = 0.0;
                continue;
            }

            /* Setup golden-section search interval for b */
            double bL = 0.1 * Tau_first;
            if (!(bL > 0.0)) bL = 1e-6; /* fallback */
            double bR = 10.0 * Tau_last;
            if (!(bR > bL)) bR = bL * 100.0;

            /* evaluate error at two interior points c1, c2 */
            double aL = bL;
            double aR = bR;
            double c1 = aR - phi*(aR - aL);
            double c2 = aL + phi*(aR - aL);

            /* Inline function: evaluate residual sum of squares for a given b */
            auto eval_err = [&](double b)->double {
                if (!(b > 0.0)) return 1e300;
                mwSize t;
                double sR = 0.0, sRR = 0.0, sY = 0.0, sRY = 0.0;
                for (t = 0; t < frames; ++t) {
                    double rVal = 1.0 / (1.0 + Tau[t] / b);
                    R[t] = rVal;
                    sR += rVal;
                    sRR += rVal * rVal;
                    sY += Y[t];
                    sRY += rVal * Y[t];
                }
                double det = (double)frames * sRR - sR * sR;
                if (fabs(det) < 1e-30) return 1e300; /* degenerate */
                double c_hat = (sRR * sY - sR * sRY) / det;
                double a_hat = ((double)frames * sRY - sR * sY) / det;
                double rss = 0.0;
                for (t = 0; t < frames; ++t) {
                    double f = a_hat * R[t] + c_hat;
                    double d = Y[t] - f;
                    rss += d*d;
                }
                return rss;
            };

            double f1 = eval_err(c1);
            double f2 = eval_err(c2);

            /* iterate golden-section search */
            int iter;
            const int MAX_ITER = 25; /* enough for double precision */
            for (iter = 0; iter < MAX_ITER; ++iter) {
                if (f1 > f2) {
                    aL = c1;
                    c1 = c2;
                    f1 = f2;
                    c2 = aL + phi*(aR - aL);
                    f2 = eval_err(c2);
                } else {
                    aR = c2;
                    c2 = c1;
                    f2 = f1;
                    c1 = aR - phi*(aR - aL);
                    f1 = eval_err(c1);
                }
            }

            double b_est = 0.5 * (aL + aR);

            /* Save Lifetime */
            outLifetime[r + c*rows] = b_est;

            /* compute D and n if requested */
            double Dval = 0.0;
            double nval = 0.0;
            if (b_est > 0.0) {
                Dval = sqrt(C) / (4.0 * b_est);
                /* n = (kB*T) / (6*pi*eta*D*1e-12) * 1e3;  matches your earlier formula */
                double D_si = Dval * 1e-12; /* convert to m^2/s ? previously you multiplied/divided: follow your formula */
                double denom = viscosity_term * D_si;
                if (denom <= 0.0) nval = 0.0;
                else nval = (kB * Tkelvin) / denom * 1e3;
            }

            if (outD) outD[r + c*rows] = Dval;
            if (outn) outn[r + c*rows] = nval;
        }
    }

    /* free temporaries */
    mxFree(Y);
    mxFree(R);
}
