#define BUILDING_FDSP 1
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fdsp.h"

typedef union _f32 {
    float f;
    uint32_t i;
} _f32;

FDSP_EXPORT void vector_fmul_c(float *dst, const float *src0, const float *src1, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++)
        dst[i] = src0[i] * src1[i];
}

FDSP_EXPORT void vector_fmac_scalar_c(float *dst, const float *src, float mul, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++)
        dst[i] += src[i] * mul;
}

FDSP_EXPORT void vector_fmul_scalar_c(float *dst, const float *src, float mul, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++)
        dst[i] = src[i] * mul;
}

FDSP_EXPORT void vector_fmul_window_c(float *dst, const float *src0, const float *src1, const float *win, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        float s0 = src0[i];
        float s1 = src1[len - i - 1];
        float wi = win[i];
        float wj = win[2*len - i - 1];
        dst[i]             = s0 * wj - s1 * wi;
        dst[2*len - i - 1] = s0 * wi + s1 * wj;
    }
}

FDSP_EXPORT void vector_fmul_copy_c(float *dst, const float *src, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++)
        dst[i] = src[i];
}

FDSP_EXPORT void vector_fmul_add_c(float *dst, const float *src0, const float *src1, const float *src2, unsigned int len){
    unsigned int i;
    for (i = 0; i < len; i++)
        dst[i] = src0[i] * src1[i] + src2[i];
}

FDSP_EXPORT void vector_fmul_reverse_c(float *dst, const float *src0, const float *src1, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++)
        dst[i] = src0[i] * src1[len - i - 1];
}

FDSP_EXPORT void butterflies_float_c(float *__restrict v1, float *__restrict v2, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        float s = v1[i] + v2[i];
        float t = v1[i] - v2[i];
        v1[i] = s;
        v2[i] = t;
    }
}

FDSP_EXPORT float scalarproduct_float_c(const float *v1, const float *v2, unsigned int len)
{
    float p = 0.0;
    unsigned int i;
    for (i = 0; i < len; i++) {
        p += v1[i] * v2[i];
    }
    return p;
}

FDSP_EXPORT float scalarproduct_symmetric_fir_float_c(const float *v1, const float *v2, unsigned int len)
{
    float p = 0.0;
    unsigned int i;
    for (i = 0; i < len; i++) {
        float t = v1[i] + v1[2*len - i];
        p += t * v2[i];
    }
    p += v1[len] * v2[len];
    return p;
}

FDSP_EXPORT void vector_clipf_c(float *dst, const float *src,
                                float min, float max, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        float f = src[i];
        if (f < min) f = min;
        if (f > max) f = max;
        dst[i] = f;
    }
}

FDSP_EXPORT void sbr_sum64x5_c(float *z)
{
    unsigned int k;
    for (k = 0; k < 64; k++) {
        float f = z[k] + z[k + 64] + z[k + 128] + z[k + 192] + z[k + 256];
        z[k] = f;
    }
}

FDSP_EXPORT void sbr_qmf_pre_shuffle_c(float *z)
{
    unsigned int k;
    float tmp1 = z[0], tmp2 = z[1];
    //z[64] = z[0];
    //z[65] = z[1];
    for (k = 0; k < 32; k++) {
        z[64+2*k  ] = -z[64 - k];
        z[64+2*k+1] =  z[ k + 1];
    }
    z[64] = tmp1;
}

FDSP_EXPORT void sbr_qmf_post_shuffle_c(FFTComplex W[32], float *z)
{
    unsigned int k;
    for (k = 0; k < 32; k++) {
        W[k].re = z[k];
        W[k].im = z[63-k];
    }
}

FDSP_EXPORT void sbr_qmf_deint_bfly_c(float *v, const float *src0, const float *src1)
{
    unsigned int i;
    for (i = 0; i < 64; i++) {
        v[      i] = src1[63 - i] + src0[i];
        v[127 - i] = src1[63 - i] - src0[i];
    }
}

FDSP_EXPORT void sbr_hf_g_filt_c(FFTComplex *Y, FFTComplex (*X_high)[40],
                     const float *g_filt, size_t m_max, size_t ixh)
{
    size_t m;
    for (m = 0; m < m_max; m++) {
        Y[m].re = X_high[m][ixh].re * g_filt[m];
        Y[m].im = X_high[m][ixh].im * g_filt[m];
    }
}

FDSP_EXPORT void sbr_hf_gen_c(FFTComplex *X_high, FFTComplex *X_low,
                              float alpha[4], unsigned int start, unsigned int end)
{
    int i;

    for (i = start; i < end; i++) {
        X_high[i].re = X_low[i].re +
            X_low[i - 2].re * alpha[0] + X_low[i - 2].im * alpha[1] +
            X_low[i - 1].re * alpha[2] + X_low[i - 1].im * alpha[3];
        X_high[i].im = X_low[i].im +
            X_low[i - 2].im * alpha[0] - X_low[i - 2].re * alpha[1] +
            X_low[i - 1].im * alpha[2] - X_low[i - 1].re * alpha[3];
    }
}

FDSP_EXPORT void sbr_qmf_synthesis_window_c(float *out, float *v, float *sbr_qmf_window, unsigned int k)
{
    const size_t _k = ((k == 32) ? 32 : 64);
    unsigned int n, mask = (k-1);
    for (n = 0; n < _k; n++) {
        out[n] = v[n]*sbr_qmf_window[n] + v[n+19*_k]*sbr_qmf_window[n+9*_k];
    }
    for (n = 0; n < 2*_k; n++) {
        float t  = v[n+ 3*_k]*sbr_qmf_window[n+  _k];
              t += v[n+ 7*_k]*sbr_qmf_window[n+3*_k];
              t += v[n+11*_k]*sbr_qmf_window[n+5*_k];
              t += v[n+15*_k]*sbr_qmf_window[n+7*_k];
        out[n&mask] += t;
    }
}

FDSP_EXPORT void sbrenc_sum128x5_c(float *z)
{
    unsigned int k;
    for (k = 0; k < 128; k++) {
        float f = z[k] + z[k + 128] + z[k + 256] + z[k + 384] + z[k + 512];
        z[k] = f;
    }
}

FDSP_EXPORT void sbrenc_qmf_deint_bfly_c(float *v, const float *src)
{
    unsigned int i;
    for (i = 0; i < 64; i++) {
        v[      i] = -src[127 - i] + src[i];
        v[127 - i] =  src[127 - i] + src[i];
    }
}

FDSP_EXPORT void aacenc_calc_expspec_c(float *expspec, float *mdct_spectrum, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        _f32 tmp;
        uint32_t ix;
        float x, y, t;

        x = mdct_spectrum[i];
        tmp.f = x;
        tmp.i &= 0x7fffffff;
        x = tmp.f;
        ix = tmp.i;

        ix  = 0x4f58cae5 - (ix >> 2);
        tmp.i = ix;
        y = tmp.f;
        t = y * y;
        y  = y * 0.25f * (5.0f - (x * t * t));   // 1st iteration
        t = y * y;
        y  = y * 0.25f * (5.0f - (x * t * t));   // 2nd iteration
        t = y * y;
        y  = y * 0.25f * (5.0f - (x * t * t));   // 3rd iteration
        y = y * x;
        expspec[i] = y;
    }
}

FDSP_EXPORT void vorbis_inverse_coupling_c(float *mag, float *ang, unsigned int blocksize)
{
    unsigned int i;
    for (i = 0;  i < blocksize; i++) {
        if (mag[i] > 0.0f) {
            if (ang[i] > 0.0f) {
                ang[i] = mag[i] - ang[i];
            } else {
                float temp = ang[i];
                ang[i]     = mag[i];
                mag[i]    += temp;
            }
        } else {
            if (ang[i] > 0.0f) {
                ang[i] = mag[i] + ang[i];
            } else {
                float temp = ang[i];
                ang[i]     = mag[i];
                mag[i]    -= temp;
            }
        }
    }
}

void sbr_qmf_deint_neg_c(float *v, const float *src)
{
    unsigned int i;
    for (i = 0; i < 32; i++) {
        v[     i] =  src[63 - 2*i    ];
        v[63 - i] = -src[63 - 2*i - 1];
    }
}

void sbr_autocorrelate_c(const FFTComplex x[40], float phi[5], unsigned int ac_len)
{
    float real_sum2 = 0.0f, imag_sum2 = 0.0f;
    float real_sum1 = 0.0f, imag_sum1 = 0.0f, real_sum0 = 0.0f;
    unsigned int i;
    for (i = 1; i < ac_len; i++) {
        real_sum0 +=  x[i].re * x[i  ].re + x[i].im * x[i  ].im;
        real_sum1 +=  x[i].re * x[i+1].re + x[i].im * x[i+1].im;
        imag_sum1 += -x[i].re * x[i+1].im + x[i].im * x[i+1].re;
        real_sum2 +=  x[i].re * x[i+2].re + x[i].im * x[i+2].im;
        imag_sum2 += -x[i].re * x[i+2].im + x[i].im * x[i+2].re;
    }
    phi[2] = real_sum2;
    phi[3] = imag_sum2;
    phi[4] = real_sum0;
    phi[0] = real_sum1;
    phi[1] = imag_sum1;
}

FDSP_EXPORT void conv_fltp_to_flt_2ch_c(float *dst, float *src[2], unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        dst[2*i] = src[0][i];
        dst[2*i+1] = src[1][i];
    }
}

FDSP_EXPORT void conv_flt_to_fltp_2ch_c(float *dst[2], float *src, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        dst[0][i] = src[2*i];
        dst[1][i] = src[2*i+1];
    }
}

FDSP_EXPORT void conv_s16p_to_s16_2ch_c(short *dst, short *src[2], unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        dst[2*i] = src[0][i];
        dst[2*i+1] = src[1][i];
    }
}

FDSP_EXPORT void conv_s16_to_s16p_2ch_c(short *dst[2], short *src, unsigned int len)
{
    unsigned int i;
    for (i = 0; i < len; i++) {
        dst[0][i] = src[2*i];
        dst[1][i] = src[2*i+1];
    }
}

