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
    unsigned int i, j;
    for (i = 0, j = len - 1; i < len; i++, j--) {
        float s0 = src0[i];
        float s1 = src1[len - i - 1];
        float wi = win[i];
        float wj = win[len + j];
        dst[i]       = s0 * wj - s1 * wi;
        dst[len + j] = s0 * wi + s1 * wj;
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
    z[64] = z[0];
    z[65] = z[1];
    for (k = 1; k < 32; k++) {
        z[64+2*k  ] = -z[64 - k];
        z[64+2*k+1] =  z[ k + 1];
    }
}

FDSP_EXPORT void sbr_qmf_post_shuffle_c(FFTComplex W[32], float *z)
{
    unsigned int k;
    for (k = 0; k < 32; k++) {
        W[k].re = -z[k];
        W[k].im =  z[63-k];
    }
}

FDSP_EXPORT void sbr_qmf_deint_bfly_c(float *v, const float *src0, const float *src1)
{
    unsigned int i;
    for (i = 0; i < 64; i++) {
        v[      i] = src1[63 - i] - src0[i];
        v[127 - i] = src1[63 - i] + src0[i];
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
        X_high[i].re =
            X_low[i - 2].re * alpha[0] -
            X_low[i - 2].im * alpha[1] +
            X_low[i - 1].re * alpha[2] -
            X_low[i - 1].im * alpha[3] +
            X_low[i].re;
        X_high[i].im =
            X_low[i - 2].im * alpha[0] +
            X_low[i - 2].re * alpha[1] +
            X_low[i - 1].im * alpha[2] +
            X_low[i - 1].re * alpha[3] +
            X_low[i].im;
    }
}

FDSP_EXPORT void sbr_qmf_synthesis_window_c(float *out, float *v, float *sbr_qmf_window, unsigned int n)
{
    for (n = 0; n < 64; n++) {
        float t  = v[n]*sbr_qmf_window[n];
              t += v[n+192]*sbr_qmf_window[n+64];
              t += v[n+256]*sbr_qmf_window[n+128];
              t += v[n+448]*sbr_qmf_window[n+192];
              t += v[n+512]*sbr_qmf_window[n+256];
              t += v[n+704]*sbr_qmf_window[n+320];
              t += v[n+768]*sbr_qmf_window[n+384];
              t += v[n+960]*sbr_qmf_window[n+448];
              t += v[n+1024]*sbr_qmf_window[n+512];
              t += v[n+1216]*sbr_qmf_window[n+576];
        out[n] = t;
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

