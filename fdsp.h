#ifndef LIBFDSP_H
#define LIBFDSP_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef FFTCOMPLEX_T_DEFINED
typedef struct FFTComplex {
    float re, im;
} FFTComplex;
#define FFTCOMPLEX_T_DEFINED
#endif

#if (defined(_WIN32) || defined(_WIN64)) && defined(FDSP_DLL)
#ifdef BUILDING_FDSP
#define FDSP_EXPORT __declspec(dllexport)
#else
#define FDSP_EXPORT __declspec(dllimport)
#endif
#else
#define FDSP_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * DSP utility functions.
 */
FDSP_EXPORT void vector_fmul_window_c(float *dst, const float *src0, const float *src1, const float *win, uint32_t len);
FDSP_EXPORT void vector_fmul_c(float *dst, const float *src0, const float *src1, uint32_t len);
FDSP_EXPORT void vector_fmul_scalar_c(float *dst, const float *src, float mul, uint32_t len);
FDSP_EXPORT void vector_fmac_scalar_c(float *dst, const float *src, float mul, uint32_t len);
FDSP_EXPORT void vector_fmul_reverse_c(float *dst, const float *src0, const float *src1, uint32_t len);
FDSP_EXPORT void vector_fmul_add_c(float *dst, const float *src0, const float *src1,
                       const float *src2, uint32_t len);
FDSP_EXPORT void vector_fmul_copy_c(float *dst, const float *src, uint32_t len);
FDSP_EXPORT float scalarproduct_float_c(const float *v1, const float *v2, uint32_t len);
FDSP_EXPORT void butterflies_float_c(float *src0, float *src1, uint32_t len);

FDSP_EXPORT void vector_fmul_window_sse(float *dst, const float *src0, const float *src1, const float *win, uint32_t len);
FDSP_EXPORT void vector_fmul_sse(float *dst, const float *src0, const float *src1, uint32_t len);
FDSP_EXPORT void vector_fmul_scalar_sse(float *dst, const float *src, float mul, uint32_t len);
FDSP_EXPORT void vector_fmac_scalar_sse(float *dst, const float *src, float mul, uint32_t len);
FDSP_EXPORT void vector_fmul_reverse_sse(float *dst, const float *src0, const float *src1, uint32_t len);
FDSP_EXPORT void vector_fmul_add_sse(float *dst, const float *src0, const float *src1,
                         const float *src2, uint32_t len);
FDSP_EXPORT void vector_fmul_copy_sse(float *dst, const float *src, uint32_t len);
FDSP_EXPORT float scalarproduct_float_sse(const float *v1, const float *v2, uint32_t len);
FDSP_EXPORT void butterflies_float_sse(float *src0, float *src1, uint32_t len);

FDSP_EXPORT void vector_fmul_window_neon(float *dst, const float *src0, const float *src1, const float *win, uint32_t len);
FDSP_EXPORT void vector_fmul_neon(float *dst, const float *src0, const float *src1, uint32_t len);
FDSP_EXPORT void vector_fmul_scalar_neon(float *dst, const float *src, float mul, uint32_t len);
FDSP_EXPORT void vector_fmac_scalar_neon(float *dst, const float *src, float mul, uint32_t len);
FDSP_EXPORT void vector_fmul_reverse_neon(float *dst, const float *src0, const float *src1, uint32_t len);
FDSP_EXPORT void vector_fmul_add_neon(float *dst, const float *src0, const float *src1,
                          const float *src2, uint32_t len);
FDSP_EXPORT void vector_fmul_copy_neon(float *dst, const float *src, uint32_t len);
FDSP_EXPORT float scalarproduct_float_neon(const float *v1, const float *v2, uint32_t len);
FDSP_EXPORT void butterflies_float_neon(float *src0, float *src1, uint32_t len);

FDSP_EXPORT void aacenc_calc_expspec_c(float *expspec, float *mdct_spectrum, unsigned int len);
FDSP_EXPORT void vorbis_inverse_coupling_c(float *mag, float *ang, unsigned int blocksize);

FDSP_EXPORT void aacenc_calc_expspec_sse(float *expspec, float *mdct_spectrum, unsigned int len);
FDSP_EXPORT void vorbis_inverse_coupling_sse(float *mag, float *ang, unsigned int blocksize);

FDSP_EXPORT void aacenc_calc_expspec_neon(float *expspec, float *mdct_spectrum, unsigned int len);
FDSP_EXPORT void vorbis_inverse_coupling_neon(float *mag, float *ang, unsigned int blocksize);

FDSP_EXPORT void sbr_sum64x5_c(float *z);
FDSP_EXPORT void sbr_qmf_pre_shuffle_c(float *z);
FDSP_EXPORT void sbr_qmf_post_shuffle_c(FFTComplex W[32], float *z);
FDSP_EXPORT void sbr_qmf_deint_bfly_c(float *v, const float *src0, const float *src1);
FDSP_EXPORT void sbr_hf_g_filt_c(FFTComplex *Y, FFTComplex (*X_high)[40],
                     const float *g_filt, size_t m_max, size_t ixh);
FDSP_EXPORT void sbrenc_sum128x5_c(float *z);
FDSP_EXPORT void sbrenc_qmf_deint_bfly_c(float *v, const float *src);
FDSP_EXPORT void sbr_hf_gen_c(FFTComplex *X_high, FFTComplex *X_low, float alpha[4], uint32_t start, uint32_t end);
FDSP_EXPORT void sbr_qmf_synthesis_window_c(float *out, float *v, float *sbr_qmf_window, uint32_t n);
FDSP_EXPORT void sbr_qmf_deint_neg_c(float *v, const float *src);
FDSP_EXPORT void sbr_autocorrelate_c(const FFTComplex x[40], float acorr_sums[5]);

FDSP_EXPORT void sbr_sum64x5_sse(float *z);
FDSP_EXPORT void sbr_qmf_pre_shuffle_sse(float *z);
FDSP_EXPORT void sbr_qmf_post_shuffle_sse(FFTComplex W[32], float *z);
FDSP_EXPORT void sbr_qmf_deint_bfly_sse(float *v, const float *src0, const float *src1);
FDSP_EXPORT void sbr_hf_g_filt_sse(FFTComplex *Y, FFTComplex (*X_high)[40],
                       const float *g_filt, size_t m_max, size_t ixh);
FDSP_EXPORT void sbrenc_sum128x5_sse(float *z);
FDSP_EXPORT void sbrenc_qmf_deint_bfly_sse(float *v, const float *src);
FDSP_EXPORT void sbr_hf_gen_sse(FFTComplex *X_high, FFTComplex *X_low, float alpha[4], uint32_t start, uint32_t end);
FDSP_EXPORT void sbr_qmf_synthesis_window_sse(float *out, float *v, float *sbr_qmf_window, uint32_t n);
FDSP_EXPORT void sbr_qmf_deint_neg_sse(float *v, const float *src);
FDSP_EXPORT void sbr_autocorrelate_sse(const FFTComplex x[40], float acorr_sums[5]);

FDSP_EXPORT void sbr_sum64x5_neon(float *z);
FDSP_EXPORT void sbr_qmf_pre_shuffle_neon(float *z);
FDSP_EXPORT void sbr_qmf_post_shuffle_neon(FFTComplex W[32], float *z);
FDSP_EXPORT void sbr_qmf_deint_bfly_neon(float *v, const float *src0, const float *src1);
FDSP_EXPORT void sbr_hf_g_filt_neon(FFTComplex *Y, FFTComplex (*X_high)[40],
                        const float *g_filt, size_t m_max, size_t ixh);
FDSP_EXPORT void sbrenc_sum128x5_neon(float *z);
FDSP_EXPORT void sbrenc_qmf_deint_bfly_neon(float *v, const float *src);
FDSP_EXPORT void sbr_hf_gen_neon(FFTComplex *X_high, FFTComplex *X_low, float alpha[4], uint32_t start, uint32_t end);
FDSP_EXPORT void sbr_qmf_synthesis_window_neon(float *out, float *v, float *sbr_qmf_window, uint32_t n);
FDSP_EXPORT void sbr_qmf_deint_neg_neon(float *v, const float *src);
FDSP_EXPORT void sbr_autocorrelate_neon(const FFTComplex x[40], float acorr_sums[5]);

FDSP_EXPORT void conv_fltp_to_flt_2ch_c(float *dst, float *src[2], unsigned int len);
FDSP_EXPORT void conv_flt_to_fltp_2ch_c(float *dst[2], float *src, unsigned int len);
FDSP_EXPORT void conv_s16p_to_s16_2ch_c(short *dst, short *src[2], unsigned int len);
FDSP_EXPORT void conv_s16_to_s16p_2ch_c(short *dst[2], short *src, unsigned int len);

FDSP_EXPORT void conv_fltp_to_flt_2ch_sse(float *dst, float *src[2], unsigned int len);
FDSP_EXPORT void conv_flt_to_fltp_2ch_sse(float *dst[2], float *src, unsigned int len);
FDSP_EXPORT void conv_s16p_to_s16_2ch_sse(short *dst, short *src[2], unsigned int len);
FDSP_EXPORT void conv_s16_to_s16p_2ch_sse(short *dst[2], short *src, unsigned int len);

FDSP_EXPORT void conv_fltp_to_flt_2ch_neon(float *dst, float *src[2], unsigned int len);
FDSP_EXPORT void conv_flt_to_fltp_2ch_neon(float *dst[2], float *src, unsigned int len);
FDSP_EXPORT void conv_s16p_to_s16_2ch_neon(short *dst, short *src[2], unsigned int len);
FDSP_EXPORT void conv_s16_to_s16p_2ch_neon(short *dst[2], short *src, unsigned int len);

#ifdef __cplusplus
};
#endif

#endif /* LIBFDSP_H */

