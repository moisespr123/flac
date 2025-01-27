#include "private/cpu.h"

#ifndef FLAC__INTEGER_ONLY_LIBRARY
#ifndef FLAC__NO_ASM
#if defined FLAC__CPU_AARCH64 && FLAC__HAS_A64NEONINTRIN
#include "private/lpc.h"
#include "FLAC/assert.h"
#include "FLAC/format.h"
#include <arm_neon.h>

void FLAC__lpc_compute_autocorrelation_intrin_neon_lag_8(const FLAC__real data[], uint32_t data_len, uint32_t lag, double autoc[])
{
    // This function calculates autocorrelation with NEON
    // vector functions up to a lag of 8 (or max LPC order of 7)
    int i;
    float64x2_t sum0 = vdupq_n_f64(0.0f);
    float64x2_t sum1 = vdupq_n_f64(0.0f);
    float64x2_t sum2 = vdupq_n_f64(0.0f);
    float64x2_t sum3 = vdupq_n_f64(0.0f);
	float64x2_t d0 = vdupq_n_f64(0.0f);
	float64x2_t d1 = vdupq_n_f64(0.0f);
	float64x2_t d2 = vdupq_n_f64(0.0f);
	float64x2_t d3 = vdupq_n_f64(0.0f);
	float64x2_t d;

    (void)lag;
    FLAC__ASSERT(lag <= 8);

	// Loop backwards through samples from data_len to 0
	for (i = data_len - 1; i >= 0; i--)
	{
		d = vdupq_n_f64(data[i]); // Create vector with 2 entries data[i]

		// The next lines of code right-shift the elements through the vectors d0..d3.
		// The last line adds the newly loaded element to d0. This works like a stack, where
		// data[i] is pushed onto the stack every time and the last element falls off
		d3 = vextq_f64(d2,d3,1);
		d2 = vextq_f64(d1,d2,1);
		d1 = vextq_f64(d0,d1,1);
		d0 = vextq_f64(d,d0,1);

		// Fused multiply-add sum += d * d0..d3
		sum0 = vfmaq_f64(sum0, d, d0);
		sum1 = vfmaq_f64(sum1, d, d1);
		sum2 = vfmaq_f64(sum2, d, d2);
		sum3 = vfmaq_f64(sum3, d, d3);
	}

    // Store sum0..sum4 in autoc[0..8]
    vst1q_f64(autoc, sum0);
    vst1q_f64(autoc + 2, sum1);
    vst1q_f64(autoc + 4, sum2);
    vst1q_f64(autoc + 6, sum3);
}

void FLAC__lpc_compute_autocorrelation_intrin_neon_lag_10(const FLAC__real data[], uint32_t data_len, uint32_t lag, double autoc[])
{
    // This function calculates autocorrelation with NEON
    // vector functions up to a lag of 10 (or max LPC order of 9)
    int i;
    float64x2_t sum0 = vdupq_n_f64(0.0f);
    float64x2_t sum1 = vdupq_n_f64(0.0f);
    float64x2_t sum2 = vdupq_n_f64(0.0f);
    float64x2_t sum3 = vdupq_n_f64(0.0f);
    float64x2_t sum4 = vdupq_n_f64(0.0f);
	float64x2_t d0 = vdupq_n_f64(0.0f);
	float64x2_t d1 = vdupq_n_f64(0.0f);
	float64x2_t d2 = vdupq_n_f64(0.0f);
	float64x2_t d3 = vdupq_n_f64(0.0f);
	float64x2_t d4 = vdupq_n_f64(0.0f);
	float64x2_t d;

    (void)lag;
    FLAC__ASSERT(lag <= 10);

	// Loop backwards through samples from data_len to 0
	for (i = data_len - 1; i >= 0; i--)
	{
		d = vdupq_n_f64(data[i]); // Create vector with 2 entries data[i]

		// The next lines of code right-shift the elements through the vectors d0..d4.
		// The last line adds the newly loaded element to d0. This works like a stack, where
		// data[i] is pushed onto the stack every time and the last element falls off
		d4 = vextq_f64(d3,d4,1);
		d3 = vextq_f64(d2,d3,1);
		d2 = vextq_f64(d1,d2,1);
		d1 = vextq_f64(d0,d1,1);
		d0 = vextq_f64(d,d0,1);

		// Fused multiply-add sum += d * d0..d4
		sum0 = vfmaq_f64(sum0, d, d0);
		sum1 = vfmaq_f64(sum1, d, d1);
		sum2 = vfmaq_f64(sum2, d, d2);
		sum3 = vfmaq_f64(sum3, d, d3);
		sum4 = vfmaq_f64(sum4, d, d4);
	}

    // Store sum0..sum4 in autoc[0..10]
    vst1q_f64(autoc, sum0);
    vst1q_f64(autoc + 2, sum1);
    vst1q_f64(autoc + 4, sum2);
    vst1q_f64(autoc + 6, sum3);
    vst1q_f64(autoc + 8, sum4);
}

void FLAC__lpc_compute_autocorrelation_intrin_neon_lag_14(const FLAC__real data[], uint32_t data_len, uint32_t lag, double autoc[])
{
    // This function calculates autocorrelation with NEON
    // vector functions up to a lag of 14 (or max LPC order of 13)
    int i;
    float64x2_t sum0 = vdupq_n_f64(0.0f);
    float64x2_t sum1 = vdupq_n_f64(0.0f);
    float64x2_t sum2 = vdupq_n_f64(0.0f);
    float64x2_t sum3 = vdupq_n_f64(0.0f);
    float64x2_t sum4 = vdupq_n_f64(0.0f);
    float64x2_t sum5 = vdupq_n_f64(0.0f);
    float64x2_t sum6 = vdupq_n_f64(0.0f);
	float64x2_t d0 = vdupq_n_f64(0.0f);
	float64x2_t d1 = vdupq_n_f64(0.0f);
	float64x2_t d2 = vdupq_n_f64(0.0f);
	float64x2_t d3 = vdupq_n_f64(0.0f);
	float64x2_t d4 = vdupq_n_f64(0.0f);
	float64x2_t d5 = vdupq_n_f64(0.0f);
	float64x2_t d6 = vdupq_n_f64(0.0f);
	float64x2_t d;

    (void)lag;
    FLAC__ASSERT(lag <= 14);

	// Loop backwards through samples from data_len to 0
	for (i = data_len - 1; i >= 0; i--)
	{
		d = vdupq_n_f64(data[i]); // Create vector with 2 entries data[i]

		// The next 6 lines of code right-shift the elements through the 7 vectors d0..d6.
		// The 7th line adds the newly loaded element to d0. This works like a stack, where
		// data[i] is pushed onto the stack every time and the 9th element falls off
		d6 = vextq_f64(d5,d6,1);
		d5 = vextq_f64(d4,d5,1);
		d4 = vextq_f64(d3,d4,1);
		d3 = vextq_f64(d2,d3,1);
		d2 = vextq_f64(d1,d2,1);
		d1 = vextq_f64(d0,d1,1);
		d0 = vextq_f64(d,d0,1);

		// Fused multiply-add sum += d * d0..d6
		sum0 = vfmaq_f64(sum0, d, d0);
		sum1 = vfmaq_f64(sum1, d, d1);
		sum2 = vfmaq_f64(sum2, d, d2);
		sum3 = vfmaq_f64(sum3, d, d3);
		sum4 = vfmaq_f64(sum4, d, d4);
		sum5 = vfmaq_f64(sum5, d, d5);
		sum6 = vfmaq_f64(sum6, d, d6);
	}

    // Store sum0..sum6 in autoc[0..14]
    vst1q_f64(autoc, sum0);
    vst1q_f64(autoc + 2, sum1);
    vst1q_f64(autoc + 4, sum2);
    vst1q_f64(autoc + 6, sum3);
    vst1q_f64(autoc + 8, sum4);
    vst1q_f64(autoc + 10, sum5);
    vst1q_f64(autoc + 12, sum6);
}
#endif /* FLAC__CPU_AARCH64 && FLAC__HAS_A64NEONINTRIN */
#endif /* FLAC__NO_ASM */
#endif /* FLAC__INTEGER_ONLY_LIBRARY */
