#include <cmath>
#include <vector>
#include <complex>
#include <immintrin.h>
#include <iostream>
#include <numeric>

class FFT_AVX2
{
public:
    /**
     * @brief Cooley-Tukey FFT Algorithm Implementation
     */
    static void fft(std::complex<float>* data, int n)
    {
        if (n <= 1)
            return;

        // Bit-reversal permutation
        bitReverse(data, n);

        // FFT computation
        for (int i = 1; i <= std::log2(n); ++i)
        {
            int m = 1 << i;
            int m2 = m >> 1;
            
            // Twiddle factor base for this stage
            constexpr double DEGREES = -2.0 * M_PI;
            std::complex<float> wm(std::cos(DEGREES / m), std::sin(DEGREES / m));

            // Process all blocks of size m
            for (int j = 0; j < n; j += m)
            {
                std::complex<float> w(1.0, 0.0);

                // Process butterflies within this block
                for (int k = 0; k < m2;)
                {
                    // The parallel operation sums at the same time 4 complex numbers,
                    // thus we are analyzing the code 4 elements at a time
                    if (k + 4 <= m2)
                    {
                        butterflyParallel(data, j, k, m2, w, wm);
                        k += 4;
                    }
                    else
                    {
                        butterfly(data[j + k], data[j + k + m2], w);
                        w *= wm;
                        k += 1;
                    }
                }
            }
        }
    }

    /**
     * @brief Bit-reversal permutation
     */
    static void bitReverse(std::complex<float>* data, int n)
    {
        int j = 0;
        for (int i = 0; i < n - 1; ++i)
        {
            if (i < j)
                std::swap(data[i], data[j]);

            int k = n >> 1;
            while (k <= j)
            {
                j -= k;
                k >>= 1;
            }
            j += k;
        }
    }

    /**
     * @brief Scalar butterfly operation
     */
    static void butterfly(std::complex<float>& a, std::complex<float>& b, const std::complex<float>& w)
    {
        std::complex<float> t = w * b;
        b = a - t;
        a = a + t;
    }

    /**
     * @brief AVX2-optimized butterfly parallel operation for 4 complex numbers
     * This is basically a parallel implementation of the complex multiplication (a+bi)(c+di) = (ac-bd) + (ad+bc)i
     * We compute 4 elements at a time because each complex is composed by 2 floats made of 32 bits each.
     * This should take a theoretica speedup of 4x
     */
    static void butterflyParallel(std::complex<float>* data, const int j, const int k, const int m2, std::complex<float>& w, const std::complex<float>& wm)
    {
        // Loads 4 complex numbers (8 floats) in a 256 bit register -> 8 floats * 32 bits
        // Returning a 256 vector containing all the 8 floats in a consecutive manner, line (r0, i0) etc.
        // where that couple is a complex number
        // it is unaligned, meaning that floats can be not consecutive
        __m256 upper = _mm256_loadu_ps(reinterpret_cast<float*>(&data[j + k]));
        __m256 lower = _mm256_loadu_ps(reinterpret_cast<float*>(&data[j + k + m2]));
        
        // Compute twiddle factors

        // We align the struct to 32 bits in order to be sure that the opration _mm256_load_ps works correctly,
        // because it requires that the elements are 32-bit aligned.
        // This way we ensurre that twiddles start at an index that ends with 00
        alignas(32) float twiddles[8];
        std::complex<float> wCurr = w;
        for (int i = 0; i < 4; ++i)
        {
            twiddles[i * 2] = wCurr.real();
            twiddles[i * 2 + 1] = wCurr.imag();
            wCurr *= wm;
        }
        w = wCurr;
        
        __m256 twiddle = _mm256_load_ps(twiddles);
        
        // Duplicate real and imaginary parts
        // High are ODD indexes, Lox are EVEN indexes
        // Given the complex-numbers multiplication formula
        // (a+bi)(c+di) = (ac-bd) + (ad+bc)i, we duplicate real an imaignary parts
        // in order to ease the various multiplications
        // Suppose we start with lower = [r0, i0, r1, i1, r2, i2, r3, i3], we obtain:
        __m256 lowerReal = _mm256_moveldup_ps(lower);  // [r0, r0, r1, r1, r2, r2, r3, r3]
        __m256 lowerImag = _mm256_movehdup_ps(lower);  // [i0, i0, i1, i1, i2, i2, i3, i3]
        
        // Elements-wise operation for the 8 real floats
        __m256 acbd = _mm256_mul_ps(lowerReal, twiddle);
        
        // Rearranged bits with each 128-bit lane, based on 0xB1 permutation controller, where
        // B1 stands for 10 11 00 01, representing the final position of the permutation.
        // Supposed we have as a starting position 11 10 01 00, we note that the permutation is
        // structured to swap adjacent elements in a couple. Supposed we have a 128-bits line
        // of type ((A B) (C D)), we have ((B A)(D C))
        __m256 twiddle_swap = _mm256_permute_ps(twiddle, 0xB1);
        
        // Elements-wise operation for the 8 imaginary floats
        __m256 adbc = _mm256_mul_ps(lowerImag, twiddle_swap);
        
        // For each element-wise couple, executes addition or subtraction
        __m256 t = _mm256_addsub_ps(acbd, adbc);
        
        // Element-wise addtions and subtractions
        __m256 resultUpper = _mm256_add_ps(upper, t);
        __m256 resultLower = _mm256_sub_ps(upper, t);
        
        // Stores back
        _mm256_storeu_ps(reinterpret_cast<float*>(&data[j + k]), resultUpper);
        _mm256_storeu_ps(reinterpret_cast<float*>(&data[j + k + m2]), resultLower);
    }

    static void ifft(std::complex<float>* data, int n)
    {
        // Computes the conjugate for the provided data
        // We use this because, since the fft and ifft formula differ only for the 
        // singn of the exponent of the single element, we optimize it by computing the
        // conjugate of the complet number, which passed from a + bi to a - bi.
        // Then, we can use the same formula as fft
        for (int i{0}; i < n; ++i)
            data[i] = std::conj(data[i]);

        fft(data, n);

        // Computes the inverse conjugate and scales data
        for (int i{0}; i < n; ++i)
            data[i] = std::conj(data[i] / static_cast<float>(n));
    }

    static void suppressFeedback(std::complex<float>* data, int n, float suppressionFactor = 0.3f)
    {
        // Calculate average magnitude across all bins, in order to get all the bins that excess it
        // Using reduce to be potentially parallel
        const float totalMagnitude = std::reduce(data, data + n, 0.f,
            [](float sum, const std::complex<float>& val) { return sum + std::abs(val); });
        
            const float avgMagnitude = totalMagnitude / n;

        const float threshold = avgMagnitude * 3.0f;  // Suppress frequencies > 3x average magnitude
        for (int i = 0; i < n; ++i)
        {
            // If the current bin is over the limit, applies suppression
            if (std::abs(data[i]) > threshold)
                data[i] *= suppressionFactor;
        }
    }

};
