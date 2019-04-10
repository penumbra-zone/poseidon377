#ifndef HASHMIMC_F_2_63_OPS_H
#define HASHMIMC_F_2_63_OPS_H

#include <cstdint>
#include <smmintrin.h>
#include <wmmintrin.h>

namespace hashMIMC {
    namespace F_2_63 {

        // X^63 + X + 1
        static const uint64_t IRREDUCIBLE_POLY = 0x8000000000000003ULL;

        // TODO: this is the performance bottleneck, can we do faster?
        inline static uint64_t reduce(__m128i in) {
            uint64_t res = _mm_extract_epi64(in, 0);
            uint64_t b = (res >> 63) & 0x1;
            b = b ^ _mm_extract_epi64(_mm_slli_epi64(in, 1), 1);
            b = b ^ (b >> 62);
            b = b ^ (b << 1);
            return (res ^ b) & 0x7fffffffffffffffULL;
        }

        // variants of reduce that work on the result of pclmul(b, a<<1),
        // so the boundary of the two halves of the result is exactly at 64-bit
        inline static uint64_t reduceOPT1(__m128i in) {
            uint64_t res = _mm_extract_epi64(in, 0);
            __m128i a = _mm_xor_si128(_mm_srli_epi64(in, 62), in);
            a = _mm_xor_si128(_mm_slli_epi64(a, 1), a);
            uint64_t b = _mm_extract_epi64(a, 1);
            return ((res>>1) ^ b) & 0x7fffffffffffffffULL;
        }

        inline static uint64_t reduceOPT2(__m128i in) {
            uint64_t res = _mm_extract_epi64(in, 0);
            uint64_t b = _mm_extract_epi64(in, 1);
            b = b ^ (b >> 62);
            b = b ^ (b << 1);
            return ((res>>1) ^ b) & 0x7fffffffffffffffULL;
        }

        inline static __m128i reduceOPT3(__m128i in) {
            __m128i a = _mm_xor_si128(_mm_srli_epi64(in, 62), in);
            a = _mm_xor_si128(_mm_slli_epi64(a, 1), a);
            //in = _mm_srli_epi64(in, 1);
            a = _mm_srli_si128(a, 8);
            a = _mm_slli_epi64(a, 1);
            return _mm_xor_si128(a,in);
        }

        // TODO is there a more efficient way than above to get result in half of SSE register?
//        inline static __m128i reduceToSSE(__m128i in) {
//        }

        inline static uint64_t mulReduce(uint64_t a, uint64_t b) {
            __m128i tmp = _mm_set_epi64x(b, a<<1);
            __m128i result = _mm_clmulepi64_si128(tmp, tmp, 0x01);
            return reduceOPT1(result);
        }

        inline static __m128i mulReduce2(uint64_t a, uint64_t b) {
            __m128i tmp = _mm_set_epi64x(b, a<<1);
            __m128i result = _mm_clmulepi64_si128(tmp, tmp, 0x01);
            return reduceOPT3(result);
        }

    }

}

#endif //HASHMIMC_F_2_63_OPS_H
