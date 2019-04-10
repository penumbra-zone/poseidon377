
#include <utils/F_2_63_ops.h>
#include <ios>
#include <iostream>

int testReduction() {
    __m128i unreduced = _mm_set_epi64x(0x5d961c716063c2aULL, 0x1e6765bc2ed6aab8ULL);
    uint64_t reduced = hashMIMC::F_2_63::reduce(unreduced);
    return reduced != 0x2b0212e5ac22244ULL;
}

int testMulReduction() {
    uint64_t a =  0x4008010020020082ULL;
    uint64_t b = 0x6400018400064004;
    uint64_t c = hashMIMC::F_2_63::mulReduce(a,b);
    return  c != 0x1b62cd1a00e64aaeULL;
}

int testSquareReduction() {
    uint64_t a =  0x4008010020020082ULL;
    uint64_t b = hashMIMC::F_2_63::mulReduce(a,a);
    return  b != 0x6400018400064004ULL;
}

int main() {
    int ret = 0;
    ret |= testReduction();
    ret |= testSquareReduction();
    ret |= testMulReduction();
    return ret;

}