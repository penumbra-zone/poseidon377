
#include <HadesMIMC.h>
#include <ios>
#include <iostream>
#include <array>

int test25_12_34_test1() {
    std::cout << "Testing HashMIMC_25_12_34...\n";
    std::cout << "initial state = 25x 0x1 ...\n";

    hashMIMC::HadesMIMCPermutation_63_25_12_34 perm;
    std::array<uint64_t, 25> state;
    for(size_t i = 0; i < 25; i++) {
        state[i] = 0x1;
    }
    perm.permute(state.data());
    for(size_t i = 0; i < 25; i++) {
        std::cout << std::hex << state[i] << ",";
    }
    std::cout << "\n\n";
    return 0;
}

int test25_12_34_test2() {
    std::cout << "Testing HashMIMC_25_12_34...\n";
    std::cout << "initial state = 25x 0x7FFFFFFFFFFFFFFF ...\n";

    hashMIMC::HadesMIMCPermutation_63_25_12_34 perm;
    std::array<uint64_t, 25> state;
    for(size_t i = 0; i < 25; i++) {
        state[i] = 0x7FFFFFFFFFFFFFFFULL;
    }
    perm.permute(state.data());
    for(size_t i = 0; i < 25; i++) {
        std::cout << std::hex << state[i] << ",";
    }
    std::cout << "\n\n";
    return 0;
}

int test25_12_34_test3() {
    std::cout << "Testing HashMIMC_25_12_34...\n";
    std::cout << "initial state = 25x 0x0 ...\n";

    hashMIMC::HadesMIMCPermutation_63_25_12_34 perm;
    std::array<uint64_t, 25> state;
    for(size_t i = 0; i < 25; i++) {
        state[i] = 0x0;
    }
    perm.permute(state.data());
    for(size_t i = 0; i < 25; i++) {
        std::cout << std::hex << state[i] << ",";
    }
    std::cout << "\n\n";
    return 0;
}

int test25_12_34_test4() {
    std::cout << "Testing HashMIMC_25_12_34...\n";
    std::cout << "initial state = 25x i ...\n";

    hashMIMC::HadesMIMCPermutation_63_25_12_34 perm;
    std::array<uint64_t, 25> state;
    for(size_t i = 0; i < 25; i++) {
        state[i] = i;
    }
    perm.permute(state.data());
    for(size_t i = 0; i < 25; i++) {
        std::cout << std::hex << state[i] << ",";
    }
    std::cout << "\n\n";
    return 0;
}

int main() {
    int ret = 0;
    ret |= test25_12_34_test1();
    ret |= test25_12_34_test2();
    ret |= test25_12_34_test3();
    ret |= test25_12_34_test4();
    return ret;

}