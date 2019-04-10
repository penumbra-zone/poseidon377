#ifndef HASHMIMC_HADESMIMC_H
#define HASHMIMC_HADESMIMC_H


#include <cstddef>
#include <cstdint>

namespace hashMIMC {

class HadesMIMCPermutation {
public:
    virtual ~HadesMIMCPermutation() = default;
    virtual void permute (uint64_t* buffer) = 0;
    virtual size_t get_num_sboxes () = 0;

    virtual void SBOX (uint64_t* state) = 0;
    virtual void SBOX_PART (uint64_t* state) = 0;
    virtual void LIN (uint64_t* state) = 0;
    virtual void ARK (uint64_t* state, size_t round) = 0;
};

//class HadesMIMCPermutation_63_24_10 : public HadesMIMCPermutation{
//public:
//    static constexpr size_t sbox_width = 63;
//    static constexpr size_t num_sboxes = 24;
//    static constexpr size_t state_size = (sbox_width) * num_sboxes;
//    static constexpr size_t state_size_bytes = (state_size + 7) / 8;
//    static constexpr size_t num_rounds = 10;
//
//    virtual size_t get_num_sboxes () { return num_sboxes; }
//    void toState (uint64_t* state, const uint8_t* buffer);
//    void fromState (uint8_t* buffer, const uint64_t* state);
//    void permute (uint64_t* buffer);
//
//    void SBOX (uint64_t* state);
//    void SBOX_PART (uint64_t* state);
//    //void SBOX2 (uint64_t* state);
//    void LIN (uint64_t* state);
//    void ARK (uint64_t* state, size_t round);
//private:
//};

class HadesMIMCPermutation_63_25_12_34 : public HadesMIMCPermutation {
public:
    static constexpr size_t sbox_width = 63;
    static constexpr size_t num_sboxes = 25;
    static constexpr size_t state_size = (sbox_width) * num_sboxes;
    static constexpr size_t num_full_rounds = 12;
    static constexpr size_t num_part_rounds = 34;


    size_t get_num_sboxes () override { return num_sboxes; }
    void toState (uint64_t* state, const uint8_t* buffer);
    void fromState (uint8_t* buffer, const uint64_t* state);
    void permute (uint64_t* buffer) override;

    void SBOX (uint64_t* state) override;
    void SBOX_PART (uint64_t* state) override;
    void SBOX2 (uint64_t* state);
    void LIN (uint64_t* state) override;
    void LIN2 (uint64_t* state);
    void ARK (uint64_t* state, size_t round) override;
private:
};

}

#endif //HASHMIMC_HADESMIMC_H
