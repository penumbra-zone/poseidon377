from sage.rings.polynomial.polynomial_gf2x import GF2X_BuildIrred_list

if len(sys.argv) < 4:
    print "Usage: <script> <field_size> <num_constants> <mode> (<prime_number_hex>)"
    print "mode = 0 for GF(2^n), mode = 1 for GF(p)"
    exit()

FIELD_SIZE = int(sys.argv[1]) # size (in bits) of each constant
NUM_CONSTANTS = int(sys.argv[2]) # total number of constants needed
MODE = int(sys.argv[3]) # 0 .. GF(2^n), 1 .. GF(p)

PRIME_NUMBER = 0
if MODE == 1 and len(sys.argv) != 5:
    print "Please specify a prime number (in hex format)!"
    exit()
elif MODE == 1 and len(sys.argv) == 5:
    PRIME_NUMBER = int(sys.argv[4], 16) # e.g. 0xa7, 0xFFFFFFFFFFFFFEFF, 0xa1a42c3efd6dbfe08daa6041b36322ef

def grain_sr_generator():
    bit_sequence = [1] * 80
    for _ in range(0, 160):
        new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)
        
    while True:
        new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)
        while new_bit == 0:
            new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
            bit_sequence.pop(0)
            bit_sequence.append(new_bit)
            new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
            bit_sequence.pop(0)
            bit_sequence.append(new_bit)
        new_bit = bit_sequence[62] ^^ bit_sequence[51] ^^ bit_sequence[38] ^^ bit_sequence[23] ^^ bit_sequence[13] ^^ bit_sequence[0]
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)
        yield new_bit
grain_gen = grain_sr_generator()
        
def grain_random_bits(num_bits):
    random_bits = [grain_gen.next() for i in range(0, num_bits)]
    random_int = int("".join(str(i) for i in random_bits), 2)
    return random_int

def generate_constants(field_size, num_constants, mode):
    round_constants = []
    if mode == 0:
        for i in range(0, num_constants):
            random_int = grain_random_bits(field_size)
            round_constants.append(random_int)
        print "Round constants for GF(2^n):"
    elif mode == 1:
        for i in range(0, num_constants):
            random_int = grain_random_bits(field_size)
            while random_int >= PRIME_NUMBER:
                # print "[Info] Round constant is not in prime field! Taking next one."
                random_int = grain_random_bits(field_size)
            round_constants.append(random_int)
        print "Round constants for GF(p):"
    hex_length = int(ceil(float(field_size) / 4)) + 2 # +2 for "0x"
    print ["{0:#0{1}x}".format(entry, hex_length) for entry in round_constants]


generate_constants(FIELD_SIZE, NUM_CONSTANTS, MODE)