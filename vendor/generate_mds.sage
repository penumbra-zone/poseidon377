# This is a notably edited copy of STARKWare's code for generating MDS.
# New Fields, optimized MDS matrix generation for low state size,
# Poseidon parameterization, and generating Cpp/Rust code has been added here.
# Copying their license here
# Copyright 2019 StarkWare Industries Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License").
# You may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# https://www.starkware.co/open-source-license/
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions
# and limitations under the License.
# language enum
cpp = 0
rust = 1
# Prime fields.
F61 = GF(2 ** 61 + 20 * 2 ** 32 + 1)
F81 = GF(2 ** 81 + 80 * 2 ** 64 + 1)
F91 = GF(2 ** 91 + 5 * 2 ** 64 + 1)
F125 = GF(2 ** 125 + 266 * 2 ** 64 + 1)
F161 = GF(2 ** 161 + 23 * 2 ** 128 + 1)
F253 = GF(2 ** 253 + 2 ** 199 + 1)
AltBn = GF(
    21888242871839275222246405745257275088548364400416034343698204186575808495617
)
bw6 = GF(
    258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
)
# Binary fields.
X = GF(2)["X"].gen()
Bin63 = GF(2 ** 63, name="a", modulus=X ** 63 + X + 1)
Bin81 = GF(2 ** 81, name="a", modulus=X ** 81 + X ** 4 + 1)
Bin91 = GF(2 ** 91, name="a", modulus=X ** 91 + X ** 8 + X ** 5 + X + 1)
Bin127 = GF(2 ** 127, name="a", modulus=X ** 127 + X + 1)
Bin161 = GF(2 ** 161, name="a", modulus=X ** 161 + X ** 18 + 1)
Bin255 = GF(2 ** 255, name="a", modulus=X ** 255 + X ** 5 + X ** 3 + X ** 2 + 1)


def sponge(permutation_func, inputs, params):
    """
    Applies the sponge construction to permutation_func.
    inputs should be a vector of field elements whose size is divisible by
    params.r.
    permutation_func should be a function which gets (state, params) where state
    is a vector of params.m field elements, and returns a vector of params.m
    field elements.
    """
    assert parent(inputs) == VectorSpace(
        params.field, len(inputs)
    ), "inputs must be a vector of field elements. Found: %r" % parent(inputs)
    assert (
        len(inputs) % params.r == 0
    ), "Number of field elements must be divisible by %s. Found: %s" % (
        params.r,
        len(inputs),
    )
    state = vector([params.field(0)] * params.m)
    for i in range(0, len(inputs), params.r):
        state[: params.r] += inputs[i : i + params.r]
        state = permutation_func(state, params)
    # We do not support more than r output elements, since this requires
    # additional invocations of permutation_func.
    assert params.output_size <= params.r
    return state[: params.output_size]


def generate_round_constant(fn_name, field, idx):
    """
    Returns a field element based on the result of sha256.
    The input to sha256 is the concatenation of the name of the hash function
    and an index.
    For example, the first element for MiMC will be computed using the value
    of sha256('MiMC0').
    """
    from hashlib import sha256

    val = int(sha256(("%s%d" % (fn_name, idx)).encode("utf8")).hexdigest(), 16)
    if field.is_prime_field():
        return field(val)
    else:
        return int2field(field, val % field.order())


def int2field(field, val):
    """
    Converts val to an element of a binary field according to the binary
    representation of val.
    For example, 11=0b1011 is converted to 1*a^3 + 0*a^2 + 1*a + 1.
    """
    assert field.characteristic() == 2
    assert (
        0 <= val < field.order()
    ), "Value %d out of range. Expected 0 <= val < %d." % (val, field.order())
    res = field(map(int, bin(val)[2:][::-1]))
    assert res.integer_representation() == val
    return res


def binary_vector(field, values):
    """
    Converts a list of integers to field elements using int2field.
    """
    return vector(field, [int2field(field, val) for val in values])


def binary_matrix(field, values):
    """
    Converts a list of lists of integers to field elements using int2field.
    """
    return matrix(field, [[int2field(field, val) for val in row] for row in values])


def generate_mds_matrix(name, field, m, optimize_mds=True):
    """
    Generates an MDS matrix of size m x m over the given field, with no
    eigenvalues in the field.
    Given two disjoint sets of size m: {x_1, ..., x_m}, {y_1, ..., y_m} we set
    A_{ij} = 1 / (x_i - y_j).
    """
    for attempt in range(100):
        x_values = [
            generate_round_constant(name + "x", field, attempt * m + i)
            for i in range(m)
        ]
        y_values = [
            generate_round_constant(name + "y", field, attempt * m + i)
            for i in range(m)
        ]
        # Make sure the values are distinct.
        assert (
            len(set(x_values + y_values)) == 2 * m
        ), "The values of x_values and y_values are not distinct"
        mds_proto = [
            [1 / (x_values[i] - y_values[j]) for j in range(m)] for i in range(m)
        ]
        if optimize_mds:
            # massive reduction in constraint complexity, and computation time
            # These are near-MDS matrices
            if m == 3:
                print('we are using m=3 optimization')
                # [[1, 0, 1],
                #  [1, 1, 0],
                #  [0, 1, 1]]
                mds_proto = [
                    [field(1), field(0), field(1)],
                    [field(1), field(1), field(0)],
                    [field(0), field(1), field(1)],
                ]
            elif m == 4:
                mds_proto = [
                    [field(1) if x != y else field(0) for x in range(m)]
                    for y in range(m)
                ]
            else:
                raise RuntimeError
            mds = matrix(mds_proto)
            return mds
        mds = matrix(mds_proto)
        assert mds.determinant() != 0
        if not optimize_mds:
            # Sanity check: check the determinant of the matrix.
            x_prod = product(
                [x_values[i] - x_values[j] for i in range(m) for j in range(i)]
            )
            y_prod = product(
                [y_values[i] - y_values[j] for i in range(m) for j in range(i)]
            )
            xy_prod = product(
                [x_values[i] - y_values[j] for i in range(m) for j in range(m)]
            )
            expected_det = (1 if m % 4 < 2 else -1) * x_prod * y_prod / xy_prod
            det = mds.determinant()
            assert det != 0
            assert det == expected_det, "Expected determinant %s. Found %s" % (
                expected_det,
                det,
            )
        if len(mds.characteristic_polynomial().roots()) == 0:
            # There are no eigenvalues in the field.
            return mds
        #print(mds.characteristic_polynomial().roots())
    raise Exception("No good MDS found")


def generate_mds_cpp(mds):
    s = ["std::vector<std::vector<FieldT>> mds_matrix;"]
    # bigint<FieldT::num_limbs>("1234")
    for r in mds:
        line = "mds_matrix.push_back(std::vector<FieldT>({"
        for c in r:
            line += 'FieldT(bigint<FieldT::num_limbs>("' + str(c) + '")),'
        line = line[:-1]
        line += "}));"
        s += [line]
    print("\n".join(s))


def generate_mds_rs(mds):
    # let mds = vec![vec![F::one(),  F::zero(), F::one()],
    #                vec![F::one(),  F::one(),  F::zero()],
    #                vec![F::zero(), F::one(),  F::one()]];
    # but F::from_str(\"" + str(c) + "\").map_err(|_| ()).unwrap()
    s = ["let mds = vec!["]
    for r in mds:
        line = "vec!["
        for c in r:
            line += 'F::from_str("' + str(c) + '").map_err(|_| ()).unwrap(),'
        line = line[:-1]
        line += "],"
        s += [line]
    s[-1] += "];"
    print("\n".join(s))


def generate_ark(
    hash_name, field, state_size, num_rounds, optimize_ark=False, R_p=0, mds=None
):
    ark = [
        vector(
            generate_round_constant("Hades", field, state_size * i + j)
            for j in range(state_size)
        )
        for i in range(num_rounds)
    ]
    s = ["std::vector<std::vector<FieldT>> ark_matrix;"]
    # bigint<FieldT::num_limbs>("1234")
    for r in ark:
        line = "ark_matrix.push_back(std::vector<FieldT>({"
        for c in r:
            line += 'FieldT(bigint<FieldT::num_limbs>("' + str(c) + '")),'
        line = line[:-1]
        line += "}));"
        s += [line]
    print("\n".join(s))
    if optimize_ark:
        # Remove Rf
        ark = ark[4:-4]
        # ark = ark[::-1]
        mds_inv = mds.inverse()
        mds_pow = mds_inv
        linear_constants = ark[0]
        linear_constants[0] = field(0)
        non_linear_constants = [ark[0][0]]
        for i in range(1, R_p):
            cur = mds_pow * ark[i]
            non_linear_constants = non_linear_constants + [cur[0]]
            cur[0] = field(0)
            linear_constants += cur
            mds_pow *= mds_inv
        s = "std::vector<FieldT> rp_linear_ark_constants = std::vector<FieldT>({"
        for c in linear_constants:
            s += 'FieldT(bigint<FieldT::num_limbs>("' + str(c) + '")),'
        s += "});"
        print(s)
        s = "std::vector<FieldT> rp_non-linear_ark_constants = std::vector<FieldT>({"


def generate_ark_rs(
    hash_name, field, alpha, state_size, num_rounds, optimize_ark=False, R_p=0, mds=None
):
    ark = [
        vector(
            generate_round_constant("Hades", field, state_size * i + j)
            for j in range(state_size)
        )
        for i in range(num_rounds)
    ]
    s = ["let ark = vec!["]
    # bigint<FieldT::num_limbs>("1234")
    for r in ark:
        line = "vec!["
        for c in r:
            line += 'F::from_str("' + str(c) + '").map_err(|_| ()).unwrap(),'
        line = line[:-1]
        line += "],"
        s += [line]
    print("\n".join(s) + "];")
    # if optimize_ark:
    #     # Remove Rf
    #     ark = ark[4:-4]
    #     # ark = ark[::-1]
    #     mds_inv = mds.inverse()
    #     mds_pow = mds_inv
    #     linear_constants = ark[0]
    #     linear_constants[0] = field(0)
    #     non_linear_constants = [ark[0][0]]
    #     for i in range(1, R_p):
    #         cur = mds_pow * ark[i]
    #         non_linear_constants = non_linear_constants + [cur[0]]
    #         cur[0] = field(0)
    #         linear_constants += cur
    #         mds_pow *= mds_inv
    #     s = "std::vector<FieldT> rp_linear_ark_constants = std::vector<FieldT>({"
    #     for c in linear_constants:
    #         s += "FieldT(bigint<FieldT::num_limbs>(\"" + str(c) + "\")),"
    #     s += "});"
    #     print(s)
    #     s = "std::vector<FieldT> rp_non-linear_ark_constants = std::vector<FieldT>({"


def generate_mds_code(hash_name, field, state_size, optimize_mds, lang=cpp):
    mds = generate_mds_matrix(hash_name + "MDS", field, state_size, optimize_mds)
    if lang == cpp:
        generate_mds_cpp(mds)
    else:
        generate_mds_rs(mds)
    return mds


def generate_poseidon_param_code(
    hash_name,
    field,
    alpha,
    state_size,
    num_rounds,
    optimize_mds=True,
    optimize_ark=False,
):
    # Assume that the capacity is 1
    rate = state_size - 1
    R_f = 8
    R_p = num_rounds - R_f

    opening = f"""
/// Parameters for the rate-{rate} instance of Poseidon.
///
/// Note: `F` must be the BLS12-377 scalar field.
pub fn rate_{rate}<F: PrimeField>() -> PoseidonParameters<F> {{
"""

    closing = f"""
    PoseidonParameters {{
        full_rounds: {num_rounds - R_p},
        partial_rounds: {R_p},
        alpha: {alpha},
        ark,
        mds,
        rate: {rate},
        capacity: 1,
    }}
}}
"""
   
    print(opening)

    mds = generate_mds_code(hash_name, field, state_size, optimize_mds, lang=rust)
    generate_ark_rs(
        hash_name, field, alpha, state_size, num_rounds, optimize_ark, R_p, mds
    )
    print(closing)


# This just calculates the minimum number rounds per the dominating constraint
# We do a full calculation against all relevant equations in the CPP logic, this is just for sanity checks.
# The dominating constraint is the defense against interpolation attack.
# TODO: Include analysis not in the paper of capacity > 1 improving num rounds
def calculate_num_poseidon_rounds(field, sec, alpha, num_capacity_elems, state_size, p):
    assert num_capacity_elems * len(bin(field.order())[2:]) > sec

    # From Eqn 3 section 5.5.2 on p.10, where M is the security level
    min_num_rounds = ceil(log(2, alpha) * min(sec, log(p, state_size))) + ceil(log(state_size, 2))

    full_rounds = 6
    partial_rounds = min_num_rounds - full_rounds

    # Section 5.4 of the paper mentions that in addition to the 7.5% security margin,
    # two additional rounds with full S-box layers are also applied.
    partial_rounds = int(ceil(partial_rounds * 1.075))
    full_rounds += 2

    return partial_rounds + full_rounds


poseidon_hash_name = "Hades"
# Penumbra
alpha=17
capacity_size = 1
sec = 128
p = 8444461749428370424248824938781546531375899335154063827935233455917409239041
bls377 = GF(p)
# You can set optimize_mds to True if you want to take a heuristic on using near-MDS matrices
# The paper authors were of the opinion that this worked (with an update to the differential analysis given)
# which will be satisfied for any large field

print("""//! Parameters for various Poseidon instances over the BLS12-377 scalar field.
//!
//! All parameter instances have capacity 1, targeting the 128-bit security
//! level.

// Generated with `generate_mds.sage`. Do not edit manually.
// Regenerate with: `sage vendor/generate_mds.sage > src/params.rs`

use ark_ff::PrimeField;
use ark_sponge::poseidon::PoseidonParameters;

""")

arity = 1
state_size = arity + capacity_size
num_rounds = calculate_num_poseidon_rounds(bls377, sec, alpha, capacity_size, state_size, p)
generate_poseidon_param_code(
    poseidon_hash_name, bls377, alpha, state_size, num_rounds, optimize_mds=False,
)

arity = 2
state_size = arity + capacity_size
num_rounds = calculate_num_poseidon_rounds(bls377, sec, alpha, capacity_size, state_size, p)
generate_poseidon_param_code(
    poseidon_hash_name, bls377, alpha, state_size, num_rounds, optimize_mds=False
)

arity = 4
state_size = arity + capacity_size
num_rounds = calculate_num_poseidon_rounds(bls377, sec, alpha, capacity_size, state_size, p)
generate_poseidon_param_code(
    poseidon_hash_name, bls377, alpha, state_size, num_rounds, optimize_mds=False,
)

arity = 5
state_size = arity + capacity_size
num_rounds = calculate_num_poseidon_rounds(bls377, sec, alpha, capacity_size, state_size, p)
generate_poseidon_param_code(
    poseidon_hash_name, bls377, alpha, state_size, num_rounds, optimize_mds=False,
)