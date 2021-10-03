from math import *
import sys
import Crypto.Util.number

def sat_inequiv_alpha(p, t, R_F, R_P, alpha, M):
    n = ceil(log(p, 2))
    N = int(n * t)
    if alpha > 0:
        R_F_1 = 6 if M <= ((floor(log(p, 2) - ((alpha-1)/2.0))) * (t + 1)) else 10 # Statistical
        R_F_2 = 1 + ceil(log(2, alpha) * min(M, n)) + ceil(log(t, alpha)) - R_P # Interpolation
        #R_F_3 = ceil(min(n, M) / float(3*log(alpha, 2))) - R_P # Groebner 1
        #R_F_3 = ((log(2, alpha) / float(2)) * min(n, M)) - R_P # Groebner 1
        R_F_3 = 1 + (log(2, alpha) * min(M/float(3), log(p, 2)/float(2))) - R_P # Groebner 1
        R_F_4 = t - 1 + min((log(2, alpha) * M) / float(t+1), ((log(2, alpha)*log(p, 2)) / float(2))) - R_P # Groebner 2
        #R_F_5 = ((1.0/(2*log((alpha**alpha)/float((alpha-1)**(alpha-1)), 2))) * min(n, M) + t - 2 - R_P) / float(t - 1) # Groebner 3
        R_F_max = max(ceil(R_F_1), ceil(R_F_2), ceil(R_F_3), ceil(R_F_4))
        return (R_F >= R_F_max)
    elif alpha == (-1):
        R_F_1 = 6 if M <= ((floor(log(p, 2) - 2)) * (t + 1)) else 10 # Statistical
        R_P_1 = 1 + ceil(0.5 * min(M, n)) + ceil(log(t, 2)) - floor(R_F * log(t, 2)) # Interpolation
        R_P_2 = 1 + ceil(0.5 * min(M, n)) + ceil(log(t, 2)) - floor(R_F * log(t, 2))
        R_P_3 = t - 1 + ceil(log(t, 2)) + min(ceil(M / float(t+1)), ceil(0.5*log(p, 2))) - floor(R_F * log(t, 2)) # Groebner 2
        R_F_max = ceil(R_F_1)
        R_P_max = max(ceil(R_P_1), ceil(R_P_2), ceil(R_P_3))
        return (R_F >= R_F_max and R_P >= R_P_max)
    else:
        print("Invalid value for alpha!")
        exit(1)

# Hypothesis: Multiplying the number of Sboxes by 5 (depth of addition chain) should yield number of constraints.
def get_sbox_cost(R_F, R_P, N, t):
    return int(t * R_F + R_P)

def get_size_cost(R_F, R_P, N, t):
    n = ceil(float(N) / t)
    return int((N * R_F) + (n * R_P))

def get_depth_cost(R_F, R_P, N, t):
    return int(R_F + R_P)

def find_FD_round_numbers(p, t, alpha, M, cost_function, security_margin):
    n = ceil(log(p, 2))
    N = int(n * t)

    sat_inequiv = sat_inequiv_alpha
    
    R_P = 0
    R_F = 0
    min_cost = float("inf")
    max_cost_rf = 0
    # Brute-force approach
    for R_P_t in range(1, 500):
        for R_F_t in range(4, 100):
            if R_F_t % 2 == 0:
                if (sat_inequiv(p, t, R_F_t, R_P_t, alpha, M) == True):
                    # From paper in section 5.4: arbitrarily add two full rounds and 7.5% more partial rounds.
                    if security_margin == True:
                        R_F_t += 2
                        R_P_t = int(ceil(float(R_P_t) * 1.075))
                    cost = cost_function(R_F_t, R_P_t, N, t)
                    if (cost < min_cost) or ((cost == min_cost) and (R_F_t < max_cost_rf)):
                        R_P = ceil(R_P_t)
                        R_F = ceil(R_F_t)
                        min_cost = cost
                        max_cost_rf = R_F
    return (int(R_F), int(R_P))

def calc_final_numbers_fixed(p, t, alpha, M, security_margin):
    # [Min. S-boxes] Find best possible for t and N
    n = ceil(log(p, 2))
    N = int(n * t)
    cost_function = get_sbox_cost
    ret_list = []
    (R_F, R_P) = find_FD_round_numbers(p, t, alpha, M, cost_function, security_margin)
    min_sbox_cost = cost_function(R_F, R_P, N, t)
    ret_list.append(R_F)
    ret_list.append(R_P)
    ret_list.append(min_sbox_cost)

    # [Min. Size] Find best possible for t and N
    # Minimum number of S-boxes for fixed n results in minimum size also (round numbers are the same)!
    min_size_cost = get_size_cost(R_F, R_P, N, t)
    ret_list.append(min_size_cost)

    return ret_list # [R_F, R_P, min_sbox_cost, min_size_cost]

def print_latex_table_combinations(combinations, alpha, security_margin):
    for comb in combinations:
        N = comb[0]
        t = comb[1]
        M = comb[2]
        n = int(N / t)
        prime = Crypto.Util.number.getPrime(n)
        ret = calc_final_numbers_fixed(prime, t, alpha, M, security_margin)
        field_string = "\mathbb F_{p}"
        sbox_string = "x^{" + str(alpha) + "}"
        print("$" + str(M) + "$ & $" + str(N) + "$ & $" + str(n) + "$ & $" + str(t) + "$ & $" + str(ret[0]) + "$ & $" + str(ret[1]) + "$ & $" + field_string + "$ & $" + str(ret[2]) + "$ & $" + str(ret[3]) + "$ \\\\")

# Single tests
# print calc_final_numbers_fixed(Crypto.Util.number.getPrime(64), 24, 3, 128, True)
# print calc_final_numbers_fixed(Crypto.Util.number.getPrime(253), 6, -1, 128, True)
print(calc_final_numbers_fixed(Crypto.Util.number.getPrime(255), 3, 5, 128, True))
print(calc_final_numbers_fixed(Crypto.Util.number.getPrime(255), 6, 5, 128, True))
print(calc_final_numbers_fixed(Crypto.Util.number.getPrime(254), 3, 5, 128, True))
print(calc_final_numbers_fixed(Crypto.Util.number.getPrime(254), 6, 5, 128, True))
print(calc_final_numbers_fixed(Crypto.Util.number.getPrime(64), 24, 3, 128, True))

# x^5 (254-bit prime number)
#prime = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
x_5_combinations = [
    [1536, 2, 128], [1536, 4, 128], [1536, 6, 128], [1536, 8, 128], [1536, 16, 128],
    [1536, 2, 256], [1536, 4, 256], [1536, 6, 256], [1536, 8, 256], [1536, 16, 256]
]

# With security margin
print("--- Table x^5 WITH security margin ---")
print_latex_table_combinations(x_5_combinations, 5, True)

# Without security margin
print("--- Table x^5 WITHOUT security margin ---")
print_latex_table_combinations(x_5_combinations, 5, False)

x_3_combinations = [
    [1536, 2, 128], [1536, 4, 128], [1536, 6, 128], [1536, 8, 128], [1536, 16, 128],
    [1536, 2, 256], [1536, 4, 256], [1536, 6, 256], [1536, 8, 256], [1536, 16, 256]
]

# With security margin
print("--- Table x^3 WITH security margin ---")
print_latex_table_combinations(x_3_combinations, 3, True)

# Without security margin
print("--- Table x^3 WITHOUT security margin ---")
print_latex_table_combinations(x_3_combinations, 3, False)

x_inv_combinations = [
    [1536, 2, 128], [1536, 4, 128], [1536, 6, 128], [1536, 8, 128], [1536, 16, 128],
    [1536, 2, 256], [1536, 4, 256], [1536, 6, 256], [1536, 8, 256], [1536, 16, 256]
]

# With security margin
print("--- Table x^(-1) WITH security margin ---")
print_latex_table_combinations(x_inv_combinations, -1, True)

# Without security margin
print("--- Table x^(-1) WITHOUT security margin ---")
print_latex_table_combinations(x_inv_combinations, -1, False)

bls377_prime = 8444461749428370424248824938781546531375899335154063827935233455917409239041

M = 128  # Bits of security
security_margin = True

addition_chain_depth = 5
alpha_vals = [11, 17]
t_vals = [2, 3, 5]
for alpha in alpha_vals:
    for t in t_vals:
        print(f'for alpha {alpha} and t {t}')
        R_F, R_P, min_sbox_cost, min_size_cost = calc_final_numbers_fixed(bls377_prime, t, alpha, M, security_margin)
        print(f'results: R_F {R_F}, R_P {R_P}, min_sbox_cost {min_sbox_cost}, min_size_cost {min_size_cost}')
        num_constraints = min_sbox_cost * addition_chain_depth
        print(f'num_constraints for Groth16: {num_constraints}')

# We select alpha=17 as it has less partial rounds than 11. 
