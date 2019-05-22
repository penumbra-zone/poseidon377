# For GF(2^n)
from sage.rings.polynomial.polynomial_gf2x import GF2X_BuildIrred_list

if len(sys.argv) < 3:
    print "Usage: <script> <N> <t>"
    exit()

N = int(sys.argv[1])
t = int(sys.argv[2])
n = int(N / t)
irred = GF(2)['x'](GF2X_BuildIrred_list(n))
F.<x> = GF(2**n, name='x', modulus = irred)

def print_matrix_format(M_int, n, t):
    print "n:", n
    print "t:", t
    print "N:", (n * t)
    hex_length = int(ceil(float(n) / 4)) + 2 # +2 for "0x"
    print "Irreducible polynomial:", irred
    print "MDS matrix (rows):"
    for i in range(0, t):
        print ["{0:#0{1}x}".format(entry, hex_length) for entry in M_int[i]]

def matrix_entries_to_int(M, t):
    M_int = []
    for i in range(0, t):
        M_int.append([])
        for j in range(0, t):
            M_int[i].append(M[i, j].integer_representation())
    return M_int

def create_mds(n, t, start):
    M = matrix(F, t, t)
    xs = []
    ys = []
    
    for i in range(0, t):
        xs.append(F.fetch_int(start + i))
        ys.append(F.fetch_int(start + t + i))
    
    for i in range(0, t):
        for j in range(0, t):
            entry = (xs[i] + ys[j])^(-1)
            M[i, j] = entry
    return M
    
mds_matrix = create_mds(n, t, 0)
mds_matrix_int = matrix_entries_to_int(mds_matrix, t)
print_matrix_format(mds_matrix_int, n, t)