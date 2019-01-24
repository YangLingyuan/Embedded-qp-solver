import sys
import numpy as np
import qpsolvers as qp

def eval_qp(P, q, x):
    return .5 * x.T @ P @ x + q @ x

def read_qp(filename):
    import struct
    import numpy as np
    try:
        fp = open(filename, 'rb')
    except OSError as msg:
        print(msg)
        return None
    bin_data = fp.read()
    fp.close()
    num_bytes = len(bin_data)
    double_size = 8
    fmt = 'd' * (num_bytes // double_size)
    unpacked = struct.unpack(fmt, bin_data)
    n = int(unpacked[0])
    p = np.zeros((n, n))
    q = np.zeros((n, 1))
    for i, j in ((i, j) for i in range(n) for j in range(n)):
        p[i, j] = unpacked[1 + i * n + j]
    for i in range(n):
        q[i] = unpacked[1 + n * n + i]

    return n, p, q

FILENAME = sys.argv[1]
n, p, q = read_qp(FILENAME)
q = q.reshape((n, ))
x = qp.solve_qp(p, q, G=np.zeros((n, n)), h=np.zeros((n,)))
print("reference implementation gave:")
print("f(x_min) = %e" % eval_qp(p, q, x))
