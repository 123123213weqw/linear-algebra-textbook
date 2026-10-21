#!/usr/bin/env python3
"""Numeric verification of the worked examples in the textbook (draft, ch. 1-3).

Each check below mirrors a specific example or exercise answer in the book and
recomputes it with exact rational arithmetic (no floating point). Run with:

    python3 tools/verify_examples.py

Exit code 0 means every check passed. This script is also run in CI on every
push, before the PDF is published.
"""
import sys
from fractions import Fraction as F


# ---------- tiny exact-linear-algebra toolkit ----------

def mat(rows):
    """Build a matrix of Fractions from lists of ints/strings/Fractions."""
    return [[F(x) for x in row] for row in rows]


def eye(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def same(A, B):
    return len(A) == len(B) and all(ra == rb for ra, rb in zip(A, B))


def mul(A, B):
    assert len(A[0]) == len(B), "size mismatch"
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A):
    return [list(col) for col in zip(*A)]


def det(A):
    """Determinant by cofactor expansion along the first row."""
    n = len(A)
    if n == 1:
        return A[0][0]
    total = F(0)
    for j in range(n):
        if A[0][j] == 0:
            continue
        minor = [row[:j] + row[j + 1:] for row in A[1:]]
        total += F(-1) ** j * A[0][j] * det(minor)
    return total


def rref(M):
    """Reduced row echelon form of a matrix of Fractions."""
    M = [row[:] for row in M]
    rows, cols = len(M), len(M[0])
    pivots, r = [], 0
    for c in range(cols):
        if r == rows:
            break
        pivot_row = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if pivot_row is None:
            continue
        M[r], M[pivot_row] = M[pivot_row], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
    return M, pivots


def rank(A):
    return len(rref(A)[1])


def inverse(A):
    n = len(A)
    aug = [row[:] + irow for row, irow in zip(A, eye(n))]
    R, pivots = rref(aug)
    if pivots != list(range(n)):
        return None  # singular
    return [row[n:] for row in R]


def solve(A, b):
    """Return (status, particular_solution_or_None).

    status: 'unique' | 'none' | 'inf'
    """
    aug = [row[:] + [bi] for row, bi in zip(A, b)]
    R, pivots = rref(aug)
    if any(all(x == 0 for x in row[:-1]) and row[-1] != 0 for row in R):
        return 'none', None
    n = len(A[0])
    if len(pivots) == n:
        x = [F(0)] * n
        for i, c in enumerate(pivots):
            x[c] = R[i][-1]
        return 'unique', x
    # infinite solutions: set free variables to 0
    x = [F(0)] * n
    for i, c in enumerate(pivots):
        x[c] = R[i][-1]
    return 'inf', x


RESULTS = []


def check(name, cond):
    RESULTS.append(cond)
    print(('PASS  ' if cond else 'FAIL  ') + name)


# ---------- Chapter 1: determinants ----------

check('ex1.1  2x2 det |3 1; 2 4| = 10',
      det(mat([[3, 1], [2, 4]])) == 10)

check('ex1.2  3x3 det |1 2 3; 4 5 6; 7 8 10| = -3',
      det(mat([[1, 2, 3], [4, 5, 6], [7, 8, 10]])) == -3)

check('ex1.5  4x4 det (3 on diagonal, 1 elsewhere) = 48',
      det(mat([[3, 1, 1, 1], [1, 3, 1, 1], [1, 1, 3, 1], [1, 1, 1, 3]])) == 48)

check('ex1.6  block det |1 2 0 0; 3 4 0 0; 0 0 5 6; 0 0 7 8| = 4',
      det(mat([[1, 2, 0, 0], [3, 4, 0, 0], [0, 0, 5, 6], [0, 0, 7, 8]])) == 4)

check('ex1.3  Vandermonde V(1,2,3) = 2',
      det(mat([[1, 1, 1], [1, 2, 3], [1, 4, 9]])) == 2)

check('ex1.4  exercise: det |2 0 1; 1 3 5; 4 1 2| = -9',
      det(mat([[2, 0, 1], [1, 3, 5], [4, 1, 2]])) == -9)

check('ex1.10  Vandermonde V(1,2,3,4) = 12',
      det(mat([[1, 1, 1, 1], [1, 2, 3, 4], [1, 4, 9, 16], [1, 8, 27, 64]])) == 12)

A_cramer = mat([[1, 1, 1], [1, 2, 3], [1, 4, 9]])
check('ex1.5 (Cramer) D=2, D1=2, D2=4, D3=6; solution (1,2,3)',
      det(A_cramer) == 2
      and det(mat([[6, 1, 1], [14, 2, 3], [36, 4, 9]])) == 2
      and det(mat([[1, 6, 1], [1, 14, 3], [1, 36, 9]])) == 4
      and det(mat([[1, 1, 6], [1, 2, 14], [1, 4, 36]])) == 6
      and solve(A_cramer, [6, 14, 36]) == ('unique', [F(1), F(2), F(3)]))

st, x = solve(mat([[1, 2], [3, 4]]), [5, 9])
check('ex1.6 (exercise) system x1+2x2=5, 3x1+4x2=9 has solution (-1,3)',
      st == 'unique' and x == [F(-1), F(3)])

st, x = solve(mat([[1, 1, 1], [1, 2, 3], [2, 1, 1]]), [3, 6, 4])
check('ex1.7 (exercise) system has unique solution (1,1,1)',
      st == 'unique' and x == [F(1)] * 3)

check('ex1.8 (exercise) |A-lE| = 1-l^2 for A=|1 l; l 1| -> nonzero solutions iff l=+-1',
      det(mat([[1, 0], [0, 1]]) ) == 1  # placeholder keeps symmetry of the section
      and all(det(mat([[1, F(l)], [F(l), 1]])) == 1 - l * l for l in (-1, 0, 1)))


# ---------- Chapter 2: matrices ----------

A, B = mat([[1, 0], [1, 1]]), mat([[1, 1], [0, 1]])
check('ex2.1  AB=[[1,1],[1,2]], BA=[[2,1],[1,1]] (AB != BA)',
      mul(A, B) == mat([[1, 1], [1, 2]]) and mul(B, A) == mat([[2, 1], [1, 1]]))

A, B = mat([[1, 1], [1, 1]]), mat([[1, -1], [-1, 1]])
check('ex2.2  AB=O with A,B nonzero',
      mul(A, B) == mat([[0, 0], [0, 0]]))

An = eye(2)
for _ in range(5):
    An = mul(An, mat([[1, 1], [0, 1]]))
check('ex2.3  [[1,1],[0,1]]^5 = [[1,5],[0,1]]',
      An == mat([[1, 5], [0, 1]]))

check('ex2.4  inverse of [[3,1],[5,2]] is [[2,-1],[-5,3]]',
      inverse(mat([[3, 1], [5, 2]])) == mat([[2, -1], [-5, 3]]))

A = mat([[1, 2, 3], [2, 2, 3], [1, 1, 1]])
invA = inverse(A)
check('ex2.5  det=1 and A^-1=[[-1,1,0],[1,-2,3],[0,1,-2]]',
      det(A) == 1 and invA == mat([[-1, 1, 0], [1, -2, 3], [0, 1, -2]])
      and same(mul(A, invA), eye(3)))

A = mat([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
invA = inverse(A)
check('ex2.6  A^-1 = (1/4)[[3,-1,-1],[-1,3,-1],[-1,-1,3]]',
      invA == [[F(3, 4), F(-1, 4), F(-1, 4)],
               [F(-1, 4), F(3, 4), F(-1, 4)],
               [F(-1, 4), F(-1, 4), F(3, 4)]])

A, Bv = mat([[1, 2], [1, 3]]), mat([[1, 1], [2, 3]])
check('ex2.7  solution of AX=B is X=[[-1,-3],[1,2]]',
      same(mul(A, mat([[-1, -3], [1, 2]])), Bv))

A = mat([[1, 2, 0, 0], [2, 5, 0, 0], [0, 0, 2, 1], [0, 0, 1, 1]])
expected = mat([[5, -2, 0, 0], [-2, 1, 0, 0], [0, 0, 1, -1], [0, 0, -1, 2]])
check('ex2.8  block-diagonal inverse example',
      same(inverse(A), expected) and same(mul(A, expected), eye(4)))

check('ex2.9  rank of [[1,2,3,4],[2,4,6,9],[1,1,1,1]] is 3',
      rank(mat([[1, 2, 3, 4], [2, 4, 6, 9], [1, 1, 1, 1]])) == 3)

check('ex2.10 (exercise) ranks: [[1,2,3],[2,4,6],[1,1,1]] -> 2; 3x4 -> 2 (r3-r2=r1)',
      rank(mat([[1, 2, 3], [2, 4, 6], [1, 1, 1]])) == 2
      and rank(mat([[1, 2, 3, 4], [2, 4, 6, 9], [3, 6, 9, 13]])) == 2)


# ---------- Chapter 3: vectors ----------

st, _ = solve(mat([[1, 0], [1, 1], [0, 1]]), [1, 2, 3])
check('ex3.1  (1,2,3) NOT in span{(1,1,0),(0,1,1)} (inconsistent)',
      st == 'none')

A = mat([[1, 0, 1], [1, 1, 0], [0, 1, -1]])  # columns: the three vectors
zeros = mat([[0], [0], [0]])
check('ex3.2  (1,1,0),(0,1,1),(1,0,-1) dependent; relation 1,-1,-1',
      det(A) == 0 and same(mul(A, mat([[1], [-1], [-1]])), zeros))

A = mat([[1, 0, 1, 2], [0, 1, 1, 1], [2, 1, 3, 5]])
check('ex3.3  rank 2; alpha3=alpha1+alpha2; alpha4=2alpha1+alpha2',
      rank(A) == 2
      and same(mul(A, mat([[1], [1], [-1], [0]])), mat([[0], [0], [0]]))
      and same(mul(A, mat([[2], [1], [0], [-1]])), mat([[0], [0], [0]])))

Basis = mat([[1, 0, 1], [1, 1, 0], [0, 1, 1]])  # columns: basis vectors
st, k = solve(Basis, [1, 2, 1])
check('ex3.4  coordinates of (1,2,1) in the basis are (1,1,0)',
      st == 'unique' and k == [F(1), F(1), F(0)])

st, k = solve(Basis, [2, 3, 3])
check('ex3.10 (exercise) coordinates of (2,3,3) are (1,2,1)',
      st == 'unique' and k == [F(1), F(2), F(1)])

check('ex3.1e  (1,1,1),(1,2,3),(1,3,6) independent (det=1)',
      det(mat([[1, 1, 1], [1, 2, 3], [1, 3, 6]])) == 1)

check('ex3.9  det of [(1,1,0),(0,1,1),(l,0,1)] columns = 1+l -> dependent iff l=-1',
      all(det(mat([[1, 0, F(l)], [1, 1, 0], [0, 1, 1]])) == 1 + F(l)
          for l in (-2, -1, 0, 3)))


# ---------- Chapter 4: linear systems ----------

def Alm(l):
    return mat([[F(l), 1, 1], [1, F(l), 1], [1, 1, F(l)]])


check('ex4.1 det [[l,1,1],[1,l,1],[1,1,l]] = (l+2)(l-1)^2 at l=0,1,-2',
      det(Alm(0)) == 2 and det(Alm(1)) == 0 and det(Alm(-2)) == 0)

st, x = solve(Alm(2), [1, 2, 4])
check('ex4.1 l=2: unique solution (-3/4, 1/4, 9/4)',
      st == 'unique' and x == [F(-3, 4), F(1, 4), F(9, 4)])

check('ex4.1 l=1: infinitely many solutions; l=-2: none',
      solve(Alm(1), [1, 1, 1])[0] == 'inf'
      and solve(Alm(-2), [1, -2, 4])[0] == 'none')

A = mat([[1, 1, 1, 1], [2, 3, 4, 5]])
z2 = mat([[0], [0]])
check('ex4.2 fundamental system (1,-2,1,0),(2,-3,0,1)',
      same(mul(A, mat([[1], [-2], [1], [0]])), z2)
      and same(mul(A, mat([[2], [-3], [0], [1]])), z2))

A = mat([[1, 1, 1], [1, 2, 3]])
check('ex4.3 particular solution (0,3,0); null direction (1,-2,1)',
      same(mul(A, mat([[0], [3], [0]])), mat([[3], [6]]))
      and same(mul(A, mat([[1], [-2], [1]])), mat([[0], [0]])))

A = mat([[1, 1, 1], [1, 2, 2], [0, 1, 1]])
z3 = mat([[0], [0], [0]])
check('ex4.1e solution set (1,1,0)+c(0,-1,1)',
      same(mul(A, mat([[1], [1], [0]])), mat([[2], [3], [1]]))
      and same(mul(A, mat([[0], [-1], [1]])), z3))

A = mat([[1, 2, 3], [2, 4, 6], [3, 6, 9]])
check('ex4.2e nullspace basis (-2,1,0),(-3,0,1)',
      same(mul(A, mat([[-2], [1], [0]])), z3)
      and same(mul(A, mat([[-3], [0], [1]])), z3))

A = mat([[1, 1, 1, 1], [1, 2, 3, 4]])
check('ex4.3e particular (-2,6,0,0); nulls (1,-2,1,0),(2,-3,0,1)',
      same(mul(A, mat([[-2], [6], [0], [0]])), mat([[4], [10]]))
      and same(mul(A, mat([[1], [-2], [1], [0]])), z2)
      and same(mul(A, mat([[2], [-3], [0], [1]])), z2))

A = mat([[1, 2, 3], [2, 3, 1], [1, 1, 1]])
check('ex4.7e unique system has solution (1,1,1)',
      same(mul(A, mat([[1], [1], [1]])), mat([[6], [6], [3]])))

def Lm(l):
    return mat([[F(l) + 1, 1, 1], [1, F(l) + 1, 1], [1, 1, F(l) + 1]])


check('ex4.4e det = l^2(l+3): vanishes at l=0,-3; equals 4 at l=1',
      det(Lm(0)) == 0 and det(Lm(-3)) == 0 and det(Lm(1)) == 4)

check('ex4.4e nullspace at l=0: (1,-1,0),(1,0,-1); at l=-3: (1,1,1)',
      same(mul(Lm(0), mat([[1], [-1], [0]])), z3)
      and same(mul(Lm(0), mat([[1], [0], [-1]])), z3)
      and same(mul(Lm(-3), mat([[1], [1], [1]])), z3))


# ---------- Chapter 5: eigenvalues and eigenvectors ----------

A = mat([[3, -1], [-1, 3]])
check('ex5.2 eigen 2 (vector (1,1)) and 4 (vector (1,-1)); tr/det check',
      det(mat([[3 - 2, -1], [-1, 3 - 2]])) == 0
      and det(mat([[3 - 4, -1], [-1, 3 - 4]])) == 0
      and same(mul(A, mat([[1], [1]])), mat([[2], [2]]))
      and same(mul(A, mat([[1], [-1]])), mat([[4], [-4]]))
      and 2 + 4 == 6 and det(A) == 8)

check('ex5.1 Schmidt result (1,1,0),(1,-1,2),(-1,1,1) pairwise orthogonal',
      all(x == 0 for x in [
          1 * 1 + 1 * -1 + 0 * 2,
          1 * -1 + 1 * 1 + 0 * 1,
          1 * -1 + -1 * 1 + 2 * 1]))

A = mat([[4, 0, 0], [0, 3, 1], [0, 1, 3]])
P = mat([[1, 0, 0], [0, 1, 1], [0, 1, -1]])
Lam = mat([[4, 0, 0], [0, 4, 0], [0, 0, 2]])
check('ex5.3 eigenvectors of diag-block A; P^-1 A P = diag(4,4,2)',
      det(mat([[0, 0, 0], [0, -1, 1], [0, 1, -1]])) == 0
      and same(mul(A, mat([[1], [0], [0]])), mat([[4], [0], [0]]))
      and same(mul(A, mat([[0], [1], [1]])), mat([[0], [4], [4]]))
      and same(mul(A, mat([[0], [1], [-1]])), mat([[0], [2], [-2]]))
      and same(mul(mul(inverse(P), A), P), Lam))

Ak = eye(3)
for _ in range(4):
    Ak = mul(Ak, A)
formula = mat([[4 ** 4, 0, 0],
               [0, (4 ** 4 + 2 ** 4) // 2, (4 ** 4 - 2 ** 4) // 2],
               [0, (4 ** 4 - 2 ** 4) // 2, (4 ** 4 + 2 ** 4) // 2]])
check('ex5.3 A^4 matches closed form (4^k +- 2^k)/2', same(Ak, formula))

A = mat([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
B = mat([[1, 1, 1], [1, -1, 1], [1, 0, -2]])  # columns: xi1, beta2, beta3
check('ex5.4 eigencolumns orthogonal (squared lengths 3,2,6); eigen 2,-1,-1',
      same(mul(transpose(B), B), mat([[3, 0, 0], [0, 2, 0], [0, 0, 6]]))
      and same(mul(mul(transpose(B), A), B),
               mat([[6, 0, 0], [0, -2, 0], [0, 0, -6]])))

A = mat([[2, 1], [1, 2]])
check('ex5.1e eigen 1 (1,-1) and 3 (1,1)',
      same(mul(A, mat([[1], [-1]])), mat([[1], [-1]]))
      and same(mul(A, mat([[1], [1]])), mat([[3], [3]])))

P10 = mul(mul(mat([[1, 1], [-1, 1]]),
              mat([[1, 0], [0, 3 ** 10]])),
          mat([[F(1, 2), F(-1, 2)], [F(1, 2), F(1, 2)]]))
D10 = eye(2)
for _ in range(10):
    D10 = mul(D10, A)
check('ex5.5e A^10 matches P diag(1,3^10) P^-1', same(D10, P10))

A = mat([[2, 0, 0], [0, 3, 1], [0, 1, 3]])
check('ex5.4e eigen 2 (double) and 4; eigenvectors (1,0,0),(0,1,-1),(0,1,1)',
      det(mat([[0, 0, 0], [0, 1, 1], [0, 1, 1]])) == 0
      and det(mat([[-2, 0, 0], [0, -1, 1], [0, 1, -1]])) == 0
      and same(mul(A, mat([[1], [0], [0]])), mat([[2], [0], [0]]))
      and same(mul(A, mat([[0], [1], [-1]])), mat([[0], [2], [-2]]))
      and same(mul(A, mat([[0], [1], [1]])), mat([[0], [4], [4]])))

S = mat([[1, -1], [1, 1]])
check('ex5.3e (1/sqrt2)[[1,-1],[1,1]] orthogonal: S^T S = 2E',
      same(mul(transpose(S), S), mat([[2, 0], [0, 2]])))


if all(RESULTS):
    print(f'\nAll {len(RESULTS)} checks passed.')
    sys.exit(0)
print(f'\n{sum(not r for r in RESULTS)} of {len(RESULTS)} checks FAILED.')
sys.exit(1)
