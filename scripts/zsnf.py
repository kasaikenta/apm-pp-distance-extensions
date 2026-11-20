"""zsnf -- exact integer Smith normal form and COMPLETE solution modules mod m.

Motivation: apx_toy.integer_kernel returns a Q-kernel basis (primitive integer
vectors).  Reduced mod m those vectors generate only a subgroup of the true
solution module {x in Z_m^n : A x == 0 mod m}; this file computes the module
completely and exactly.

API
---
smith(A) -> (D, Vinv_cols)
    Decomposition A = U * D * V with U (rows x rows) and V (cols x cols)
    unimodular over Z, D diagonal, d_1 | d_2 | ... | d_r >= 1, then zeros.
    Returns:
      D          list of ALL min(rows, cols) diagonal entries (nonnegative,
                 divisibility chain, trailing zeros kept).
      Vinv_cols  V^{-1} as a LIST OF COLUMNS:  Vinv_cols[i][r] = (V^{-1})[r][i].
                 With the substitution y = V x the inverse map is
                 x = sum_i y_i * Vinv_cols[i].
    U is not returned (never needed for kernel work).
    Certificate checked exactly on EVERY call (AssertionError on violation):
      * final working matrix is diagonal with the divisibility chain,
      * P * A * Q == D_matrix where P (row-op product) and Q = V^{-1}
        (column-op product) were tracked alongside,
      * P * Pinv == I and Q * Qinv == I exactly, where Pinv/Qinv were built
        by applying the INVERSE of every elementary op -- an integer inverse
        is a complete unimodularity proof (det = +-1).

solve_mod(A, m) -> (gens, module_size)
    COMPLETE generating set of {x in Z_m^n : A x == 0 mod m} plus the exact
    number of solutions (python int).
    Method: y = V x.  Pivot i (d_i != 0): d_i y_i == 0 mod m forces
    y_i in (m // gcd(d_i, m)) * Z_m, i.e. gcd(d_i, m) choices; generator
    x = (m // gcd(d_i, m)) * Vinv_cols[i] mod m (skipped when gcd == 1, the
    generator would be 0).  Free coordinate i: m choices, generator
    Vinv_cols[i] mod m.  module_size = prod_pivots gcd(d_i, m) * m^(#free).

sample_uniform(gens, m, rng[, n]) -> uniform random element of the module
    Sum of rng.randrange(m) multiples of the generators mod m.  Uniform
    because Z_m^len(gens) -> module is a surjective group homomorphism.

Self-tests: run  python3 zsnf.py  -- brute-force verification, prints
PASS/FAIL per case, exit code 1 on any failure.
"""

from __future__ import annotations

from math import gcd


# ------------------------------------------------------------------------- SNF

def smith(A, verify=True):
    """Integer Smith normal form.  Returns (D, Vinv_cols); see module docstring."""
    rows = len(A)
    assert rows > 0, "empty matrix"
    cols = len(A[0])
    assert all(len(r) == cols for r in A), "ragged matrix"
    W = [[int(v) for v in r] for r in A]
    A0 = [r[:] for r in W]
    # P accumulates row ops (list of rows), Qc accumulates column ops as
    # COLUMNS: Qc[c] is column c of Q.  Invariant: W == P * A0 * Q.
    # Pinv/Qinv accumulate the INVERSE elementary ops so that P*Pinv == I and
    # Q*Qinv == I exactly (checked when verify=True): unimodularity proof.
    # Pinv is stored as COLUMNS (a row op on P is a column op on Pinv),
    # Qinv as ROWS (a column op on Q is a row op on Qinv).
    P = [[int(i == j) for j in range(rows)] for i in range(rows)]
    Pinvc = [[int(i == j) for i in range(rows)] for j in range(rows)]
    Qc = [[int(i == j) for i in range(cols)] for j in range(cols)]
    Qinv = [[int(i == j) for j in range(cols)] for i in range(cols)]

    def rowop_add(src, dst, c):          # row dst += c * row src
        W[dst] = [x + c * y for x, y in zip(W[dst], W[src])]
        P[dst] = [x + c * y for x, y in zip(P[dst], P[src])]
        # Pinv -> Pinv * E^{-1},  E = I + c e_dst e_src^T:  col src -= c*col dst
        Pinvc[src] = [x - c * y for x, y in zip(Pinvc[src], Pinvc[dst])]

    def rowop_swap(a, b):
        W[a], W[b] = W[b], W[a]
        P[a], P[b] = P[b], P[a]
        Pinvc[a], Pinvc[b] = Pinvc[b], Pinvc[a]

    def rowop_neg(a):
        W[a] = [-x for x in W[a]]
        P[a] = [-x for x in P[a]]
        Pinvc[a] = [-x for x in Pinvc[a]]

    def colop_add(src, dst, c):          # col dst += c * col src
        for i in range(rows):
            W[i][dst] += c * W[i][src]
        Qc[dst] = [x + c * y for x, y in zip(Qc[dst], Qc[src])]
        # Qinv -> F^{-1} * Qinv,  F = I + c e_src e_dst^T:  row src -= c*row dst
        Qinv[src] = [x - c * y for x, y in zip(Qinv[src], Qinv[dst])]

    def colop_swap(a, b):
        for i in range(rows):
            W[i][a], W[i][b] = W[i][b], W[i][a]
        Qc[a], Qc[b] = Qc[b], Qc[a]
        Qinv[a], Qinv[b] = Qinv[b], Qinv[a]

    def bq(a, b):
        """balanced quotient, b > 0: remainder a - q*b lies in [-b//2, b/2]"""
        return (a + (b >> 1)) // b

    mn = min(rows, cols)
    t = 0
    while t < mn:
        # pivot: nonzero entry of MINIMAL |value| (anti-swell), early out on 1
        pos = None
        best = None
        for i in range(t, rows):
            for j in range(t, cols):
                w = W[i][j]
                if w:
                    if best is None or abs(w) < best:
                        best = abs(w)
                        pos = (i, j)
                        if best == 1:
                            break
            if best == 1:
                break
        if pos is None:
            break                        # remaining submatrix is zero
        if pos[0] != t:
            rowop_swap(t, pos[0])
        if pos[1] != t:
            colop_swap(t, pos[1])
        if W[t][t] < 0:
            rowop_neg(t)
        while True:
            # clear column t below the pivot (Euclidean row steps)
            for i in range(t + 1, rows):
                while W[i][t]:
                    q = bq(W[i][t], W[t][t])
                    rowop_add(t, i, -q)
                    if W[i][t]:          # |remainder| <= pivot/2: swap it in
                        rowop_swap(t, i)
                        if W[t][t] < 0:
                            rowop_neg(t)
            # clear row t right of the pivot (Euclidean column steps);
            # a column swap here may refill column t below the pivot
            for j in range(t + 1, cols):
                while W[t][j]:
                    q = bq(W[t][j], W[t][t])
                    colop_add(t, j, -q)
                    if W[t][j]:
                        colop_swap(t, j)
                        if W[t][t] < 0:
                            rowop_neg(t)
            if any(W[i][t] for i in range(t + 1, rows)):
                continue
            if any(W[t][j] for j in range(t + 1, cols)):
                continue
            # divisibility fix-up: d_t must divide the whole remaining block
            d = W[t][t]
            bad = None
            for i in range(t + 1, rows):
                for j in range(t + 1, cols):
                    if W[i][j] % d:
                        bad = i
                        break
                if bad is not None:
                    break
            if bad is None:
                break
            rowop_add(bad, t, 1)         # pull the offending row up, re-clear
        if W[t][t] < 0:
            rowop_neg(t)
        t += 1

    D = [W[i][i] for i in range(mn)]

    if verify:
        for i in range(rows):
            for j in range(cols):
                assert (i == j) or W[i][j] == 0, "not diagonal"
        assert all(x >= 0 for x in D), "negative divisor"
        fz = next((k for k, x in enumerate(D) if x == 0), mn)
        assert all(x == 0 for x in D[fz:]), "zero before nonzero on diagonal"
        nz = [x for x in D if x]
        for a, b in zip(nz, nz[1:]):
            assert b % a == 0, "divisibility chain broken"
        # exact decomposition:  P * A0 * Q == W
        PA = [[sum(P[i][k] * A0[k][j] for k in range(rows)) for j in range(cols)]
              for i in range(rows)]
        for i in range(rows):
            for j in range(cols):
                v = sum(PA[i][k] * Qc[j][k] for k in range(cols))
                assert v == W[i][j], "P*A*Q != D"
        # unimodularity: exact integer inverses exist
        for i in range(rows):
            Pi = P[i]
            for j in range(rows):
                Pj = Pinvc[j]
                v = sum(Pi[k] * Pj[k] for k in range(rows))
                assert v == int(i == j), "P * Pinv != I"
        for i in range(cols):
            for j in range(cols):
                v = sum(Qc[k][i] * Qinv[k][j] for k in range(cols))
                assert v == int(i == j), "Q * Qinv != I"

    return D, Qc


# ------------------------------------------------------------------- solve mod

def solve_mod(A, m):
    """Complete generating set + exact size of {x in Z_m^n : A x == 0 mod m}."""
    cols = len(A[0])
    D, Qc = smith(A)
    gens = []
    size = 1
    for i in range(cols):
        d = D[i] if i < len(D) else 0
        if d:
            g = gcd(d, m)
            if g == 1:
                continue                 # y_i forced to 0: no generator
            size *= g
            s = m // g
            gens.append([(s * x) % m for x in Qc[i]])
        else:
            size *= m
            gens.append([x % m for x in Qc[i]])
    return gens, size


def sample_uniform(gens, m, rng, n=None):
    """Uniform random element of the module generated by gens in Z_m^len."""
    if not gens:
        return [0] * (n if n is not None else 0)
    v = [0] * len(gens[0])
    for g in gens:
        c = rng.randrange(m)
        if c:
            for k in range(len(v)):
                v[k] = (v[k] + c * g[k]) % m
    return v


# ------------------------------------------------------------------ self-tests

def _dp_count(A, m, n):
    """Exact solution count of A x == 0 mod m by syndrome dynamic programming
    (exhaustive over Z_m^n, organised column by column)."""
    rows = len(A)
    zero = (0,) * rows
    states = {zero: 1}
    for k in range(n):
        tab = [tuple((A[i][k] * v) % m for i in range(rows)) for v in range(m)]
        new = {}
        for s, c in states.items():
            for tv in tab:
                key = tuple((a + b) % m for a, b in zip(s, tv))
                new[key] = new.get(key, 0) + c
        states = new
    return states.get(zero, 0)


def _enum_count(A, m, n):
    """Plain product() enumeration (cross-check of the DP itself)."""
    from itertools import product
    cnt = 0
    for x in product(range(m), repeat=n):
        if all(sum(a * xi for a, xi in zip(row, x)) % m == 0 for row in A):
            cnt += 1
    return cnt


def _closure_size(gens, m, n):
    """BFS closure of the subgroup generated by gens in Z_m^n."""
    seen = {(0,) * n}
    stack = [(0,) * n]
    gt = [tuple(g) for g in gens]
    while stack:
        x = stack.pop()
        for g in gt:
            y = tuple((a + b) % m for a, b in zip(x, g))
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen)


def _selftest():
    import random
    import sys

    rng = random.Random(20260813)
    NMAX = {4: 8, 6: 8, 8: 7, 9: 6, 12: 6, 16: 5}   # m^n <= 3e6, cols <= 8
    MS = [4, 6, 8, 9, 12, 16]
    N_CASES = 30
    passed = 0
    failed = 0
    resampled = 0

    # fixed deterministic cases first (known answers by hand)
    fixed = [
        ([[2, 4], [4, 8]], 4),           # SNF diag(2,0)
        ([[6, 0], [0, 6]], 4),           # 6x==0 mod 4 -> x in {0,2}: 2*2 sols
        ([[1, 0], [0, 1]], 9),           # only 0
        ([[0, 0, 3]], 9),                # 3z==0 mod 9: 3 * 9^2
    ]
    cases = []
    for A, m in fixed:
        cases.append((A, m, len(A[0])))
    ci = 0
    while len(cases) < N_CASES + len(fixed):
        m = MS[ci % len(MS)]
        ci += 1
        n = rng.randint(2, NMAX[m])
        nr = rng.randint(1, 6)
        A = [[rng.randint(-9, 9) for _ in range(n)] for _ in range(nr)]
        if all(v == 0 for row in A for v in row):
            continue
        _, size = solve_mod(A, m)
        if size > 1_000_000:             # keep BFS closure tractable
            resampled += 1
            continue
        cases.append((A, m, n))

    for idx, (A, m, n) in enumerate(cases):
        gens, size = solve_mod(A, m)     # smith() self-verifies internally
        ok = True
        why = []
        # (a) exact count
        cnt = _dp_count(A, m, n)
        if cnt != size:
            ok = False
            why.append(f"count {cnt} != module_size {size}")
        if m ** n <= 40000:              # validate the DP itself on small cases
            cnt2 = _enum_count(A, m, n)
            if cnt2 != cnt:
                ok = False
                why.append(f"enum {cnt2} != dp {cnt}")
        # (b) every generator solves
        for g in gens:
            if any(sum(a * x for a, x in zip(row, g)) % m for row in A):
                ok = False
                why.append("generator does not solve")
                break
        # (c) BFS closure has exactly module_size elements
        cl = _closure_size(gens, m, n)
        if cl != size:
            ok = False
            why.append(f"closure {cl} != module_size {size}")
        # (d) uniform samples solve
        for _ in range(5):
            v = sample_uniform(gens, m, rng, n=n)
            if any(sum(a * x for a, x in zip(row, v)) % m for row in A):
                ok = False
                why.append("sample does not solve")
                break
        tag = "PASS" if ok else "FAIL"
        if ok:
            passed += 1
        else:
            failed += 1
        print(f"[{idx:02d}] m={m:2d} shape={len(A)}x{n} size={size:<9d} {tag}"
              + ("" if ok else "  <-- " + "; ".join(why)))

    print(f"SELFTEST: {passed}/{passed + failed} PASS"
          f" (random resamples due to >1e6 module: {resampled})")
    if failed:
        sys.exit(1)
    return passed, passed + failed


if __name__ == "__main__":
    _selftest()
