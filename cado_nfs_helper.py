from sage.all import *
import subprocess
from sympy.ntheory.modular import crt
from math import gcd
from collections import defaultdict
cado_path = "/home/connor/Documents/cado-nfs/cado-nfs.py"

def n_order(a, n):
    if gcd(a, n) != 1:
        raise ValueError("The two numbers should be relatively prime")
    factors = defaultdict(int)
    for px, kx in dict(factor(n)).items():
        if kx > 1:
            factors[px] += kx - 1
        fpx = dict(factor(px - 1))
        for py, ky in fpx.items():
            factors[py] += ky
    group_order = 1
    for px, kx in factors.items():
        group_order *= px**kx
    order = 1
    ff = {}
    if a > n:
        a = a % n
    for p, e in factors.items():
        exponent = group_order
        for f in range(e + 1):
            if pow(a, exponent, n) != 1:
                order *= p ** (e - f + 1)
                ff[p] = e - f + 1
                break
            exponent = exponent // p
    return order, ff
    
def cado_dlog_unknown_base(p, ell, target):
    out = subprocess.run(["python", f"{cado_path}", "-dlp", 
                          f"{ell=}", f"{target=}", f"{p}"], 
                         capture_output=True)
    for line in out.stderr.decode().splitlines():
        logbase_ = line.split("logbase = ")
        logtarget_ = line.split("log(target) = ")
        if len(logbase_) > 1:
            logbase = int(logbase_[-1])
        if len(logtarget_) > 1:
            logtarget = int(logtarget_[-1].split()[0])

    assert 1 == pow(pow(logbase, logtarget, p) * pow(target, -1, p), (p-1)//ell, p)
    return logbase, logtarget

def cado_dlog(p, ell, a, base):
    lb1, log_a = cado_dlog_unknown_base(p, ell, a)
    lb2, log_b = cado_dlog_unknown_base(p, ell, base)
    assert lb1 == lb2
    dlog = (log_a * pow(log_b, -1, ell)) % ell
    return dlog

def pohlig_hellman(n, a, b, order=None):
    if order is None:
        order, f = n_order(b, n)
    l = [0] * len(f)

    for i, (pi, ri) in enumerate(f.items()):
        print(f'factor {pi}:')
        for j in range(ri):
            gj = pow(b, l[i], n)
            aj = pow(a * pow(gj, -1, n), order // pi**(j + 1), n)
            bj = pow(b, order // pi, n)
            if pi < 2**58:
                if is_prime(n):
                    cj = GF(n, proof=False)(aj).log(bj)
                else:
                    cj = Zmod(n)(aj).log(bj)
            else:
                cj = cado_dlog(n, pi, aj, bj)
            assert 1 == pow(pow(bj, cj, n) * pow(aj, -1, n), order//pi, n)
            l[i] += cj * pi**j

    d, _ = crt([pi**ri for pi, ri in f.items()], l)
    return d
