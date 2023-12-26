# cado-nfs-discrete-log

<br>

<https://gitlab.inria.fr/cado-nfs/cado-nfs>
<br>

<https://gitlab.inria.fr/cado-nfs/cado-nfs/-/blob/master/README.dlp>
<br>
This is a helper script to compute discrete logarithms using the cado-nfs sieve software. 


<br>

An experiment to test it's speed for differet max sized factors of the order for cado-nfs vs sage's baby step giant step:

```python
import subprocess
from Crypto.Util.number import *
from random import randint
from math import prod, isqrt, gcd
from sympy.ntheory.modular import crt
from collections import defaultdict
import time
cado_path = "/home/connor/Documents/cado-nfs/cado-nfs.py"

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
    assert 2**91 < p < 2**330
    lb1, log_a = cado_dlog_unknown_base(p, ell, a)
    lb2, log_b = cado_dlog_unknown_base(p, ell, base)
    assert lb1 == lb2
    dlog = (log_a * pow(log_b, -1, ell)) % ell
    return dlog


def n_order(a, n, isprime=False):
    if gcd(a, n) != 1:
        raise ValueError("The two numbers should be relatively prime")
    factors = defaultdict(int)
    if isprime:
        t1 = {n: 1}
    else:
        t1 = dict(factor(n)).items()
    for px, kx in t1.items():
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

def pohlig_hellman(n, a, b, into, order=None, f=None):
    if order is None:
        order, ff = n_order(b, n)
    if f is None:
        f = dict(factor(order))
    l = [0] * len(f)

    for i, (pi, ri) in enumerate(f.items()):
        #print(f'factor {pi}:')
        for j in range(ri):
            gj = pow(b, l[i], n)
            aj = pow(a * pow(gj, -1, n), order // pi**(j + 1), n)
            bj = pow(b, order // pi, n)
            if into == "bsgs":
                cj = GF(n)(aj).log(bj)
            elif into == "cado":
                cj = cado_dlog(n, pi, aj, bj)
            assert 1 == pow(pow(bj, cj, n) * pow(aj, -1, n), order//pi, n)
            l[i] += cj * pi**j

    d, _ = crt([pi**ri for pi, ri in f.items()], l)
    return d


def get_prime(mxsize):
    while True:
        factors = {getPrime(i):1 for i in range(mxsize-200//mxsize, mxsize)}
        factors[2] = 1
        n = prod([(a**b) for a, b in factors.items()]) + 1
        if isPrime(n):
            print(factors)
            print(n)
            return n


for s in range(50, 60):
    print(f"\n\n{s} bit factor:")
    n = get_prime(s)
    b = 3
    order, f = n_order(b, n, isprime=True)
    d = randint(1, order)
    a = pow(b, d, n)

    start = time.time()
    assert d == pohlig_hellman(n, a, b, "cado", order, f)
    end = time.time()
    print(f"cado: {(end - start)/60:0.1f} minutes")

    start = time.time()
    assert d == pohlig_hellman(n, a, b, "bsgs", order, f)
    end = time.time()
    print(f"bsgs: {(end - start)/60:0.1f} minutes")

```

<br>

```
40 bit factor:
cado: 11.8 minutes
bsgs: 0.1 minutes


45 bit factor:
cado: 8.1 minutes
bsgs: 0.2 minutes


50 bit factor:
cado: 9.5 minutes
bsgs: 1.0 minutes


51 bit factor:
cado: 4.1 minutes
bsgs: 0.9 minutes


52 bit factor:
cado: 4.6 minutes
bsgs: 1.6 minutes


53 bit factor:
cado: 4.5 minutes
bsgs: 2.1 minutes


54 bit factor:
cado: 5.3 minutes
bsgs: 3.2 minutes


55 bit factor:
cado: 4.4 minutes
bsgs: 2.6 minutes


56 bit factor:
cado: 5.4 minutes
bsgs: 6.5 minutes


57 bit factor:
cado: 6.3 minutes
bsgs: 5.5 minutes


58 bit factor:
cado: 6.1 minutes
bsgs: 8.2 minutes


59 bit factor:
cado: 6.4 minutes
bsgs: 22.7 minutes


60 bit factor:
cado: 8.4 minutes
bsgs: 44.0 minutes
```
