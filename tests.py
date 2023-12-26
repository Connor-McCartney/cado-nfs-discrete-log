from cado_nfs_helper import n_order, pohlig_hellman
from random import randint

def test1():
    # ~15 mins
    p = 33184772290615481426295675425316668758122179640330548849957081783509
    g = 5
    order, f = n_order(g, p, isprime=True)
    d = randint(0, order)
    a = pow(g, d, p)
    print(d == pohlig_hellman(p, a, g, order, f))

test1()
