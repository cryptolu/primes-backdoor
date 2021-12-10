from math import log2

def mulZsqrtD(v1,v2,D):
    a1,b1 = v1
    a2,b2 = v2
    return (a1*a2 + b1*b2*D, a2*b1+a1*b2)

def gen_smooth(bits,FB,D):
    P = list(FB.keys())
    t,v = (D % 2) + abs(D), (D % 2, 1)
    while True:
        pi = random.choice(P)
        t *= pi
        v = mulZsqrtD(v,FB[pi],D)
        if (bits-log2(t)) > 0 and (bits-log2(t)) < 1:
            break 
        if log2(t) > bits:
            t,v = (D % 2) + abs(D), (D % 2, 1)
    a,b = v
    v = (abs(a),abs(b))
    return t,v

def genFB(bound, D):
    FB = {a**2 + abs(D)*b**2 : (a,b) for a in range(bound) for b in range(bound) if (is_prime(a**2 + abs(D)*b**2)) and ((a**2 + abs(D)*b**2)%4 == 1 if (D == -4) else 1)}
    FB[(D % 2) + abs(D)] = (D % 2, 1)
    return FB

def check_prime(t,v,D, safe_prime=False):
    a,b = v
    if D != -4:
        p = ((a+1)**2 + abs(D)*b**2)
        if is_prime(p) and (is_prime((p-1)//2) if safe_prime else 1):
            return p, (a+1,b)
        p = ((a-1)**2 + abs(D)*b**2)
        if is_prime(p) and (is_prime((p-1)//2) if safe_prime else 1):
            return p, (a-1,b)
    if D == -4:
        #we absorb |D| in b, so that t = a^2+b^2
        b = 2*b
        p = ((a+1)**2 + b**2)
        if gp.isprime(p) and (gp.isprime((p-1)//2) if safe_prime else 1):
            return p, (a+1,b)
        p = ((a-1)**2 + b**2)
        if gp.isprime(p) and (gp.isprime((p-1)//2) if safe_prime else 1):
            return p, (a-1,b)
        p = (a**2 + (b+1)**2)
        if gp.isprime(p) and (gp.isprime((p-1)//2) if safe_prime else 1):
            return p, (a,b+1)
        p = (a**2 + (b-1)**2)
        if gp.isprime(p) and (gp.isprime((p-1)//2) if safe_prime else 1):
            return p, (a,b-1)
    return None 

def gen_prime(bits,FB,D,safe_prime=False):
    while True:
        t,v = gen_smooth(bits,FB,D)
        a,b = v
        if ((t % 2) == 0):
            res = check_prime(t,v,D,safe_prime=safe_prime)
            if res != None:
                break
    p,v = res
    return res

############################

D = -4
assert (D % 4 == 0) or (D % 4  == 1)
bound = 2**8
bits = 128

FB = genFB(bound,D)

p,v = gen_prime(bits,FB,D,safe_prime=False)
a,b = v

if D != -4:
    assert (p == a**2 + abs(D)*b**2)
else:
    #we rewrite a,b as in Washington's Theorem
    if a % 2 == 0:
        a,b = b,a
    if (a+b) % 4 == 3:
        a = -a
    assert(((a+b) % 4 == 1) and (b%2 == 0))
    assert(p == a**2 + b**2)

p,a,b