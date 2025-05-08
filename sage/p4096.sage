import time

def cornacchiax(d,m): # finds the x-coordinate of a primitive solution to x^2 + dy^2 = m
    if not m > 0:
        return 0
    if not mod(-d,m).is_square():
        return 0
    r = int(mod(-d,m).sqrt())
    a = m
    b = r
    while  a^2 >= m:
        c = a%b
        a = b
        b = c
    s2 = (m-a^2)/d # if s2 = s^2 is a square, then a^2 + ds^2 = m is a primitive solution
    if s2.is_square():
        return a
    return 0

B = ell = 3

ells = [3]

while n(log(B,2)) <= 4110:
    ell = ell.next_prime()
    ells.append(ell)
    B *= ell

ells.remove(1033)
ells.remove(1187)
number_of_odd_primes = len(ells)
M = 4*prod(ells)
p = M-1
Fp = GF(p)

print('the value of p is',p)
print('the bitsize of p is',n(log(p,2)))
print('the number of odd ell_i is',number_of_odd_primes)
print('the largest ell_i is',max(ells))

# We want to find an element of End(E) of norm prod_i l_i^2.
# Specifically, we look for an element sigma of the form
# sigma = a + b*i + c*(i+pi)/2 + d*(1+i*pi)/2,
# where x := tr(sigma) = 2*a+d is as small as possible.
# We do this by solving a quaternion norm equation; this happens in the outer while loop.
# Specifically, the norm equation is (2a+d)^2 + (2b+c)^2 + p(c^2+d^2) = 4M^2.
# Substituting x = 2a+d and y = 2b+c, we obtain the equation x^2 + y^2 + p(c^2+d^2) = 4M^2.
# This may be solved by trying successive (small) values for x and using Cornacchia's algorithm.
# Since p%4 = 3, we may ensure that x is odd, y is even, c is even, and d is odd.
# The values of a and b are then obtained as a = (x-d)/2 and b = (y-c)/2.

# Additionally, we would like the discriminant of sigma not to contain small prime factors.
# This is ensured by the inner while loop.
# Note that Disc(sigma) = tr(sigma)^2 - 4*N(sigma) = -(2*M-x)(2*M+x).
# We will thus restrict ourselves to the case that 2*M-x and 2*M+x are both prime.
# Since M is divisible by each ell_i, we know that x must be coprime to ell_i for every i.
# We simply try successive prime numbers x > max_i(ell_i).

target_norm = M^2

x = 366593496 #precomputed value

while True:
    while True:
        x = x.next_prime()
        factor_1 = 2*M + x
        factor_2 = 2*M - x
        if factor_1.is_prime() and factor_2.is_prime():
            break
    
    y_squared_mod_p = Fp(4*target_norm-x^2)
    if y_squared_mod_p.is_square():
        y_mod_p = Fp(y_squared_mod_p).sqrt()
        y = ZZ(y_mod_p)
        if is_odd(y):
            y = p-y
        m = (4*target_norm - (x^2 + y^2)) // p
        if m%4 == 1 and m.is_prime():
            start = time.time()
            c = cornacchiax(1,m)
            end = time.time()
            print('cornacchias algorithm took time', end-start)
            if not c == 0:
                d = ZZ((m-c^2).sqrt())
                if is_odd(c):
                    c,d = d,c
                break
                
a = ZZ((x-d)/2)
b = ZZ((y-c)/2)

sigma = [a,b,c,d]
sigma_dual = [a+d,-b,-c,-d]

print('Orientation found.')

print('sigma =', sigma)

print('the trace of sigma is',x)

Disc_sigma = (2*M + x)*(2*M - x)

print('The bitsize of the discriminant of sigma is',n(log(Disc_sigma,2)))

print('We will now construct the public parameters of the scheme')

print('That is, we will construct a basis for the endomorphism given by sigma on the base curve...')


R.<X> = Fp[]
Fq.<i> = Fp.extension(X^2+1)
E = EllipticCurve([Fq(1),0])

# Functions to evaluate endomorphisms from Z[1,i,pi,i*pi] on elements of E.

def aut_i(P):
    E = P.curve()
    if P == E(0):
        return P
    else:
        xP = P[0]
        yP = P[1]
        xQ = -xP
        yQ = i*yP
        return E(xQ,yQ)

def frob_p(P):
    E = P.curve()
    if P == E(0):
        return P
    else:
        xP = P[0]
        yP = P[1]
        xQ = xP^p
        yQ = yP^p
        return E(xQ,yQ)

def endo_basic(vec,P):
    
    # Input: vec = [a,b,c,d] in Z^4, P in E.
    # Output: the evaluation of a+b*i+c*pi+d*i*pi at P.
    
    [a,b,c,d] = vec
    iP = aut_i(P)
    jP = frob_p(P)
    kP = aut_i(jP)
    Q = a*P + b*iP + c*jP + d*kP
    return Q

# Evaluate endomorphisms from the maximal order Z[1,i,(i+pi)/2,(1+i*pi)/2].

def endo_max(vec,P):
    
    # Input: vec = [a,b,c,d] in Z^4, and P in E.
    # Output: the evaluation of a+b*i+c*(i+pi)/2+d*(1+i*pi)/2 in End(E) at P.
    
    E = P.curve()
    [a,b,c,d] = vec
    N = order(P)
    if is_odd(N): # In this case, dividing by two is easy.
        inverse_of_2_mod_N = mod(2,N)^(-1)
        P_divided_by_two = ZZ(inverse_of_2_mod_N) * P
    else:
        divide_by_two_field.<X> = P.division_points(2,poly_only=True).splitting_field()
        Eq2 = E.base_extend(divide_by_two_field)
        P_divided_by_two = (Eq2(P).division_points(2))[0]
    Q = endo_basic([2*a+d,2*b+c,c,d], P_divided_by_two)
    return E(Q)

# Find a basis P,Q of E[ell] such that P in ker(sigma) and Q in ker(sigma_dual)

def find_kernel_basis_max(sigma, ell):

    [a,b,c,d] = sigma
    
    sigma_dual = [a+d,-b,-c,-d]

    while True:
        R = E.random_point()
        cofactor = (p+1) // ell
        R *= cofactor
        P = endo_max(sigma_dual, R)
        if order(P) == ell:
            break
    
    while True:
        R = E.random_point()
        cofactor = (p+1) // ell
        R *= cofactor
        Q = endo_max(sigma, R)
        if order(Q) == ell:
            break

    return P,Q

# We construct a basis P,Q for E[M] (where M = 4 * prod_i ell_i),
# such that P in ker(sigma) and Q in ker(sigma_dual).

P = E(0)
Q = E(0)

for ell in ells + [4]:
    P0,Q0 = find_kernel_basis_max(sigma, ell)
    P += P0
    Q += Q0

assert order(P) == M
assert order(Q) == M
assert endo_max(sigma, P) == E(0)
assert endo_max(sigma_dual, Q) == E(0)

print('public parameters found!')

print('executing an example key exchange...')

# Using Kummer arithmetic

import sys

sys.path.append('KummerIsogeny')

from kummer_line import KummerLine
from kummer_isogeny import KummerLineIsogeny

# Public parameters

K = KummerLine(E)
xP = K(P)
xQ = K(Q)

# Compute the "companion" elliptic curve; i.e. [1,1,1,...,1]E,
# together with a generator for the dual isogeny.
# This is also part of the public parameters.

def compute_companion_curve(K,xP,xQ,codomain=False):
    
    cofactor = M
    
    for i in range(number_of_odd_primes):
        ell = ells[i]
        cofactor = cofactor // ell
        kernel_point = cofactor * xP
        phi = KummerLineIsogeny(K, kernel_point, ell)
        K = phi.codomain()
        xP = phi(xP)
        xQ = phi(xQ)

    ell = 4
    cofactor = cofactor // ell
    E = K.curve()
    kernel_point = (cofactor * xP).curve_point()
    if not codomain:
        phi = E.isogeny(kernel_point, model="montgomery")
    else:
        phi = E.isogeny(kernel_point, codomain=codomain.curve()) 
    K = KummerLine(phi.codomain())
    xQ = K(phi(xQ.curve_point()))
    
    return K,xQ

K_companion,xQ_companion = compute_companion_curve(K,xP,xQ)
K_companion_left,xP_companion = compute_companion_curve(K,xQ,xP,codomain=K_companion)

assert K_companion.a() == K_companion_left.a()

# We compute the group action for a "square-free" set of exponents vec
# That is, the entries of vec are either 0 or 1

def group_action_square_free(K,K_companion,xP,xQ,xP_companion,xQ_companion,vec):
    
    K_left = K_companion
    K_right = K

    left_kernel = xQ_companion
    right_kernel = xP
    
    cofactor = M
    
    for i in range(number_of_odd_primes):
        ell = ells[i]
        exponent = vec[i]
        cofactor = cofactor // ell
        if exponent == 0:
            
            # In this case, we have to walk in the left direction
            # The kernel of the isogeny is generated by kernel_point, which is of order ell:
            
            kernel_point = cofactor * left_kernel
            phi = KummerLineIsogeny(K_left, kernel_point, ell)
            K_left = phi.codomain()
            
            # We push through the point of order M that walks in the right (i.e. opposite) direction
            
            xP_companion = phi(xP_companion)
            
            # We compute the new points walking left and right
            # Their order has been reduced by a factor of ell
            
            left_kernel = phi(left_kernel)
            right_kernel *= ell
            
        if exponent == 1:
            
            # In this case, we have to walk in the right direction
            # The kernel of the isogeny is generated by kernel_point, which is of order ell:
            
            kernel_point = cofactor * right_kernel
            phi = KummerLineIsogeny(K_right, kernel_point, ell)
            K_right = phi.codomain()
            
            # We push through the point of order M that walks in the left (i.e. opposite) direction
            
            xQ = phi(xQ)
            
            # We compute the new points walking right and left
            # Their order has been reduced by a factor of ell
            
            right_kernel = phi(right_kernel)
            left_kernel *= ell

    ell = 4
    cofactor = cofactor // ell
    
    if exponent == 0:

        E = K_left.curve()
        kernel_point = (cofactor * left_kernel).curve_point()
        phi = E.isogeny(kernel_point, codomain=K_right.curve())
        K = KummerLine(phi.codomain())
        xP = K(phi(xP_companion.curve_point()))
        xQ = xQ
        
    if exponent == 1:
        
        E = K_right.curve()
        kernel_point = (cofactor * right_kernel).curve_point()
        phi = E.isogeny(kernel_point, codomain=K_left.curve())
        K = KummerLine(phi.codomain())
        xQ = K(phi(xQ.curve_point()))
        xP = xP_companion
            
    return K, xP, xQ

# Now also outputs the companion curve:

def group_action_square_free_companion(K,K_companion,xP,xQ,xP_companion,xQ_companion,vec):
    
    (K, xP, xQ),\
    (K_companion, xP_companion, xQ_companion) = \
    (group_action_square_free(K,K_companion,xP,xQ,xP_companion,xQ_companion,vec)), \
    (group_action_square_free(K_companion,K,xP_companion,xQ_companion,xP,xQ,vec))

    return K, K_companion, xP, xQ, xP_companion, xQ_companion

# Computing the full group action by breaking it up into "square-free" steps.

def group_action(K,K_companion,xP,xQ,xP_companion,xQ_companion,vec):

    while sum(vec) > 0:
        vec_square_free = [min(1,x) for x in vec]
        K,K_companion,xP,xQ,xP_companion,xQ_companion = \
        group_action_square_free_companion(K,K_companion,xP,xQ,xP_companion,xQ_companion,vec_square_free)
        vec = [vec[i] - vec_square_free[i] for i in range(len(vec))]

    return K,K_companion,xP,xQ,xP_companion,xQ_companion

secret_alice = [randint(0,2) for ell in ells+[4]]
secret_bob = [randint(0,2) for ell in ells+[4]]

print("Alice's secret key is:", secret_alice)
print("Bob's secret key is:", secret_bob)

start = time.time()

public_key_alice = group_action(K,K_companion,xP,xQ,xP_companion,xQ_companion,secret_alice)
[KA, KA_companion, xPA, xQA, xPA_companion, xQA_companion] = public_key_alice

end = time.time()

print("Alice's public key is:", KA.a(),KA_companion.a(),xPA.XZ(),xQA.XZ(),xPA_companion.XZ(),xQA_companion.XZ())

print("This took time:", end-start)

start = time.time()

public_key_bob = group_action(K,K_companion,xP,xQ,xP_companion,xQ_companion,secret_bob)
[KB, KB_companion, xPB, xQB, xPB_companion, xQB_companion] = public_key_bob

end = time.time()

print("Bob's public key is:", KB.a(),KB_companion.a(),xPB.XZ(),xQB.XZ(),xPB_companion.XZ(),xQB_companion.XZ())
print("This took time:", end-start)

start = time.time()

shared_secret_alice = group_action(KB, KB_companion, xPB, xQB, xPB_companion, xQB_companion, secret_alice)

end = time.time()

print("Alice's shared secret is:", shared_secret_alice[0].a())
print("This took time:", end-start)

start = time.time()

shared_secret_bob = group_action(KA, KA_companion, xPA, xQA, xPA_companion, xQA_companion, secret_bob)

end = time.time()

print("Bob's shared secret is:", shared_secret_bob[0].a())
print("This took time:", end-start)

print("The shared secrets are equal:", shared_secret_alice[0].a() == shared_secret_bob[0].a())

