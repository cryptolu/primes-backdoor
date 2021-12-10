#XZ-arithmetic using Montgomery Ladder
#Alg. 3, A.3, A.5 from "Fast Elliptic Curve Multiplications Resistant against Side Channel Attacks" by Tetsuya IZU and Tsuyoshi TAKAGI
def xECDBL(P,EC):
	A,B = EC
	R1,R2 = P
	#########
	R3 = R1^2
	R4 = R2^2
	R5 = R4*A
	R6 = R4*B
	R4 = R6*R4
	R2 = R1*R2
	R1 = R3-R5
	R5 = R3+R5
	R1 = R1^2
	R3 = R6*R2
	R3 = R3+R3
	R3 = R3+R3
	R3 = R3+R3
	R1 = R1-R3
	R5 = R5*R2
	R4 = R5+R4
	R4 = R4+R4
	R4 = R4+R4
	return (R1,R4)

def xECADD(P,Q,x,EC):
	A,B = EC
	R1,R2 = P
	R3,R4 = Q
	#########
	R5 = R1*R3
	R6 = R2*R4
	R1 = R1*R4
	R2 = R2*R3
	R3 = A*R6
	R4 = R5-R3
	R4 = R4^2
	R5 = B*R6
	R5 = R5+R5
	R5 = R5+R5
	R6 = R1+R2
	R6 = R5*R6
	R6 = R4-R6
	R5 = R1-R2
	R5 = R5^2
	R5 = R5*x
	return (R6,R5)

def ladder(k,P,EC):
	d = k.bits()
	n = len(d)
	Q = [0,0,0]
	Q[0] = P
	Q[1] = xECDBL(P,EC)
	for i in range(n-2,-1,-1):
		Q[2] = xECDBL(Q[d[i]],EC)
		Q[1] = xECADD(Q[0], Q[1], P[0], EC)
		Q[0] = Q[2-d[i]]
		Q[1] = Q[1+d[i]]
	return Q[0]

#Inversion in ZN(j)
def inv(v):
	R = parent(v)
	P.<x> = Integers()[]
	g,a,b = xgcd(P(list(v)),P(R.modulus()))	
	try:
		assert(g.degree() == 0)
	except:
		print(v, "is a divisor of zero. A factor of HD(x) is ", g)
	return R(a*inverse_mod(ZZ(g),N))

def get_curve(N,D):
	Zn = Integers(N)
	if D == -3:
		c = Zn.random_element()
		A = 0
		B = c**3
		print("c : ", c)
		return ([A,B],Zn)
	if D == -4:
		c = Zn.random_element()
		A = -c
		B = 0
		print("c : ", c)
		return ([A,B],Zn)
	else:
		PZn.<x> = PolynomialRing(Zn)
		HD = hilbert_class_polynomial(D)
		R.<j> = QuotientRing(PZn, PZn.ideal(HD))
		k = j*inv(j-1728)
		c = Zn.random_element()
		A = -3*k*c**2
		B = 2*k*c**3
		return ([A,B],R)

def random_twist(N,EC):
	A,B = EC
	while True:
		c = Integers(N).random_element()
		if c != 0:
			break
	if D == -4:
		A = -c
		B = 0
	else:
		A *= c**2
		B *= c**3
	return ([A,B])

def genFB(bound, D):
    FB = {a**2 + abs(D)*b**2 : (a,b) for a in range(bound) for b in range(bound) if (is_prime(a**2 + abs(D)*b**2)) and ((a**2 + abs(D)*b**2)%4 == 1 if (D == -4) else 1)}
    FB[(D % 2) + abs(D)] = (D % 2, 1)
    return FB

####################

#p = 10069891168272853289682414533444101158961971160721810960425299727500040856155399910408097910806556095082529154204912221872656086515020034916524535013055607
#q = 9163585010630240914315854502079840960495506763984904313022218586074134346260658498746192771596180175665775150505169330955024838846917582563246839347729083
#N = p*q
N = 108218055286220658892305686966796450005656093853561880573303757779870852502715750896717380601416332282635416751505870710179577464401193642824027551121223839066820321283484303809645541513752756237536860056248506557928305788452287016998812668348785242282044946059872302424866654255328887383184108204782014745221

D = -107
bound = 2^8

FB = genFB(bound,D)
ell = 2
k = prod(FB.keys())**ell

Zn = Integers(N)
PZn.<x> = PolynomialRing(Zn)
HD = PZn(hilbert_class_polynomial(D))

EC,R = get_curve(N,D)
while True:
	P = (Zn.random_element(), Zn(1))
	Q = ladder(k,P,EC)
	Qz = PZn(list(R(Q[1])))
	g = gcd(Qz.resultant(HD),N) 
	if g not in [1,N]:
		break
	else:
		EC = random_twist(N,EC)

p,q = ZZ(g),N//ZZ(g)
assert N == p*q

print("p = ",p)
print("q = ",q)