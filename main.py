from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair
from openpyxl import Workbook
import numpy as np
from functools import reduce
class GS():
    def __init__(self, groupObj):
        global util, group
        group = groupObj

    def Trusted_Setup(self,pp):
    # To sample four different scalars from \Z_p.
        rho, zeta, sigma, omega = [group.random() for _ in range(4)]
        vv1 = [pp['G1']**zeta, pp['G1']]
        vv2 = [pp['G2']**omega, pp['G2']]
        ww1 = [pp['G1']**(rho*zeta), pp['G1']**rho]
        ww2 = [pp['G2']**(sigma*omega), pp['G2']*sigma]
        Zeta = [-zeta**(-1), 1]
        Omega = [-omega**(-1), 1]
        crs = {'vv1':vv1, 'vv2':vv2, 'ww1':ww1, 'ww2':ww2}
        # trapdoors: one can use public randomness techniques to avoid them.
        tpd = {'crs':crs, 'Zeta':Zeta, 'Omega': Omega}
        return crs, tpd

    def Transpatent_Setup(self,pp):
        vv1 = [group.random(G1), pp['G1']]
        vv2 = [group.random(G2), pp['G2']]
        ww1 = [group.random(G1), pp['G1']]
        ww2 = [group.random(G2), pp['G2']]
        crs = {'vv1':vv1, 'vv2':vv2, 'ww1':ww1, 'ww2':ww2}
        tpd = {'empty'}
        return crs, tpd
    
    def commit(self, crs , X, Y, C_x, C_y):
        com_x = {}; com_y = {}
        N = len(C_x) #= len(C_y)
        # First sample fresh randomness for the witnesses and assign 0 for public variables.
        r = list(map(lambda i: [0,0] if C_x[i] == 0 \
                    else [group.random(),group.random()], range(N)))
        s = list(map(lambda i: [0,0] if C_y[i] == 0 \
                    else [group.random(),group.random()], range(N)))
        # After forming vectors r and s, compute the commitments to both X and Y elements.
        com_x = [[(crs['vv1'][0]**r[i][0])*(crs['ww1'][0]**s[i][0]),\
                X[i]*(crs['vv1'][1]**r[i][0])*(crs['ww1'][1]**s[i][0])] for i in range(N)]
        com_y = [[(crs['vv2'][0]**r[i][1])*(crs['ww2'][0]**s[i][1]),\
                Y[i]*(crs['vv2'][1]**r[i][1])*(crs['ww2'][1]**s[i][1])] for i in range(N)]
        return com_x, com_y, r, s
    def prove(self, crs, X, r, s, com_y):
        # Define prod function that multiplies elements of any given list
        def prod(list):
            result = 1
            for x in list:
                result *= x
            return result
        # Generate random values alpha, beta, gamma, and delta
        alpha, beta, gamma, delta = [group.random() for _ in range(4)]
        # Compute proof components; pi_v1, pi_v2, pi_w1, and pi_w2
        pi_v1 = [prod([com_y[i][0]**r[i][0] for i in range(len(X))]) * crs['vv2'][0]**alpha * crs['ww2'][0]**beta,\
                prod([com_y[i][1]**r[i][0] for i in range(len(X))]) * crs['vv2'][1]**alpha * crs['ww2'][1]**beta]
        pi_w1 = [prod([com_y[i][0]**s[i][0] for i in range(len(X))]) * crs['vv2'][0]**gamma * crs['ww2'][0]**delta,\
                prod([com_y[i][1]**s[i][0] for i in range(len(X))]) * crs['vv2'][1]**gamma * crs['ww2'][1]**delta]
        pi_v2 = [crs['vv1'][0]**-alpha * crs['ww1'][0]**(-gamma),\
                prod([X[i]**r[i][1] for i in range(len(X))]) * crs['vv1'][1]**-alpha * crs['ww1'][1]**-gamma]
        pi_w2 = [crs['vv1'][0]**-beta * crs['ww1'][0]**(-delta),\
                prod([X[i]**s[i][1] for i in range(len(X))]) * crs['vv1'][1]**-beta * crs['ww1'][1]**-delta]
        # Return the proof values as a dictionary with 8 group elements
        return {'pi_v1': pi_v1, 'pi_w1': pi_w1, 'pi_v2': pi_v2, 'pi_w2': pi_w2}
    def verify(self, pp, crs, Pi, com_x, com_y):
        # Initialize dictionaries and LHS
        p1 = {}; p2 = {}; LHS = 1
        # Set N to the length of com_x and the lengh of com_y
        m = len(com_x) #= len(com_y)
        # Compute an extended bilinear pairing on the received valus
        for vv1 in [0, 1]:
            for vv2 in [0, 1]:
                p1.update({i: com_x[i][vv1] if i < m else \
                        (crs['vv1'][vv1]**-1 if i == m else crs['ww1'][vv1]**-1 if i == m+1 \
                            else Pi['pi_v2'][vv1] if i == m+2 else Pi['pi_w2'][vv1]) for i in range(m+4)})
                p2.update({i: com_y[i][vv2] if i < n else \
                        (Pi['pi_v1'][vv2] if i == m else Pi['pi_w1'][vv2] if i == m+1 \
                            else crs['vv2'][vv2]**-1 if i == m+2 else crs['ww2'][vv2]**-1) for i in range(m+4)})
                # Compute the pairing of each element in p1 and p2, and multiply them all and keep them in LHS
                LHS = reduce(lambda x, y: x * y, [pair(p1[k], p2[k]) for k in range(m+4)])
                print(LHS)
        # Checrs if LHS is equal to the identity value in GT, i.e. pp['GT']**0, and return the result
        return LHS == pp['GT']**0
    # The batched verification algorithm reduces the number of pairings to N+4
    def Batched_verify(self, pp, crs, Pi, com_x, com_y):
        # Initialize dictionaries and LHS
        p1 = {}; p2 = {}; LHS = 1;
        # Set m to the length of com_x and n to the lengh of com_y
        m = len(com_x) #= len(com_y)
        P1 = {}; P2 = {}
        S = [group.random(), group.random()]
        R = [group.random(), group.random()]
        # Loop over all possible combinations of vv1 and vv2
        for vv1 in [0, 1]:
            p1.update({i: com_x[i][vv1] if i < m else \
                        (crs['vv1'][vv1]**-1 if i == m else crs['ww1'][vv1]**-1 if i == m+1 \
                        else Pi['pi_v2'][vv1] if i == m+2 else Pi['pi_w2'][vv1]) for i in range(m+4)})
            P1[vv1] = p1
            p2.update({i: com_y[i][vv1] if i < n else \
                        (Pi['pi_v1'][vv1] if i == m else Pi['pi_w1'][vv1] if i == m+1 \
                        else crs['vv2'][vv1]**-1 if i == m+2 else crs['ww2'][vv1]**-1) for i in range(m+4)})
            P2[vv1] = p2
            # Compute the pairing of each element in p1 and p2, and multiply them all and keep them in LHS
        P1 = [(P1[0][i]**S[0])*(P1[1][i]**S[1]) for i in range(len(P1[0]))]
        P2 = [(P2[0][i]**R[0])*(P2[1][i]**R[1]) for i in range(len(P2[0]))]
        LHS = reduce(lambda x, y: x * y, [pair(P1[k], P2[k]) for k in range(m+4)])
        # Checrs if LHS is equal to the identity value in GT, i.e. pp['GT']**0, and return the result
        return LHS == pp['GT'] ** 0

group=PairingGroup('BN254')
GS= GS(group)

def start_bench(group):
    group.InitBenchmark()
    group.StartBenchmark(["RealTime", "Pair"])

def end_bench(group):
    group.EndBenchmark()
    benchmarks = group.GetGeneralBenchmarks()
    real_time = benchmarks['RealTime'], benchmarks["Pair"]
    return real_time

g1,g2=group.random(G1),group.random(G2)
pp={'G1':g1,'G2':g2,'GT':pair(g1,g2)} 

def example():
    c_x = [10, 4]
    c_y = [2, group.init(ZR,-5)]

    # The values a,b only hide 2 and 5: "None" means "hide this values under commitment".
    # Non-hidden values must be the same as in X, Y.
    # So essentially it makes sure that \exist W1 W2 s.t.
    #   e(10[G],W1)e(2[G],(-1*)W2) = 1
    c_a = [10,None]
    c_b = [2,None]
    crs,td = GS.sampleParams(pp)
    r,s = GS.ParamGen(c_a,c_b)
    x = [pp['G1']**val for val in c_x]
    y = [pp['G2']**val for val in c_y]
    com_x, com_y = GS.commit(crs,x,y,r,s)
    pi = GS.prove(crs, x, y, r, s, com_x, com_y)
    start_bench(group)
    out = GS.Bverify(pp,crs,pi,com_x,com_y)
    verify_time, verify_pair = end_bench(group)
    print(verify_pair)
    print(out)
#example()

def main(n):
    result = [n]
    l_1 = [group.random(ZR) for _ in range(n-1)]
    l_2 = [group.random(ZR) for _ in range(n-1)]
    p = group.order()
    l_1.append(p-(np.sum([x * y for x, y in zip(l_1,l_2)])))
    l_2.append(1)
    print('IP(l_1,l_2)={}'.format(np.sum([x * y for x, y in zip(l_1, l_2)])))
    x = [pp['G1']**val for val in l_1]
    y = [pp['G2']**val for val in l_2]
    c_a = [None]*(n-1); c_a.append(x[n-1])
    c_b = [None]*(n-1); c_b.append(y[n-2])

    crs,td = GS.Transpatent_Setup(pp)
    #r,s = GS.ParamGen(c_a,c_b)
    commit_time=0
    start_bench(group)
    com_x, com_y, r, s = GS.commit(crs,x,y,c_a,c_b)
    commit_time, commit_pair=end_bench(group)
    result.append(commit_time)
    prove_time=0
    start_bench(group)
    pi = GS.prove(crs, x, r, s, com_y)
    prove_time, prove_pair = end_bench(group)
    result.append(prove_time)
    verify_time=0
    start_bench(group)
    out = GS.verify(pp,crs,pi,com_x,com_y)
    verify_time, verify_pair = end_bench(group)
    result.append(verify_time); result.append(verify_pair)
    verify_time=0
    start_bench(group)
    out = GS.Batched_verify(pp,crs,pi,com_x,com_y)
    verify_time, verify_pair = end_bench(group)
    result.append(verify_time); result.append(verify_pair)
    print("The verification returned {} in n={}".format(out,n))
    return result

book = Workbook()
data = book.active
title = ["n", "commit_time", "prove_time", "verify_time", "#Pairings in Verify", "Batched-verify_time", "#Pairings in BVerify"]
data.append(title)
for n in range(10,101,5):
    data.append(main(n))
book.save("result.xlsx")
