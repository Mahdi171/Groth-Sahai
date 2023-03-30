from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair
from openpyxl import Workbook
import numpy as np
from functools import reduce
class GS():
    def __init__(self, groupObj):
        global util, group      
        group = groupObj

    def sampleParams(self,pp):
        rho, zeta, sigma, omega = [group.random() for _ in range(4)]
        vv1 = [pp['G1']**zeta, pp['G1']]
        vv2 = [pp['G2']**omega, pp['G2']]
        ww1 = [pp['G1']**(rho*zeta), pp['G1']**rho]
        ww2 = [pp['G2']**(sigma*omega), pp['G2']*sigma]
        uu1 = [0, ww1[1]*pp['G1']]
        uu2 = [0, ww2[1]*pp['G2']]
        Zeta = [-zeta**(-1), 1]
        Omega = [-omega**(-1), 1]
        ck = {'uu1':uu1, 'vv1':vv1, 'ww1':ww1, 'uu2':uu2, 'vv2':vv2, 'ww2':ww2}
        xk = {'ck':ck, 'Zeta':Zeta, 'Omega': Omega}
        return ck, xk

    def ParamGen(self,c_a,c_b):
        m = len(c_a); n = len(c_b)
        r = list(map(lambda i: [0,0] if c_a[i] != None \
                      else [group.random(),group.random()], range(m)))
        s = list(map(lambda i: [0,0] if c_b[i] != None \
                     else [group.random(),group.random()], range(n)))
        return r,s
    
    def commit(self, ck ,X, Y, r, s):
        cc_x={}; cc_y={}
        cc_x = [[(ck['vv1'][0]**r[i][0])*(ck['ww1'][0]**s[i][0]),\
                  X[i]*(ck['vv1'][1]**r[i][0])*(ck['ww1'][1]**s[i][0])] for i in range(len(X))]
        cc_y = [[(ck['vv2'][0]**r[i][1])*(ck['ww2'][0]**s[i][1]),\
                 Y[i]*(ck['vv2'][1]**r[i][1])*(ck['ww2'][1]**s[i][1])] for i in range(len(Y))]
        return cc_x, cc_y


    def prove(self, ck, X, Y, r, s, cc_x, cc_y):
        # Define prod function that multiplies elements of any given list
        def prod(list):
            result = 1
            for x in list:
                result *= x
            return result
        # Generate random values alpha, beta, gamma, and delta
        alpha, beta, gamma, delta = [group.random() for _ in range(4)]
        # Compute pi_v1, pi_v2, pi_w1, and pi_w2
        pi_v1 = [prod([cc_y[i][0]**r[i][0] for i in range(len(X))]) * ck['vv2'][0]**alpha * ck['ww2'][0]**beta,\
                prod([cc_y[i][1]**r[i][0] for i in range(len(X))]) * ck['vv2'][1]**alpha * ck['ww2'][1]**beta]
        pi_w1 = [prod([cc_y[i][0]**s[i][0] for i in range(len(X))]) * ck['vv2'][0]**gamma * ck['ww2'][0]**delta,\
                prod([cc_y[i][1]**s[i][0] for i in range(len(X))]) * ck['vv2'][1]**gamma * ck['ww2'][1]**delta]
        pi_v2 = [ck['vv1'][0]**-alpha * ck['ww1'][0]**-gamma,\
                prod([X[i]**r[i][1] for i in range(len(X))]) * ck['vv1'][1]**-alpha * ck['ww1'][1]**-gamma]
        pi_w2 = [ck['vv1'][0]**-beta * ck['ww1'][0]**-delta,\
                prod([X[i]**s[i][1] for i in range(len(X))]) * ck['vv1'][1]**-beta * ck['ww1'][1]**-delta]
        # Return the proof values as a dictionary
        return {'pi_v1': pi_v1, 'pi_w1': pi_w1, 'pi_v2': pi_v2, 'pi_w2': pi_w2}
    
    def verify(self, pp, ck, Pi, cc_x, cc_y):
        # Initialize dictionaries and LHS        
        p1 = {}; p2 = {}; LHS = 1
        # Set m to the length of cc_x and n to the lengh of cc_y
        m = len(cc_x); n = len(cc_y)
        # Loop over all possible combinations of vv1 and vv2
        for vv1 in [0, 1]:
            for vv2 in [0, 1]:
                p1.update({i: cc_x[i][vv1] if i < m else \
                           (ck['vv1'][vv1]**-1 if i == m else ck['ww1'][vv1]**-1 if i == m+1 \
                            else Pi['pi_v2'][vv1] if i == m+2 else Pi['pi_w2'][vv1]) for i in range(m+4)})
                p2.update({i: cc_y[i][vv2] if i < n else \
                           (Pi['pi_v1'][vv2] if i == m else Pi['pi_w1'][vv2] if i == m+1 \
                            else ck['vv2'][vv2]**-1 if i == m+2 else ck['ww2'][vv2]**-1) for i in range(m+4)})
                # Compute the pairing of each element in p1 and p2, and multiply them all and keep them in LHS
                LHS = reduce(lambda x, y: x * y, [pair(p1[k], p2[k]) for k in range(m+4)])
                print(LHS)
        # Check if LHS is equal to the identity value in GT, i.e. pp['GT']**0, and return the result
        return LHS == pp['GT']**0

    # The batched verification algorithm reduces the number of pairings to min(m,n)+4
    def Bverify(self, pp, ck, Pi, cc_x, cc_y):
        # Initialize dictionaries and LHS        
        p1 = {}; p2 = {}; LHS = 1; 
        # Set m to the length of cc_x and n to the lengh of cc_y
        m = len(cc_x); n = len(cc_y)
        P1={}; P2={}
        #P1= [0 for _ in range(m+4)]; P2= [0 for _ in range(n+4)]
        S = [group.random(), group.random()]
        R = [group.random(), group.random()]
        # Loop over all possible combinations of vv1 and vv2
        for vv1 in [0, 1]:
            p1.update({i: cc_x[i][vv1] if i < m else \
                        (ck['vv1'][vv1]**-1 if i == m else ck['ww1'][vv1]**-1 if i == m+1 \
                        else Pi['pi_v2'][vv1] if i == m+2 else Pi['pi_w2'][vv1]) for i in range(m+4)})
            P1[vv1] = p1
            p2.update({i: cc_y[i][vv1] if i < n else \
                        (Pi['pi_v1'][vv1] if i == m else Pi['pi_w1'][vv1] if i == m+1 \
                        else ck['vv2'][vv1]**-1 if i == m+2 else ck['ww2'][vv1]**-1) for i in range(m+4)})
            P2[vv1] = p2
            # Compute the pairing of each element in p1 and p2, and multiply them all and keep them in LHS
        P1 = [(P1[0][i]**S[0])*(P1[1][i]**S[1]) for i in range(len(P1[0]))]
        P2 = [(P2[0][i]**R[0])*(P2[1][i]**R[1]) for i in range(len(P2[0]))]
        LHS = reduce(lambda x, y: x * y, [pair(P1[k], P2[k]) for k in range(m+4)])
        print(LHS)
        # Check if LHS is equal to the identity value in GT, i.e. pp['GT']**0, and return the result
        return LHS == pp['GT']**0

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
    ck,td = GS.sampleParams(pp)

    r,s = GS.ParamGen(c_a,c_b)
    commit_time=0
    start_bench(group)
    cc_x, cc_y = GS.commit(ck,x,y,r,s)
    commit_time, commit_pair=end_bench(group)
    result.append(commit_time)
    prove_time=0
    start_bench(group)
    pi = GS.prove(ck, x, y, r, s, cc_x, cc_y)
    prove_time, prove_pair = end_bench(group)
    result.append(prove_time)
    verify_time=0
    start_bench(group)
    out = GS.Bverify(pp,ck,pi,cc_x,cc_y)
    verify_time, verify_pair = end_bench(group)
    result.append(verify_time); result.append(verify_pair)
    print("The verification returned {} in n={}".format(out,n))
    return result

book = Workbook()
data = book.active
title = ["n", "commit_time", "prove_time", "verify_time", "#Pairings in Verify"]
data.append(title)
for n in range(20,31):
    data.append(main(n))
book.save("result.xlsx")
