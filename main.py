from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair
from openpyxl import Workbook
import numpy as np
from functools import reduce
import time
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
        com_x = []
        com_y = []
        n = len(X); m = len(Y)
        r = [([group.init(ZR,0),group.init(ZR,0)] if C_x[i] != "None" else [group.random(),group.random()]) for i in range(len(C_x))]
        s = [([group.init(ZR,0),group.init(ZR,0)] if C_y[i] != "None" else [group.random(),group.random()]) for i in range(len(C_y))]
        for i in range(n):
            com_x.append([(crs['vv1'][0]**r[i][0])*(crs['ww1'][0]**r[i][1]),\
                           X[i]*(crs['vv1'][1]**r[i][0])*(crs['ww1'][1]**r[i][1])])
        for i in range(m):
            com_y.append([(crs['vv2'][0]**s[i][0])*(crs['ww2'][0]**s[i][1]),\
                           Y[i]*(crs['vv2'][1]**s[i][0])*(crs['ww2'][1]**s[i][1])])
        return com_x, com_y, r, s
    
    def prove(self, crs, X, Y, r, s, com_y, GammaT):
        proof = {}
        n = len(X); m = len(Y)
        for ii in range(len(GammaT)):
            gammaT = GammaT[ii]; Com_y = {}; Xp={}
            alpha, beta, gamma, delta = [group.random() for _ in range(4)]
            for j in range(n):
                aux1 = 1; aux2 = 1; aux3 = 1
                for k in range(m):
                    aux1 *= com_y[k][0]**gammaT[k][j]
                    aux2 *= com_y[k][1]**gammaT[k][j]
                    aux3 *= X[k]**gammaT[j][k]
                Com_y[j] = [aux1, aux2]
                Xp[j] = aux3
            
            pi_v1 = [reduce(lambda x, y: x * y, [Com_y[i][0]**r[i][0] for i in range(m)]) * crs['vv2'][0]**alpha * crs['ww2'][0]**beta,\
                    reduce(lambda x, y: x * y, [Com_y[i][1]**r[i][0] for i in range(m)]) * crs['vv2'][1]**alpha * crs['ww2'][1]**beta]
            pi_w1 = [reduce(lambda x, y: x * y, [Com_y[i][0]**s[i][0] for i in range(m)]) * crs['vv2'][0]**gamma * crs['ww2'][0]**delta,\
                    reduce(lambda x, y: x * y, [Com_y[i][1]**s[i][0] for i in range(m)]) * crs['vv2'][1]**gamma * crs['ww2'][1]**delta]
            pi_v2 = [crs['vv1'][0]**-alpha * crs['ww1'][0]**(-gamma),\
                    reduce(lambda x, y: x * y, [Xp[i]**r[i][1] for i in range(n)]) * crs['vv1'][1]**-alpha * crs['ww1'][1]**-gamma]
            pi_w2 = [crs['vv1'][0]**-beta * crs['ww1'][0]**(-delta),\
                    reduce(lambda x, y: x * y, [Xp[i]**s[i][1] for i in range(n)]) * crs['vv1'][1]**-beta * crs['ww1'][1]**-delta]
            proof[ii]={'pi_v1': pi_v1, 'pi_w1': pi_w1, 'pi_v2': pi_v2, 'pi_w2': pi_w2}
        return proof
    
    def verify(self, pp, crs, proof, com_x, com_y, GammaT):
        # Initialize dictionaries and LHS
        p1 = {}; p2 = {}; LHS = 1
        # Set N to the length of com_x and the lengh of com_y
        n = len(com_x); m = len(com_y)
        # Compute an extended bilinear pairing on the received valus
        for ii in range(len(GammaT)):
            gammaT = GammaT[ii]
            Pi = proof[ii]
            for vv1 in [0, 1]:
                for vv2 in [0, 1]:
                    for i in range(n):
                        p1[i] = com_x[i][vv1]
                        p2[i] = 1
                        for j in range(m):
                            p2[i] *= com_y[j][vv2]**gammaT[j][i]
                    
                    p1.update({i:crs['vv1'][vv1]**-1 if i == m else crs['ww1'][vv1]**-1 if i == m+1 \
                                else Pi['pi_v2'][vv1] if i == m+2 else Pi['pi_w2'][vv1] for i in range(m,m+4)})
                    p2.update({i:Pi['pi_v1'][vv2] if i == m else Pi['pi_w1'][vv2] if i == m+1 \
                                else crs['vv2'][vv2]**-1 if i == m+2 else crs['ww2'][vv2]**-1 for i in range(m,m+4)})
                    # Compute the pairing of each element in p1 and p2, and multiply them all and keep them in LHS
                    LHS = reduce(lambda x, y: x * y, [pair(p1[k], p2[k]) for k in range(m+4)])
                    if LHS != pp['GT']**0:
                        return False
            # Checrs if LHS is equal to the identity value in GT, i.e. pp['GT']**0, and return the result
        return True
    # The batched verification algorithm reduces the number of pairings to N+4
    def Batched_verify(self, pp, crs, proof, com_x, com_y, GammaT):
        # Initialize dictionaries and LHS
        p1 = {}; p2 = {}; LHS = 1;
        # Set m to the length of com_x and n to the lengh of com_y
        m = len(com_x) #= len(com_y)
        P1 = {}; P2 = {}
        m = len(com_x); n= len(com_y)
        S = [group.random(), group.random()]
        R = [group.random(), group.random()]
        # Loop over all possible combinations of vv1 and vv2
        for ii in range(len(GammaT)):
            gammaT = GammaT[ii]
            Pi = proof[ii]
            for vv1 in [0, 1]:
                for vv2 in [0, 1]:
                    for i in range(m):
                        p1[i] = com_x[i][vv1]
                        p2[i] = 1
                        for j in range(n):
                            p2[i] *= com_y[j][vv2]**gammaT[j][i]
        for vv1 in [0, 1]:
            p1.update({i:(crs['vv1'][vv1]**-1 if i == m else crs['ww1'][vv1]**-1 if i == m+1 \
                        else Pi['pi_v2'][vv1] if i == m+2 else Pi['pi_w2'][vv1]) for i in range(m,m+4)})
            P1[vv1] = p1
            p2.update({i: (Pi['pi_v1'][vv1] if i == m else Pi['pi_w1'][vv1] if i == m+1 \
                        else crs['vv2'][vv1]**-1 if i == m+2 else crs['ww2'][vv1]**-1) for i in range(m,m+4)})
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


###########################################################
# A simple Test & main function for arbitrary number of PPE
###########################################################



def example1():
    c_x = [10, 4]
    c_y = [2, 5]
    # This example is the same as one of examples in this repo: volhovm/groth-sahai-python

    # The values a,b only hide 2 and 5: "None" means "hide this values under commitment".
    # Non-hidden values must be the same as in X, Y.
    # So essentially it makes sure that \exist W1 W2 s.t.
    #   e([10]_1,W1)e([2]_2,(-1*)W2) = 1
    #([com_y[k][0]**gammaT[k][i] for k in range(len(com_y))])
    c_a = [10,None]
    c_b = [2,None]
    gammaT = [[1,0],[0,-1]]
    crs,td = GS.Transpatent_Setup(pp)
    x = [pp['G1']**val for val in c_x]
    y = [pp['G2']**val for val in c_y]
    com_x, com_y, r, s = GS.commit(crs,x,y,c_a,c_b)
    pi = GS.prove(crs, x, y, r, s, com_y,[gammaT])
    start_bench(group)
    t0=time.time()
    out = GS.verify(pp,crs,pi,com_x,com_y,[gammaT])
    print(time.time()-t0)
    verify_time, verify_pair = end_bench(group)
    print(verify_pair)
    #print(verify_pair)
    t0=time.time()
    start_bench(group)
    out = GS.Batched_verify(pp,crs,pi,com_x,com_y,[gammaT])
    verify_time, verify_pair = end_bench(group)
    print(verify_pair)
    print(time.time()-t0)
    print(out)

example1()

def example2():
    msg = group.init(ZR,1) # or 0
    r = 14352345
    sk = 36534152
    ct1 = r
    ct2 = sk * r + msg

    w1 = r
    w2 = msg
    w3 = msg

    c_x = [w2, ct1, ct2, sk, group.init(ZR,1)]
    c_a = ["None", ct1, ct2, sk, group.init(ZR,1)]

    c_y = [w1, w3, group.init(ZR,1)]
    c_b = ["None", "None", group.init(ZR,1)]
    print(c_b)
    x = [pp['G1']**val for val in c_x]
    y = [pp['G2']**val for val in c_y]
    gammaT_E1 = [[0,0,0,0,-1],[0,0,0,0,0],[0,1,0,0,0]]
    gammaT_E2 = [[0,0,0,-1,0],[0,0,0,0,0],[-1,0,1,0,0]]
    gammaT_E3 = [[0,0,0,0,0],[0,0,0,0,-1],[1,0,0,0,0]]
    gammaT_E4 = [[0,0,0,0,0],[1,0,0,0,0],[-1,0,0,0,0]]
    gammaT = [gammaT_E1,gammaT_E2,gammaT_E3,gammaT_E4]
    
    crs,td = GS.Transpatent_Setup(pp)
    com_x, com_y, r, s = GS.commit(crs,x,y,c_a,c_b)
    pi = GS.prove(crs, x, y, r, s, com_y,[gammaT])
    start_bench(group)
    out = GS.verify(pp,crs,pi,com_x,com_y,[gammaT])
    verify_time, verify_pair = end_bench(group)
    print(out)

#example2()

def example3():
    c_x = [10, 5, 3, 9]
    c_y = [2, 4, 6, 2]
    # This example is the same as one of examples in this repo: volhovm/groth-sahai-python

    # The values a,b only hide 2 and 5: "None" means "hide this values under commitment".
    # Non-hidden values must be the same as in X, Y.
    # So essentially it makes sure that \exist W1 W2 s.t.
    #   e([10]_1,W1)e([2]_2,(-1*)W2) = 1
    #([com_y[k][0]**gammaT[k][i] for k in range(len(com_y))])
    c_a = [10,None,3,9]
    c_b = [2,4,6,None]
    gammaT = [[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]]
    crs,td = GS.Transpatent_Setup(pp)
    x = [pp['G1']**val for val in c_x]
    y = [pp['G2']**val for val in c_y]
    com_x, com_y, r, s = GS.commit(crs,x,y,c_a,c_b)
    pi = GS.prove(crs, x, y, r, s, com_y,[gammaT])
    start_bench(group)
    out = GS.verify(pp,crs,pi,com_x,com_y,[gammaT])
    verify_time, verify_pair = end_bench(group)
    #print(verify_pair)
    print(out)

example3()

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
    pi = GS.prove(crs, x, y, r, s, com_y)
    prove_time, prove_pair = end_bench(group)
    result.append(prove_time)
    verify_time=0
    start_bench(group)
    out = GS.verify(pp,crs,pi,com_x,com_y)
    verify_time, verify_pair = end_bench(group)
    result.append(verify_time); result.append(verify_pair)
    verify_time=0
    start_bench(group)
    #out = GS.Batched_verify(pp,crs,pi,com_x,com_y)
    #verify_time, verify_pair = end_bench(group)
    #result.append(verify_time); result.append(verify_pair)
    print("The verification returned {} in n={}".format(out,n))
    return result
'''
book = Workbook()
data = book.active
title = ["n", "commit_time", "prove_time", "verify_time", "#Pairings in Verify", "Batched-verify_time", "#Pairings in BVerify"]
data.append(title)
for n in range(10,101,5):
    data.append(main(n))
book.save("result.xlsx")
'''