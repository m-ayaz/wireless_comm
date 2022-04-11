# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 12:19:03 2022

@author: M
"""

import pulp
from pulp import lpSum
# import numpy as np
#from pulp.const import LpBinary

#LpBinary=pulp.const.LpBinary

def distance(x1,y1,x2,y2):
    return ((x2-x1)**2+(y2-y1)**2)**0.5

if __name__=='__main__':
    
    """Table 1: data set"""
    
    C={
       "CO1": (0,0)
       }
    
    M={}
    
    O={}
    
    B={}
    
    """Table 2: network parameters"""
    
    n_c=len(C) # Verified
    
    n_m=len(M) # Verified
    
    n_o=len(O) # Verified
    
    n_b=len(B) # Verified
    
    df={}
    
    dd={}
    
    t={}
    for k in O:
        for l in B:
            if 1:
                t[k,l]=1
            else:
                t[k,l]=0
                
    
    d_max=0
    
    n_s=1
    
    n_p=16 # Verified
    
    n_r=3 # Verified
    
    p_c=0.99 # verified
    
    c_h=25 # MHz, verified
    
    c_r=1000 # MHz, verified
    
    bigM=100000
    
    """Table 3: cost components"""
    
    alpha_ff=41.6 # relative, verified
    
    alpha_fd=10.4 # relative, verified
    
    alpha_ft=93.7 # relative, verified
    
    alpha_fr=1 # relative, verified
    
    alpha_oc=118.75 # relative, verified
    
    alpha_lc=13.3 # relative, verified
    
    alpha_ot=3.6 # relative, verified
    
    alpha_fc=3.1 # relative, verified
    
    alpha_si=3.3 # relative, verified
    
    alpha_on=4.7 # relative, verified
    
    alpha_rr=20.8 # relative, verified
    
    alpha_bb=729.1 # relative, verified
    
    '''Temp'''
    alpha_fi=alpha_fb=alpha_bi=alpha_ss=alpha_ri=0
    
    '''Variables'''
    f={} # if an OLT at the ith CO is connected to a splitter at the jth FAP
    for i in C:
        for j in M:
            f[i,j]=pulp.LpVariable(f"f_{i}_{j}",cat="Binary")
    
    f_tilde={} # number of fiber connections between the ith CO and the jth FAP
    for i in C:
        for j in M:
            f_tilde[i,j]=pulp.LpVariable(f"f_tilde_{i}_{j}",cat="Integer",lowBound=0)
    
    d={} # if a splitter at the jth FAP is connected to an ONU at the kth node
    for j in M:
        for k in O:
            d[j,k]=pulp.LpVariable(f"d_{j}_{k}",cat="Binary")
    
    d_tilde={} # number of fiber connections between the jth FAP and the kth node
    for j in M:
        for k in O:
            d_tilde[j,k]=pulp.LpVariable(f"d_tilde_{j}_{k}",cat="Integer",lowBound=0)
    
    s={} # if jth FAP is selected to place splitter/s.
    for j in M:
        s[j]=pulp.LpVariable(f"s_{j}",cat="Binary")
    
    s_tilde={} # number of splitters placed at the jth FAP
    for j in M:
        s_tilde[j]=pulp.LpVariable(f"s_tilde_{j}",cat="Integer",lowBound=0)
    
    r={} # if kth node is selected to place RRH/s
    for k in O:
        r[k]=pulp.LpVariable(f"r_{k}",cat="Binary",upBound=n_r)
    
    r_tilde={} # number of RRHs placed at the kth node
    for k in O:
        r_tilde[k]=pulp.LpVariable(f"r_tilde_{k}",cat="Integer",lowBound=0)
    
    x={} # if the ith CO is selected to place a BBU
    for i in C:
        x[i]=pulp.LpVariable(f"x_{i}",cat="Binary")
    
    y_tilde={} # number of line cards in the ith CO.
    for i in C:
        y_tilde[i]=pulp.LpVariable(f"r_{i}",cat="Binary",upBound=n_r)
    
    h={} # if the lth household is covered
    for l in B:
        h[l]=pulp.LpVariable(f"h_{l}",cat="Binary")
    
    z={} # if lth household is served by an RRH placed at the kth node
    for l in B:
        for k in O:
            z[l,k]=pulp.LpVariable(f"z_{l}_{k}",cat="Binary")
    
    """Set up a problem"""
    prob=pulp.LpProblem(name="Cost_Optimization",sense=pulp.LpMinimize)
    
    """Constraints"""
    
    """Coverage Requirement"""
    """eq2""" """Typo Warning"""
    terms_eq2=[]
    for l in B:
        terms_eq2.append(h[l])
    prob+=pulp.lpSum(terms_eq2)>=p_c*n_b/100
    """eq3"""
    for l in B:
        terms_eq3=[]
        for k in O:
            terms_eq3.append(t[k,l]*r[k])
        prob+=h[l]<=lpSum(terms_eq3)
    """eq4"""
    for l in B:
        terms_eq4=[]
        for k in O:
            terms_eq4.append(t[k,l]*r[k]/bigM)
        prob+=h[l]>=lpSum(terms_eq4)

    """Capacity Requirement"""
    """eq5"""
    for k in O:
        terms_eq5=[]
        for l in B:
            terms_eq5.append(z[k,l])
        prob+=r_tilde[k]*c_r>=c_h*lpSum(terms_eq5)
    """eq6"""
    for l in B:
        terms_eq6=[]
        for k in O:
            terms_eq6.append(z[k,l])
        prob+=lpSum(terms_eq6)>=h[l]
    """eq7"""
    for k in O:
        for l in B:
            prob+=z[k,l]<=t[k,l]
    
    """Point to Multipoint Deployment"""
    """eq8"""
    for j in M:
        terms_eq8=[]
        for i in C:
            terms_eq8.append(f_tilde[i,j])
        prob+=lpSum(terms_eq8)==s_tilde[j]
    """eq9"""
    for k in O:
        terms_eq9=[]
        for j in M:
            terms_eq9.append(d_tilde[j,k])
        prob+=lpSum(terms_eq9)==r_tilde[k]
    """eq10"""
    for k in O:
        terms_eq10=[]
        for j in M:
            terms_eq10.append(d[j,k])
        prob+=lpSum(terms_eq10)==r[k]
    """eq11"""
    for j in M:
        for k in O:
            prob+=d[j,k]<=s[j]
    """eq12"""
    for j in M:
        terms_eq12=[]
        for k in O:
            terms_eq12.append(d_tilde[j,k])
        prob+=lpSum(terms_eq12)<=n_s*s_tilde[j]
    
    """Network Span"""
    """eq13"""
    for i in C:
        for j in M:
            for k in O:
                prob+=f[i,j]*df[i,j]+d[i,k]*dd[j,k]<=d_max
    
    """Line Card Selection at the CO"""
    """eq14""" """Typo Warning"""
    for i in C:
        terms_eq14=[]
        for j in M:
            # print("sadasd/asdasd")
            terms_eq14.append(f_tilde[i,j]/n_p)
        prob+=y_tilde[i]<=lpSum(terms_eq14)+1
    """eq15""" """Typo Warning"""
    for i in C:
        terms_eq15=[]
        for j in M:
            terms_eq15.append(f_tilde[i,j]/n_p)
        prob+=y_tilde[i]>=lpSum(terms_eq15)
    
    """Placement of BBU"""
    """eqs 16 and 17""" """Typo Warning"""
    for i in C:
        prob+=y_tilde[i]*1/bigM<=x[i]
        prob+=y_tilde[i]>=x[i]
    
    """Nonlinear Relationships Between Variables"""
    """eqs 18 and 19"""
    for i in C:
        for j in M:
            prob+=f_tilde[i,j]>=f[i,j]
            prob+=f_tilde[i,j]/bigM<=f[i,j]
    """eqs 20 and 21"""
    for j in M:
        for k in O:
            prob+=d_tilde[j,k]>=d[j,k]
            prob+=d_tilde[j,k]/bigM<=d[j,k]
    """eqs 22 and 23"""
    for j in M:
        prob+=s_tilde[j]>=s[j]
        prob+=s_tilde[j]/bigM<=s[j]
    """eqs 24 and 25"""
    for k in O:
        prob+=r_tilde[j]>=r[j]
        prob+=r_tilde[j]/bigM<=r[j]
    
    """Bounds of Integer Decision Variables"""
    """These constraints are previously fulfilled in variable declarations."""
    
    """Objective function"""
    terms_obj=[]
    for i in C:
        for j in M:
            terms_obj.append(alpha_fr*f_tilde[i,j])
            terms_obj.append(alpha_ff*f_tilde[i,j]*df[i,j])
    for j in M:
        for k in O:
            terms_obj.append((alpha_ft+alpha_fi)*d[j,k]*dd[j,k])
            terms_obj.append(alpha_fb*d_tilde[j,k]*dd[j,k])
    for i in C:
        terms_obj.append(alpha_lc*y_tilde[i])
        terms_obj.append((alpha_bb+alpha_bi+alpha_oc)*x[i])
    for j in M:
        terms_obj.append((alpha_ot+alpha_fc)*s_tilde[j])
        terms_obj.append(alpha_si*s[j])
        terms_obj.append(alpha_ss*s_tilde[j])
    for k in O:
        terms_obj.append((alpha_on+alpha_rr)*r_tilde[k])
        terms_obj.append(alpha_ri*r[k])
    
    prob+=lpSum(terms_obj)
    
    """Solve the problem"""
    
    
    
    prob.solve()
        
    # print(prob)
    
    # prob.value()
    
    
    """Print results"""
    print("\f")
    print("Optimal value =",prob.objective.value())