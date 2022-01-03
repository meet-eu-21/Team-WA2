#!/usr/bin/python
import sys
import pandas as pd
from math import sqrt

def overlap(a,b):
  c=[0,0,0] #overlaping fragment
  if (a[1]<b[1] and a[2]>b[2]):
    c[1]=b[1]
    c[2]=b[2]
  elif (a[1]>=b[1] and a[1]<=b[2]):
    c[1]=a[1]
    c[2]=min(a[2],b[2])
  elif (a[2]>=b[1] and a[2]<=b[2]):
    c[2]=a[2]
    c[1]=max(a[1],b[1])
  else:
    return 0

  return c[2]-c[1]

def compute_MoC(P,Q): # returns unnormalized sum
  # and normalizing factor (MOC*sqrt(P.shape[0]*Q.shape[0])-1) and sqrt(P.shape[0]*Q.shape[0])-1))
  sum=0
  for i, p in P.iterrows(): #index, row
    for j, q in Q.iterrows():
      if p[0] == q[0]: #same chromosomes
        F = overlap(p,q)
        if F:
          sum+= F*F / ((p[2]-p[1])*(q[2]-q[1]))
  return (sum,(sqrt(P.shape[0]*Q.shape[0])-1))  #without -1 returns MoC=1.0 for equals TAD sets

file1= sys.argv[1]
file2= sys.argv[2]
P= pd.read_csv(file1,sep=',')
P= P[P["tag"]=="domain"]
Q= pd.read_csv(file2,sep=',')
Q= Q[Q["tag"]=="domain"]

big_sum=0
MoC=[]
chroms = list(set(P["chr"]))

def chroms_val(e):
  e=e.replace("chr","")
  if(e == 'X'):
    return 100
  elif (e == 'Y'):
    return 101
  return int(e)

chroms.sort(key=chroms_val)

for chr in chroms:
  sum,norm= compute_MoC(P[P["chr"]==chr],Q[Q["chr"]==chr])
  big_sum+=sum
  MoC.append(sum/norm)
  print(sum/norm)

chroms.append("Whole genome")
res=pd.DataFrame(chroms,columns=["Chr"])
MoC.append(big_sum/ (sqrt(P.shape[0]*Q.shape[0])-1))
res["MoC"]=MoC

res.to_csv(str(file1).replace(".csv",'')+".x."+ str(file2).replace(".csv",'')+".MoC.scores.csv",index=False)
