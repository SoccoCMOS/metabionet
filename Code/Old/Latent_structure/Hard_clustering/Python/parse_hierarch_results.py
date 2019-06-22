import numpy as np
import pandas as pd

filename='groupsW-10.txt'

with open(filename) as f:
    clusters = f.readlines()

lt=list()
lc=list()
ltg=list()
ltgd=list()

cpt=-1
for c in clusters:
    cpt+=1
    taxac=c[:-1].split("\t")
    for t in taxac:
        lt.append(t)
        lc.append(cpt)

ltf=pd.Series(lt)
lcf=pd.Series(lc)
df=pd.concat([ltf,lcf],axis=1)
df.columns=["retained_tax","class"]

df.to_csv("groupsW-10.txt")
        
    
    


