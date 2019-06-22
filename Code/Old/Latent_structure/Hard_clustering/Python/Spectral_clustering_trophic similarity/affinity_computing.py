import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import binarize

folder="../../data/input/taxo_metabar/adj/none/"
taxo_scales=["species","genus","family"]

for taxo in taxo_scales:
    path=folder+taxo+".csv"
    adj=pd.read_csv(path,sep=";",decimal=".",index_col=0)
    adj["esp"]=adj.index
    adj=adj.set_index("esp")

    #adj=pd.DataFrame(data=binarize(adjw,threshold=0.0,copy=False),columns=list(adjw.columns))

    sumcols=[adj.loc[:,x].sum() for x in adj.columns]
    sumrows=[adj.loc[x,:].sum() for x in adj.columns]

    ### Sum_row=0 => all rows for x are set to zero => no matter transfer from x to any other specie => top predator
    toppreds=[adj.columns[i] for i in range(0,len(sumrows)) if (sumrows[i]==0 and sumcols[i]>0)]
    nopreds=[adj.columns[i] for i in range(0,len(sumrows)) if sumrows[i]==0]

    ### Inverse logic => no matter transfer from any other specie to x => resource
    res=[adj.columns[i] for i in range(0,len(sumcols)) if (sumcols[i]==0 and sumrows[i]>0)]
    nores=[adj.columns[i] for i in range(0,len(sumcols)) if sumcols[i]==0]

    ### Integrity test, isolated nodes removal
    isol=set(nopreds).intersection(set(nores))
    filt=adj.drop(isol,axis=1)
    filt=filt.drop(isol,axis=0)

    ##### Pairwise cosine similarity
    ## Prey/Resource similarity
    Ri=cosine_similarity(filt,filt)
    dfi=pd.DataFrame(data=Ri,columns=filt.columns)

    #dfi.to_csv("../../data/in_simil.csv",sep=";",decimal=".")

    ## Predator/Consumer similarity
    Ro=cosine_similarity(filt.transpose(),filt.transpose())
    dfo=pd.DataFrame(data=Ro,columns=filt.columns)

    #dfo.to_csv("../../data/out_simil.csv",sep=";",decimal=".")

    df=(dfi+dfo)/2  ### Aggregating in/out trophic similarities => average trophic similarity 
    #df.to_csv("../data/simil.csv",sep=";",decimal=".")

    ### Forcing intra-spec similarity to 1
    for i in range(0,len(df)):
        df.iloc[i,i]=1

    df["esp"]=pd.Series(df.columns)
    df=df.set_index("esp")
    df.to_csv(folder+taxo+"-level_affinity.csv",sep=";",decimal=".",index=True)
