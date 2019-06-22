import numpy as np
import scipy
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from matplotlib import pyplot as plt

folder="../../data/input/taxo_metabar/adj/none/"
outfolder="../../data/input/taxo_metabar/adj/none/out/"
taxo_scales=["family","genus","species"]

n=10

for taxo in taxo_scales:
    ### Read affinity matrice
    path=folder+taxo+"-level_affinity.csv"
    aff=pd.read_csv(path,sep=";",decimal=".",index_col=0)

    ### Agglomerative clustering
    Z = linkage(aff, 'complete')
    fc=np.transpose([fcluster(Z, height, criterion='distance') for height in range(30)])

    out=pd.DataFrame(data=fc)

    uniques=[len(np.unique(fc[:,i])) for i in range(0,fc.shape[1])]
    out.loc[out.shape[0]]=uniques

    col=list(aff.columns)
    col.append("counts")
    out["esp"]=pd.Series(col)
    out=out.set_index("esp")
    out.to_csv(outfolder+taxo+"_hier.csv",sep=";",decimal=".",index=True)

    dn= dendrogram(Z)
    plt.show()
    #pd.DataFrame(dn.get("icoord")).to_csv(outfolder+taxo+"_dn.csv",sep=";")
    
    #plt.savefig(outfolder+taxo+".png")
    
##    data=[pd.Series(aff.columns,name="SPEC"),pd.Series(y,name="LABEL")]
##    comm_assign=pd.concat(data,axis=1)
##
##    comm_assign.to_csv(folder+taxo+"_communities.csv",sep=";",decimal=".",index=False)
    




