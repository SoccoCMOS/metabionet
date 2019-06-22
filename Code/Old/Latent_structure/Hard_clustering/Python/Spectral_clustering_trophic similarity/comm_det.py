import numpy as np
import scipy
import pandas as pd
from sklearn.metrics import jaccard_similarity_score
from sklearn.cluster import SpectralClustering

### Read affinity matrix

aff_mat=pd.read_csv("../../data/simil_forced.csv",sep=";",decimal=".",index_col=0)

n=15 ### Number of clusters

spectral=SpectralClustering(n, eigen_solver='arpack',
        affinity="precomputed",random_state=0)
##
spectral.fit(aff_mat)
labels=spectral.fit_predict(aff_mat)

data=[pd.Series(aff_mat.columns,name="SPEC"),pd.Series(labels,name="LABEL")]
comm_assign=pd.concat(data,axis=1)

comm_assign.to_csv("../../data/"+str(n)+"_communities.csv",sep=";",decimal=".",index=False)
