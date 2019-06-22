import numpy as np
import pandas as pd
from graph_tool.all import *

adj=pd.read_csv("../../data/adjacence_esp_bin.csv",sep=";",decimal=".",index_col=0)

