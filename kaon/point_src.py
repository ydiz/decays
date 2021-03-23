import os
import numpy as np
import pandas as pd

path = '/home/G-parity/ydzhao/24ID/point_s/2300'
files = os.listdir(path)

pts = [np.array(f.split(','), dtype=int) for f in files]


dist = [u - v for u in pts for v in pts]
dist = pd.Series(dist)

print(dist.value_counts()[0:100])
