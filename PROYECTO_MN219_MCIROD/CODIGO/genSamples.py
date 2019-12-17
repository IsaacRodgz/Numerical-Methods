import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

'''
data = np.random.normal(2, 0.5, 25)
data = np.concatenate((data, np.random.normal(4, 0.5, 25)))
data = np.concatenate((data, np.random.normal(6, 0.5, 25)))
data = np.concatenate((data, np.random.normal(8, 0.5, 25)))
'''

#_ = plt.hist(data, bins=40)

means = 2, 4, 6, 8
stdevs = 0.3, 0.3, 0.3, 0.3
dist = pd.DataFrame(np.random.normal(loc=means, scale=stdevs, size=(50, 4)), columns=['mu = 2', 'mu = 4', 'mu = 6', 'mu = 8'])

'''
fig, ax = plt.subplots()
dist.plot.kde(ax=ax, legend=False)
dist.plot.hist(ax=ax, bins=50, alpha=0.5)
plt.show()
'''

s1 = dist['mu = 2'].tolist()
s2 = dist['mu = 4'].tolist()
s3 = dist['mu = 6'].tolist()
s4 = dist['mu = 8'].tolist()

s = s1 + s2 + s3 + s4

with open("data.txt", 'w') as f:

    for n in s:
        f.write(str(n)+"\n")
