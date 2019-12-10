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
dist = pd.DataFrame(np.random.normal(loc=means, scale=stdevs, size=(50, 4)), columns=['a', 'b', 'c', 'd'])

#fig, ax = plt.subplots()
#dist.plot.kde(ax=ax, legend=False)
#dist.plot.hist(ax=ax, bins=50, alpha=0.5)
#plt.show()

s1 = dist['a'].tolist()
s2 = dist['b'].tolist()
s3 = dist['c'].tolist()
s4 = dist['d'].tolist()

s = s1 + s2 + s3 + s4

with open("data.txt", 'w') as f:

    for n in s:
        f.write(str(n)+"\n")