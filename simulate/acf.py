#!/usr/bin/env python
import matplotlib; matplotlib.use("Qt4Agg")
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf

data = pd.read_table("/home/luozhihui/Desktop/accurity/rc_ratio_window_count_smoothed.tsv")
vc = data.window_count_smoothed
lag_acf = acf(vc, nlags=1000)
plt.plot(lag_acf)
plt.show()