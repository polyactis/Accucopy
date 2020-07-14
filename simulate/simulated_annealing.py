import os, sys
import pandas as pd
data = pd.read_table("/home/luozhihui/Desktop/accurity/rc_ratio_window_count_smoothed.tsv", header=0, sep="\t")
windows =  data.read_count_ratio

