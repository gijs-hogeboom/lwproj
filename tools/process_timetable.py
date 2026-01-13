import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

path = '/home/gijs/mclw/data_output/time_table.dat'

df = pd.read_csv(path)
print(df)