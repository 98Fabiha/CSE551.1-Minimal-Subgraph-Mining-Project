# -*- coding: utf-8 -*-
"""Minimal Subgraph Mining.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1Mp95hPoo1xoQoDp44CX-SXBLSMHJ2NFp
"""

from google.colab import drive
drive.mount('/content/drive')

!python --version

!pip install numpy

import tensorflow as tf

import pandas as pd

import numpy as np

import seaborn as sns

!pip install rpy2

!install.packages("sdcMicro")

import pandas as pd
from google.colab import files

uploaded = files.upload()
df = pd.read_csv('panda_LAML.csv')

print(df.head())

!install.packages(sdcMicro)

!pip install rpy2

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

pandas2ri.activate()

!r_df = pandas2ri.py2rpy(df)

!ro.r('''
library(sdcMicro)

!data <- as.data.frame(r_df)

pip install suda