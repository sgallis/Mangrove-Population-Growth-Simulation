import numpy as np
from functions import Forest, H
from constants import L, years, pop
import matplotlib.pyplot as plt



def solve():
    mangrove = Forest(seed=1, n_cold=4)
    mangrove.initialize(nb_trees=pop)
    mangrove.update(years)




solve()
