from halton_seq import *
from gp_tools import *
from tokamak import *
from tqdm import tqdm
from pygmo import *
import numpy as np
import itertools
import pickle
import json
import multiprocessing
import glob
import os

def init(outer, major, coil, halton, m_batches, m_error, neutrons):
    number_of_datapoints=halton
    lower_x=coil
    upper_x=coil+outer-0.01
    lower_y=coil+0.00001
    upper_y=coil+outer-0.001

    points_to_search=[]
    leak=[]
    leak_err=[]
    points=[]
    points_to_search_double_list = halton_sequence(number_of_datapoints, 2)

    for x,y in zip(points_to_search_double_list[0],points_to_search_double_list[1]):
        new_x = rescale(x, 0.0, 1.0, lower_x, upper_x)
        new_y = rescale(y, 0.0, 1.0, lower_y, upper_y)
        points_to_search.append([new_x,new_y])

    j=0
    for i in range(len(points_to_search)):
       if points_to_search[i][1]>points_to_search[i][0]:
           points.append(points_to_search[i])
           j+=1

    print('Sampling...', multiprocessing.cpu_count(),'cores.')
    pbar = tqdm(total=j)

    for i in range(len(points)):
        leakage, leakage_error = shield(points[i], major, coil, False, outer, m_batches, m_error, neutrons)
        leak.append(leakage)
        leak_err.append(leakage_error)
        pbar.update(1)

    points = np.array(points)
    pbar.close()
    print('Sampling finished.')
    coords = list(zip(points[:,0],points[:,1]))
    print('GP model working...', multiprocessing.cpu_count(),'cores.')
    GP = GpRegressor(coords, leak, leak_err)

    return GP, points, leak, leak_err
