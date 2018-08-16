from halton_seq import *
from gp_tools import *
from model import *
from tqdm import tqdm
from pygmo import *
import numpy as np
import itertools
import pickle
from init import *
import matplotlib.pyplot as plt
import multiprocessing
import glob
import os
from terminal_color import bcolors
import time
start = time.time()


"""
outer : sets the outer radius of the shield (i.e., total size) in cm
halton : sets the number of points used to build the first GP function.
m_batches : sets maximum number of batches for openmc
m_error : sets the maximum relative error admissible for openmc in %
neutrons : sets number of neutrons for openmc

leakage and error returned is in %

go to model.py to change the materials (default : tungsten and water)

BE CAREFUL with the directory in model.py
"""

def optimize_with_size(outer, halton, m_batches, m_error, neutrons):

    gp, points, leak, leak_err = init(outer, halton, m_batches, m_error, neutrons)

    print('GP model calculated')
    print('Optimizing...', multiprocessing.cpu_count(),'cores.')

    if __name__ == "__main__":

        class my_problem():
            global outer
            def fitness(self, x):

                leak, error = gp([(x[0], x[1])])
                ci1 = -(x[1]-x[0])

                return (leak[0], ci1)

            def get_nic(self):
                return 1

            def get_nobj(self):
                return 1

            def get_bounds(self):
                lb = [1, 1]
                ub = [outer-0.01, outer-0.001]
                return (lb,ub)

        gen = 10
        size = 1000
        algo = algorithm(sga(gen = 1, mutation='uniform', crossover='sbx', cr=0.50, m=0.10))
        prob = unconstrain(prob = my_problem(), method = 'death penalty')

        gp_leak = []

        pbar = tqdm(total=len(points))

        for k in range(len(points)):

            pop = population(prob,size)

            for i in range(gen):
                pop = algo.evolve(pop)

            leakage, leakage_error = shield(pop.champion_x, False, outer, m_batches, m_error, neutrons)

            points[k] = pop.champion_x
            leak[k] = leakage
            leak_err[k] = leakage_error
            gp_leak.append(pop.champion_f[0])

            gp = GpRegressor(points, leak, leak_err)

            pbar.update(1)

    pbar.close()
    opt = leak.index(min(leak))

    print(bcolors.BOLD +'---------------')
    print('ITERATION TERMINATED')
    print('BEST SOLUTION FOUND')
    print('////////////////////')
    print('R1:', points[opt][0], 'R2:', points[opt][1])
    print('Leakage (GP) %:', gp_leak[opt])
    print('Leakage (MC) %:', leak[opt])
    print('Leakage error (MC) %:', leak_err[opt])
    print('Relative error (MC) %:', leak_err[opt]/leak[opt]*100)
    print('////////////////////'+ bcolors.ENDC)

optimize_with_size(10, 21, 500, 1, 1000)

end = time.time()
print(end - start, 'seconds.')
