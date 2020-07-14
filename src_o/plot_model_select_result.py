#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats.distributions import norm
import argparse
import h5py
import os
import sys

def result_plot(data,result,filename,text):

    plt.subplot(121)
    plt.hist(data,100,alpha=0.5,density=True,color="gray")
    x_grid = np.linspace(min(data),max(data),1000)
    kde = gaussian_kde(data)
    plt.plot(x_grid,kde.evaluate(x_grid),alpha=0.5,label='kde',linewidth=2,
            color='red')
    i = 0
    ss = ""
    for minor,major in result["model"]:
        alpha1,mu,var,convergence = result["arg"][i]
        label_txt = "%d:%d" % (minor,major)
        plt.plot(x_grid,alpha1*norm(mu,math.sqrt(var)).pdf(x_grid)
                 + (1-alpha1)*norm(-mu,math.sqrt(var)).pdf(x_grid), label= label_txt,
                 alpha=0.5, linewidth=2)
        ss += "Model(%d:%d): logL: %.6f\n" % (minor,major,
               result["logL"][i])
        i += 1
    ss += "Selection: %d-%d\n" % (result["selection"][0],result["selection"][1])
    plt.legend()
    plt.subplot(122)
    plt.axis([0,10,0,10])
    
    plt.text(0.1,0.1,text+ss)
    plt.savefig(filename,format='png',dpi=300)
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file",help="The HDF5 file containing the data")
    parser.add_argument("-o","--output", help="The output directory")
    args = parser.parse_args()
    
    if os.path.isfile(args.output):
        sys.stderr.write("Output dir %s is a file. Remove it.\n" % args.output)
        os.remove(args.output)
        os.mkdir(args.output)
    elif os.path.isdir(args.output):
        sys.stderr.write("%s aleady exists.\n" % args.output)
    else:
        os.mkdir(args.output)
    f = h5py.File(args.file, 'r')
    for each_group in f.keys():
        data = None
        result = dict()
        title = None
        group = f[each_group]
        for name, value in group.iteritems():
            if name == 'model_list':
                result["model"] = list(value)
            if name == 'arg_list':
                result["arg"] = list(value)
            if name == 'logL_list':
                result["logL"] = list(value)
            if name == 'selection':
                result['selection'] = list(value)
            if name == 'data':
                data = list(value)
            if name == 'title':
                title = list(value)[0]
        outputfile = os.path.join(args.output, title + ".png")
        result_plot(data,result,outputfile,title+"\n")