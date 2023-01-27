import numpy as np
import matplotlib.pyplot as plt
import os, sys
import scipy
from scipy import interpolate

def main():
    PDD(100)
    plt.show()

    CR(100,15, "OFF")
    CR(100,15, "ON")
    plt.show()

    Diag(15,"OFF")
    Diag(15,"ON")
    plt.show()

def PDD(fs):
    filename = os.path.join(sys.path[0], "6X_FFF_PDD.csv")
    dose_pdd = np.array([], dtype='f')
    depth_pdd = np.array([], dtype='f')

    f = open(filename, 'r')
   
    for line in f:
        if line.startswith('# FLSZ'):
            spliti = line.split(" ")
            FS = spliti[2]
            XxX = int(FS.split("*")[1])
            if (XxX == fs):
                header = [f.readline().strip() for i in range(3)]
                pnts = int((header[0]).split(" ")[2])
                pdd_data = [f.readline().strip() for i in range(3,pnts+3,1)]
                for j in range(len(pdd_data)):
                    pdd_strip = pdd_data[j].split(" ")
                    depth_pdd = np.append(depth_pdd,float(pdd_strip[0]))
                    dose_pdd = np.append(dose_pdd,float(pdd_strip[1]))

    f.close()
    plt.plot(depth_pdd, dose_pdd, '.-') 
    plt.xlabel('Depth (mm)')
    plt.ylabel('Dose (%)')
    plt.gcf().set_size_inches(15,10)
    plt.grid(True) 


def CR(fs,depth,conv):
    filename = os.path.join(sys.path[0], "6X_FFF_CR_all.csv")
    dose_dp = np.array([], dtype='f')
    depth_pdd = np.array([], dtype='f')

    f = open(filename, 'r')
   
    for line in f:
        if line.startswith('# FLSZ'):
            spliti = line.split(" ")
            FS = spliti[2]
            XxX = int(FS.split("*")[1])
            if (XxX == fs):
                header = [f.readline().strip() for i in range(5)]
                dep = int((header[2]).split(" ")[2])
                convo = str((header[3]).split(" ")[2])
                if (dep == depth) & (convo == conv):
                    pnts = int((header[0]).split(" ")[2])
                    dp_data = [f.readline().strip() for i in range(5,pnts+5,1)]
                    for j in range(len(dp_data)):
                        dp_split = dp_data[j].split(" ")
                        depth_pdd = np.append(depth_pdd,float(dp_split[0]))
                        dose_dp = np.append(dose_dp,float(dp_split[1]))

    f.close()
    plt.plot(depth_pdd, dose_dp, '.-') 
    plt.xlabel('Off-axis distance (mm)')
    plt.ylabel('Dose (%)')
    plt.gcf().set_size_inches(15,10)
    plt.grid(True) 


def Diag(depth,conv):

    filename = os.path.join(sys.path[0], "6X_FFF_Diag.csv")
    dose_dp = np.array([], dtype='f')
    depth_pdd = np.array([], dtype='f')
    fs = 400

    f = open(filename, 'r')
   
    for line in f:
        if line.startswith('# FLSZ'):
            spliti = line.split(" ")
            FS = spliti[2]
            XxX = int(FS.split("*")[1])
            if (XxX == fs):
                header = [f.readline().strip() for i in range(5)]
                dep = int((header[2]).split(" ")[2])
                convo = str((header[3]).split(" ")[2])
                if (dep == depth) & (convo == conv):
                    pnts = int((header[0]).split(" ")[2])
                    dp_data = [f.readline().strip() for i in range(5,pnts+5,1)]
                    for j in range(len(dp_data)):
                        dp_split = dp_data[j].split(" ")
                        depth_pdd = np.append(depth_pdd,float(dp_split[0]))
                        dose_dp = np.append(dose_dp,float(dp_split[1]))

    f.close()
    plt.plot(depth_pdd, dose_dp, '.-') 
    plt.xlabel('Diagonal (mm)')
    plt.ylabel('Dose (%)')
    plt.gcf().set_size_inches(15,10)
    plt.grid(True) 

if __name__ == '__main__':
    main()
