import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import pandas as pd
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams.update({'font.size': 12})

def cutregion_AgTri(eField,zlow=10.1e-9,zhigh=11.1e-9,xlow=-1e-9,xhigh=90e-9,ylow=-43e-9,yhigh=43e-9):
    x = eField[:,0].astype(float)
    y = eField[:,1].astype(float)
    z = eField[:,2].astype(float)
    Ez = eField[:,3]
    Ey = eField[:,4]
    Ex = eField[:,5]
    xlist = []
    ylist = []
    zlist = []
    Ezlist = []
    Eylist = []
    Exlist = []
    for k in range(len(z)):
        if z[k] < zhigh and z[k] > zlow:
            if x[k] < xhigh and x[k] > xlow:
                if y[k] < yhigh and y[k] > ylow:
                    xlist.append(x[k])
                    ylist.append(y[k])
                    zlist.append(z[k])
                    Ezlist.append(Ez[k])
                    Eylist.append(Ey[k])
                    Exlist.append(Ex[k])

    xlist = np.asarray(xlist)*1e9 #express in nanometers
    ylist = np.asarray(ylist)*1e9
    zlist = np.asarray(zlist)*1e9
    Ezlist = np.complex_(np.asarray(Ezlist))*1e3 #express in millivolts per meter
    Eylist = np.complex_(np.asarray(Eylist))*1e3
    Exlist = np.complex_(np.asarray(Exlist))*1e3
    return xlist,ylist,zlist,Ezlist,Eylist,Exlist


def cutregion_ceria(eField,zlow=60e-9,zhigh=100e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9):
    x = eField[:,0].astype(float)
    y = eField[:,1].astype(float)
    z = eField[:,2].astype(float)
    Ez = eField[:,3]
    Ey = eField[:,4]
    Ex = eField[:,5]
    xlist = []
    ylist = []
    zlist = []
    Ezlist = []
    Eylist = []
    Exlist = []
    for k in range(len(z)):
        if z[k] < zhigh and z[k] > zlow:
            if x[k] < xhigh and x[k] > xlow:
                if y[k] < yhigh and y[k] > ylow:
                    xlist.append(x[k])
                    ylist.append(y[k])
                    zlist.append(z[k])
                    Ezlist.append(Ez[k])
                    Eylist.append(Ey[k])
                    Exlist.append(Ex[k])

    xlist = np.asarray(xlist)*1e9 #express in nanometers
    ylist = np.asarray(ylist)*1e9
    zlist = np.asarray(zlist)*1e9
    Ezlist = np.complex_(np.asarray(Ezlist))*1e3 #express in millivolts per meter
    Eylist = np.complex_(np.asarray(Eylist))*1e3
    Exlist = np.complex_(np.asarray(Exlist))*1e3
    return xlist,ylist,zlist,Ezlist,Eylist,Exlist


###### 1_,2_,3_,4_,5_ correspond to 3.1883e-19,4.2538e-19,4.6944e-19,5.1911e-19,and 5.7198e-19 J modes ######
###### _1 = corner impact (A) 4 nm away, _2 = side impact (B) 4 nm away, _3 = bulk center (C) impact, _4 = ceria cube with Ag triangle, _5 = just ceria cube, _6 = side impact 164 nm away ######


def useless():
    # eField2_1 = np.asarray(pd.read_csv('p_AgTriangle_4.2538_Aparam_electricField.txt',sep='\t',dtype=str))
    # x2_1,y2_1,z2_1,Ez2_1,Ey2_1,Ex2_1 = cutregion_AgTri(eField2_1)
    # eField2_3 = np.asarray(pd.read_csv('p_AgTriangle_4.2538_Cparam_electricField.txt',sep='\t',dtype=str))
    # x2_3,y2_3,z2_3,Ez2_3,Ey2_3,Ex2_3 = cutregion_AgTri(eField2_3)
    # eField2_4 = np.asarray(pd.read_csv('p_ceria160nmAndAgTriangle78nm_4.2538mode_Efield.txt',sep='\t',dtype=str))
    # x2_4,y2_4,z2_4,Ez2_4,Ey2_4,Ex2_4 = cutregion_AgTri(eField2_4,zlow=90.5e-9,zhigh=91.5e-9,xlow=0e-9,xhigh=250e-9,ylow=-85e-9,yhigh=85e-9)
    # eField2_5 = np.asarray(pd.read_csv('p_ceria160nm_4.2538mode_electricField.txt',sep='\t',dtype=str))
    # x2_5,y2_5,z2_5,Ez2_5,Ey2_5,Ex2_5 = cutregion_AgTri(eField2_5,zlow=90.5e-9,zhigh=91.5e-9,xlow=0e-9,xhigh=250e-9,ylow=-85e-9,yhigh=85e-9)
    # eField3_4 = np.asarray(pd.read_csv('p_ceria160nmAndAgTriangle78nm_4.6944mode_Efield.txt',sep='\t',dtype=str))
    # x3_4,y3_4,z3_4,Ez3_4,Ey3_4,Ex3_4 = cutregion_AgTri(eField3_4,zlow=90.5e-9,zhigh=91.5e-9,xlow=0e-9,xhigh=250e-9,ylow=-85e-9,yhigh=85e-9)
    # eField3_5 = np.asarray(pd.read_csv('p_ceria160nm_4.6944mode_electricField.txt',sep='\t',dtype=str))
    # x3_5,y3_5,z3_5,Ez3_5,Ey3_5,Ex3_5 = cutregion_AgTri(eField3_5,zlow=90.5e-9,zhigh=91.5e-9,xlow=0e-9,xhigh=250e-9,ylow=-85e-9,yhigh=85e-9)
    # eField4_2 = np.asarray(pd.read_csv('p_AgTriangle_5.1911_Bparam_electricField.txt',sep='\t',dtype=str))
    # x4_2,y4_2,z4_2,Ez4_2,Ey4_2,Ex4_2 = cutregion_AgTri(eField4_2,zlow=85e-9,zhigh=115e-9,xlow=-160e-9,xhigh=90e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_2 = np.sqrt((np.real(Ex4_2))**2+(np.real(Ey4_2))**2+(np.real(Ez4_2))**2)
    # x4_5_2,y4_5_2,z4_5_2,Ez4_5_2,Ey4_5_2,Ex4_5_2 = cutregion_ceria(eField4_5,zlow=60e-9,zhigh=100e-9,xlow=20e-9,xhigh=40e-9)
    # energyDensity4_5_2 = np.sqrt((np.real(Ex4_5_2))**2+(np.real(Ey4_5_2))**2+(np.real(Ez4_5_2))**2)
    # x4_5_4,y4_5_4,z4_5_4,Ez4_5_4,Ey4_5_4,Ex4_5_4 = cutregion_ceria(eField4_5,zlow=60e-9,zhigh=100e-9,xlow=60e-9,xhigh=80e-9)
    # energyDensity4_5_4 = np.sqrt((np.real(Ex4_5_4))**2+(np.real(Ey4_5_4))**2+(np.real(Ez4_5_4))**2)
    # x4_5_6,y4_5_6,z4_5_6,Ez4_5_6,Ey4_5_6,Ex4_5_6 = cutregion_ceria(eField4_5,zlow=60e-9,zhigh=100e-9,xlow=100e-9,xhigh=120e-9)
    # energyDensity4_5_6 = np.sqrt((np.real(Ex4_5_6))**2+(np.real(Ey4_5_6))**2+(np.real(Ez4_5_6))**2)
    # x4_5_8,y4_5_8,z4_5_8,Ez4_5_8,Ey4_5_8,Ex4_5_8 = cutregion_ceria(eField4_5,zlow=60e-9,zhigh=100e-9,xlow=140e-9,xhigh=160e-9)
    # energyDensity4_5_8 = np.sqrt((np.real(Ex4_5_8))**2+(np.real(Ey4_5_8))**2+(np.real(Ez4_5_8))**2)
    # x4_5_9,y4_5_9,z4_5_9,Ez4_5_9,Ey4_5_9,Ex4_5_9 = cutregion_ceria(eField4_5,xlow=240e-9,xhigh=260e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_5_9 = np.sqrt((np.real(Ex4_5_9))**2+(np.real(Ey4_5_9))**2+(np.real(Ez4_5_9))**2)
    # x4_5_10,y4_5_10,z4_5_10,Ez4_5_10,Ey4_5_10,Ex4_5_10 = cutregion_ceria(eField4_5,xlow=260e-9,xhigh=280e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_5_10 = np.sqrt((np.real(Ex4_5_10))**2+(np.real(Ey4_5_10))**2+(np.real(Ez4_5_10))**2)
    # x4_5_11,y4_5_11,z4_5_11,Ez4_5_11,Ey4_5_11,Ex4_5_11 = cutregion_ceria(eField4_5,xlow=280e-9,xhigh=300e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_5_11 = np.sqrt((np.real(Ex4_5_11))**2+(np.real(Ey4_5_11))**2+(np.real(Ez4_5_11))**2)
    # x4_5_12,y4_5_12,z4_5_12,Ez4_5_12,Ey4_5_12,Ex4_5_12 = cutregion_ceria(eField4_5,xlow=300e-9,xhigh=320e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_5_12 = np.sqrt((np.real(Ex4_5_12))**2+(np.real(Ey4_5_12))**2+(np.real(Ez4_5_12))**2)
    # x4_5_13,y4_5_13,z4_5_13,Ez4_5_13,Ey4_5_13,Ex4_5_13 = cutregion_ceria(eField4_5,xlow=320e-9,xhigh=340e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_5_13 = np.sqrt((np.real(Ex4_5_13))**2+(np.real(Ey4_5_13))**2+(np.real(Ez4_5_13))**2)
    # x4_5_14,y4_5_14,z4_5_14,Ez4_5_14,Ey4_5_14,Ex4_5_14 = cutregion_ceria(eField4_5,xlow=340e-9,xhigh=360e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_5_14 = np.sqrt((np.real(Ex4_5_14))**2+(np.real(Ey4_5_14))**2+(np.real(Ez4_5_14))**2)
    # x4_5_15,y4_5_15,z4_5_15,Ez4_5_15,Ey4_5_15,Ex4_5_15 = cutregion_ceria(eField4_5,xlow=360e-9,xhigh=380e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_5_15 = np.sqrt((np.real(Ex4_5_15))**2+(np.real(Ey4_5_15))**2+(np.real(Ez4_5_15))**2)
    # x4_5_16,y4_5_16,z4_5_16,Ez4_5_16,Ey4_5_16,Ex4_5_16 = cutregion_ceria(eField4_5,xlow=380e-9,xhigh=400e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_5_16 = np.sqrt((np.real(Ex4_5_16))**2+(np.real(Ey4_5_16))**2+(np.real(Ez4_5_16))**2)
    # eField4_7 = np.asarray(pd.read_csv('p_ceria160nm_5.3833mode_electricField.txt',sep='\t',dtype=str))
    # x4_7,y4_7,z4_7,Ez4_7,Ey4_7,Ex4_7 = cutregion_AgTri(eField4_7,zlow=60e-9,zhigh=100e-9,xlow=20e-9,xhigh=150e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_7 = np.sqrt((np.real(Ex4_7))**2+(np.real(Ey4_7))**2+(np.real(Ez4_7))**2)
    # x4_8_2,y4_8_2,z4_8_2,Ez4_8_2,Ey4_8_2,Ex4_8_2 = cutregion_ceria(eField4_8,xlow=20e-9,xhigh=30e-9)
    # energyDensity4_8_2 = np.sqrt((np.real(Ex4_8_2))**2+(np.real(Ey4_8_2))**2+(np.real(Ez4_8_2))**2)
    # x4_8_4,y4_8_4,z4_8_4,Ez4_8_4,Ey4_8_4,Ex4_8_4 = cutregion_ceria(eField4_8,xlow=60e-9,xhigh=70e-9)
    # energyDensity4_8_4 = np.sqrt((np.real(Ex4_8_4))**2+(np.real(Ey4_8_4))**2+(np.real(Ez4_8_4))**2)
    # x4_8_6,y4_8_6,z4_8_6,Ez4_8_6,Ey4_8_6,Ex4_8_6 = cutregion_ceria(eField4_8,xlow=100e-9,xhigh=110e-9)
    # energyDensity4_8_6 = np.sqrt((np.real(Ex4_8_6))**2+(np.real(Ey4_8_6))**2+(np.real(Ez4_8_6))**2)
    # x4_8_8,y4_8_8,z4_8_8,Ez4_8_8,Ey4_8_8,Ex4_8_8 = cutregion_ceria(eField4_8,xlow=140e-9,xhigh=150e-9)
    # energyDensity4_8_8 = np.sqrt((np.real(Ex4_8_8))**2+(np.real(Ey4_8_8))**2+(np.real(Ez4_8_8))**2)
    # x4_8_9,y4_8_9,z4_8_9,Ez4_8_9,Ey4_8_9,Ex4_8_9 = cutregion_ceria(eField4_8,xlow=240e-9,xhigh=250e-9)
    # energyDensity4_8_9 = np.sqrt((np.real(Ex4_8_9))**2+(np.real(Ey4_8_9))**2+(np.real(Ez4_8_9))**2)
    # x4_8_11,y4_8_11,z4_8_11,Ez4_8_11,Ey4_8_11,Ex4_8_11 = cutregion_ceria(eField4_8,xlow=280e-9,xhigh=290e-9)
    # energyDensity4_8_11 = np.sqrt((np.real(Ex4_8_11))**2+(np.real(Ey4_8_11))**2+(np.real(Ez4_8_11))**2)
    # x4_8_13,y4_8_13,z4_8_13,Ez4_8_13,Ey4_8_13,Ex4_8_13 = cutregion_ceria(eField4_8,xlow=320e-9,xhigh=330e-9)
    # energyDensity4_8_13 = np.sqrt((np.real(Ex4_8_13))**2+(np.real(Ey4_8_13))**2+(np.real(Ez4_8_13))**2)
    # x4_8_15,y4_8_15,z4_8_15,Ez4_8_15,Ey4_8_15,Ex4_8_15 = cutregion_ceria(eField4_8,xlow=360e-9,xhigh=370e-9)
    # energyDensity4_8_15 = np.sqrt((np.real(Ex4_8_15))**2+(np.real(Ey4_8_15))**2+(np.real(Ez4_8_15))**2)
    # eField4_9 = np.asarray(pd.read_csv('p_ceriaDimer160nm_5.1991mode_electricField.txt',sep='\t',dtype=str))
    # x4_9,y4_9,z4_9,Ez4_9,Ey4_9,Ex4_9 = cutregion_AgTri(eField4_9,zlow=60e-9,zhigh=100e-9,xlow=0e-9,xhigh=400e-9)
    # energyDensity4_9 = np.sqrt((np.real(Ex4_9))**2+(np.real(Ey4_9))**2+(np.real(Ez4_9))**2)
    # x4_10_2,y4_10_2,z4_10_2,Ez4_10_2,Ey4_10_2,Ex4_10_2 = cutregion_ceria(eField4_10,xlow=20e-9,xhigh=30e-9)
    # energyDensity4_10_2 = np.sqrt((np.real(Ex4_10_2))**2+(np.real(Ey4_10_2))**2+(np.real(Ez4_10_2))**2)
    # x4_10_4,y4_10_4,z4_10_4,Ez4_10_4,Ey4_10_4,Ex4_10_4 = cutregion_ceria(eField4_10,xlow=60e-9,xhigh=70e-9)
    # energyDensity4_10_4 = np.sqrt((np.real(Ex4_10_4))**2+(np.real(Ey4_10_4))**2+(np.real(Ez4_10_4))**2)
    # x4_10_6,y4_10_6,z4_10_6,Ez4_10_6,Ey4_10_6,Ex4_10_6 = cutregion_ceria(eField4_10,xlow=100e-9,xhigh=110e-9)
    # energyDensity4_10_6 = np.sqrt((np.real(Ex4_10_6))**2+(np.real(Ey4_10_6))**2+(np.real(Ez4_10_6))**2)
    # x4_10_8,y4_10_8,z4_10_8,Ez4_10_8,Ey4_10_8,Ex4_10_8 = cutregion_ceria(eField4_10,xlow=140e-9,xhigh=150e-9)
    # energyDensity4_10_8 = np.sqrt((np.real(Ex4_10_8))**2+(np.real(Ey4_10_8))**2+(np.real(Ez4_10_8))**2)
    # x4_10_9,y4_10_9,z4_10_9,Ez4_10_9,Ey4_10_9,Ex4_10_9 = cutregion_ceria(eField4_10,xlow=240e-9,xhigh=250e-9)
    # energyDensity4_10_9 = np.sqrt((np.real(Ex4_10_9))**2+(np.real(Ey4_10_9))**2+(np.real(Ez4_10_9))**2)
    # x4_10_11,y4_10_11,z4_10_11,Ez4_10_11,Ey4_10_11,Ex4_10_11 = cutregion_ceria(eField4_10,xlow=280e-9,xhigh=290e-9)
    # energyDensity4_10_11 = np.sqrt((np.real(Ex4_10_11))**2+(np.real(Ey4_10_11))**2+(np.real(Ez4_10_11))**2)
    # x4_10_13,y4_10_13,z4_10_13,Ez4_10_13,Ey4_10_13,Ex4_10_13 = cutregion_ceria(eField4_10,xlow=320e-9,xhigh=330e-9)
    # energyDensity4_10_13 = np.sqrt((np.real(Ex4_10_13))**2+(np.real(Ey4_10_13))**2+(np.real(Ez4_10_13))**2)
    # x4_10_16,y4_10_16,z4_10_16,Ez4_10_16,Ey4_10_16,Ex4_10_16 = cutregion_ceria(eField4_10,xlow=380e-9,xhigh=390e-9)
    # energyDensity4_10_16 = np.sqrt((np.real(Ex4_10_16))**2+(np.real(Ey4_10_16))**2+(np.real(Ez4_10_16))**2)
    # eField4_11 = np.asarray(pd.read_csv('p_ceria250nm_5.135mode_electricField.txt',sep='\t',dtype=str))
    # eField4_12 = np.asarray(pd.read_csv('p_ceria250nm_5.1911mode_electricField.txt',sep='\t',dtype=str))
    # x4_12,y4_12,z4_12,Ez4_12,Ey4_12,Ex4_12 = cutregion_AgTri(eField4_12,zlow=0e-9,zhigh=250e-9,xlow=0e-9,xhigh=250e-9,ylow=-125e-9,yhigh=125e-9)
    # energyDensity4_12 = np.sqrt((np.real(Ex4_12))**2+(np.real(Ey4_12))**2+(np.real(Ez4_12))**2)
    # eField4_13 = np.asarray(pd.read_csv('p_ceria250nm_5.3192mode_electricField.txt',sep='\t',dtype=str))
    # eField4_14 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_5.1991mode_Efield.txt',sep='\t',dtype=str))
    # x4_14,y4_14,z4_14,Ez4_14,Ey4_14,Ex4_14 = cutregion_AgTri(eField4_14,zlow=85e-9,zhigh=115e-9,xlow=10e-9,xhigh=250e-9,ylow=-80e-9,yhigh=80e-9)
    # energyDensity4_14 = np.sqrt((np.real(Ex4_14))**2+(np.real(Ey4_14))**2+(np.real(Ez4_14))**2)
    # eField5_4 = np.asarray(pd.read_csv('p_ceria160nmAndAgTriangle78nm_5.7198mode_Efield.txt',sep='\t',dtype=str))
    # x5_4,y5_4,z5_4,Ez5_4,Ey5_4,Ex5_4 = cutregion_AgTri(eField5_4,zlow=90.5e-9,zhigh=91.5e-9,xlow=0e-9,xhigh=250e-9,ylow=-85e-9,yhigh=85e-9)
    # eField5_5 = np.asarray(pd.read_csv('p_ceria160nm_5.7198mode_electricField.txt',sep='\t',dtype=str))
    # x5_5,y5_5,z5_5,Ez5_5,Ey5_5,Ex5_5 = cutregion_AgTri(eField5_5,zlow=90.5e-9,zhigh=91.5e-9,xlow=0e-9,xhigh=250e-9,ylow=-85e-9,yhigh=85e-9)
    def ampl(Ez):
        amp = np.max(np.real(Ez)) - np.min(np.real(Ez))
        return amp
        amp1_1 = ampl(Ez1_1)
        amp1_2 = ampl(Ez1_2)
        amp1_3 = ampl(Ez1_3)

        print('amplitude of A for 1.99eV: ', amp1_1)
        print('amplitude of B for 1.99eV: ', amp1_2)
        print('amplitude of C for 1.99eV: ', amp1_3)

        amp2_1 = ampl(Ez2_1)
        amp2_2 = ampl(Ez2_2)
        amp2_3 = ampl(Ez2_3)

        print('amplitude of A for 2.66eV: ', amp2_1)
        print('amplitude of B for 2.66eV: ', amp2_2)
        print('amplitude of C for 2.66eV: ', amp2_3)

        amp3_1 = ampl(Ez3_1)
        amp3_2 = ampl(Ez3_2)
        amp3_3 = ampl(Ez3_3)

        print('amplitude of A for 2.93eV: ', amp3_1)
        print('amplitude of B for 2.93eV: ', amp3_2)
        print('amplitude of C for 2.93eV: ', amp3_3)

        amp5_1 = ampl(Ez5_1)
        amp5_2 = ampl(Ez5_2)
        amp5_3 = ampl(Ez5_3)

        print('amplitude of A for 3.57eV: ', amp5_1)
        print('amplitude of B for 3.57eV: ', amp5_2)
        print('amplitude of C for 3.57eV: ', amp5_3)

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, projection='3d')
        ax1.scatter(xs=x16,ys=y16,zs=np.real(Ez16))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')
        ax2.scatter(xs=x2,ys=y2,zs=np.real(Ez2))
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111, projection='3d')
        ax3.scatter(xs=x3,ys=y3,zs=np.real(Ez3))
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111, projection='3d')
        ax4.scatter(xs=x4,ys=y4,zs=np.real(Ez4))
        fig5 = plt.figure()
        ax5 = fig5.add_subplot(111, projection='3d')
        ax5.scatter(xs=x5,ys=y5,zs=np.real(Ez5))
        fig6 = plt.figure()
        ax6 = fig6.add_subplot(111, projection='3d')
        ax6.scatter(xs=x6,ys=y6,zs=np.real(Ez6))

set0 = 0
set1 = 0
set2 = 0
set3 = 1
sum4 = 0
set4 = 0
set5 = 0


if set0 == 1:
    eField1_2 = np.asarray(pd.read_csv('p_AgTriangle_3.1883_Bparam_electricField.txt',sep='\t',dtype=str))
    x1_2,y1_2,z1_2,Ez1_2,Ey1_2,Ex1_2 = cutregion_AgTri(eField1_2,xlow=20e-9,xhigh=60e-9)
    energyDensity1_2 = np.sqrt((np.real(Ex1_2))**2+(np.real(Ey1_2))**2+(np.real(Ez1_2))**2)
    eField1_6 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_3.1883mode_Efield.txt',sep='\t',dtype=str))
    x1_6,y1_6,z1_6,Ez1_6,Ey1_6,Ex1_6 = cutregion_AgTri(eField1_6,zlow=90.5e-9,zhigh=91.5e-9,xlow=180e-9,xhigh=220e-9,ylow=-43e-9,yhigh=43e-9)
    energyDensity1_6 = np.sqrt((np.real(Ex1_6))**2+(np.real(Ey1_6))**2+(np.real(Ez1_6))**2)
    # eField2_2 = np.asarray(pd.read_csv('p_AgTriangle_4.2538_Bparam_electricField.txt',sep='\t',dtype=str))
    # x2_2,y2_2,z2_2,Ez2_2,Ey2_2,Ex2_2 = cutregion_AgTri(eField2_2,xlow=20e-9,xhigh=60e-9)
    # eField2_6 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_4.2538mode_Efield.txt',sep='\t',dtype=str))
    # x2_6,y2_6,z2_6,Ez2_6,Ey2_6,Ex2_6 = cutregion_AgTri(eField2_6,zlow=90.5e-9,zhigh=91.5e-9,xlow=180e-9,xhigh=220e-9,ylow=-43e-9,yhigh=43e-9)
    eField3_2 = np.asarray(pd.read_csv('p_AgTriangle_4.6944_Bparam_electricField.txt',sep='\t',dtype=str))
    x3_2,y3_2,z3_2,Ez3_2,Ey3_2,Ex3_2 = cutregion_AgTri(eField3_2,xlow=20e-9,xhigh=60e-9)
    eField3_6 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_4.6944mode_Efield.txt',sep='\t',dtype=str))
    x3_6,y3_6,z3_6,Ez3_6,Ey3_6,Ex3_6 = cutregion_AgTri(eField3_6,zlow=90.5e-9,zhigh=91.5e-9,xlow=180e-9,xhigh=220e-9,ylow=-43e-9,yhigh=43e-9)
    # eField5_2 = np.asarray(pd.read_csv('p_AgTriangle_5.7198_Bparam_electricField.txt',sep='\t',dtype=str))
    # x5_2,y5_2,z5_2,Ez5_2,Ey5_2,Ex5_2 = cutregion_AgTri(eField5_2,zlow=4.8e-9,zhigh=5.2e-9,)
    # eField5_6 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_5.7198mode_Efield.txt',sep='\t',dtype=str))
    # x5_6,y5_6,z5_6,Ez5_6,Ey5_6,Ex5_6 = cutregion_AgTri(eField5_6,zlow=90.5e-9,zhigh=91.5e-9,xlow=180e-9,xhigh=220e-9,ylow=-43e-9,yhigh=43e-9)
    ymax0_1 = np.max(np.real(Ez1_6))
    ymin0_1 = np.min(np.real(Ez1_6))
    ymax0_2 = np.max(np.real(Ez3_6))
    ymin0_2 = np.min(np.real(Ez3_6))
    # ymax0_3 =  np.max(np.real(Ez5_6))
    # ymin0_3 = np.min(np.real(Ez5_6))

    eels1 = np.asarray(pd.read_csv('p_AgtriangularNP_Aparam.txt',sep='\t'))
    e1 = eels1[:,0]/1.602e-19
    s1 = eels1[:,1]
    eels2 = np.asarray(pd.read_csv('p_AgtriangularNP_Bparam.txt',sep='\t'))
    e2 = eels2[:,0]/1.602e-19
    s2 = eels2[:,1]
    eels3 = np.asarray(pd.read_csv('p_AgtriangularNP_Cparam.txt',sep='\t'))
    e3 = eels3[:,0]/1.602e-19
    s3 = eels3[:,1]


def plotset0():

    plt.figure(0)
    plt.plot(e1,s1/np.max(s1),color='black')
    plt.plot(e2,s2/np.max(s1),color='darkslateblue')
    plt.plot(e3,s3/np.max(s1),color='red')
    plt.ylim(0,1.05)
    plt.yticks([],[])
    plt.show()
    q

    xn = 2
    yn = 4
    plt.figure(1)
    # plt.suptitle('Electric Field Z-Component Distributions for 4nm and 164nm Ag Nanoparticle',fontsize=20)
    plt.subplot2grid((xn,yn),(0,0))
    # plt.subplots_adjust(hspace=0.3)
    # plt.ylabel(r'$E_z$ (mV/m)',fontsize=18)
    plt.title('2.0 eV mode',fontsize=18)
    plt.plot(-x1_2+78,np.real(Ez1_2),'o',ms=5,color='coral')
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.subplot2grid((xn,yn),(0,1))
    plt.title('2.9 eV mode',fontsize=18)
    plt.plot(-x3_2+78,np.real(Ez3_2),'o',ms=5,color='coral')
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.subplot2grid((xn,yn),(1,0))
    plt.plot(x1_6-160,np.real(Ez1_6),'o',ms=5,color='deepskyblue')
    plt.ylim(1.3*ymin0_1,1.3*ymax0_1)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.subplot2grid((xn,yn),(1,1))
    plt.plot(x3_6-160,np.real(Ez3_6),'o',ms=5,color='deepskyblue')
    plt.ylim(1.3*ymin0_2,1.3*ymax0_2)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.subplot2grid((xn,yn),(0,2))
    plt.title('2.0 eV mode',fontsize=18)
    plt.plot(y1_2,np.real(Ez1_2),'o',ms=5,color='firebrick')
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.subplot2grid((xn,yn),(0,3))
    plt.title('2.9 eV mode',fontsize=18)
    plt.plot(y3_2,np.real(Ez3_2),'o',ms=5,color='firebrick')
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.subplot2grid((xn,yn),(1,2))
    plt.plot(y1_6,np.real(Ez1_6),'o',ms=5,color='cornflowerblue')
    plt.ylim(1.3*ymin0_1,1.3*ymax0_1)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.subplot2grid((xn,yn),(1,3))
    plt.plot(y3_6,np.real(Ez3_6),'o',ms=5,color='cornflowerblue')
    plt.ylim(1.3*ymin0_2,1.3*ymax0_2)
    plt.tick_params(axis='x',which='minor',length=0)


if set1 == 1:
    eField1_1 = np.asarray(pd.read_csv('p_AgTriangle_3.1883_Aparam_electricField.txt',sep='\t',dtype=str))
    x1_1x,y1_1x,z1_1x,Ez1_1x,Ey1_1x,Ex1_1x = cutregion_AgTri(eField1_1)
    x1_1y,y1_1y,z1_1y,Ez1_1y,Ey1_1y,Ex1_1y = cutregion_AgTri(eField1_1)
    # energyDensity1_1 = np.sqrt((np.real(Ex1_1))**2+(np.real(Ey1_1))**2+(np.real(Ez1_1))**2)
    eField1_3 = np.asarray(pd.read_csv('p_AgTriangle_3.1883_Cparam_electricField.txt',sep='\t',dtype=str))
    x1_3x,y1_3x,z1_3x,Ez1_3x,Ey1_3x,Ex1_3x = cutregion_AgTri(eField1_3)
    x1_3y,y1_3y,z1_3y,Ez1_3y,Ey1_3y,Ex1_3y = cutregion_AgTri(eField1_3)
    # energyDensity1_3 = np.sqrt((np.real(Ex1_3))**2+(np.real(Ey1_3))**2+(np.real(Ez1_3))**2)


    eField3_1 = np.asarray(pd.read_csv('p_AgTriangle_4.6944_Aparam_electricField.txt',sep='\t',dtype=str))
    x3_1x,y3_1x,z3_1x,Ez3_1x,Ey3_1x,Ex3_1x = cutregion_AgTri(eField3_1)
    x3_1y,y3_1y,z3_1y,Ez3_1y,Ey3_1y,Ex3_1y = cutregion_AgTri(eField3_1)
    eField3_3 = np.asarray(pd.read_csv('p_AgTriangle_4.6944_Cparam_electricField.txt',sep='\t',dtype=str))
    x3_3x,y3_3x,z3_3x,Ez3_3x,Ey3_3x,Ex3_3x = cutregion_AgTri(eField3_3)
    x3_3y,y3_3y,z3_3y,Ez3_3y,Ey3_3y,Ex3_3y = cutregion_AgTri(eField3_3)


    # eField5_1 = np.asarray(pd.read_csv('p_AgTriangle_5.7198_Aparam_electricField.txt',sep='\t',dtype=str))
    # x5_1,y5_1,z5_1,Ez5_1,Ey5_1,Ex5_1 = cutregion_AgTri(eField5_1,zlow=4.8e-9,zhigh=5.2e-9)
    # eField5_3 = np.asarray(pd.read_csv('p_AgTriangle_5.7198_Cparam_electricField.txt',sep='\t',dtype=str))
    # x5_3,y5_3,z5_3,Ez5_3,Ey5_3,Ex5_3 = cutregion_AgTri(eField5_3,zlow=4.8e-9,zhigh=5.2e-9)

    if set0 != 1:
        # eField5_2 = np.asarray(pd.read_csv('p_AgTriangle_5.7198_Bparam_electricField.txt',sep='\t',dtype=str))
        # x5_2,y5_2,z5_2,Ez5_2,Ey5_2,Ex5_2 = cutregion_AgTri(eField5_2,zlow=4.8e-9,zhigh=5.2e-9)
        eField3_2 = np.asarray(pd.read_csv('p_AgTriangle_4.6944_Bparam_electricField.txt',sep='\t',dtype=str))
        x3_2x,y3_2x,z3_2x,Ez3_2x,Ey3_2x,Ex3_2x = cutregion_AgTri(eField3_2)
        x3_2y,y3_2y,z3_2y,Ez3_2y,Ey3_2y,Ex3_2y = cutregion_AgTri(eField3_2)
        eField1_2 = np.asarray(pd.read_csv('p_AgTriangle_3.1883_Bparam_electricField.txt',sep='\t',dtype=str))
        x1_2x,y1_2x,z1_2x,Ez1_2x,Ey1_2x,Ex1_2x = cutregion_AgTri(eField1_2)
        x1_2y,y1_2y,z1_2y,Ez1_2y,Ey1_2y,Ex1_2y = cutregion_AgTri(eField1_2)
        # energyDensity1_2 = np.sqrt((np.real(Ex1_2))**2+(np.real(Ey1_2))**2+(np.real(Ez1_2))**2)

def plotset1():

    plt.rcParams.update({'font.size': 12})

    xn = 3
    yn = 4
    all1 = np.r_[Ez1_1x,Ez1_2x,Ez1_3x,Ez1_1y,Ez1_2y,Ez1_3y,Ez3_1x,Ez3_2x,Ez3_3x,Ez3_1y,Ez3_2y,Ez3_3y]
    ymin = np.min(np.real(all1))
    ymax = np.max(np.real(all1))
    # all3 = np.r_[Ez3_1x,Ez3_2x,Ez3_3x,Ez3_1y,Ez3_2y,Ez3_3y]
    # ymin3 = np.min(np.real(all3))
    # ymax3 = np.max(np.real(all3))
    # all5 = np.r_[Ez5_1,Ez5_2,Ez5_3]
    # ymin5 = np.min(np.real(all5))
    # ymax5 = np.max(np.real(all5))

    plt.figure(11)
    # plt.suptitle(r'Plasmonic Mode $E_z$ Distributions for Corner, Edge, and Bulk Beam Placements',fontsize=22)
    plt.subplot2grid((xn,yn),(0,0))
    # plt.subplots_adjust(hspace=0.3)
    plt.title('2.0 eV mode',fontsize=18)
    # plt.ylabel(r'$E_z$ (mV/m)',fontsize=18)
    plt.plot(x1_1x,np.real(Ez1_1x),'o',ms=7.5,color='blueviolet')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='x',labelsize=0)

    plt.subplot2grid((xn,yn),(0,1))
    plt.title('2.9 eV mode',fontsize=18)
    plt.plot(x3_1x,np.real(Ez3_1x),'o',ms=7.5,color='blueviolet')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='both',labelsize=0)

    plt.subplot2grid((xn,yn),(1,0))
    # plt.ylabel('Edge')
    plt.plot(-x1_2x+78,np.real(Ez1_2x),'o',ms=7.5,color='coral')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='x',labelsize=0)

    plt.subplot2grid((xn,yn),(1,1))
    plt.plot(-x3_2x+78,np.real(Ez3_2x),'o',ms=7.5,color='coral')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='both',labelsize=0)

    plt.subplot2grid((xn,yn),(2,0))
    # plt.ylabel('Bulk')
    plt.plot(x1_3x,np.real(Ez1_3x),'o',ms=7.5,color='springgreen')
    plt.ylim(1.4*ymin,1.4*ymax)

    plt.subplot2grid((xn,yn),(2,1))
    plt.plot(x3_3x,np.real(Ez3_3x),'o',ms=7.5,color='springgreen')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='y',labelsize=0)

    plt.subplot2grid((xn,yn),(0,2))
    plt.title('2.0 eV mode',fontsize=18)
    plt.plot(y1_1y,np.real(Ez1_1y),'o',ms=7.5,color='darkslateblue')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='both',labelsize=0)

    plt.subplot2grid((xn,yn),(0,3))
    plt.title('2.9 eV mode',fontsize=18)
    plt.plot(y3_1y,np.real(Ez3_1y),'o',ms=7.5,color='darkslateblue')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='both',labelsize=0)

    plt.subplot2grid((xn,yn),(1,2))
    plt.plot(y1_2y,np.real(Ez1_2y),'o',ms=7.5,color='firebrick')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='both',labelsize=0)

    plt.subplot2grid((xn,yn),(1,3))
    plt.plot(y3_2y,np.real(Ez3_2y),'o',ms=7.5,color='firebrick')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='both',labelsize=0)

    plt.subplot2grid((xn,yn),(2,2))
    plt.plot(y1_3y,np.real(Ez1_3y),'o',ms=7.5,color='darkolivegreen')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='y',labelsize=0)
    # plt.xlabel('X coordinate (nm)',fontsize=18)

    plt.subplot2grid((xn,yn),(2,3))
    plt.plot(y3_3y,np.real(Ez3_3y),'o',ms=7.5,color='darkolivegreen')
    plt.ylim(1.4*ymin,1.4*ymax)
    plt.tick_params(axis='y',labelsize=0)
    # plt.xlabel('Y coordinate (nm)',fontsize=18)


if set2 == 1:
    eField4_5 = np.asarray(pd.read_csv('p_ceria160nm_5.175mode_electricField.txt',sep='\t',dtype=str))
    x4_5_1,y4_5_1,z4_5_1,Ez4_5_1,Ey4_5_1,Ex4_5_1 = cutregion_ceria(eField4_5,zlow=60e-9,zhigh=100e-9,xlow=0e-9,xhigh=30e-9)
    energyDensity4_5_1 = np.sqrt((np.real(Ex4_5_1))**2+(np.real(Ey4_5_1))**2+(np.real(Ez4_5_1))**2)
    x4_5_3,y4_5_3,z4_5_3,Ez4_5_3,Ey4_5_3,Ex4_5_3 = cutregion_ceria(eField4_5,zlow=60e-9,zhigh=100e-9,xlow=40e-9,xhigh=60e-9)
    energyDensity4_5_3 = np.sqrt((np.real(Ex4_5_3))**2+(np.real(Ey4_5_3))**2+(np.real(Ez4_5_3))**2)
    x4_5_5,y4_5_5,z4_5_5,Ez4_5_5,Ey4_5_5,Ex4_5_5 = cutregion_ceria(eField4_5,zlow=60e-9,zhigh=100e-9,xlow=80e-9,xhigh=100e-9)
    energyDensity4_5_5 = np.sqrt((np.real(Ex4_5_5))**2+(np.real(Ey4_5_5))**2+(np.real(Ez4_5_5))**2)
    x4_5_7,y4_5_7,z4_5_7,Ez4_5_7,Ey4_5_7,Ex4_5_7 = cutregion_ceria(eField4_5,zlow=60e-9,zhigh=100e-9,xlow=120e-9,xhigh=140e-9)
    energyDensity4_5_7 = np.sqrt((np.real(Ex4_5_7))**2+(np.real(Ey4_5_7))**2+(np.real(Ez4_5_7))**2)

    eField4_8 = np.asarray(pd.read_csv('p_ceriaDimer160nm_5.175mode_electricField.txt',sep='\t',dtype=str))
    x4_8_1,y4_8_1,z4_8_1,Ez4_8_1,Ey4_8_1,Ex4_8_1 = cutregion_ceria(eField4_8,xlow=0e-9,xhigh=30e-9)
    energyDensity4_8_1 = np.sqrt((np.real(Ex4_8_1))**2+(np.real(Ey4_8_1))**2+(np.real(Ez4_8_1))**2)
    x4_8_3,y4_8_3,z4_8_3,Ez4_8_3,Ey4_8_3,Ex4_8_3 = cutregion_ceria(eField4_8,xlow=40e-9,xhigh=50e-9)
    energyDensity4_8_3 = np.sqrt((np.real(Ex4_8_3))**2+(np.real(Ey4_8_3))**2+(np.real(Ez4_8_3))**2)
    x4_8_5,y4_8_5,z4_8_5,Ez4_8_5,Ey4_8_5,Ex4_8_5 = cutregion_ceria(eField4_8,xlow=80e-9,xhigh=90e-9)
    energyDensity4_8_5 = np.sqrt((np.real(Ex4_8_5))**2+(np.real(Ey4_8_5))**2+(np.real(Ez4_8_5))**2)
    x4_8_7,y4_8_7,z4_8_7,Ez4_8_7,Ey4_8_7,Ex4_8_7 = cutregion_ceria(eField4_8,xlow=120e-9,xhigh=130e-9)
    energyDensity4_8_7 = np.sqrt((np.real(Ex4_8_7))**2+(np.real(Ey4_8_7))**2+(np.real(Ez4_8_7))**2)
    x4_8_10,y4_8_10,z4_8_10,Ez4_8_10,Ey4_8_10,Ex4_8_10 = cutregion_ceria(eField4_8,xlow=260e-9,xhigh=290e-9)
    energyDensity4_8_10 = np.sqrt((np.real(Ex4_8_10))**2+(np.real(Ey4_8_10))**2+(np.real(Ez4_8_10))**2)
    x4_8_12,y4_8_12,z4_8_12,Ez4_8_12,Ey4_8_12,Ex4_8_12 = cutregion_ceria(eField4_8,xlow=300e-9,xhigh=310e-9)
    energyDensity4_8_12 = np.sqrt((np.real(Ex4_8_12))**2+(np.real(Ey4_8_12))**2+(np.real(Ez4_8_12))**2)
    x4_8_14,y4_8_14,z4_8_14,Ez4_8_14,Ey4_8_14,Ex4_8_14 = cutregion_ceria(eField4_8,xlow=340e-9,xhigh=350e-9)
    energyDensity4_8_14 = np.sqrt((np.real(Ex4_8_14))**2+(np.real(Ey4_8_14))**2+(np.real(Ez4_8_14))**2)
    x4_8_16,y4_8_16,z4_8_16,Ez4_8_16,Ey4_8_16,Ex4_8_16 = cutregion_ceria(eField4_8,xlow=380e-9,xhigh=390e-9)
    energyDensity4_8_16 = np.sqrt((np.real(Ex4_8_16))**2+(np.real(Ey4_8_16))**2+(np.real(Ez4_8_16))**2)

    eField4_10 = np.asarray(pd.read_csv('p_ceriaDimer160nm_5.3833mode_electricField.txt',sep='\t',dtype=str))
    x4_10_1,y4_10_1,z4_10_1,Ez4_10_1,Ey4_10_1,Ex4_10_1 = cutregion_ceria(eField4_10,xlow=0e-9,xhigh=30e-9)
    energyDensity4_10_1 = np.sqrt((np.real(Ex4_10_1))**2+(np.real(Ey4_10_1))**2+(np.real(Ez4_10_1))**2)
    x4_10_3,y4_10_3,z4_10_3,Ez4_10_3,Ey4_10_3,Ex4_10_3 = cutregion_ceria(eField4_10,xlow=40e-9,xhigh=50e-9)
    energyDensity4_10_3 = np.sqrt((np.real(Ex4_10_3))**2+(np.real(Ey4_10_3))**2+(np.real(Ez4_10_3))**2)
    x4_10_5,y4_10_5,z4_10_5,Ez4_10_5,Ey4_10_5,Ex4_10_5 = cutregion_ceria(eField4_10,xlow=80e-9,xhigh=90e-9)
    energyDensity4_10_5 = np.sqrt((np.real(Ex4_10_5))**2+(np.real(Ey4_10_5))**2+(np.real(Ez4_10_5))**2)
    x4_10_7,y4_10_7,z4_10_7,Ez4_10_7,Ey4_10_7,Ex4_10_7 = cutregion_ceria(eField4_10,xlow=120e-9,xhigh=130e-9)
    energyDensity4_10_7 = np.sqrt((np.real(Ex4_10_7))**2+(np.real(Ey4_10_7))**2+(np.real(Ez4_10_7))**2)
    x4_10_10,y4_10_10,z4_10_10,Ez4_10_10,Ey4_10_10,Ex4_10_10 = cutregion_ceria(eField4_10,xlow=260e-9,xhigh=290e-9)
    energyDensity4_10_10 = np.sqrt((np.real(Ex4_10_10))**2+(np.real(Ey4_10_10))**2+(np.real(Ez4_10_10))**2)
    x4_10_12,y4_10_12,z4_10_12,Ez4_10_12,Ey4_10_12,Ex4_10_12 = cutregion_ceria(eField4_10,xlow=300e-9,xhigh=310e-9)
    energyDensity4_10_12 = np.sqrt((np.real(Ex4_10_12))**2+(np.real(Ey4_10_12))**2+(np.real(Ez4_10_12))**2)
    x4_10_14,y4_10_14,z4_10_14,Ez4_10_14,Ey4_10_14,Ex4_10_14 = cutregion_ceria(eField4_10,xlow=340e-9,xhigh=350e-9)
    energyDensity4_10_14 = np.sqrt((np.real(Ex4_10_14))**2+(np.real(Ey4_10_14))**2+(np.real(Ez4_10_14))**2)
    x4_10_15,y4_10_15,z4_10_15,Ez4_10_15,Ey4_10_15,Ex4_10_15 = cutregion_ceria(eField4_10,xlow=360e-9,xhigh=370e-9)
    energyDensity4_10_15 = np.sqrt((np.real(Ex4_10_15))**2+(np.real(Ey4_10_15))**2+(np.real(Ez4_10_15))**2)

    eels4_5 = np.asarray(pd.read_csv('p_ceria160nm_EELspectrum.txt',sep='\t'))
    e4_5 = eels4_5[:,0]/1.602e-19
    s4_5 = eels4_5[:,1]
    eels4_8 = np.asarray(pd.read_csv('p_ceriaDimer160nm_EELspectrum.txt',sep='\t'))
    e4_8 = eels4_8[:,0]/1.602e-19
    s4_8 = eels4_8[:,1]


def plotset2():


    all2 = np.r_[Ez4_5_1,Ez4_8_1,Ez4_8_10,Ez4_10_1,Ez4_10_10]
    ymin2 = np.min(np.real(all2))
    ymax2 = np.max(np.real(all2))

    xn = 3
    yn = 2

    plt.rcParams.update({'font.size': 10})

    plt.figure(21)
    plt.subplot2grid((xn,yn),(0,0))
    plt.tick_params(axis='x',labelsize=0)
    plt.ylim(1.2*ymin2,1.2*ymax2)
    plt.plot(y4_5_1,np.real(Ez4_5_1),'o',ms=4.5,color='indianred')
    plt.subplot2grid((xn,yn),(1,0))
    plt.ylim(1.2*ymin2,1.2*ymax2)
    plt.plot(y4_8_1,np.real(Ez4_8_1),'o',ms=4.5,color='royalblue')
    plt.tick_params(axis='x',labelsize=0)
    plt.subplot2grid((xn,yn),(2,0))
    plt.ylim(1.2*ymin2,1.2*ymax2)
    plt.plot(y4_8_10,np.real(Ez4_8_10),'o',ms=4.5,color='skyblue')
    plt.subplot2grid((xn,yn),(1,1))
    plt.tick_params(axis='both',labelsize=0)
    plt.ylim(1.2*ymin2,1.2*ymax2)
    plt.plot(y4_10_1,np.real(Ez4_10_1),'o',ms=4.5,color='darkolivegreen')
    plt.subplot2grid((xn,yn),(2,1))
    plt.tick_params(axis='y',labelsize=0)
    plt.ylim(1.2*ymin2,1.2*ymax2)
    plt.plot(y4_10_10,np.real(Ez4_10_10),'o',ms=4.5,color='darkseagreen')

    plt.figure(22)
    plt.plot(e4_5[:400],s4_5[:400]/np.max(s4_5[:400]))
    plt.plot(e4_8[:200],s4_8[:200]/np.max(s4_5[:400]))


if set3 == 1:

    eField1_4 = np.asarray(pd.read_csv('p_ceria160nmAndAgTriangle78nm_3.1883mode_Efield.txt',sep='\t',dtype=str))
    x1_4_1,y1_4_1,z1_4_1,Ez1_4_1,Ey1_4_1,Ex1_4_1 = cutregion_AgTri(eField1_4,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x1_4_2,y1_4_2,z1_4_2,Ez1_4_2,Ey1_4_2,Ex1_4_2 = cutregion_AgTri(eField1_4,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x1_4 = np.r_[x1_4_1,x1_4_2]
    y1_4 = np.r_[y1_4_1,y1_4_2]
    z1_4 = np.r_[z1_4_1,z1_4_2]
    Ez1_4 = np.r_[Ez1_4_1,Ez1_4_2]
    Ex1_4 = np.r_[Ex1_4_1,Ex1_4_2]
    Ey1_4 = np.r_[Ey1_4_1,Ey1_4_2]
    energyDensity1_4 = np.sqrt((np.real(Ex1_4))**2+(np.real(Ey1_4))**2+(np.real(Ez1_4))**2)

    eField1_5 = np.asarray(pd.read_csv('p_ceria160nm_3.1883mode_electricField.txt',sep='\t',dtype=str))
    x1_5_1,y1_5_1,z1_5_1,Ez1_5_1,Ey1_5_1,Ex1_5_1 = cutregion_AgTri(eField1_5,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x1_5_2,y1_5_2,z1_5_2,Ez1_5_2,Ey1_5_2,Ex1_5_2 = cutregion_AgTri(eField1_5,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x1_5 = np.r_[x1_5_1,x1_5_2]
    y1_5 = np.r_[y1_5_1,y1_5_2]
    z1_5 = np.r_[z1_5_1,z1_5_2]
    Ez1_5 = np.r_[Ez1_5_1,Ez1_5_2]
    Ex1_5 = np.r_[Ex1_5_1,Ex1_5_2]
    Ey1_5 = np.r_[Ey1_5_1,Ey1_5_2]
    energyDensity1_5 = np.sqrt((np.real(Ex1_5))**2+(np.real(Ey1_5))**2+(np.real(Ez1_5))**2)

    eField1_6 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_3.1883mode_Efield.txt',sep='\t',dtype=str))
    x1_6_1,y1_6_1,z1_6_1,Ez1_6_1,Ey1_6_1,Ex1_6_1 = cutregion_AgTri(eField1_6,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x1_6_2,y1_6_2,z1_6_2,Ez1_6_2,Ey1_6_2,Ex1_6_2 = cutregion_AgTri(eField1_6,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x1_6 = np.r_[x1_6_1,x1_6_2]
    y1_6 = np.r_[y1_6_1,y1_6_2]
    z1_6 = np.r_[z1_6_1,z1_6_2]
    Ez1_6 = np.r_[Ez1_6_1,Ez1_6_2]
    Ex1_6 = np.r_[Ex1_6_1,Ex1_6_2]
    Ey1_6 = np.r_[Ey1_6_1,Ey1_6_2]
    energyDensity1_6 = np.sqrt((np.real(Ex1_6))**2+(np.real(Ey1_6))**2+(np.real(Ez1_6))**2)

    eField3_4 = np.asarray(pd.read_csv('p_ceria160nmAndAgTriangle78nm_4.9507mode_Efield.txt',sep='\t',dtype=str))
    x3_4_1,y3_4_1,z3_4_1,Ez3_4_1,Ey3_4_1,Ex3_4_1 = cutregion_AgTri(eField3_4,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x3_4_2,y3_4_2,z3_4_2,Ez3_4_2,Ey3_4_2,Ex3_4_2 = cutregion_AgTri(eField3_4,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x3_4 = np.r_[x3_4_1,x3_4_2]
    y3_4 = np.r_[y3_4_1,y3_4_2]
    z3_4 = np.r_[z3_4_1,z3_4_2]
    Ez3_4 = np.r_[Ez3_4_1,Ez3_4_2]
    Ex3_4 = np.r_[Ex3_4_1,Ex3_4_2]
    Ey3_4 = np.r_[Ey3_4_1,Ey3_4_2]
    energyDensity3_4 = np.sqrt((np.real(Ex3_4))**2+(np.real(Ey3_4))**2+(np.real(Ez3_4))**2)

    eField3_5 = np.asarray(pd.read_csv('p_ceria160nm_4.9507mode_electricField.txt',sep='\t',dtype=str))
    x3_5_1,y3_5_1,z3_5_1,Ez3_5_1,Ey3_5_1,Ex3_5_1 = cutregion_AgTri(eField3_5,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x3_5_2,y3_5_2,z3_5_2,Ez3_5_2,Ey3_5_2,Ex3_5_2 = cutregion_AgTri(eField3_5,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x3_5 = np.r_[x3_5_1,x3_5_2]
    y3_5 = np.r_[y3_5_1,y3_5_2]
    z3_5 = np.r_[z3_5_1,z3_5_2]
    Ez3_5 = np.r_[Ez3_5_1,Ez3_5_2]
    Ex3_5 = np.r_[Ex3_5_1,Ex3_5_2]
    Ey3_5 = np.r_[Ey3_5_1,Ey3_5_2]
    energyDensity3_5 = np.sqrt((np.real(Ex3_5))**2+(np.real(Ey3_5))**2+(np.real(Ez3_5))**2)

    eField3_6 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_4.9507mode_Efield.txt',sep='\t',dtype=str))
    x3_6_1,y3_6_1,z3_6_1,Ez3_6_1,Ey3_6_1,Ex3_6_1 = cutregion_AgTri(eField3_6,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x3_6_2,y3_6_2,z3_6_2,Ez3_6_2,Ey3_6_2,Ex3_6_2 = cutregion_AgTri(eField3_6,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x3_6 = np.r_[x3_6_1,x3_6_2]
    y3_6 = np.r_[y3_6_1,y3_6_2]
    z3_6 = np.r_[z3_6_1,z3_6_2]
    Ez3_6 = np.r_[Ez3_6_1,Ez3_6_2]
    Ex3_6 = np.r_[Ex3_6_1,Ex3_6_2]
    Ey3_6 = np.r_[Ey3_6_1,Ey3_6_2]
    energyDensity3_6 = np.sqrt((np.real(Ex3_6))**2+(np.real(Ey3_6))**2+(np.real(Ez3_6))**2)

    eField4_4 = np.asarray(pd.read_csv('p_ceria160nmAndAgTriangle78nm_5.1991mode_Efield.txt',sep='\t',dtype=str)) ### this is not a mistake -- this one is at 5.1991 instead. I should fix it... ###
    x4_4_1,y4_4_1,z4_4_1,Ez4_4_1,Ey4_4_1,Ex4_4_1 = cutregion_AgTri(eField4_4,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x4_4_2,y4_4_2,z4_4_2,Ez4_4_2,Ey4_4_2,Ex4_4_2 = cutregion_AgTri(eField4_4,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x4_4 = np.r_[x4_4_1,x4_4_2]
    y4_4 = np.r_[y4_4_1,y4_4_2]
    z4_4 = np.r_[z4_4_1,z4_4_2]
    Ez4_4 = np.r_[Ez4_4_1,Ez4_4_2]
    Ex4_4 = np.r_[Ex4_4_1,Ex4_4_2]
    Ey4_4 = np.r_[Ey4_4_1,Ey4_4_2]
    energyDensity4_4 = np.sqrt((np.real(Ex4_4))**2+(np.real(Ey4_4))**2+(np.real(Ez4_4))**2)

    eField4_5 = np.asarray(pd.read_csv('p_ceria160nm_5.1911mode_electricField.txt',sep='\t',dtype=str))
    x4_5_1,y4_5_1,z4_5_1,Ez4_5_1,Ey4_5_1,Ex4_5_1 = cutregion_AgTri(eField4_5,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x4_5_2,y4_5_2,z4_5_2,Ez4_5_2,Ey4_5_2,Ex4_5_2 = cutregion_AgTri(eField4_5,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x4_5 = np.r_[x4_5_1,x4_5_2]
    y4_5 = np.r_[y4_5_1,y4_5_2]
    z4_5 = np.r_[z4_5_1,z4_5_2]
    Ez4_5 = np.r_[Ez4_5_1,Ez4_5_2]
    Ex4_5 = np.r_[Ex4_5_1,Ex4_5_2]
    Ey4_5 = np.r_[Ey4_5_1,Ey4_5_2]
    energyDensity4_5 = np.sqrt((np.real(Ex4_5))**2+(np.real(Ey4_5))**2+(np.real(Ez4_5))**2)
    # x4_5,y4_5,z4_5,Ez4_5,Ey4_5,Ex4_5 = cutregion_AgTri(eField4_5,zlow=90.5e-9,zhigh=91.5e-9,xlow=0e-9,xhigh=250e-9,ylow=-85e-9,yhigh=85e-9)

    eField4_6 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_5.1991mode_Efield.txt',sep='\t',dtype=str))
    x4_6_1,y4_6_1,z4_6_1,Ez4_6_1,Ey4_6_1,Ex4_6_1 = cutregion_AgTri(eField4_6,zlow=85e-9,zhigh=115e-9,xlow=0e-9,xhigh=160e-9,ylow=-80e-9,yhigh=80e-9)
    x4_6_2,y4_6_2,z4_6_2,Ez4_6_2,Ey4_6_2,Ex4_6_2 = cutregion_AgTri(eField4_6,zlow=90.5e-9,zhigh=91.5e-9,xlow=160e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    x4_6 = np.r_[x4_6_1,x4_6_2]
    y4_6 = np.r_[y4_6_1,y4_6_2]
    z4_6 = np.r_[z4_6_1,z4_6_2]
    Ez4_6 = np.r_[Ez4_6_1,Ez4_6_2]
    Ex4_6 = np.r_[Ex4_6_1,Ex4_6_2]
    Ey4_6 = np.r_[Ey4_6_1,Ey4_6_2]
    energyDensity4_6 = np.sqrt((np.real(Ex4_6))**2+(np.real(Ey4_6))**2+(np.real(Ez4_6))**2)
    # x4_6,y4_6,z4_6,Ez4_6,Ey4_6,Ex4_6 = cutregion_AgTri(eField4_6,zlow=60e-9,zhigh=100e-9,xlow=0e-9,xhigh=250e-9,ylow=-40e-9,yhigh=40e-9)
    # energyDensity4_6 = np.sqrt((np.real(Ex4_6))**2+(np.real(Ey4_6))**2+(np.real(Ez4_6))**2)

    eels4 = np.asarray(pd.read_csv('p_NOceria_AgTriangle78nm_EELspectrum.txt',sep='\t'))
    e4 = eels4[:,0]/1.602e-19
    s4 = eels4[:,1]
    eels5 = np.asarray(pd.read_csv('p_ceria160nm_EELspectrum.txt',sep='\t'))
    e5 = eels5[:,0]/1.602e-19
    s5 = eels5[:,1]
    eels6 = np.asarray(pd.read_csv('p_ceria160nmAndAgTriangle78nm_EELspectrum.txt',sep='\t'))
    e6 = eels6[:,0]/1.602e-19
    s6 = eels6[:,1]


def plotset3():
    plt.rcParams.update({'font.size': 9})
    ymin1 = np.min(np.real(Ez1_6))
    ymin2 = np.min(np.real(Ez4_4))
    ymax1 = np.max(np.real(Ez1_6))
    ymax2 = np.max(np.real(Ez4_4))
    xn = 3
    yn = 3
    plt.figure(31)
    plt.subplot2grid((xn,yn),(0,0))
    # plt.subplots_adjust(wspace=0.4)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.title('Ag Plasmonic Mode, 2.0 eV', fontsize=18)
    plt.plot(x1_6,np.real(Ez1_6),'o',ms=4.5,color='coral')
    plt.ylim(-35,35)
    plt.xlim(0,250)

    plt.subplot2grid((xn,yn),(0,1))
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.title('Ag Plasmonic Mode, 2.9 eV', fontsize=18)
    plt.plot(x3_6,np.real(Ez3_6),'o',ms=4.5,color='coral')
    plt.ylim(-35,35)
    plt.xlim(0,250)

    plt.subplot2grid((xn,yn),(0,2))
    plt.ylim(-35,35)
    plt.xlim(0,250)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.title('Ceria Cavity Mode, 3.2 eV', fontsize=18)
    plt.plot(x4_6,np.real(Ez4_6),'o',ms=4.5,color='coral')

    plt.subplot2grid((xn,yn),(1,0))
    plt.ylim(-35,35)
    plt.xlim(0,250)
    plt.ylabel(r'$E_z$ (mV/m)',fontsize=18)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.plot(x1_5,np.real(Ez1_5),'o',ms=4.5,color='darkviolet')

    plt.subplot2grid((xn,yn),(1,1))
    plt.ylim(-35,35)
    plt.xlim(0,250)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.plot(x3_5,np.real(Ez3_5),'o',ms=4.5,color='darkviolet')

    plt.subplot2grid((xn,yn),(1,2))
    plt.ylim(-35,35)
    plt.xlim(0,250)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.plot(x4_5,np.real(Ez4_5),'o',ms=4.5,color='darkviolet')

    plt.subplot2grid((xn,yn),(2,0))
    plt.ylim(-35,35)
    plt.xlim(0,250)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.tick_params(axis='x',labelsize=12)
    plt.plot(x1_4,np.real(Ez1_4),'o',ms=4.5,color='deepskyblue')

    plt.subplot2grid((xn,yn),(2,1))
    plt.ylim(-35,35)
    plt.xlim(0,250)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.tick_params(axis='x',labelsize=12)
    plt.plot(x3_4,np.real(Ez3_4),'o',ms=4.5,color='deepskyblue')

    plt.subplot2grid((xn,yn),(2,2))
    plt.ylim(-35,35)
    plt.xlim(0,250)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.tick_params(axis='x',labelsize=12)
    plt.plot(x4_4,np.real(Ez4_4),'o',ms=4.5,color='deepskyblue')

    plt.figure(32)
    plt.subplot2grid((xn,yn),(0,0))
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.title('Ag Plasmonic Mode, 2.0 eV', fontsize=18)
    plt.plot(x1_6,np.real(energyDensity1_6)**2,'o',ms=4.5,color='coral')
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.yscale('log')

    plt.subplot2grid((xn,yn),(0,1))
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.title('Ag Plasmonic Mode, 2.9 eV', fontsize=18)
    plt.plot(x3_6,np.real(energyDensity3_6)**2,'o',ms=4.5,color='coral')
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.yscale('log')

    plt.subplot2grid((xn,yn),(0,2))
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.title('Ceria Cavity Mode, 3.2 eV', fontsize=18)
    plt.plot(x4_6,np.real(energyDensity4_6)**2,'o',ms=4.5,color='coral')
    plt.yscale('log')

    plt.subplot2grid((xn,yn),(1,0))
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.ylabel(r'$|\mathbf{E}|^2$ (mV$^2$/m$^2$)',fontsize=18)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.plot(x1_5,np.real(energyDensity1_5)**2,'o',ms=4.5,color='darkviolet')
    plt.yscale('log')

    plt.subplot2grid((xn,yn),(1,1))
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.plot(x3_5,np.real(energyDensity3_5)**2,'o',ms=4.5,color='darkviolet')
    plt.yscale('log')

    plt.subplot2grid((xn,yn),(1,2))
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.tick_params(axis='x',labelsize=0)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.plot(x4_5,np.real(energyDensity4_5)**2,'o',ms=4.5,color='darkviolet')
    plt.yscale('log')

    plt.subplot2grid((xn,yn),(2,0))
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.tick_params(axis='x',labelsize=12)
    plt.plot(x1_4,np.real(energyDensity1_4)**2,'o',ms=4.5,color='deepskyblue')
    plt.yscale('log')

    plt.subplot2grid((xn,yn),(2,1))
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.tick_params(axis='x',labelsize=12)
    plt.plot(x3_4,np.real(energyDensity3_4)**2,'o',ms=4.5,color='deepskyblue')
    plt.yscale('log')

    plt.subplot2grid((xn,yn),(2,2))
    plt.ylim(1e-4,4000)
    plt.xlim(0,250)
    plt.tick_params(axis='x',which='minor',length=0)
    plt.tick_params(axis='x',labelsize=12)
    plt.plot(x4_4,np.real(energyDensity4_4)**2,'o',ms=4.5,color='deepskyblue')
    plt.yscale('log')

    plt.figure(33)
    plt.plot(e4,s4/np.max(s5),color='coral',linewidth=3)
    plt.plot(e5,s5/np.max(s5),color='darkviolet',linewidth=3)
    plt.plot(e6,s6/np.max(s5),color='deepskyblue',linewidth=3)

def plotset4():

    plt.figure(14)
    plt.suptitle('Ag Triangle, Ceria Cube, and Both for 2 eV Mode')
    plt.subplot(321)
    plt.plot(x1_5,np.real(Ez1_5),'o',ms=3.5)
    plt.title('Ceria Cube Only')
    plt.ylabel('Ez Field (mV/m)')
    plt.ylim(-ylim,ylim)
    plt.subplot(322)
    plt.plot(y1_5,np.real(Ez1_5),'o',ms=3.5)
    plt.title('Ceria Cube Only')
    plt.ylim(-ylim,ylim)
    plt.subplot(323)
    plt.plot(x1_4,np.real(Ez1_4),'o',ms=3.5)
    plt.ylabel('Ez Field (mV/m)')
    plt.title('Ceria Cube with Ag Triangle')
    plt.ylim(-ylim,ylim)
    plt.subplot(324)
    plt.plot(y1_4,np.real(Ez1_4),'o',ms=3.5)
    plt.title('Ceria Cube with Ag Triangle')
    plt.ylim(-ylim,ylim)
    plt.subplot(325)
    plt.plot(x1_6,np.real(Ez1_6),'o',ms=3.5)
    plt.ylabel('Ez Field (mV/m)')
    plt.title('Ag Triangle, 164 nm from e-beam')
    plt.ylim(-ylim,ylim)
    plt.subplot(326)
    plt.plot(y1_6,np.real(Ez1_6),'o',ms=3.5)
    plt.title('Ag Triangle, 164 nm from e-beam')
    plt.ylim(-ylim,ylim)

def plotset5():
    ylim1 = np.max(abs(np.real(Ez1_2)))

    plt.figure(15)
    plt.suptitle('Ag Triangle, Ceria Cube, and Both for 2 eV Mode')
    plt.subplot(321)
    plt.plot(x1_5,np.real(energyDensity1_5),'o',ms=3.5)
    plt.title('Ceria Cube Only')
    plt.ylabel('Energy Density')
    plt.ylim(-ylim1,ylim1)
    plt.subplot(322)
    plt.plot(y1_5,np.real(energyDensity1_5),'o',ms=3.5)
    plt.title('Ceria Cube Only')
    plt.ylim(-ylim1,ylim1)
    plt.subplot(323)
    plt.plot(x1_4,np.real(energyDensity1_4),'o',ms=3.5)
    plt.ylabel('Energy Density')
    plt.title('Ceria Cube with Ag Triangle')
    plt.ylim(-ylim1,ylim1)
    plt.subplot(324)
    plt.plot(y1_4,np.real(energyDensity1_4),'o',ms=3.5)
    plt.title('Ceria Cube with Ag Triangle')
    plt.ylim(-ylim1,ylim1)
    plt.subplot(325)
    plt.plot(x1_6,np.real(energyDensity1_6),'o',ms=3.5)
    plt.ylabel('Energy Density')
    plt.title('Ag Triangle, 164 nm from e-beam')
    plt.ylim(-ylim1,ylim1)
    plt.subplot(326)
    plt.plot(y1_6,np.real(energyDensity1_6),'o',ms=3.5)
    plt.title('Ag Triangle, 164 nm from e-beam')
    plt.ylim(-ylim1,ylim1)


if set0 == 1:
    plotset0()
if set1 == 1:
    plotset1()
if set2 == 1:
    plotset2()
if set3 == 1:
    plotset3()
if set4 == 1:
    plotset4()
if set5 == 1:
    plotset5()


plt.show()
