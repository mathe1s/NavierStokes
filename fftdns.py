#! /usr/bin/python 
#FFT of DNS code to find the 
#spatial eigenvalues 
# Create by: Jhonatan
# date: 06-03-2018
#Module for symbolic albegra 
import numpy as np
import math  as mt
#from math import exp, expm1
import matplotlib.pyplot as plt
#from sci to make fft
import scipy.fftpack
################################
# Letura de aquivos para obter as transformadas

file1='./ekx1_' 
f = open('taxafft.dat','w+')

# number of frequency files 
N       = 11
# frequency  
wa      = range(0,N,1)

pi      = mt.pi
################################
#step of point in time 
dt      = 0.01
################################
#step of point in x direction 
dx      = 0.05
################################
x_ini   = 10.0
################################

#mean of the alpha using fortwar differences
alpham1  =   np.zeros(N)
#mean of the alpha using central differences
alpham3  =   np.zeros(N)
#mean of the alpha using four orde scheme 
alpham4  =   np.zeros(N)

for i in range(2,N):
#for i in range(1,2):

    #To read the first line of the files, which is the frequency 
    #infile    = open('./energias/ekx1_%d.dat'%i,'rU')
    infile    = open('%s%d.dat'%(file1,i),'rU')
    w         = infile.readline()
    wa[i]     = float(w) 
    print(wa[i])

    #Number of columns of the Matrix,Matrix[t ekx1 ekx1+dx ekx1+2dx ekx1+3dx...........]
    #################################
    #M1  =   np.loadtxt('./energias/ekx1_%d.dat'%i,unpack=True,skiprows=1) 
    M1  =   np.loadtxt('%s%d.dat'%(file1,i),unpack=True,skiprows=1) 
    N1  =   M1.shape[1]
    N2  =   M1.shape[0]
    x1  =   M1[0:N1][0]
    f1  =   M1[0:N1][1]
    k1  =   range(0,N1,1)

    x   =   range(0,N2,1)
    x1  =   x_ini+np.asfarray(x,dtype='float')*dx
    #print(x1)

    H        =   np.zeros((N1,N2))
    alpha1   =   np.zeros((N,N2))
    alpha3   =   np.zeros((N,N2))
    alpha4   =   np.zeros((N,N2))

    ######################################
    #fft of the temporal Marix M
    for l in range(1,N2):
        H[:,l]   =   scipy.fftpack.fft(M1[0:N1][l])
        H        =   H.real


    #####################################
    #Created the frequency vetor 
    dw   =    pi/(dt*N1)
    #dw   =    1/(dt*N1*pi)
    print(dw,'dw')
    k    =    range(0,N1,1)
    kf   =    np.asfarray(k,dtype='float')
    kf   =    kf*dw

    #plt.plot(kf,H[:,3],'*b')
    #plt.plot(kf,H[:,4],'-r')
    #plt.grid() 
    #plt.show()

    
    ######################################

    #for l in range(1,N2-1):
    #for l in range(2,N2-4):

    #for l in range(2,N-2-4): 
    for l in range(5,6): 

            for j in range(0,N1): 

                #if (kf[j]>(wa[i]-dw-0.01)) and (kf[j]<(wa[i]+dw-0.001)):
                diff=np.abs(wa[i]-kf[j])

                if (diff<dw):

                    anterior    =np.abs(wa[i]-kf[j-1])
                    posterior   =np.abs(wa[i]-kf[j+1])

                    if ((anterior<diff)or(posterior<diff)):

                        if (anterior<posterior):
                      
                            kk=j-1

                        else:
                            kk=j+1


                    else: 

                       kk=j 

           

            print(l,kf[kk],wa[i],'kf',H[kk][l+1],H[kk][l],H[kk][l+1]/H[kk][l]) 

            #gain         = (H[kk][l+1]+H[0][l+1])/(H[kk][l]+H[0][l])
            gain         = (H[0][l+1])/(H[0][l])
            alpha1[i][l] = -1.0/2.0*(mt.log(float(gain)))*1.0/dx
            #gain3           = (H[kk][l+1]+H[0][l+1])/(H[kk][l-1]+H[0][l-1])
            gain3           = (H[0][l+1])/(H[0][l-1])
            alpha3[i][l] = -1.0/2.0*(mt.log(float(gain3)))*1.0/(2.0*dx)
            #print(alpha3[l],'alpha3') 
            #alphai3s[i] = -1/2*float(math.log(H3[j]/H1[j]))/(2*dx)
            #gain4        = H[j][l+3]**8*H[j][l]/(H[j][l+4]*H[j][l+1]**8)
            gain4        = H[0][l+3]**8*H[0][l]/(H[0][l+4]*H[0][l+1]**8)
            alpha4[i][l] = -1.0/2.0*(mt.log(float(gain4)))*1.0/(12.0*dx)
            #print(alpha4[l],'alpha4') 




            #print(alpha1[i][l],'alpha1') 

    
    #alpham1[i]=np.mean(alpha1[i][2:N2-4],axis=0)
    #alpham3[i]=np.mean(alpha3[i][2:N2-4],axis=0)
    #alpham4[i]=np.mean(alpha4[i][2:N2-4],axis=0)

    alpham1[i]=alpha1[i][l]
    alpham3[i]=alpha3[i][l]
    alpham4[i]=alpha4[i][l]



print(i,alpham1[i],alpham3[i],alpham4,'media') 
plt.plot(wa,-alpham1,'*b',label='alpha1, Forwart')
plt.plot(wa,-alpham3,'*y',label='alpha3, central')
plt.plot(wa,-alpham4,'*m',label='alpha4, Four Order')

for i in range(0,N):

    f.write("%f\t%f\t%f\t%f\n"%(wa[i],-alpham1[i],-alpham3[i],-alpham4[i]))

f.close

#plt.plot(x1[2:N2-4],alpha1[2:N2-4],'-*b',label='alpha1,forwart')
#plt.plot(x1[2:N2-4],alpha3[2:N2-4],'-+r',label='alpha3,central')
#plt.plot(x1[2:N2-4],alpha4[2:N2-4],'-.m',label='alpha4,4 order')
plt.grid()
#plt.axis((0.0,2.0,0.0,0.6))
#plt.legend(handles=[line_up, line_down])
plt.legend()
plt.show()

#plt.plot(x1[2:N2-4],alpha1[2:N2-4],'-*b',label='alpha1,forwart')
#plt.plot(x1[2:N2-4],alpha3[2:N2-4],'-+r',label='alpha3,central')
#plt.plot(x1[2:N2-4],alpha4[2:N2-4],'-.m',label='alpha4,4 order')
#plt.grid()
##plt.axis((0.01, 10, 0, 1000))
##plt.legend(handles=[line_up, line_down])
#plt.legend()
#plt.show()


