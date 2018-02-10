#!/usr/bin/env python
from numpy import *
import numpy as np
import sys

#NOW FOLLOWING REV.MOD.PHYS.81.109(2009)
Nkx=51
Nky=51

ts=1.0
tsp=0.0
phi=0.0
Mh=0.0

a=1.0
a0=sqrt(3.0)
pauli_0 = array([[1, 0],[ 0, 1]])
pauli_x = array([[0, 1],[ 1, 0]])
pauli_y = array([[0, -1j],[1j, 0]])
pauli_z = array([[1, 0],[0, -1]])


#Lattice basis (a=1; a0=sqrt3*a) is:
#\a_1 = a0 [ sqrt3/2 , 1/2 ]
#\a_2 = a0 [ sqrt3/2 ,-1/2 ]
#
#
#nearest neighbor: A-->B, B-->A
d1= a*np.array([  1.0/2.0 , sqrt(3.0)/2.0 ])
d2= a*np.array([  1.0/2.0 ,-sqrt(3.0)/2.0 ])
d3= a*np.array([ -1.0     , 0.0           ])
#
#next nearest-neighbor displacements: A-->A, B-->B \== \nu_1,\nu_2, \nu_3=\nu_1-\nu_2
a1=a0*np.array([ sqrt(3.0)/2.0, 1.0/2.0])
a2=a0*np.array([ sqrt(3.0)/2.0,-1.0/2.0])
a3=a2-a1


def hk_graphene_model(kpoint):
    kdotd = np.zeros(3)
    kdota = np.zeros(3)
    kdotd[0] = np.inner(kpoint,d1)
    kdotd[1] = np.inner(kpoint,d2)
    kdotd[2] = np.inner(kpoint,d3)
    kdota[0] = np.inner(kpoint,a1)
    kdota[1] = np.inner(kpoint,a2)
    kdota[2] = np.inner(kpoint,a3)
    h0 = 2*tsp*cos(phi)*np.sum( cos(kdota[:]) )
    hx =-ts*np.sum( cos(kdotd[:]) )
    hy =-ts*np.sum( sin(kdotd[:]) )
    hz = 2*tsp*sin(phi)*np.sum( sin(kdota[:]) ) + Mh
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
    # hk=np.zeros((2,2),dtype=complex)
    # hk[0,1] = -ts*sum( exp(-kdotd[:]*1j) )
    # hk[1,0] = -ts*sum( exp(kdotd[:]*1j) )
    return hk

#RECIPROCAL LATTICE VECTORS:
blen=4*pi/3
b1=blen*np.array([ 1/2.0 ,  sqrt(3.0)/2.0 ])
b2=blen*np.array([ 1/2.0 , -sqrt(3.0)/2.0 ])

K = np.array([2*pi/3,2*pi/3/sqrt(3)])
Kp= np.array([2*pi/3,-2*pi/3/sqrt(3)])



p = np.zeros((2,Nkx*Nky))
hk = np.zeros((2,2,Nkx*Nky),dtype=complex)

#SIZE OF THE K-GRID:
KYlen=4*pi/3/sqrt(3)
KXlen=2*pi/3
i=-1
file1=open("BZhex.dat","w")
for iy in range(Nky):
    ky = KYlen*(-1.0 + 2.0/Nky*iy)
    for ix in range(Nkx):
        kx = KXlen*(-1.0 + 2.0/Nkx*ix)
        #
        kvec = np.array( [kx,ky] )
        klimit1 = np.inner(kvec,b1/blen)
        klimit2 = np.inner(kvec,b2/blen)
        klimit3 = np.inner(kvec,np.array([1,0]))
        Kbound = blen/2
        if ( abs(klimit1)<=Kbound and abs(klimit2)<=Kbound and abs(klimit3)<=Kbound ):
            i+=1
            p[:,i] = kvec
            print >> file1, kvec[0],kvec[1]
            hk[:,:,i] = hk_graphene_model(p[:,i])


file2=open("BZpoints.dat","w")
print >> file2, 0,0
print >> file2, b1[0]/blen,b1[1]/blen
print >> file2
print >> file2, 0,0
print >> file2, b2[0]/blen,b2[1]/blen
print >> file2
print >> file2, 0,0
print >> file2, 1,0
print >> file2 
print >> file2, K[0],K[1]
print >> file2 
print >> file2, Kp[0],Kp[1]


# area_hex = area_rect * #points_inside/#points_generated
area_rect = (2*KXlen)*(2*KYlen)
points_in=i
points_tot=Nkx*Nky
points_ratio=1.0*points_in/points_tot
area_hex = area_rect*points_ratio
print "Points ratio:", points_ratio
print "Area HexBZ, Err(%)=",area_hex,str(abs(area_hex-8.0*pi**2/3.0/sqrt(3.0))*100)+"%"

print sum(hk,axis=2)/points_in


file1=open("BZrombo.dat","w")
i=-1
for iy in range(Nky):
    ky = float(iy)/Nky
    for ix in range(Nkx):
        kx = float(ix)/Nkx
        kvec = np.array( [kx,ky] )
        i+=1
        p[:,i] = kx*b1 + ky*b2
        print >> file1, p[0,i],p[1,i]
        hk[:,:,i] = hk_graphene_model(p[:,i])

print sum(hk,axis=2)/Nkx/Nky
