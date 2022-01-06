#!/usr/bin/env python


import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c

def B(ku,m,k,r):
    return(k**2.0*(n.cross(n.cross(ku,m),ku))*n.exp(1j*k*r)/r+(3*ku*n.dot(ku,m)-m)*(1.0/(r**3.0) - 1j*k/(r**2.0))*n.exp(1j*k*r))

    
def dipole_B(k=2.0*n.pi,r=10**(n.linspace(-3,3,num=1000))):
    m=n.array([0,1.0,0])
    x=n.array([1.0,0,0])
    y=n.array([0.0, 1.0, 0])
    z=n.array([0,0,1.0])
    Bx=n.zeros(len(r))
    By=n.zeros(len(r))
    Bz=n.zeros(len(r))
    for i in range(len(r)):
        Bnow=B(x,m,k,r[i])
        Bx[i]=n.sqrt(n.dot(Bnow,n.conj(Bnow)))
        Bnow=B(z,m,k,r[i])
        Bz[i]=n.sqrt(n.dot(Bnow,n.conj(Bnow)))
        Bnow=B(y,m,k,r[i])
        By[i]=n.sqrt(n.dot(Bnow,n.conj(Bnow)))
        
    plt.loglog(r,n.abs(Bz),label="z")
    plt.loglog(r,n.abs(Bx),label="x")
    plt.loglog(r,n.abs(By),label="y")
    plt.title("$m=[0,1,0]^T$")
    plt.legend()
    plt.xlabel("Range (wavelengths)")
    plt.ylabel("Magnetic flux density $|B|$")
    plt.show()



def dipole_B_2d(k=2.0*n.pi,r=10**(n.linspace(-3,3,num=100))):
    m=n.array([0,1.0,0])

    Np=len(r)
    Bmag=n.zeros([Np,Np])
    
    for i in range(Np):
        for j in range(Np):
            ku=n.array([r[i],r[j],0.0])
            R=n.sqrt(n.dot(ku,ku))
            ku=ku/R
            Bnow=B(ku,m,k,R)
            
            Bmag[i,j]=n.sqrt(n.dot(Bnow,n.conj(Bnow)))

    dlog=n.log10(r)
    plt.pcolormesh(dlog,dlog,n.log10(Bmag))
    plt.xlabel("$\log_{10}(y/\lambda)$")
    plt.ylabel("$\log_{10}(x/\lambda)$")    
    cb=plt.colorbar()
    cb.set_label("Magnetic flux density $\log_{10}(|B|)$")
    plt.show()

    

def nearfield_B(m,ks,r=1.0):
    # Jackson, Classical Electrodynamics, 2nd Ed. p. 398
    B=n.zeros([N_angs,3])
    dots=3.0*(k[:,0]*m[0]+k[:,1]*m[1]+k[:,2]*m[2])
    B[:,0]=k[:,0]*dots
    B[:,1]=k[:,1]*dots
    B[:,2]=k[:,2]*dots
    
    B[:,0]=(B[:,0]-m[0])/r**3.0
    B[:,1]=(B[:,1]-m[1])/r**3.0
    B[:,2]=(B[:,2]-m[2])/r**3.0

    Bmag=n.sqrt(B[:,0]**2.0+B[:,1]**2.0+B[:,2]**2.0)
    return(B,Bmag)

dipole_B()
dipole_B_2d()
# zero angle aligned with magnetic dipole moment
angs=n.linspace(-n.pi,n.pi,num=1000)
N_angs=len(angs)
# magnetic dipole moment
m = n.array([0,0,1])

# y,z plane
k = n.zeros([N_angs,3])
k[:,0]=0.0
k[:,1]=n.sin(angs)
k[:,2]=n.cos(angs)
B,Bmag=nearfield_B(m,k)

# x,y plane
k = n.zeros([N_angs,3])
k[:,0]=n.sin(angs)
k[:,1]=n.cos(angs)
k[:,2]=0
Bxy,Bmagxy=nearfield_B(m,k)

# x,z plane
k = n.zeros([N_angs,3])
k[:,0]=n.sin(angs)
k[:,1]=0
k[:,2]=n.cos(angs)
Bxz,Bmagxz=nearfield_B(m,k)

plt.plot(angs,Bmag,label="y,z & x,z plane")
plt.plot(angs,Bmagxy,label="x,y plane")
#plt.plot(angs,Bmagxz,label="x,z plane")
plt.xlabel("Angle (rad)")
plt.title("$m=[0,0,1]^T$")
plt.ylabel("Magnetic flux density $|B|$")
plt.legend()
plt.show()



