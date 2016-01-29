#Program which solves the one-particle Schrodinger equation 
#for a potential specified in function
#potential(). This example is for the harmonic oscillator in 3d

from  matplotlib import pyplot as plt
import numpy as np
#Function for initialization of parameters
def initialize():
    RMin = 0.0
    RMax = 10.0
    lOrbital = 0
    Dim = 400
    return RMin, RMax, lOrbital, Dim
# Here we set up the harmonic oscillator potential
def potential(r):
    return r*r

#set up Woods-Saxon potential
def WSpotential(r):
    v0 = 50
    A = 100
    a = 0.5
    r0 = 1.25
    R = r0 * (A**(0.3333))
    return -v0 / (1.0 + np.exp((r-R) / a))

#Get the boundary, orbital momentum and number of integration points
RMin, RMax, lOrbital, Dim = initialize()

#Initialize constants
#renormalize energy into base unit of hbar^2 / 2m(nucleons)*A, and using length unit in fm, say we're calculating in a system with A = 100 nucleons
alp = 197.0**2 / (938.0 * 2 * 100)
Step    = RMax/(Dim+1)
DiagConst = 2.0 * alp / (Step*Step)
NondiagConst =  -1.0 * alp / (Step*Step)
OrbitalFactor = alp * lOrbital * (lOrbital + 1.0)

#Calculate array of potential values
v = np.zeros(Dim)
r = np.linspace(RMin,RMax,Dim)
for i in xrange(Dim):
    r[i] = RMin + (i+1) * Step;
    v[i] = WSpotential(r[i]) + OrbitalFactor/(r[i]*r[i]);

#Setting up tridiagonal matrix and find eigenvectors and eigenvalues
Hamiltonian = np.zeros((Dim,Dim))
Hamiltonian[0,0] = DiagConst + v[0];
Hamiltonian[0,1] = NondiagConst;
for i in xrange(1,Dim-1):
    Hamiltonian[i,i-1]  = NondiagConst;
    Hamiltonian[i,i]    = DiagConst + v[i];
    Hamiltonian[i,i+1]  = NondiagConst;
Hamiltonian[Dim-1,Dim-2] = NondiagConst;
Hamiltonian[Dim-1,Dim-1] = DiagConst + v[Dim-1];
# diagonalize and obtain eigenvalues, not necessarily sorted
EigValues, EigVectors = np.linalg.eig(Hamiltonian)
# sort eigenvectors and eigenvalues
permute = EigValues.argsort()
EigValues = EigValues[permute]
EigVectors = EigVectors[:,permute]
# now plot the results for the three lowest lying eigenstates
for i in xrange(5):
    print EigValues[i]

#check normalization
sum1 = 0.0
sum2 = 0.0
sum3 = 0.0
count = 0
for i in xrange(Dim):
    sum1 += EigVectors[i,1]**2
    sum2 += EigVectors[i,2]**2
    sum3 += EigVectors[i,3]**2
    count += 1
    print sum1, sum2, sum3, count

FirstEigvector = EigVectors[:,0]
SecondEigvector = EigVectors[:,1]
ThirdEigvector = EigVectors[:,2]
plt.plot(r, FirstEigvector**2,'b-', label = 'n = 0, l = 0')
plt.plot(r, SecondEigvector**2,'g-', label = 'n = 1, l = 0')
plt.plot(r, ThirdEigvector**2,'r-', label = 'n = 2, l = 0')

plt.axis([0,4.6,0.0, 0.025])
plt.grid([0.0])
plt.xlabel(r'$r \left(fm\right)$')
plt.ylabel(r'Radial probability $r^2|R(r)|^2$')
plt.title(r'Radial probability distributions for three lowest-lying states l = 0, Woods-Saxon')
#plt.ylabel(r'Radial wfn amplitude $r|R(r)|$')
#plt.title('Eigenfunctions vs Radius for three lowest-lying states for l = 0')
plt.savefig('eigenvector.pdf')
plt.savefig('eigenvector.png')
plt.legend()
plt.show()
