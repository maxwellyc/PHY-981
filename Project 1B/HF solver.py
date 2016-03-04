import numpy as npy
from matplotlib import pyplot as plt
from copy import copy, deepcopy


#Import input data
f0 = open(str('twobody.dat'))
f1 = open(str('spdata.dat'))
lines0 = f0.readlines()
lines1 = f1.readlines()
interaction = 1 # 1 - turn on interaction, all other value turns it off
NStates = 40
threshold = 1E-6
#Iteration index
IterMax = 100
partNum = 8


type = 0    #type 1- neutron, -1- proton, 0- both, 3- Oxygen17 ?, 4- Fluorine17 ?, 5- Silicon 28, for 5 maybe we know the next s.p. energy for 0d5/2 orbit
two_iso = type

if (type == 0):
    NStates = 80
    partNum = 16  #if we just use 17, since neutron energy is lower, we'll just get binding energy of O-17
    #how to get F-17 then?, reconfigure density to sum up: 1,1,1 ... 1,(0,0,0,0,0,0),1 ?
if (type ==3 or type == 4):
    NStates = 80
    two_iso = 0
    partNum = 17

if (type == 5):
    NStates = 80
    partNum = 28
    two_iso = 0


V0 = {} #dictionary that stores interaction M.E., with key as '(a,b,c,d)'
QNum = {} #dictionary that stores all orbit number in spdata.dat
PNum = {} #Inverse structure as QNum, where keyword becomes value, value becomes keyword
Ind = {}
RInd = {}


#initialization
E = npy.zeros(NStates)
ePrev = npy.zeros(NStates)
diff = npy.zeros(NStates)
C = npy.zeros([NStates,NStates])
Halm = npy.zeros([NStates,NStates])



#Import interaction matrix element (M.E.)
for line in lines0:
    rs0 = line.split()   #raw string
    try:
        a = int(rs0[0])
        b = int(rs0[1])
        c = int(rs0[2])
        d = int(rs0[3])
        V0[(a,b,c,d)] = float(rs0[4])
    except (ValueError, IndexError):
        continue
f0.close()
count = 0
#Import quantum numbers matrix
for line in lines1:
    
    rs1 = line.split()   #raw string
    try:
        OrbNum = int(rs1[2])   #orbit number which corresponds to indices in twobody.dat
        n = int(rs1[3])
        l = int(rs1[4])
        two_j = int(rs1[5])
        two_m = int(rs1[6])
        two_tz = int(rs1[7])
        #Later We'll only use two_tz = 1 for Neutrons
        QNum[(n,l,two_j,two_m,two_tz)] = OrbNum
        PNum[(OrbNum)] = str(n) + ',' + str(l) + ',' + str(two_j) + ',' + str(two_m) + ',' + str(two_tz)
        if (type != 0):
            if (two_tz == two_iso):
                Ind[(OrbNum)] = count
                RInd[(count)] = OrbNum
                count += 1
        else:
            Ind[(OrbNum)] = count
            RInd[(count)] = OrbNum
            count += 1
    except (ValueError, IndexError):
        continue
f1.close()


#define a Kronecker delta function
def KDelta(n1,n3):
    if (n1 == n3):
        return 1.0
    else:
        return 0.0

#define density M.E.


#define function that converts quantum numbers into orbital numbers corresponding in spdata.dat
#def OrbNum(n,l,j,m,QNum):
#    for key in QNum.keys:
#        line = QNum.split()
#        if ( n == line[0] and l == line[1] and j == line[2] and m == line[3]):
#            return key[0]

def rho(gamma, delta, C):
    rho0 = 0.0
    for N in xrange(0,partNum):
        rho0 += C[gamma][N] * npy.conjugate(C[delta][N])
    if (type == 3):
        rho0 += C[gamma][16] * npy.conjugate(C[delta][16])
    elif (type == 4):
        rho0 += C[gamma][22] * npy.conjugate(C[delta][22])
    else:
        rho0 += 0

    return rho0


    
#units are in MeV, harmonic oscillator basis of hbar*omega = 10 MeV
def spHF(a, b, C, QNum):     #need to convert alpha, beta to explicit l,j,n1,n3,
    alpha = RInd[(a)]   #now alpha, beta correspond to Neutron Orbital Number
    beta = RInd[(b)]
    lineA = PNum[(alpha)].split(',')
    n1 = int(lineA[0])
    l1 = int(lineA[1])
    h0 = KDelta(a,b) * (2*n1 + l1 + 1.5) * 10
    h1 = 0.0
    N1 = alpha
    N3 = beta
    for i in xrange(0,NStates):
        N2 = RInd[(i)]
        for j in xrange(0,NStates):
            N4 = RInd[(j)]
            if ( (N1,N2,N3,N4) in V0 ):
                h1 += V0[(N1,N2,N3,N4)] * rho(i, j, C)
    if ( interaction == 1):
        return h0 + h1
    else:
        return h0

#making C0 as an Identity
for i in xrange(0,NStates):
    C[i][i] = 1.0

for i in xrange(0,NStates):
    a = RInd[(i)]
    line = PNum[(a)].split(',')
    n = int(line[0])
    l = int(line[1])
    Halm[i][i] = (2*n + l + 1.5) * 10



Iter = 0
while (Iter < IterMax):
    ePrev = deepcopy(E)
    for i in xrange(0,NStates):
        for j in xrange(0,NStates):
            Halm[i][j] = spHF(i, j, C, QNum)
    EigValues, EigVectors = npy.linalg.eig(Halm)
    permute = EigValues.argsort()
    EigValues = EigValues[permute]
    EigVectors = EigVectors[:,permute]
    #C = EigVectors.transpose()
    C = deepcopy(EigVectors)
    E = deepcopy(EigValues)
    for i in xrange(0,NStates):
        diff[i] = abs(E[i] - ePrev[i])
    diff.sort()
    Iter += 1
    print 'Iteration: ' + str(Iter) + ', diff = ' + str(diff[NStates-1]) + ' MeV\n'
    print E
    if ( diff[NStates-1] < threshold ):
        break








