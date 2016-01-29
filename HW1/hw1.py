
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#import scipy.interpolate



def PES(InFile0):
 f0 = open(str(InFile0))
 lines0 = f0.readlines()
 
 #matrix for raw data
 BEmat=[]
 #array to store plotting data, Neutron number, separartion energy (SE), SE from Liquid Drop Model, format : [neutron] ; [SE] ; [LDM SE]
 O_N = []
 O_SE = []
 O_LDM = []
 Ca_N = []
 Ca_SE = []
 Ca_LDM = []
 Ni_N = []
 Ni_SE = []
 Ni_LDM = []
 Sn_N = []
 Sn_SE = []
 Sn_LDM = []
 Pb_N = []
 Pb_SE = []
 Pb_LDM = []
 F_N = []
 F_SE = []
 F_LDM = []
 #array to store Binding Energy (BE) per nucleon for all nuclei in data base, format : [[mass, BE/A, proton]]
 AvgBEmat = []
 #array to store most stable nuclie, ie., maximum BE/A, format : [[mass, BE/A, proton]]
 AvgBE_max = []
 #array to store most stable nuclei's corresponding mass A and BE, for plotting, format: [mass] ; [BE]
 max_A = []
 max_BE = []
 #array to store LDM data, again, A is the corresponding most stable nuclei's mass #, etc. , format : [mass] ; [LDM BE/A] ...
 LDM_A = []
 LDM_avgBE0 = []
 LDM_avgBE1 = []
 LDM_avgBE12 = []
 LDM_avgBE123 = []
 #array that stores most stable nuclei's mass # A and proton # Z, we later use these nuclei as the 'most stable nuclei',
 #since we need neutron # N in LDM BE equation, format : [mass, proton]
 stableAZ = []
 #list of which elements (proton number) exists in bedata.dat file
 zList = []
 #array that stores neutron dripline from LDM calculation, format : [Proton#, Neutron#, LDM SE of next isotope]
 dripLine = []
 
 for line in lines0:
     LineRaw = line.split()
     z = int(LineRaw[0])
     n = int(LineRaw[5])
     BE = float(LineRaw[2])
     BEmat.append([z,n,BE])

 BE_len = int(len(BEmat))

 for i in range(0,BE_len):
   #Oxygen chain
   if (BEmat[i][0] == 8 and BEmat[i-1][0] == 8 and (BEmat[i][1] - BEmat[i-1][1] == 1) ):
     O_N.append(BEmat[i][1])
     O_SE.append(BEmat[i][2] - BEmat[i-1][2])
     O_LDM.append(LDM_SE(BEmat[i][1],BEmat[i][0]))
   
   #Ca chain
   if (BEmat[i][0] == 20 and BEmat[i-1][0] == 20 and (BEmat[i][1] - BEmat[i-1][1] == 1) ):
     Ca_N.append(BEmat[i][1])
     Ca_SE.append(BEmat[i][2] - BEmat[i-1][2])
     Ca_LDM.append(LDM_SE(BEmat[i][1],BEmat[i][0]))

   #Ni chain
   if (BEmat[i][0] == 28 and BEmat[i-1][0] == 28 and (BEmat[i][1] - BEmat[i-1][1] == 1) ):
     Ni_N.append(BEmat[i][1])
     Ni_SE.append(BEmat[i][2] - BEmat[i-1][2])
     Ni_LDM.append(LDM_SE(BEmat[i][1],BEmat[i][0]))

   #Sn chain
   if ( (BEmat[i][0] == 50) and (BEmat[i-1][0] == 50) and (BEmat[i][1] - BEmat[i-1][1] == 1) ):
     Sn_N.append(BEmat[i][1])
     Sn_SE.append(BEmat[i][2] - BEmat[i-1][2])
     Sn_LDM.append(LDM_SE(BEmat[i][1],BEmat[i][0]))

   #Pb chain
   if (BEmat[i][0] == 82 and BEmat[i-1][0] == 82 and (BEmat[i][1] - BEmat[i-1][1] == 1) ):
     Pb_N.append(BEmat[i][1])
     Pb_SE.append(BEmat[i][2] - BEmat[i-1][2])
     Pb_LDM.append(LDM_SE(BEmat[i][1],BEmat[i][0]))

   #F chain
   if (BEmat[i][0] == 9 and BEmat[i-1][0] == 9 and (BEmat[i][1] - BEmat[i-1][1] == 1) ):
     F_N.append(BEmat[i][1])
     F_SE.append(BEmat[i][2] - BEmat[i-1][2])
     F_LDM.append(LDM_SE(BEmat[i][1],BEmat[i][0]))
   


   
   #Entire dataset, BE/nucleon vs A
   Z = BEmat[i][0]
   A = BEmat[i][0] + BEmat[i][1]
   AvgBE = BEmat[i][2]*1.0/A
   AvgBEmat.append([A,AvgBE,Z])
   AvgBEmat_len = int(len(AvgBEmat))
   maxA = A
   #print maxA
   
 #existing Z in dataset
 for j in range(1,120):
   for i in range(0,BE_len):
     if ( BEmat[i][0] == j):
       zList.append(j)
       break
 len2 = len(zList)
 
 #find the max BE/nucleon (most stable for each element) list AvgBE_max: [A,BE/nucleon] in exp. data
 for i in range(1,len2):         #for each element, loop thru Z
   for j in range(0,AvgBEmat_len):    #loop thru entire AvgBE array
     Z = AvgBEmat[j][2]
     A = AvgBEmat[j][0]
     if (Z == i):
       AvgBE_max.append([A,0,Z])
       break
 len0 = int(len(AvgBE_max))

 for j in range(0,len0):
   for i in range(0,AvgBEmat_len):
     if (AvgBEmat[i][2] == AvgBE_max[j][2]):
       if ( (AvgBEmat[i][1] > AvgBE_max[j][1])):
         AvgBE_max[j][1] = AvgBEmat[i][1]
         AvgBE_max[j][0] = AvgBEmat[i][0]


 for i in range(0,len0):
   max_A.append(AvgBE_max[i][0])
   max_BE.append(AvgBE_max[i][1])
   stableAZ.append([AvgBE_max[i][0],AvgBE_max[i][2]])
 len1 = int(len(stableAZ))

 #print stableAZ

 #find BE/nucleon vs A using Liquid Drop Model, and separated contribution from different terms
 for i in range(0,len1):
   A = stableAZ[i][0]
   Z = stableAZ[i][1]
   LDM_A.append(A)
   LDM_avgBE1.append([LDM_BE(A,Z,0)/A])
   LDM_avgBE12.append([LDM_BE(A,Z,1)/A])
   LDM_avgBE123.append([LDM_BE(A,Z,2)/A])
   LDM_avgBE0.append([LDM_BE(A,Z,3)/A])

 #Use LDM to find neutron drip lines for Z values up to 120.
 for i in range(1,121):
   for j in range(1,301):
     Z = i
     N = j
     if (LDM_SE(N,Z) < 0):
       dripLine.append([Z,N-1,LDM_SE(N,Z)])
       break


 #write neutron dripline into file
 target = open(str('LDM_neutron_dripline.dat'),"w")
 #title line
 title = '##N is the neutron number of the last isotope to have SE > 0, ie. neutron drip line\n##all data predicted by Liquid Drop Model \n##Separation Energy colomn is the SE of next isotope \n##it should always be negative in this file \n\n'+'Z'.rjust(3)+'     '+'N'.rjust(3)+'     '+'SE (MeV)'.rjust(8)+'\n'
 target.write(title)
 #define output data string
 linew = ''
 
 for Z in range(1,121):
   linew += str(Z).rjust(3)+'     '+str(dripLine[Z-1][1]).rjust(3)+'     '+str(dripLine[Z-1][2]).rjust(8)+'\n'
 target.write(linew)
 print linew

 #plot exp. data value for BE/A vs A
 #majorLocator = MultipleLocator(2)
 #majorFormatter = FormatStrFormatter('%d')
 minorLocator0 = MultipleLocator(1)
 minorLocator1 = MultipleLocator(1)
 minorLocator2 = MultipleLocator(10)

 fig = plt.figure(figsize=(8, 6))
 gs = gridspec.GridSpec(3,1)

 ax0 = plt.subplot(gs[:2,:])
 ax1 = plt.subplot(gs[2,:])
 #yticks range and major spacing
 plt.ylim(0,16)
 plt.yticks(np.arange(0,16,1))
 plt.ylim(0,10)
 plt.yticks(np.arange(0,9,2))

 #major / minor ticks for y
 ax0.yaxis.set_minor_locator(minorLocator0)
 ax1.yaxis.set_minor_locator(minorLocator1)
 #minor ticks for x
 ax1.xaxis.set_minor_locator(minorLocator2)

 #ax1.yaxis.set_minor_locator(minorLocator)
 #subplot for LDM data
 ax0.plot(LDM_A,LDM_avgBE0)
 ax0.plot(LDM_A,LDM_avgBE1)
 ax0.plot(LDM_A,LDM_avgBE12)
 ax0.plot(LDM_A,LDM_avgBE123)
 ax0.set_ylabel('BE/A (MeV)     Liquid Drop Model')

 #subplot for experimental data
 ax1.plot(max_A,max_BE)
 ax1.set_ylabel('BE/A (MeV)')




 #shared x axis
 fig.subplots_adjust(hspace=0)
 plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
 
 plt.text(265,29.1,'A')
 plt.text(265,25.6,'B')
 plt.text(265,19.8,'C')
 plt.text(265,18.6,'D')
 plt.text(265,7.1,'exp')


 plt.xlabel('Mass Number A')
 #plt.show()




#separation energy by using liquid drop model
def LDM_SE(N,Z):
    A = N + Z
    a1 = 15.49
    a2 = 17.23
    a3 = 0.697
    a4 = 22.6
    LDM_BE0 = a1*A - a2*(A**(2.0/3)) - a3*(Z**2)/(A**(1.0/3)) - a4*((N-Z)**2)/A
    LDM_BE1 = a1*(A-1) - a2*((A-1)**(2.0/3)) - a3*(Z**2)/((A-1)**(1.0/3)) - a4*((N-1-Z)**2)/(A-1)
    return LDM_BE0 - LDM_BE1

#binding energy by using LDM
def LDM_BE(A,Z,p):
    N = A - Z
    #p = : 0: only vol. ; 1: vol+surf ; 2: vol+surf+colu ; 3: all 4 terms.
    a1 = 15.49
    a2 = 17.23*int(p == 1 or p == 2 or p == 3)
    a3 = 0.697*int(p == 2 or p == 3)
    a4 = 22.6*int(p == 3)
    LDM_BE = a1*A - a2*(A**(2.0/3)) - a3*(Z**2)/(A**(1.0/3)) - a4*((N-Z)**2)/A
    return LDM_BE

 


InFile0='bedata.dat'

PES(InFile0)


##for plotting separation energies

# plt.plot(O_N,O_SE,'b.-')
# plt.plot(O_N,O_LDM,'b+-')
# plt.text(18,-2,'O')
 
# plt.plot(Ca_N,Ca_SE,'y.-')
# plt.plot(Ca_N,Ca_LDM,'y+-')
# plt.text(33,2,'Ca')
 
# plt.plot(Ni_N,Ni_SE,'r.-')
# plt.plot(Ni_N,Ni_LDM,'r+-')
# plt.text(45,2,'Ni')
 
# plt.plot(Sn_N,Sn_SE,'m.-')
# plt.plot(Sn_N,Sn_LDM,'m+-')
# plt.text(86,2,'Sn')
 
# plt.plot(Pb_N,Pb_SE,'g.-')
# plt.plot(Pb_N,Pb_LDM,'g+-')
# plt.text(134,3,'Pb')
 
# plt.plot([120],[25],'k+')
# plt.text(120,24.7,'   Liquid Drop Model')
# plt.plot([120],[26],'k.')
# plt.text(120,25.7,'   Exp. Data')
 
# plt.xlabel('Neutron Number')
# plt.ylabel('Separation Energy (MeV)')

 #plot fluorine SE, data vs exp.
 #plt.plot(F_N,F_SE,'b.-')
 #plt.plot(F_N,F_LDM,'b+-')
 #plt.text(16,8,'F')
 #plt.title('Separation Energies of Fluorine isotopes')
 #plt.xlabel('N')
 #plt.ylabel('$S_n$ (MeV)')
 #plt.plot([16],[23],'b+')
 #plt.text(16,22.82,'   Liquid Drop Model')
 #plt.plot([16],[24],'b.')
 #plt.text(16,23.82,'   Exp. Data')
 #plt.show()