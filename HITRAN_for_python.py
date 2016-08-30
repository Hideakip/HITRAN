#Implements efficient analytical approximation for Voigt function

from numpy import *
from pylab import *
import TIP
import os as os
import csv
import time
import h5py


class HITRAN():
    
    def __init__(self):
        
        global Tl,NL,c,R,k,P0,Da, parameter_path, parameter_file,database_path,calculation_file,gas_path,gas_file,mode
        
        Tl=296                  #Reference temperature k
        NL=2.479e19             #Molecular density a 1 atm molecules/cm3   
        c=3e8                   #Velocity of light m/s
        R=8.31451               #Ryberg gas constant
        k=1.38e-23              #Boltzman's constant
        P0=1.0                  #Reference pressure in atm
        Da=1.660539e-27         #1 Dalton in kg
        parameter_path="C:\Users\Hideaki\Dropbox\Scytronix\Software\Eclipse projects\HITRAN v6\\"
        parameter_file="global_parameters.txt"
        database_path="C:\Users\Hideaki\Dropbox\Scytronix\Software\Eclipse projects\HITRAN v2\Raw data par"    
        calculation_file=parameter_path+"calcs.hdf5"
    def read_global_parameters(self):
        global T,Pr,lpath,Adet,tau,source,sourceP,wn1,wn2,res,NDstar,Dstar,DstarLamb,gas_path,gas_file,mode
        gparameters=file(parameter_path+parameter_file)
        linecount=0        
        for line in gparameters:
            linelist=line.split(",")
            if linecount==0:            
                T=float(linelist[0])        #Temperature
                print "Temperature=",T,"K"                
                linecount+=1
            elif linecount==1:
                Pr=float(linelist[0])       #Pressure
                print "Pressure=",Pr,"atm"                
                linecount+=1
            elif linecount==2:
                lpath=float(linelist[0])    #Optical path length7
                print "Optical pathe length=",lpath,"m"               
                linecount+=1
            elif linecount==3:            
                Adet=float(linelist[0])        #Detector area
                print "Detector area=",Adet,"mm2"                
                linecount+=1 
            elif linecount==4:
                tau=float(linelist[0])
                print "Integration time=",tau,"s"
                linecount+=1
            elif linecount==5:              #Source type DFB, FP, LED, Lamp, Blackbody    
                source=linelist[0]
                print "Source ",source               
                linecount+=1
            elif linecount==6:            
                sourceP=float(linelist[0])   #optical power of light source
                print "Source power=",sourceP,"mW"              
                linecount+=1
            elif linecount==7:
                lamb1=float(linelist[0])     #start wavenumber
                lamb2=float(linelist[1])      #finish wavenumber
                wn1=1e4/lamb2
                wn2=1e4/lamb1
                print "Start wavenumber=",wn1,"cm-1"                
                print "Finish wavenumber=",wn2,"cm-1"                 
                linecount+=1 
            elif linecount==8:
                res=float(linelist[0])      #resolution
                linecount+=1
                print "Resolution=",res,"cm-1"
            elif linecount==9:
                NDstar=float(linelist[0])   #Number of points in D* list
                NDstar=int(NDstar)                
                linecount+=1
   
            elif linecount==10:
                Dstar=[]
                for n in range(0,NDstar):         
                    Dstar.append(linelist[n])    #Read D* parameters into list
                linecount+=1
            elif linecount==11:
                DstarLamb=[]
                for n in range(0,NDstar):
                    DstarLamb.append(linelist[n])  #Read D*wavelengths into list
                linecount+=1
            elif linecount==12:
                gas_path=linelist[0]
                linecount+=1
            elif linecount==13:
                gas_file=linelist[0]
                print gas_path+gas_file
                linecount+=1
            elif linecount==14:
                mode=float(linelist[0])
                print "mode=",mode
                
        print "Parameters read"
        #return T,Pr,lpath,Adet,sourceP,wn1,wn2,res,Dstar,DstarLamb,gas_path,gas_file

   
    def read_gas_species(self):
        
        print parameter_path+gas_file
        gas_data=file(parameter_path+gas_file)
        gas=[]
        concentration=[]
        isotope=[]
        for n in gas_data:
            linelist=n.split(",")
            gas.append(linelist[0])
            concentration.append(float(linelist[1]))
            isotope.append(float(linelist[2]))
        print "Gases read"
        return gas,concentration,isotope
        
        
    def read_gas_file(self,gas_name,isotope,wn1,wn2):
        #Wn=Wavenmber;I=Line intensity;A=Einstien Coefficient;Abr=Air broadening;Sbr=Self broadening /n
        ##LE=Lower state energy;DTAbr=Temp dependence air width;Psh=pressure shift
        ##UVQ=Upper vibrational quanta;LVQ=lower vibrational quanta;LLQ=Lower local quanta
        #global wnpeak,S,A,Abr,Sbr,LSE,DTAbr,Psh
        filelist=["01_hit08.par", "02_hit08_f53.par", "03_hit08.par", "04_hit08.par", "05_hit08.par", "06_hit08_f53.par", "07_hit08.par", "08_hit08.par", "09_hit08.par", "10_hit08.par", "11_hit08.par", "12_hit08.par", "13_hit08.par", "14_hit08.par", "15_hit08.par", "16_hit08.par", "17_hit08.par", "18_hit08.par", "19_hit08.par", "20_hit08.par", "21_hit08.par", "22_hit08.par", "23_hit08.par", "24_hit08.par", "25_hit08.par", "26_hit08.par", "27_hit08.par", "28_hit08.par", "29_hit08.par", "30_hit08.par", "31_hit08.par", "32_hit08.par", "33_hit08.par", "34_hit08.par", "35_hit08.par", "36_hit08.par", "37_hit08.par", "38_hit08.par", "39_hit08.par", "40_hit08.par", "41_hit08.par"]
        gas_names=["H2O","CO2","O3","N2O","CO","CH4","O2","NO","SO2","NO2","NH3","HNO3","OH","HF","HCl","HBr","HI","ClO","OCS","H2CO","HOCl","N2","HCN","CH3Cl","H2O2","C2H2","C2H6","PH3","COF2","SF6","H2S","HCOOH","HO2","O","ClONO2","NO+","HOBr","C2H4"]

        index=gas_names.index(gas_name)     
        file_name=filelist[index]
        #global wnpeak,S,A,Abr,Sbr,LSE,DTAbr,Psh
        iso=[];wnpeak=[];S=[];A=[];Abr=[];Sbr=[];LSE=[];DTAbr=[];Psh=[];UVQ=[];LVQ=[];LLQ=[]
        os.chdir(database_path)
        data=file(file_name)
        i=1
        for line in data:
            #line.replace(" ","*")
            #print line.split("*")   
            #linelist=(" ".join(line.split())).split()
            #print linelist            
            #wn=float(linelist[1])
            wn=float((line[3:14]).replace(' ',''))          
            iso=float(line[2:3].replace(' ',''))

            Parameter_list=["wn","I","A","Abr","Sbr","LE","DTAbr","Psh","UVQ","LVQ","LLQ"]            
            if wn>=wn1 and iso==isotope:
                
                if wn<=wn2:
                    
                    wnpeak.append(float(wn))
                  
                    Sin=(line[15:25]).replace(' ','')
                    S.append(float(Sin))
                    Ain=(line[26:35]).replace(' ','')
                    A.append(float(Ain))
                    Abrin=(line[35:40]).replace(' ','')
                    Abr.append(float(Abrin))
                    Sbrin=(line[40:46]).replace(' ','')
                    Sbr.append(float(Sbrin))
                    LSEin=(line[46:55]).replace(' ','')
                    LSE.append(float(LSEin))
                    DTAbrin=(line[56:59]).replace(' ','')
                    DTAbr.append(float(DTAbrin))
                    Pshin=(line[59:66]).replace(' ','')
                    Psh.append(float(Pshin))
        print "Number of lines loaded: ",len(wnpeak)

        return wnpeak,S,A,Abr,Sbr,LSE,DTAbr,Psh
        
    def do_gases(self):
        t0=time.time()
        gaslist=["CH4","CO2","CO","H2O","N2O","NH3","NO", "H2S"]
        MW=[16.04246,44.0095,28.0101,18.01528,44.0128,17.03052,30.0061,34.08]           #Molecular weights g/mol      
    
        gas,concentration,isotope =self.read_gas_species()
        print "Gas species read",time.time()-t0
        ##Systematically loop through all gases and concentrations in gas file       
        wnout=[]
        dw=[]       #peak widths in wn for stick spectra
        alphaout=[]
        gasi=0
        AllGasData=[]        
        for gasn in gas:
            
            print "Processing gas ",gasn, "........"
            i=gaslist.index(gasn)           #Find index of file in filenames list     
            M=MW[i]                             #Look up associated parameters  
            #Scale molecular concentration for temperature and pressure               
            Nx=self.Nx_P_T(concentration[gasi])
            t0=time.time()
            wnpeak,S,A,Abr,Sbr,LSE,DTAbr,Psh=self.read_gas_file(gasn,isotope[gasi],wn1,wn2)
            print "Gas file read read",time.time()-t0
            wnpeakS=self.pressure_shift(wnpeak,Psh)
            print
            t0=time.time()
            gammap,gammad=self.calculate_broadening(wnpeakS,M,Abr,Sbr,DTAbr)
            print "Broadening calculated",time.time()-t0
            St=self.Temp_shift_S(gasn,S,wnpeak)           #Rescale line intensity of temperature. St list of peaks
            #AllGasData.append([wn,St,A,Abr,Sbr,LSE,DTAbr,Psh])
            ##Scale all peaks and concentration for temperature and pressure
  
            
            
            if mode==0:     #stick spectra
                 print "Stick spectra"
                 alphaL,alphaG=self.stick_spectra(St, Nx,gammap,gammad)
                 #wnt,alphaL=self.traiangle_peak(St,Nx,alphaL)
                 dwi=self.generate_width_positions(St, Nx, alphaL)
                 wnout.append(wnpeakS)

                 dw.append(dwi)
                 alphaL=alphaL*1.92      #fiddle factor
                 alphaout.append(alphaL)
                 plt.plot(wnpeak,alphaL)
                 plt.show()

            elif mode==1:   #Gauss Lorentz
                 print "Gauss/Lorentz"
                 wni=self.create_wn_list(wnpeakS)           #Make list of all wavenumbers in range wn1 to wn2 step res including wn's of peaks in line centre list
                 dvps,peakrnages,gammap,gammad=self.create_peak_ranges(wnpeakS,wni,M,Abr,Sbr,DTAbr,Psh)
                 print "Number of wavenumber samples=",len(wni)
                 alpha=self.Gauss_Lorentz(peakranges,fronts,backs)
                 alphaout.append(alpha)
                 wnout.append(wni)


            elif mode==2:           #Voigt
                 print "Voigt"
                 t0=time.time()
                 wni=self.create_wn_list(wnpeakS)           #Make list of all wavenumbers in range wn1 to wn2 step res including wn's of peaks in line centre list
                 print "wn list done",time.time()-t0
                 alphaL,alphaG=self.stick_spectra(St, Nx,gammap,gammad)
                 max_alpha=max(alphaL)
                 imax=(alphaL.tolist()).index(max_alpha)
                 t0=time.time()
                 self.create_peak_ranges(wnpeakS,wni,gammap,gammad,max_alpha,imax)   #Returns nothing because data written to hdf5 file
                 print "Peak ranges calculated",time.time()-t0
                 t0=time.time()
                 self.Voigt1(wnpeakS,gammap,gammad)
                 print "Voigt calculated",time.time()-t0
                 #gv1=self.Voigt2(dvps,fronts,backs,gammap,gammad)

                 t0=time.time()
                 self.alpha(wnpeakS,St,Nx,wni)
                 print "Alpha calculated",time.time()-t0


                 plot(1.0e-4/array(wni),alphaS)
                 show()
                 #wnA,alphaA=self.assemble_spectra(fronts,backs,peakranges,alpha)

                 #alphaout=alphaA

                 wnout.append(wni)
                 j=0


            print "Gas ",gasn, "  Complete"
            gasi+=1

        print "Gas run complete"
        return wnout,dw,alphaout


    def create_wn_list(self,wnpeakS):
        ##Interdigitates wavenumbers of peaks into the wavenumber range defined by wn1,wn2, and Red
        wni=arange(wn1,wn2,res)     #Generate list of wavenumbers at required resolution
        wni=concatenate([wni,wnpeakS])             #Append wavenumbers due to peaks
        wni=sorted(wni)                 #Sort ascending
        return wni

    def calculate_broadening(self,wnpeakS,M,Abr,Sbr,DTAbr):
        #Calculate pressure broadened half-width for each peak

        Abrm=array(Abr)
        SBrm=array(Sbr)
        DTAbrm=array(DTAbr)
        q=0.877 #concentration*Pr
        gammap=(q*Abrm+(1-q)*SBrm)*Pr*(Tl/T)**DTAbrm         #Pressure broadened half-width for each peak

        #Calculate temperature broadened half-width for each peak
        v0=array(wnpeakS)
        m=M*Da
        gammad=(v0*100/c*sqrt(2*log(2)*k*T/m))/100
        return gammap,gammad

    def create_peak_ranges(self,wnpeakS,wni,gammap,gammad,max_alpha,imax):
        fcalcs=h5py.File(calculation_file,'w')
        #Selects wave number values either side of peak position to limit size of computational matrix
        fronts=fcalcs.create_group("fronts")               #List of lists of all  the range in front of the peak where peak shape is not computed
        peakranges=fcalcs.create_group("peakranges")           #List of lists of all range over which peak shape is computed
        dvps=fcalcs.create_group("dvps")
        backs=fcalcs.create_group("backs")                #Lis of lists of all range behind peak where peak shape is not computed
        
        i=0 #Index counter for position in list of peaks wmpeakS
        N=1e-2 #Low intensity limit of calculation. Fraction of highest peak in range
        peakmin=max_alpha*N
        if gammap[imax]>gammad[imax]:
            sigma=sqrt(gammap[imax]/(pi*peakmin)-gammap[imax]**2)
        if gammad[imax]>gammap[imax]:
            sigma=sqrt(-gammad[imax]**2/log(2)*log(peakmin*gammad[imax]*sqrt(pi/log(2))))
        wni=array(wni)

        for peak in wnpeakS:
            peakname=str(peak)
            a=peak-sigma
            b=peak+sigma

            peakrange=wni[[wni>=a] and [wni<=b]]
            dvp=self.create_difference_matrix(peakrange,peak)   #Calculate displacement either side of peak instead of absolute wave number
            front=len(wni[wni<a])
            back=len(wni[wni>b])

            #add data to file
            peakranges.create_dataset(peakname,data=peakrange)
            dvps.create_dataset(peakname,data=dvp)
            fronts.create_dataset(peakname,data=front)
            backs.create_dataset(peakname,data=back)
            i+=1

        fcalcs.close()

    def create_difference_matrix(self,peakrange,peak):
        #Create matrix of difference between the peak and wavelength

        dvp=peakrange-peak
        #print size(dvp),peak
        return dvp
    
    def Nx_P_T(self,concentration):   
        Nx=concentration*NL*296.0/T*Pr       #Molecular concentration
        return Nx

    def pressure_shift(self,wnpeak,Psh):
        wnpeakS=array(wnpeak)+array(Psh)*Pr  #Pressure shift
        return wnpeakS
    
    def stick_spectra(self,St,Nx,gammap,gammad):
        #Lorentian/pressure broadening
        Gp=0.318/gammap

        #Gaussian Doppler broadening
        gammad=gammad*2.6     #fiddle factor to broaden line width
        Gd=0.469/gammad 
        alphaL=St*Gp*Nx
        alphaG=St*Gd*Nx
        return alphaL,alphaG
    
    def Gauss_Lorentz(self,wni,dvp,St,Nx,M,concentration):

        Npeaks=len(wnpeak)
        #Pressure shift
        wnpeak_corr=array(wnpeak)+array(Psh)*Pr
        Abrm=array(Abr)
        SBrm=array(Sbr)
        DTAbrm=array(DTAbr)
        Pshm=array(Psh)
        q=0.877 #concentration*Pr
        gammap=(q*Abrm+(1-q)*SBrm)*Pr*(Tl/T)**DTAbrm         #Pressure broadened half-width for each peak
        
        #Loop through each peak
        i=0
        alphaL_out=[]
        for row in range(0,Npeaks):
            dw=dvp[row,:]
            Lorentz=(gammap[row]/pi)/(dw**2+gammap[row]**2)
            alpha=St[row]*Lorentz*Nx
            alphaL_out.append(alpha)
        alphaL=array(alphaL_out).sum(axis=0)     #adds all peaks for a signal gas
        
        
        
        
        #Doppler/Gaussian broadening  
#        v0=array(wnpeak)     
#        gammad=(v0*100/c*(2*R*T*log(2)/(M/1000))**0.5)/100
#        alphaG_out=[]
#        for row in range(0,Npeaks):
#            #dw=(array(wni)-wp)**2
#            Gauss=1/gammad[row]*(log(2)/pi)**0.5*exp(-log(2)*dvp**2/gammad[row]**2)
#            alpha=St[row]*Gauss*Nx
#            alphaG_out.append(alpha)
#        alphaG=array(alphaG_out).sum(axis=0)
#        Lamb=1e4/array(wni)
#        alphaG=St*Gauss_out*Nx
#        print "35456werweewewew"
        return alphaL

    def Voigt1(self,wnpeakS,gammap,gammad):
        #Voigt profile according to Olivero and Longbothum
        fcalcs=h5py.File(calculation_file,'r+')
        gvs=fcalcs.create_group("gvs")
        #Convert HWHM to FWHM
        dvpg=fcalcs["dvps"]         #read group

        #gvP=[]      #List storing all Voigt profiles per peak
        i=0
        for peak in wnpeakS:  #Scan through wavelength range set
            peakname=str(peak)
            dvp=dvpg[peakname]
            wl=2*gammap[i]
            wd=2*gammad[i]
            wv=0.5346*wl+sqrt(0.2166*wl**2+wd**2)
            wlwvR=wl/wv
            Igvmax=1/(wv*(1.065+0.447*wlwvR)+0.058*wlwvR**2)

            DN =array(dvp)
            DNwv=(DN/wv)
            DNwvS=DNwv**2
            DNwvSF=abs(DNwv)**2.25

            gv=((1-wlwvR)*exp(-2.772*DNwvS)+(wlwvR/(1+4.0*DNwvS)+0.016*(1-wlwvR)*wlwvR*(exp(-0.4*DNwvSF)-10.0/(10.0+DNwvSF))))
            gv=Igvmax*gv


            #store peak in file
            gvs.create_dataset(peakname,data=gv)
            i+=1

        fcalcs.close()

    def Voigt2(self,dvps,fronts,backs,gammap,gammad):

        A=[-1.2150,-1.3509,-1.2150,-1.3509 ]
        B=[1.2359,0.3786,-1.2359,-0.3786]
        C=[-0.3085,0.5906 ,-0.3085,0.5906 ]
        D=[0.0210,-1.1858,-0.0210, 1.1858]

        i=0
        gv=[]
        for dvp in dvps:
            Cx=2*sqrt(log(2)/(2*gammad[i]))
            Y=(gammap[i]/gammad[i])*sqrt(log(2))
            alphaL=1/(pi*gammap[i])
            C1=gammap[i]/gammad[i]*sqrt(pi*log(2))*alphaL
            X=Cx*dvp    #(peakrange)
            S=0
            for j in range(0,4):
                Sp=(C[j]*(Y-A[j])+D[j]*(X-B[j]) )/((Y-A[j])**2+(X-B[j])**2)
                S=S+Sp
            I=C1*S


            #concatenate fronts and backs
            frontpad=[0]*len(fronts[i])
            backpad=[0]*len(backs[i])
            gvi=frontpad+I.tolist()+backpad
            gv.append(gvi)

            i+=1

        return gv

    def alpha(self,wnpeakS,St,Nx,wni):
        fcalcs=h5py.File(calculation_file,'r+')
        #create objects for parameters in file
        gvs=fcalcs["gvs"]
        fronts=fcalcs["fronts"]
        backs=fcalcs["backs"]

        #create new groups to store in file
        alphas=fcalcs.create_group("alphas")
        alphatot=fcalcs.create_group("alphatot")
        i=0

        alphatot=zeros(len(wni))
        #calculate loss coefficient for each peak
        for peak in wnpeakS:

            peakname=str(peak)
            gvi=gvs[peakname]
            alphai=St[i]*gvi*Nx*1.0e-6      #!!!!!!!!!!!!Dropped a factor of 1e-6 wrt MS reference data for H20

            #zero pad distribution
            front=fronts[str(peak)]
            back=backs[str(peak)]
            print "front=",front
            print "back=",back
            alphai=hstack([zeros(front),alphai,zeros(back)])
            alphas.create_dataset(peakname,data=alphai)
            alphatot+=alphai


        alphatot.create_dataset(peakname,data=alphatot)
        fcalcs.close()




    def assemble_spectra(self,fronts,backs,peakranges,alpha):
        #Joins uncomputed ranges either side of calculated range around peak
        alphaA=[]
        wnA=[]
        i=0
        for peak in peakranges:

            fzero=[0]*len(fronts[i])
            bzero=[0]*len(backs[i])
            aA=fzero+alpha[i].tolist()+bzero
            alphaA.append(aA)
            wA=fronts[i]+peakranges[i]+backs[i]
            wnA.append(wA)
            i+=1

        return wnA,alphaA

    def Temp_shift_S(self,gasn,S,wnpeak):

        #Run through all peaks S in line centre list and calculate temperature shift
        #in absorption cross section        
        #f filename
        #T gas temperature
        flength=len(S)
        St=[]
        Q=TIP.TIP_interp(T,gasn)
    
        Qref=TIP.TIP_interp(Tl,gasn)
        
        for i in range(0,flength):
            Stl=S[i]*(Qref/Q*exp(1.439*wnpeak[i]*(T-Tl)/(T*Tl)))
            St.append(Stl)                    
        return St                 
    
    def generate_width_positions(self,St,Nx,alpha):
        gammap=array(alpha)/(array(St)*Nx)
        gp=0.318/gammap
        #fiddle factors
 
        return gp
    
    def traiangle_peak(self,St,Nx,alpha):
        gammap=array(alpha)/(array(St)*Nx)
        gp=0.318/gammap
        #fiddle factors
        #
        gp=10.0*gp
        upper_wn=(wnpeak+gp)
        lower_wn=(wnpeak-gp)
        upper_wn=upper_wn.tolist()
        lower_wn=lower_wn.tolist()
        
        wnt=sorted(lower_wn+wnpeak+upper_wn) 
        Aout=[0]
        
        for a in alpha:
            a=a*1.92       #fiddle factor for amplitude
            Aout=Aout+[a]+[0,0]
        del Aout[-1]
        return wnt,Aout

    def write_file(self,wnout,dw,alphaout):
        print "Writting file(s)"
        gas,concentration,isotope =self.read_gas_species()
        Ngases=len(gas)   #Number of gases calculated
        foutpath="C:\Users\Hideaki\Dropbox\Scytronix\Business\Bedfont_scientific\C13C12 measurement\Spectroscopic data\\"
        print "Number of gases=",Ngases
        for i in range (0,Ngases):    
            lamb=1e4/array(wnout[i])          
            minLamb=int(round(min(lamb),3)*1000)
            maxlamb=int(round(max(lamb),3)*1000)
            
            fname=str(gas[i])+"_"+str(minLamb)+"-"+str(maxlamb)+"-"+str(isotope[i])+"_nm"+".txt"
            out=open(foutpath+fname,"w")

            for row in range(0,len(lamb)-1):
                output=str(lamb[row])+","+str(wnout[i][row])+","+","+str(alphaout[i][row])+"\n"   #str(dw[i][row])+
             
                out.write(output)

            #savetxt(foutpath+fname,output, delimiter=",")
            out.close()
            print "File written"
            
    def plot_spectra(self,wnout,alphaout):
        gas,concentration,isotope =self.read_gas_species()
        Ngases=len(gas)
        
        print size(wnout),size(alphaout)
        for i in range (0,Ngases): 
            #print len(wnout[i]),len(alphaout[i])
            lamb=1e4/array(wnout[i]) 
            print len(wnout[i]),len(alphaout[i])
            plot(lamb,alphaout[i])    
            
        show()

if __name__=="__main__":

    H=HITRAN()
    H.read_global_parameters()      #Parameters are global
    wnout,dw,alphaout=H.do_gases()
    print size(alphaout)
    #H.write_file(wnout,dw,alphaout)
    lamb=divide(1e4,wnout[0])
    #print len(wnout[0]),len(alphaout[0])
    print "Finished"
#    H.plot_spectra(wnout, alphaout)
    #plot(lamb,alphaout[0])
    
    #show()


    
    
