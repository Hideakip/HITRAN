__author__ = 'Hideaki'
import numpy as np
import array
import os as os

global database_path, gas_names
gas_names=["N2O","CO","CH4","O2","NO","SO2","NO2","NH3","HNO3","OH","HF","HCl","HBr","HI","ClO","OCS","H2CO","HOCl","N2","HCN","CH3Cl","H2O2","C2H2","C2H6","PH3","COF2","SF6","H2S","HCOOH","HO2","O","ClONO2","NO+","HOBr","C2H4"]
database_path="C:\Users\Hideaki\Dropbox\Scytronix\Software\Eclipse projects\HITRAN v2\Raw data par"

def read_gas_file(gas_name):
        #Wn=Wavenmber;I=Line intensity;A=Einstien Coefficient;Abr=Air broadening;Sbr=Self broadening /n
        ##LE=Lower state energy;DTAbr=Temp dependence air width;Psh=pressure shift
        ##UVQ=Upper vibrational quanta;LVQ=lower vibrational quanta;LLQ=Lower local quanta
        #global wnpeak,S,A,Abr,Sbr,LSE,DTAbr,Psh
        filelist=["01_hit08.par", "02_hit08_f53.par", "03_hit08.par", "04_hit08.par", "05_hit08.par", "06_hit08_f53.par", "07_hit08.par", "08_hit08.par", "09_hit08.par", "10_hit08.par", "11_hit08.par", "12_hit08.par", "13_hit08.par", "14_hit08.par", "15_hit08.par", "16_hit08.par", "17_hit08.par", "18_hit08.par", "19_hit08.par", "20_hit08.par", "21_hit08.par", "22_hit08.par", "23_hit08.par", "24_hit08.par", "25_hit08.par", "26_hit08.par", "27_hit08.par", "28_hit08.par", "29_hit08.par", "30_hit08.par", "31_hit08.par", "32_hit08.par", "33_hit08.par", "34_hit08.par", "35_hit08.par", "36_hit08.par", "37_hit08.par", "38_hit08.par", "39_hit08.par", "40_hit08.par", "41_hit08.par"]
        gas_names=["N2O","CO","CH4","O2","NO","SO2","NO2","NH3","HNO3","OH","HF","HCl","HBr","HI","ClO","OCS","H2CO","HOCl","N2","HCN","CH3Cl","H2O2","C2H2","C2H6","PH3","COF2","SF6","H2S","HCOOH","HO2","O","ClONO2","NO+","HOBr","C2H4"]

        index=gas_names.index(gas_name)
        file_name=filelist[index]

        os.chdir(database_path)
        data=file(file_name)
        i=1
        OutArray=np.array([[0,0,0,0,0,0,0,0,0]])
        for line in data:
            #line.replace(" ","*")
            #print line.split("*")
            #linelist=(" ".join(line.split())).split()
            #print linelist
            #wn=float(linelist[1])
            wn=float((line[3:14]).replace(' ',''))
            iso=float(line[2:3].replace(' ',''))

            #Parameter_list=["wn","I","A","Abr","Sbr","LE","DTAbr","Psh","UVQ","LVQ","LLQ"]

            Sin=(line[15:25]).replace(' ','')
            Ain=(line[26:35]).replace(' ','')
            Abrin=(line[35:40]).replace(' ','')
            Sbrin=(line[40:46]).replace(' ','')
            LSEin=(line[46:55]).replace(' ','')
            DTAbrin=(line[56:59]).replace(' ','')
            Pshin=(line[59:66]).replace(' ','')
            row=np.array([[float(wn),float(iso),float(Sin),float(Ain),float(Abrin),float(Sbrin),float(LSEin),float(DTAbrin),float(Pshin)]])
            #print "Row=",row
            OutArray=np.vstack((OutArray,row))

        OutArray=np.delete(OutArray,0,0)
        print "Done ",np.size(OutArray)
        np.save(gas_name+".npy",OutArray)
        return OutArray

for g in gas_names:
    print "starting gas: ",g
    read_gas_file(g)
    print "Done gas: ",g