import os, sys
from ROOT import gROOT, TCanvas, TF1, TH1F, TFile
import time

iowait_list = []
cpu_list = []
counter=0

host = sys.argv[1]
rootname = host
rootname = rootname.replace("txt","")
rootname += "root"
print host

timerange=500

c1 = TCanvas( 'IOWait', 'IOWait', 700, 500 )
c1.Divide(1,2)
hIOWait = TH1F('IOWait','IOWait',1000000,0,1000000)
hCPU = TH1F('CPU','CPU',1000000,0,1000000)
hCPU.SetLineColor(2)


var = 1
while var == 1 :  # This constructs an infinite loop
    time.sleep(5) 
    fOut = TFile(rootname,"RECREATE")
    inputfile = open(host,'r+')
    for line in inputfile:
        if "avg-cpu" in line:
            data = inputfile.next()
            iowait = data.split()[3]
            cpu = float(data.split()[5])
            #print iowait
            iowait_list.append(iowait)
            cpu_list.append(cpu)
            print iowait,counter
            hIOWait.Fill(counter,float(iowait))
            hCPU.Fill(counter,float(cpu))
            counter=counter+1
    
    if counter > (timerange+1):
        hIOWait.GetXaxis().SetRangeUser(counter-timerange-1,counter-1)
        hCPU.GetXaxis().SetRangeUser(counter-timerange-1,counter-1)
    else:
        hIOWait.GetXaxis().SetRangeUser(0,timerange)
        hCPU.GetXaxis().SetRangeUser(0,timerange)
    c1.cd(1)
    hIOWait.Draw()
    hIOWait.GetXaxis().SetTitle("Time (s)")
    hIOWait.GetYaxis().SetTitle("IO Wait (%)")
    c1.cd(2)
    hCPU.Draw()
    hCPU.GetXaxis().SetTitle("Time (s)")
    hCPU.GetYaxis().SetTitle("Percent Idle (%)")
    c1.Modified()
    c1.Update()
    print "here"
    inputfile.truncate(0)
    inputfile.close()
    fOut.cd()
    hIOWait.Write()
    hCPU.Write()
    fOut.Write()
    fOut.Close()
