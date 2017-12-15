import os, sys
from ROOT import gROOT, TCanvas, TF1, TH1F
import time

inputfile = open('iowait.txt','r+')
iowait_list = []
cpu_list = []
counter=0

host = sys.argv[1]
print host

timerange=500

c1 = TCanvas( 'IOWait', 'IOWait', 700, 500 )
c1.Divide(1,2)
hIOWait = TH1F('IOWait','IOWait',1000000,0,1000000)
hCPU = TH1F('CPU','CPU',1000000,0,1000000)
hCPU.SetLineColor(2);

var = 1
while var == 1 :  # This constructs an infinite loop
    inputfile = open(host,'r+')
    for line in inputfile:
        if "avg-cpu" in line:
            data = inputfile.next()
            iowait = data.split()[3]
            cpu = 100.0-float(data.split()[5])
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
    c1.cd(2)
    hCPU.Draw()
    c1.Modified()
    c1.Update()
    print "here"
    inputfile.truncate(0)
    inputfile.close()
    time.sleep(5) 

