### Returns variations in acceptances due to scale variation, uses Z+Jets cuts ###


from sys import argv
import ROOT
import math
import numpy as np
fileStr = argv[1]
hfile = ROOT.TFile(fileStr,"READ")
hdir = hfile.Get("demo") 
hlist = hdir.GetListOfKeys()
nJets=8

#kinematic cuts
mu1pt = (nJets+1)*[20]
mu2pt = (nJets+1)*[20]

hWeights=[]
weights = []

#get weights
for hname in hlist:
    if 'Weights' in hname.GetName():
        hWeights.append(hname.GetName())
        weights.append((hname.GetName(), hdir.Get(hname.GetName()).Integral()))

#get theoretical uncertainties
for j in range(0,nJets+1):       
        
	scale=[]
	w_scale =[]
	for hname in hlist:
	    
	    if hdir.Get(hname.GetName()).IsA() != ROOT.TH2F.Class() : continue
	    if str(j)+'j' not in hname.GetName(): continue
	    
	    if 'j1001' in hname.GetName() or 'j1005' in hname.GetName() or 'j1009' in hname.GetName(): #unaltered scale, double scale, half scale
		w_scale.append(hname.GetName())


	#scale syst
	int_scale=[]
	totalIntegral = []
	totalEntries = []


	ref_scale = 0

	for n,name in enumerate(w_scale):
	    
	    if '1001' in name : 
		w_ref_scale = hdir.Get(name)
		ref_scale = w_ref_scale.Integral(mu2pt[j]+1, 100, mu1pt[j]+1, 100)/weights[0][1]
		print str(ref_scale) + "unaltered scale"
		totalIntegral.append(w_ref_scale.Integral())
		totalEntries.append(w_ref_scale.GetEntries())
	    else:
		w_ref_scale = hdir.Get(name)
		if '1005' in name:
			print str(w_ref_scale.Integral(mu2pt[j]+1, 100, mu1pt[j]+1, 100)/weights[4][1])+ "double scale"
		        int_scale.append( w_ref_scale.Integral(mu2pt[j]+1, 100, mu1pt[j]+1, 100)/weights[4][1])
		elif '1009' in name: 
			print str(w_ref_scale.Integral(mu2pt[j]+1, 100, mu1pt[j]+1, 100)/weights[8][1]) + " half scale"
			int_scale.append( w_ref_scale.Integral(mu2pt[j]+1, 100, mu1pt[j]+1, 100)/weights[8][1])
		else:
			int_scale.append( w_ref_scale.Integral(mu2pt[j]+1, 100, mu1pt[j]+1, 100)/weights[n][1])
		totalIntegral.append(w_ref_scale.Integral())
		totalEntries.append(w_ref_scale.GetEntries())
		
	for n,value in enumerate(int_scale):
	    int_scale[n]=value

	unc= (max(int_scale)-min(int_scale))/2/ref_scale #difference in acceptance

	print 'scale %s jet: %.1f or %.1f/%.1f + Add yellow report reccomendation' %(str(j),100*unc, 100*(min(int_scale)-ref_scale)/ref_scale,100*(max(int_scale)-ref_scale)/ref_scale)

