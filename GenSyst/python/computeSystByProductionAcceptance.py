### Returns variations in acceptances due to scale variation, uses W+Jets cuts ###


from sys import argv
import ROOT
import math
import numpy as np
fileStr = argv[1]
hfile = ROOT.TFile(fileStr,"READ")
hdir = hfile.Get("demo") 
hlist = hdir.GetListOfKeys()
nJets = 8
#muon and mt cuts
mupt = (nJets+1)*[25]
mMt = (nJets+1)*[50]

hWeights=[]
weights = []

#get weights
for hname in hlist:
    if 'Weights' in hname.GetName():
        hWeights.append(hname.GetName())
        weights.append((hname.GetName(), hdir.Get(hname.GetName()).Integral()))
	print hname.GetName()

#get theoretical uncertainties
for j in range(0,nJets+1):       
        
	scale=[]
	w_scale =[]
	for hname in hlist:
	    
	    if hdir.Get(hname.GetName()).IsA() != ROOT.TH2F.Class() : continue
	    if str(j)+'j' not in hname.GetName(): continue
	    
	#    if 'j1001' in hname.GetName() or 'j1005' in hname.GetName() or 'j1009' in hname.GetName(): #unaltered scale, double scale, half scale
            if 'j100' in hname.GetName(): #all scales
		print hname
		w_scale.append(hname.GetName())


	#scale syst
	int_scale=[]
	totalIntegral = []
	totalEntries = []


	ref_scale = 0

	for n,name in enumerate(w_scale):
	    print n
	    print name
	    
	    if '1001' in name : 
		w_ref_scale = hdir.Get(name)
		ref_scale = w_ref_scale.Integral(mMt[j]+1, 100, mupt[j]+1, 100)/weights[0][1]
        #        print str(ref_scale) + "unaltered scale"
		totalIntegral.append(w_ref_scale.Integral())
		totalEntries.append(w_ref_scale.GetEntries())
	    else:
		w_ref_scale = hdir.Get(name)
		if '1005' in name:
	#		print str(w_ref_scale.Integral(mMt[j]+1, 100, mupt[j]+1, 100)/weights[n][1])+ "double scale"
			print "adding to int_scale: " + name
			int_scale.append( w_ref_scale.Integral(mMt[j]+1, 100, mupt[j]+1, 100)/weights[4][1])
		elif '1009' in name: 
	#		print str(w_ref_scale.Integral(mMt[j]+1, 100, mupt[j]+1, 100)/weights[n][1]) + " half scale"
			print "adding to int_scale: " + name
			int_scale.append( w_ref_scale.Integral(mMt[j]+1, 100, mupt[j]+1, 100)/weights[8][1])
		else:
			print "adding to int_scale: " + name
			int_scale.append( w_ref_scale.Integral(mMt[j]+1, 100, mupt[j]+1, 100)/weights[n][1])
		totalIntegral.append(w_ref_scale.Integral())
		totalEntries.append(w_ref_scale.GetEntries())
		
	#for n,value in enumerate(int_scale):
	#    int_scale[n]=value
	print ref_scale
	print int_scale
        print "jets: " + str(j) + " refscale: " + str(ref_scale)
	unc= (max(int_scale)-min(int_scale))/2/ref_scale #difference in acceptance

	print 'scale %s jet: %.1f or %.1f/%.f' %(str(j),100*unc, 100*(min(int_scale)-ref_scale)/ref_scale,100*(max(int_scale)-ref_scale)/ref_scale)

