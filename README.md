# TheoreticalUncertainties

Calculates theoretical uncertainties for W+Jets and Z+Jets

Tested in CMSSW 7_6_4

To make histograms: cmsRun ConfFile(Z)_cfg.py 

  The user may specify on the command line: doExclusive=True/False (Exclusive jet bins if true, Inclusive if false)
  
                                           jetPtCut = Minimum pt to count jets, default is 30 GeV
                                           
  W+Jets and Z+Jets files are listed in the submission directory
  
  
  
To calculate differences in acceptances due to QCD scale variations: 
  python computeSystByProductionAcceptance(Z).py [filename]
  
Where [filename] is the name of a file that has been produced from the analyzers in the plugins directory
