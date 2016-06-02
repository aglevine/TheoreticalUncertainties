// -*- C++ -*-
//
// Package:    LFVAnalysis/GenSyst
// Class:      GenSystInclusive
// 
//
// Original Author:  Silvia Taroni
// Modified by Aaron Levine for W+Jets
//         Created:  Sun, 28 Feb 2016 14:58:30 GMT
//
//

#define DEBUG 0
// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>  
#include "TH1F.h"
#include "TH2F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaPhi.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class GenSystInclusive : public edm::one::EDAnalyzer<edm::one::WatchRuns>  {
   public:
      explicit GenSystInclusive(const edm::ParameterSet&);
      ~GenSystInclusive();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual void beginRun( edm::Run const&,  edm::EventSetup const&) override;
      virtual void endRun( edm::Run const&, edm::EventSetup const &) override;

      double dPhi(math::XYZTLorentzVector&, math::XYZTLorentzVector&);

      edm::EDGetTokenT<GenEventInfoProduct>    tok_gen_;
 
      edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;
      edm::EDGetTokenT<LHEEventProduct> tok_lhe_ ;
  //      edm::EDGetTokenT<reco::GenParticleCollection> tok_genPart_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> tok_genPart_; 
      edm::EDGetTokenT<reco::GenJetCollection> tok_genJet_;
      edm::EDGetTokenT<pat::JetCollection> tok_jet_;
      int nevent_run;  
      double jetPtCut;
      bool doExclusive;

      edm::Service<TFileService> fs;      
      TH1F * wMass; 
      std::vector<TH1F *> wMass_scale;
      std::vector<TH1F *> wMass_pdf;
      TH1F * wMass_0j; 
      std::vector<TH1F *> wMass_scale_0j;
      std::vector<TH1F *> wMass_pdf_0j;
      TH1F * wMass_1j; 
      std::vector<TH1F *> wMass_scale_1j;
      std::vector<TH1F *> wMass_pdf_1j;
      TH1F * wMass_2j; 
      std::vector<TH1F *> wMass_scale_2j;
      std::vector<TH1F *> wMass_pdf_2j;
      TH1F * wMass_3j;
      std::vector<TH1F *> wMass_scale_3j;
      std::vector<TH1F *> wMass_pdf_3j;
      TH1F * wMass_4j;
      std::vector<TH1F *> wMass_scale_4j;
      std::vector<TH1F *> wMass_pdf_4j;
      TH1F * wMass_5j;
      std::vector<TH1F *> wMass_scale_5j;
      std::vector<TH1F *> wMass_pdf_5j;
      TH1F * wMass_6j;
      std::vector<TH1F *> wMass_scale_6j;
      std::vector<TH1F *> wMass_pdf_6j;
      TH1F * wMass_7j;
      std::vector<TH1F *> wMass_scale_7j;
      std::vector<TH1F *> wMass_pdf_7j;
      TH1F * wMass_8j;
      std::vector<TH1F *> wMass_scale_8j;
      std::vector<TH1F *> wMass_pdf_8j;



     std::vector<TH2F *>  wMNuPt_vs_MuPt_0j; 
     std::vector<TH2F *>  wMNuPt_vs_MuPt_1j; 
     std::vector<TH2F *>  wMNuPt_vs_MuPt_2j; 
     std::vector<TH2F *>  wMNuPt_vs_MuPt_3j;
     std::vector<TH2F *>  wMNuPt_vs_MuPt_4j;
     std::vector<TH2F *>  wMNuPt_vs_MuPt_5j;
     std::vector<TH2F *>  wMNuPt_vs_MuPt_6j;
     std::vector<TH2F *>  wMNuPt_vs_MuPt_7j;
     std::vector<TH2F *>  wMNuPt_vs_MuPt_8j;


     std::vector<TH1F *> hWeights; 
   
      TH1F * wMuPt; 
      TH1F * wMuMt;
      TH1F * wMNuPt;
      TH1F * wMNuMt;
      TH1F * wMuPt_0j; 
      TH1F * wMuMt_0j;
      TH1F * wMNuPt_0j;
      TH1F * wMNuMt_0j;
      TH1F * wMuPt_1j; 
      TH1F * wMuMt_1j;
      TH1F * wMNuPt_1j;
      TH1F * wMNuMt_1j;
      TH1F * wMuPt_2j; 
      TH1F * wMuMt_2j;
      TH1F * wMNuPt_2j;
      TH1F * wMNuMt_2j;
      TH1F * wMuPt_3j;
      TH1F * wMuMt_3j;
      TH1F * wMNuPt_3j;
      TH1F * wMNuMt_3j;
      TH1F * wMuPt_4j;
      TH1F * wMuMt_4j;
      TH1F * wMNuPt_4j;
      TH1F * wMNuMt_4j;
      TH1F * wMuPt_5j;
      TH1F * wMuMt_5j;
      TH1F * wMNuPt_5j;
      TH1F * wMNuMt_5j;
      TH1F * wMuPt_6j;
      TH1F * wMuMt_6j;
      TH1F * wMNuPt_6j;
      TH1F * wMNuMt_6j;
      TH1F * wMuPt_7j;
      TH1F * wMuMt_7j;
      TH1F * wMNuPt_7j;
      TH1F * wMNuMt_7j;
      TH1F * wMuPt_8j;
      TH1F * wMuMt_8j;
      TH1F * wMNuPt_8j;
      TH1F * wMNuMt_8j;
    
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenSystInclusive::GenSystInclusive(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   tok_gen_       = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   tok_lhe_       = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer")) ;
   lheRunInfoToken_ = consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer",""));
   //tok_genPart_ =consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
   tok_genPart_ =consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"));
   tok_genJet_ =consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets"));
   tok_jet_ =consumes<pat::JetCollection>(edm::InputTag("slimmedJets"));

   jetPtCut = iConfig.getParameter<double>("jetPt");
   doExclusive = iConfig.getParameter<bool>("doExclusive");

}


GenSystInclusive::~GenSystInclusive()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenSystInclusive::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nevent_run++;
   using namespace edm;
   using namespace reco;
   using namespace std;
   
  if (DEBUG) std::cout << __LINE__ << std::endl;
   edm::Handle<GenEventInfoProduct> genEventInfo;
   iEvent.getByToken(tok_gen_, genEventInfo);


   //std::vector<double>& evtWeights = (std::vector<double>&) genEventInfo->weights();
   double theWeight = genEventInfo->weight();

   edm::Handle<LHEEventProduct> EvtHandle ;
   iEvent.getByToken( tok_lhe_ , EvtHandle ) ;

  if (DEBUG) std::cout << __LINE__ << std::endl;

  // std::string whichWeightId = "1002";
  // for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
  //   if (EvtHandle->weights()[i].id == whichWeightId) theWeight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
  //   if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
    
  // }

  edm::Handle<std::vector<reco::GenJet>> genjets;
  iEvent.getByToken(tok_genJet_,genjets);
  

  
  
  
  edm::Handle<pat::PackedGenParticleCollection> genParticles; 
  iEvent.getByToken(tok_genPart_,genParticles);
 
  const reco::Candidate* wbos  = 0;
  const reco::Candidate* mu     = 0;
  const reco::Candidate* mnu     = 0;

  for (pat::PackedGenParticleCollection::const_iterator p = genParticles->begin();p != genParticles->end(); ++p ) {
    const reco::Candidate *mother = 0;
    if( p->numberOfMothers()!=0) mother= p->mother(0);
    //if (DEBUG) cout<< __LINE__ << " " << p->pdgId()<< " " << p->status() << " ";
    for (unsigned int im =0; im<10 ; im++){
      if (mother->numberOfMothers()!=0 ) {
	mother=&*mother->mother(0);
	//cout << mother->pdgId() << " " ;
	if (abs(mother->pdgId())==24){ // W Boson
	  if(fabs(p->pdgId())==13){ // Muon
	    mu=&*p;
	  }
	  if(fabs(p->pdgId())==14){ // Muon neutrino
            mnu=&*p;
	  }
	  wbos = &*p;
	  break;
	}
      }else{
	continue; 
      }

    }
  }
  for (unsigned int ih = 0 ; ih < 222; ih++){
    int pdfset = 1001;
    if (ih<9) {
      pdfset=1001+ih; //pdf sets for direct scale variation
    }else if (ih>=9  && ih< 111){
      pdfset=2000+ih-8; 
    }else if (ih>=111 && ih < 166){
      pdfset=3000+ih-110;
    } else if (ih>=166 && ih<222){
      pdfset=4000+ih-165;
    } else{
      cout << "ERROR: non existing pdfset" << endl;
    }
    std::string whichWeightId = std::to_string(pdfset);
    if (DEBUG) cout << __LINE__ << endl;
    double weight = theWeight;
    for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
      if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
    }
    if (DEBUG) cout << __LINE__ << endl;
    hWeights[ih] ->Fill (1, weight);
    if (DEBUG) cout << __LINE__ << endl;
    
  }
  
  int ngenjet=0;
  for (GenJetCollection::const_iterator jet = genjets->begin(); jet!=genjets->end(); jet++){
    if (jet->pt()< jetPtCut) continue; 
    if (abs(jet->eta())>2.4) continue ;
    if (DEBUG) cout<< __LINE__ << " " <<  jet->mother(0) << endl;
    if(mu!=0){
    	if ( reco::deltaR(jet->eta(),jet->phi(),mu->eta(),mu->phi()) < 0.4) continue;
    }

     ngenjet++;
  }
  //cout << "number of gen jet " << ngenjet << endl;
  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(tok_jet_,jets);
  int njet=0;
  njet=ngenjet;
  std::cout << "njets: " <<  njet << std::endl;
  
  if (DEBUG) cout << __LINE__ << endl;
  if (wbos == 0){
     return;
  }


  if (DEBUG) cout << __LINE__ << endl;

  if (DEBUG) cout << __LINE__ << endl;
  if (DEBUG) cout << __LINE__ << " " << njet << endl;
  if (mu==0 ) return;
  //if (mnu==0) return;
  if ( abs(mu->eta())> 2.4) return;  
  if (DEBUG) cout << __LINE__ << endl;

  if (DEBUG) cout << __LINE__ << endl;
  
  if (mu!=0 && mnu!=0){
    if((doExclusive && njet==0) || (!doExclusive && njet>=0)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	std::string whichWeightId = std::to_string(pdfset); //get weights that correspond to particular pdfset
	if (DEBUG) cout << __LINE__ << endl;
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
        double dPhi_0J = mu->phi()-mnu->phi();
	double mMt_0J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_0J));
        double mMt_0J = std::sqrt(mMt_0J2);
	wMNuPt_vs_MuPt_0j[ih] ->Fill (mMt_0J, mu->pt(), weight); //fill 2D histogram of kinematic variables
	if (DEBUG) cout << __LINE__ << endl;

      }
    }

    if((doExclusive && njet==1) || (!doExclusive && njet>=1)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	std::string whichWeightId = std::to_string(pdfset);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << " "<< ih<< " " << wMNuPt_vs_MuPt_1j.size() << " " << EvtHandle->weights().size() << endl;
        double dPhi_1J = mu->phi()-mnu->phi();
	double mMt_1J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_1J));
        double mMt_1J = std::sqrt(mMt_1J2);
	wMNuPt_vs_MuPt_1j[ih] ->Fill (mMt_1J, mu->pt(), weight);
	if (DEBUG) cout << __LINE__ << endl;
      }
      if (DEBUG) cout << __LINE__ << endl;

    }

    if((doExclusive && njet==2) || (!doExclusive && njet>=2)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
        int pdfset = 1001;
        if (ih<9) {
          pdfset=1001+ih;
        }else if (ih>=9  && ih< 111){
          pdfset=2000+ih-8;
        }else if (ih>=111 && ih < 166){
          pdfset=3000+ih-110;
        } else if (ih>=166 && ih<222){
          pdfset=4000+ih-165;
        } else{
          cout << "ERROR: non existing pdfset" << endl;
        }
        if (DEBUG) cout << __LINE__ << endl;

        std::string whichWeightId = std::to_string(pdfset);
        double weight = theWeight;
        for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
          if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
        }
        if (DEBUG) cout << __LINE__ << endl;
        double dPhi_2J = mu->phi()-mnu->phi();
        double mMt_2J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_2J));
        double mMt_2J = std::sqrt(mMt_2J2);
        wMNuPt_vs_MuPt_2j[ih] ->Fill (mMt_2J, mu->pt(), weight);
      }
    }

    if (DEBUG) cout << __LINE__ << endl;
    if((doExclusive && njet==3) || (!doExclusive && njet>=3)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
	int pdfset = 1001;
	if (ih<9) {
	  pdfset=1001+ih;
	}else if (ih>=9  && ih< 111){
	  pdfset=2000+ih-8; 
	}else if (ih>=111 && ih < 166){
	  pdfset=3000+ih-110;
	} else if (ih>=166 && ih<222){
	  pdfset=4000+ih-165;
	} else{
	  cout << "ERROR: non existing pdfset" << endl;
	}
	if (DEBUG) cout << __LINE__ << endl;

	std::string whichWeightId = std::to_string(pdfset);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
        double dPhi_3J = mu->phi()-mnu->phi();
	double mMt_3J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_3J));
        double mMt_3J = std::sqrt(mMt_3J2);
	wMNuPt_vs_MuPt_3j[ih] ->Fill (mMt_3J, mu->pt(), weight);
      }
    }

    if((doExclusive && njet==4) || (!doExclusive && njet>=4)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
        int pdfset = 1001;
        if (ih<9) {
          pdfset=1001+ih;
        }else if (ih>=9  && ih< 111){
          pdfset=2000+ih-8;
        }else if (ih>=111 && ih < 166){
          pdfset=3000+ih-110;
        } else if (ih>=166 && ih<222){
          pdfset=4000+ih-165;
        } else{
          cout << "ERROR: non existing pdfset" << endl;
        }
        if (DEBUG) cout << __LINE__ << endl;

        std::string whichWeightId = std::to_string(pdfset);
        double weight = theWeight;
        for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
          if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
        }
        if (DEBUG) cout << __LINE__ << endl;
        double dPhi_4J = mu->phi()-mnu->phi();
        double mMt_4J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_4J));
        double mMt_4J = std::sqrt(mMt_4J2);
        wMNuPt_vs_MuPt_4j[ih] ->Fill (mMt_4J, mu->pt(), weight);
      }
    }

    if((doExclusive && njet==5) || (!doExclusive && njet>=5)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
        int pdfset = 1001;
        if (ih<9) {
          pdfset=1001+ih;
        }else if (ih>=9  && ih< 111){
          pdfset=2000+ih-8;
        }else if (ih>=111 && ih < 166){
          pdfset=3000+ih-110;
        } else if (ih>=166 && ih<222){
          pdfset=4000+ih-165;
        } else{
          cout << "ERROR: non existing pdfset" << endl;
        }
        if (DEBUG) cout << __LINE__ << endl;

        std::string whichWeightId = std::to_string(pdfset);
        double weight = theWeight;
        for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
          if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
        }
        if (DEBUG) cout << __LINE__ << endl;
        double dPhi_5J = mu->phi()-mnu->phi();
        double mMt_5J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_5J));
        double mMt_5J = std::sqrt(mMt_5J2);
        wMNuPt_vs_MuPt_5j[ih] ->Fill (mMt_5J, mu->pt(), weight);
      }
    }

    if((doExclusive && njet==6) || (!doExclusive && njet>=6)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
        int pdfset = 1001;
        if (ih<9) {
          pdfset=1001+ih;
        }else if (ih>=9  && ih< 111){
          pdfset=2000+ih-8;
        }else if (ih>=111 && ih < 166){
          pdfset=3000+ih-110;
        } else if (ih>=166 && ih<222){
          pdfset=4000+ih-165;
        } else{
          cout << "ERROR: non existing pdfset" << endl;
        }
        if (DEBUG) cout << __LINE__ << endl;

        std::string whichWeightId = std::to_string(pdfset);
        double weight = theWeight;
        for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
          if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
        }
        if (DEBUG) cout << __LINE__ << endl;
        double dPhi_6J = mu->phi()-mnu->phi();
        double mMt_6J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_6J));
        double mMt_6J = std::sqrt(mMt_6J2);
        wMNuPt_vs_MuPt_6j[ih] ->Fill (mMt_6J, mu->pt(), weight);
      }
    }


    if((doExclusive && njet==7) || (!doExclusive && njet>=7)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
        int pdfset = 1001;
        if (ih<9) {
          pdfset=1001+ih;
        }else if (ih>=9  && ih< 111){
          pdfset=2000+ih-8;
        }else if (ih>=111 && ih < 166){
          pdfset=3000+ih-110;
        } else if (ih>=166 && ih<222){
          pdfset=4000+ih-165;
        } else{
          cout << "ERROR: non existing pdfset" << endl;
        }
        if (DEBUG) cout << __LINE__ << endl;

        std::string whichWeightId = std::to_string(pdfset);
        double weight = theWeight;
        for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
          if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
        }
        if (DEBUG) cout << __LINE__ << endl;
        double dPhi_7J = mu->phi()-mnu->phi();
        double mMt_7J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_7J));
        double mMt_7J = std::sqrt(mMt_7J2);
        wMNuPt_vs_MuPt_7j[ih] ->Fill (mMt_7J, mu->pt(), weight);
      }
    }


    if((doExclusive && njet==8) || (!doExclusive && njet>=8)) {
      for (unsigned int ih = 0 ; ih < 222; ih++){
        int pdfset = 1001;
        if (ih<9) {
          pdfset=1001+ih;
        }else if (ih>=9  && ih< 111){
          pdfset=2000+ih-8;
        }else if (ih>=111 && ih < 166){
          pdfset=3000+ih-110;
        } else if (ih>=166 && ih<222){
          pdfset=4000+ih-165;
        } else{
          cout << "ERROR: non existing pdfset" << endl;
        }
        if (DEBUG) cout << __LINE__ << endl;

        std::string whichWeightId = std::to_string(pdfset);
        double weight = theWeight;
        for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
          if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
        }
        if (DEBUG) cout << __LINE__ << endl;
        double dPhi_8J = mu->phi()-mnu->phi();
        double mMt_8J2 = 2*mu->pt()*mnu->pt()*(1-cos(dPhi_8J));
        double mMt_8J = std::sqrt(mMt_8J2);
        wMNuPt_vs_MuPt_8j[ih] ->Fill (mMt_8J, mu->pt(), weight);
      }
    }


			 
    
    if (DEBUG) cout << __LINE__ << " " << njet << endl;
    if (mu->pt() < 25) return; //pt cut
    double dPhi = mu->phi()-mnu->phi();
    if ((2*mu->pt()*mnu->pt()*(1-cos(dPhi))) < 2500 ) return; //Mt cut
    //if (tau->mt() > 50 ) return; 
    wMuPt->Fill(mu->pt());
    wMuMt->Fill(mu->mt());
    wMNuPt->Fill(mnu->pt());
    wMNuMt->Fill(mnu->mt());
    wMass->Fill((mu->p4()+mnu->p4()).M());
    for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
      std::string whichWeightId = std::to_string(1001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	//if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
      }
      if (DEBUG) cout << __LINE__ << endl;
      wMass_scale[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); //fill with invariant mass, should be W Mass
    }
    for (unsigned int ih = 0 ; ih < 56; ih++){
      std::string whichWeightId = std::to_string(4002+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      wMass_pdf[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
    }
    for (unsigned int ih = 0 ; ih < 55; ih++){
      std::string whichWeightId = std::to_string(3001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      wMass_pdf[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
    }
    for (unsigned int ih = 0 ; ih < 102; ih++){
      std::string whichWeightId = std::to_string(2001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      wMass_pdf[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
    }

    if ((doExclusive && njet==0) || (!doExclusive && njet>=0)){
      //if (mnu->mt() > 50 ) return;
      //if (reco::deltaPhi(mnu->phi() , mu->phi())<2.7)return;
      wMuPt_0j->Fill(mu->pt());
      wMuMt_0j->Fill(mu->mt());
      wMNuPt_0j->Fill(mnu->pt());
      wMNuMt_0j->Fill(mnu->mt());
      wMass_0j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_0j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_0j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_0j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_0j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==0) || (!doExclusive && njet>=0)
    if ((doExclusive && njet==1) || (!doExclusive && njet>=1)){
      wMuPt_1j->Fill(mu->pt());
      wMuMt_1j->Fill(mu->mt());
      wMNuPt_1j->Fill(mnu->pt());
      wMNuMt_1j->Fill(mnu->mt());
      wMass_1j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_1j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_1j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_1j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_1j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==1) || (!doExclusive && njet>=1)

    if ((doExclusive && njet==2) || (!doExclusive && njet>=2)){
      wMuPt_2j->Fill(mu->pt());
      wMuMt_2j->Fill(mu->mt());
      wMNuPt_2j->Fill(mnu->pt());
      wMNuMt_2j->Fill(mnu->mt());
      wMass_2j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_2j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_2j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_2j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_2j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==2) || (!doExclusive && njet>=2)


    if ((doExclusive && njet==3) || (!doExclusive && njet>=3)){
      wMuPt_3j->Fill(mu->pt());
      wMuMt_3j->Fill(mu->mt());
      wMNuPt_3j->Fill(mnu->pt());
      wMNuMt_3j->Fill(mnu->mt());
      wMass_3j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_3j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_3j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_3j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_3j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==3) || (!doExclusive && njet>=3)

    if ((doExclusive && njet==4) || (!doExclusive && njet>=4)){
      wMuPt_4j->Fill(mu->pt());
      wMuMt_4j->Fill(mu->mt());
      wMNuPt_4j->Fill(mnu->pt());
      wMNuMt_4j->Fill(mnu->mt());
      wMass_4j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_4j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_4j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_4j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_4j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==4) || (!doExclusive && njet>=4)

    if ((doExclusive && njet==5) || (!doExclusive && njet>=5)){
      wMuPt_5j->Fill(mu->pt());
      wMuMt_5j->Fill(mu->mt());
      wMNuPt_5j->Fill(mnu->pt());
      wMNuMt_5j->Fill(mnu->mt());
      wMass_5j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_5j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_5j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_5j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_5j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==5) || (!doExclusive && njet>=5)

    if ((doExclusive && njet==6) || (!doExclusive && njet>=6)){
      wMuPt_6j->Fill(mu->pt());
      wMuMt_6j->Fill(mu->mt());
      wMNuPt_6j->Fill(mnu->pt());
      wMNuMt_6j->Fill(mnu->mt());
      wMass_6j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_6j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_6j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_6j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_6j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==6) || (!doExclusive && njet>=6)


    if ((doExclusive && njet==7) || (!doExclusive && njet>=7)){
      wMuPt_7j->Fill(mu->pt());
      wMuMt_7j->Fill(mu->mt());
      wMNuPt_7j->Fill(mnu->pt());
      wMNuMt_7j->Fill(mnu->mt());
      wMass_7j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_7j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_7j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_7j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_7j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==7) || (!doExclusive && njet>=7)


    if ((doExclusive && njet==8) || (!doExclusive && njet>=8)){
      wMuPt_8j->Fill(mu->pt());
      wMuMt_8j->Fill(mu->mt());
      wMNuPt_8j->Fill(mnu->pt());
      wMNuMt_8j->Fill(mnu->mt());
      wMass_8j->Fill((mu->p4()+mnu->p4()).M());
      for (unsigned int ih = 0 ; ih < wMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	wMass_scale_8j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_8j[ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_8j[56+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	wMass_pdf_8j[111+ih]->Fill((mu->p4()+mnu->p4()).M(), weight); 
      }
    }//(doExclusive && njet==8) || (!doExclusive && njet>=8)
    

  }//hadDecay
  if (DEBUG) cout << __LINE__ << endl;

  
  



   
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenSystInclusive::beginJob()
{

  std::stringstream name;
  
  
  wMass=fs->make<TH1F> ("WMass",  "WMass", 180, 0, 180); 
  wMass_0j=fs->make<TH1F> ("WMass_0j",  "WMass_0j", 180, 0, 180); 
  wMass_1j=fs->make<TH1F> ("WMass_1j",  "WMass_1j", 180, 0, 180); 
  wMass_2j=fs->make<TH1F> ("WMass_2j",  "WMass_2j", 180, 0, 180); 
  wMass_3j=fs->make<TH1F> ("WMass_3j",  "WMass_3j", 180, 0, 180);
  wMass_4j=fs->make<TH1F> ("WMass_4j",  "WMass_4j", 180, 0, 180);
  wMass_5j=fs->make<TH1F> ("WMass_5j",  "WMass_5j", 180, 0, 180);
  wMass_6j=fs->make<TH1F> ("WMass_6j",  "WMass_6j", 180, 0, 180);
  wMass_7j=fs->make<TH1F> ("WMass_7j",  "WMass_7j", 180, 0, 180);
  wMass_8j=fs->make<TH1F> ("WMass_8j",  "WMass_8j", 180, 0, 180);
  for (unsigned int ih = 0; ih < 9 ; ih++) {
    name.str("");
    name << "WMass_scale_" << 1001+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale.push_back(histo);
    name.str("");
    name << "WMass_scale_0j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_0j.push_back(histo);
    name.str("");
    name << "WMass_scale_1j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_1j.push_back(histo);
    name.str("");
    name << "WMass_scale_2j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_2j.push_back(histo);
    name.str("");
    name << "WMass_scale_3j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_3j.push_back(histo);
    name.str("");
    name << "WMass_scale_4j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_4j.push_back(histo);
    name.str("");
    name << "WMass_scale_5j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_5j.push_back(histo);
    name.str("");
    name << "WMass_scale_6j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_6j.push_back(histo);
    name.str("");
    name << "WMass_scale_7j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_7j.push_back(histo);
    name.str("");
    name << "WMass_scale_8j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_scale_8j.push_back(histo);
    name.str("");


    name.str("");
    name<< "hWeights_"<<1001+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    name.str("");
    name << "wMNuPt_vs_MuPt_0j" << 1001+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_0j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_1j" << 1001+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_1j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_2j" << 1001+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_2j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_3j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_3j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_4j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_4j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_5j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_5j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_6j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_6j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_7j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_7j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_8j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_8j.push_back(histo2D);
    
   
  }
  for (unsigned int ih = 1; ih < 57 ; ih++) {
    name.str("");
    name << "WMass_pdf_" << 4000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf.push_back(histo);
    name.str("");
    name << "WMass_pdf_0j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_0j.push_back(histo);
    name.str("");
    name << "WMass_pdf_1j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_1j.push_back(histo);
    name.str("");
    name << "WMass_pdf_2j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_2j.push_back(histo);
    name.str("");
    name << "WMass_pdf_3j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_3j.push_back(histo);
    name.str("");
    name << "WMass_pdf_4j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_4j.push_back(histo);
    name.str("");
    name << "WMass_pdf_5j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_5j.push_back(histo);
    name.str("");
    name << "WMass_pdf_6j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_6j.push_back(histo);
    name.str("");
    name << "WMass_pdf_7j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_7j.push_back(histo);
    name.str("");
    name << "WMass_pdf_8j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_8j.push_back(histo);
    name.str("");



    name.str("");
    name<< "hWeights_"<<4000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    name.str("");
    name << "wMNuPt_vs_MuPt_0j" << 4000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_0j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_1j" << 4000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_1j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_2j" << 4000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_2j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_3j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_3j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_4j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_4j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_5j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_5j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_6j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_6j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_7j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_7j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_8j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_8j.push_back(histo2D);
    name.str("");
   
  }
  for (unsigned int ih = 1; ih < 56 ; ih++) {
    name.str("");
    name << "WMass_pdf_" <<3000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf.push_back(histo);
    name.str("");
    name << "WMass_pdf_0j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_0j.push_back(histo);
    name.str("");
    name << "WMass_pdf_1j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_1j.push_back(histo);
    name.str("");
    name << "WMass_pdf_2j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_2j.push_back(histo);
    name.str("");
    name << "WMass_pdf_3j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_3j.push_back(histo);
    name.str("");
    name << "WMass_pdf_4j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_4j.push_back(histo);
    name.str("");
    name << "WMass_pdf_5j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_5j.push_back(histo);
    name.str("");
    name << "WMass_pdf_6j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_6j.push_back(histo);
    name.str("");
    name << "WMass_pdf_7j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_7j.push_back(histo);
    name.str("");
    name << "WMass_pdf_8j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_8j.push_back(histo);
    name.str("");


    name.str("");
    name<< "hWeights_"<<3000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    name.str("");
    name << "wMNuPt_vs_MuPt_0j" << 3000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_0j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_1j" << 3000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_1j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_2j" << 3000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_2j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_3j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_3j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_4j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_4j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_5j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_5j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_6j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_6j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_7j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_7j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_8j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_8j.push_back(histo2D);
    name.str("");
   
  }
  for (unsigned int ih = 1; ih < 103 ; ih++) {
    name.str("");
    name << "WMass_pdf_" <<2000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf.push_back(histo);
    name.str("");
    name << "WMass_pdf_0j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_0j.push_back(histo);
    name.str("");
    name << "WMass_pdf_1j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_1j.push_back(histo);
    name.str("");
    name << "WMass_pdf_2j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_2j.push_back(histo);
    name.str("");
    name << "WMass_pdf_3j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_3j.push_back(histo);
    name.str("");
    name << "WMass_pdf_4j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_4j.push_back(histo);
    name.str("");
    name << "WMass_pdf_5j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_5j.push_back(histo);
    name.str("");
    name << "WMass_pdf_6j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_6j.push_back(histo);
    name.str("");
    name << "WMass_pdf_7j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_7j.push_back(histo);
    name.str("");
    name << "WMass_pdf_8j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    wMass_pdf_8j.push_back(histo);
    name.str("");


    name.str("");
    name<< "hWeights_"<<2000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    name.str("");

    name.str("");
    name << "wMNuPt_vs_MuPt_0j" << 2000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_0j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_1j" << 2000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_1j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_2j" << 2000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_2j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_3j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_3j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_4j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_4j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_5j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_5j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_6j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_6j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_7j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_7j.push_back(histo2D);
    name.str("");
    name << "wMNuPt_vs_MuPt_8j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    wMNuPt_vs_MuPt_8j.push_back(histo2D);
   
  }
  wMuPt = fs->make<TH1F> ("wMuPt", "wMuPt", 100, 0, 100 ) ; 
  wMuMt = fs->make<TH1F> ("wMuMt", "wMuMt", 100, 0, 100 ) ; 
  wMNuPt= fs->make<TH1F> ("wMNuPt", "wMNuPt", 100, 0, 100 ) ; 
  wMNuMt= fs->make<TH1F> ("wMNuMt", "wMNuMt", 100, 0, 100 ) ; 
  wMuPt_0j = fs->make<TH1F> ("wMuPt_0j", "wMuPt_0j", 100, 0, 100 ) ; 
  wMuMt_0j = fs->make<TH1F> ("wMuMt_0j", "wMuMt_0j", 100, 0, 100 ) ; 
  wMNuPt_0j= fs->make<TH1F> ("wMNuPt_0j", "wMNuPt_0j", 100, 0, 100 ) ; 
  wMNuMt_0j= fs->make<TH1F> ("wMNuMt_0j", "wMNuMt_0j", 100, 0, 100 ) ; 
  wMuPt_1j = fs->make<TH1F> ("wMuPt_1j", "wMuPt_1j", 100, 0, 100 ) ; 
  wMuMt_1j = fs->make<TH1F> ("wMuMt_1j", "wMuMt_1j", 100, 0, 100 ) ; 
  wMNuPt_1j= fs->make<TH1F> ("wMNuPt_1j", "wMNuPt_1j", 100, 0, 100 ) ; 
  wMNuMt_1j= fs->make<TH1F> ("wMNuMt_1j", "wMNuMt_1j", 100, 0, 100 ) ; 
  wMuPt_2j = fs->make<TH1F> ("wMuPt_2j", "wMuPt_2j", 100, 0, 100 ) ; 
  wMuMt_2j = fs->make<TH1F> ("wMuMt_2j", "wMuMt_2j", 100, 0, 100 ) ; 
  wMNuPt_2j= fs->make<TH1F> ("wMNuPt_2j", "wMNuPt_2j", 100, 0, 100 ) ; 
  wMNuMt_2j= fs->make<TH1F> ("wMNuMt_2j", "wMNuMt_2j", 100, 0, 100 ) ; 
  wMuPt_3j = fs->make<TH1F> ("wMuPt_3j", "wMuPt_3j", 100, 0, 100 ) ;
  wMuMt_3j = fs->make<TH1F> ("wMuMt_3j", "wMuMt_3j", 100, 0, 100 ) ;
  wMNuPt_3j= fs->make<TH1F> ("wMNuPt_3j", "wMNuPt_3j", 100, 0, 100 ) ;
  wMNuMt_3j= fs->make<TH1F> ("wMNuMt_3j", "wMNuMt_3j", 100, 0, 100 ) ;
  wMuPt_4j = fs->make<TH1F> ("wMuPt_4j", "wMuPt_4j", 100, 0, 100 ) ;
  wMuMt_4j = fs->make<TH1F> ("wMuMt_4j", "wMuMt_4j", 100, 0, 100 ) ;
  wMNuPt_4j= fs->make<TH1F> ("wMNuPt_4j", "wMNuPt_4j", 100, 0, 100 ) ;
  wMNuMt_4j= fs->make<TH1F> ("wMNuMt_4j", "wMNuMt_4j", 100, 0, 100 ) ;
  wMuPt_5j = fs->make<TH1F> ("wMuPt_5j", "wMuPt_5j", 100, 0, 100 ) ;
  wMuMt_5j = fs->make<TH1F> ("wMuMt_5j", "wMuMt_5j", 100, 0, 100 ) ;
  wMNuPt_5j= fs->make<TH1F> ("wMNuPt_5j", "wMNuPt_5j", 100, 0, 100 ) ;
  wMNuMt_5j= fs->make<TH1F> ("wMNuMt_5j", "wMNuMt_5j", 100, 0, 100 ) ;
  wMuPt_6j = fs->make<TH1F> ("wMuPt_6j", "wMuPt_6j", 100, 0, 100 ) ;
  wMuMt_6j = fs->make<TH1F> ("wMuMt_6j", "wMuMt_6j", 100, 0, 100 ) ;
  wMNuPt_6j= fs->make<TH1F> ("wMNuPt_6j", "wMNuPt_6j", 100, 0, 100 ) ;
  wMNuMt_6j= fs->make<TH1F> ("wMNuMt_6j", "wMNuMt_6j", 100, 0, 100 ) ;
  wMuPt_7j = fs->make<TH1F> ("wMuPt_7j", "wMuPt_7j", 100, 0, 100 ) ;
  wMuMt_7j = fs->make<TH1F> ("wMuMt_7j", "wMuMt_7j", 100, 0, 100 ) ;
  wMNuPt_7j= fs->make<TH1F> ("wMNuPt_7j", "wMNuPt_7j", 100, 0, 100 ) ;
  wMNuMt_7j= fs->make<TH1F> ("wMNuMt_7j", "wMNuMt_7j", 100, 0, 100 ) ;
  wMuPt_8j = fs->make<TH1F> ("wMuPt_8j", "wMuPt_8j", 100, 0, 100 ) ;
  wMuMt_8j = fs->make<TH1F> ("wMuMt_8j", "wMuMt_8j", 100, 0, 100 ) ;
  wMNuPt_8j= fs->make<TH1F> ("wMNuPt_8j", "wMNuPt_8j", 100, 0, 100 ) ;
  wMNuMt_8j= fs->make<TH1F> ("wMNuMt_8j", "wMNuMt_8j", 100, 0, 100 ) ;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenSystInclusive::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenSystInclusive::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
GenSystInclusive::beginRun( edm::Run const & iRun,  edm::EventSetup const& iSetup) {
  nevent_run = 0;
}
void
GenSystInclusive::endRun( const edm::Run& iRun, const edm::EventSetup& iSetup) {
  std::cout <<" Runnumber "<<iRun.run()<<" Nevents  "<<nevent_run;
  if (DEBUG) std::cout << __LINE__ << std::endl;
     edm::Handle<LHERunInfoProduct> run; 
     typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
     
  if (DEBUG) std::cout << __LINE__ << std::endl;
  iRun.getByToken(lheRunInfoToken_, run );
     LHERunInfoProduct myLHERunInfoProduct = *(run.product());
     
    if (DEBUG) std::cout << __LINE__ << std::endl;
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
       std::cout << iter->tag() << std::endl;
       std::vector<std::string> lines = iter->lines();
       for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
  	 std::cout << lines.at(iLine);
       }
     }
}
 
double GenSystInclusive::dPhi(math::XYZTLorentzVector& vec1, math::XYZTLorentzVector& vec2) {
  return reco::deltaPhi(vec1.phi(),vec2.phi());
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenSystInclusive);
