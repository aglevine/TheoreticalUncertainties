// -*- C++ -*-
//
// Package:    LFVAnalysis/GenSystInclusiveZ
// Class:      GenSystInclusiveZ
// 
// Original Author:  Silvia Taroni
// Modified by Aaron Levine for Z+Jets
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

class GenSystInclusiveZ : public edm::one::EDAnalyzer<edm::one::WatchRuns>  {
   public:
      explicit GenSystInclusiveZ(const edm::ParameterSet&);
      ~GenSystInclusiveZ();

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
      double jetPtCut;
      bool doExclusive;
       
      int nevent_run;  

      edm::Service<TFileService> fs;      
      TH1F * zMass; 
      std::vector<TH1F *> zMass_scale;
      std::vector<TH1F *> zMass_pdf;
      TH1F * zMass_0j; 
      std::vector<TH1F *> zMass_scale_0j;
      std::vector<TH1F *> zMass_pdf_0j;
      TH1F * zMass_1j; 
      std::vector<TH1F *> zMass_scale_1j;
      std::vector<TH1F *> zMass_pdf_1j;
      TH1F * zMass_2j; 
      std::vector<TH1F *> zMass_scale_2j;
      std::vector<TH1F *> zMass_pdf_2j;
      TH1F * zMass_3j;
      std::vector<TH1F *> zMass_scale_3j;
      std::vector<TH1F *> zMass_pdf_3j;
      TH1F * zMass_4j;
      std::vector<TH1F *> zMass_scale_4j;
      std::vector<TH1F *> zMass_pdf_4j;
      TH1F * zMass_5j;
      std::vector<TH1F *> zMass_scale_5j;
      std::vector<TH1F *> zMass_pdf_5j;
      TH1F * zMass_6j;
      std::vector<TH1F *> zMass_scale_6j;
      std::vector<TH1F *> zMass_pdf_6j;
      TH1F * zMass_7j;
      std::vector<TH1F *> zMass_scale_7j;
      std::vector<TH1F *> zMass_pdf_7j;
      TH1F * zMass_8j;
      std::vector<TH1F *> zMass_scale_8j;
      std::vector<TH1F *> zMass_pdf_8j;



     std::vector<TH2F *>  zMuPt1_vs_MuPt2_0j; 
     std::vector<TH2F *>  zMuPt1_vs_MuPt2_1j; 
     std::vector<TH2F *>  zMuPt1_vs_MuPt2_2j; 
     std::vector<TH2F *>  zMuPt1_vs_MuPt2_3j;
     std::vector<TH2F *>  zMuPt1_vs_MuPt2_4j;
     std::vector<TH2F *>  zMuPt1_vs_MuPt2_5j;
     std::vector<TH2F *>  zMuPt1_vs_MuPt2_6j;
     std::vector<TH2F *>  zMuPt1_vs_MuPt2_7j;
     std::vector<TH2F *>  zMuPt1_vs_MuPt2_8j;

     std::vector<TH1F *> hWeights; 
   
      TH1F * zMuPt1; 
      TH1F * zMu1Mt;
      TH1F * zMuPt2;
      TH1F * zMu2Mt;
      TH1F * zMuPt1_0j; 
      TH1F * zMu1Mt_0j;
      TH1F * zMuPt2_0j;
      TH1F * zMu2Mt_0j;
      TH1F * zMuPt1_1j; 
      TH1F * zMu1Mt_1j;
      TH1F * zMuPt2_1j;
      TH1F * zMu2Mt_1j;
      TH1F * zMuPt1_2j; 
      TH1F * zMu1Mt_2j;
      TH1F * zMuPt2_2j;
      TH1F * zMu2Mt_2j;
      TH1F * zMuPt1_3j;
      TH1F * zMu1Mt_3j;
      TH1F * zMuPt2_3j;
      TH1F * zMu2Mt_3j;
      TH1F * zMuPt1_4j;
      TH1F * zMu1Mt_4j;
      TH1F * zMuPt2_4j;
      TH1F * zMu2Mt_4j;
      TH1F * zMuPt1_5j;
      TH1F * zMu1Mt_5j;
      TH1F * zMuPt2_5j;
      TH1F * zMu2Mt_5j;
      TH1F * zMuPt1_6j;
      TH1F * zMu1Mt_6j;
      TH1F * zMuPt2_6j;
      TH1F * zMu2Mt_6j;
      TH1F * zMuPt1_7j;
      TH1F * zMu1Mt_7j;
      TH1F * zMuPt2_7j;
      TH1F * zMu2Mt_7j;
      TH1F * zMuPt1_8j;
      TH1F * zMu1Mt_8j;
      TH1F * zMuPt2_8j;
      TH1F * zMu2Mt_8j;

    
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
GenSystInclusiveZ::GenSystInclusiveZ(const edm::ParameterSet& iConfig)

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


GenSystInclusiveZ::~GenSystInclusiveZ()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenSystInclusiveZ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
 
  const reco::Candidate* zbos  = 0;
  const reco::Candidate* mu1     = 0;
  const reco::Candidate* mu2     = 0;

  for (pat::PackedGenParticleCollection::const_iterator p = genParticles->begin();p != genParticles->end(); ++p ) {
    const reco::Candidate *mother = 0;
    if( p->numberOfMothers()!=0) mother= p->mother(0);
    //if (DEBUG) cout<< __LINE__ << " " << p->pdgId()<< " " << p->status() << " ";
    for (unsigned int im =0; im<10 ; im++){
      if (mother->numberOfMothers()!=0 ) {
	mother=&*mother->mother(0);
	//cout << mother->pdgId() << " " ;
	if (abs(mother->pdgId())==23){ //Select Z->mumu
	  if(p->pdgId()==13){
	    mu1=&*p;
	  }
	  if(p->pdgId()==-13){
            mu2=&*p;
	  }
	  zbos = &*p;
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
  int njet=0;
  
  for (GenJetCollection::const_iterator jet = genjets->begin(); jet!=genjets->end(); jet++){
    if (jet->pt() < jetPtCut) continue;
    if (abs(jet->eta())>2.4) continue;
    if (DEBUG) cout<< __LINE__ << " " <<  jet->mother(0) << endl;
    if(mu1!=0){
    	if ( reco::deltaR(jet->eta(),jet->phi(),mu1->eta(),mu1->phi()) < 0.4) continue;
    }
    if(mu2!=0){
        if ( reco::deltaR(jet->eta(),jet->phi(),mu2->eta(),mu2->phi()) < 0.4) continue;
    }

    
    ngenjet++;
    
   
    
    
  }
  njet = ngenjet;
  std::cout << "njet: " << njet << std::endl;
  //cout << "number of gen jet " << ngenjet << endl;
  
  if (DEBUG) cout << __LINE__ << endl;
  if (zbos == 0){
     return;
  }


  if (DEBUG) cout << __LINE__ << endl;

  if (DEBUG) cout << __LINE__ << endl;
  if (DEBUG) cout << __LINE__ << " " << njet << endl;
  if (mu1==0 ) return;
  if (mu2==0) return;
  if ( abs(mu1->eta())> 2.4) return;  
  if ( abs(mu2->eta())> 2.4) return;
  if (DEBUG) cout << __LINE__ << endl;

  if (DEBUG) cout << __LINE__ << endl;
  
  if (mu1!=0 && mu2!=0){

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
	std::string whichWeightId = std::to_string(pdfset);
	if (DEBUG) cout << __LINE__ << endl;
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
        double dPhi_0J = mu1->phi()-mu2->phi();
	double mMt_0J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_0J));
        double mMt_0J = std::sqrt(mMt_0J2);
	zMuPt1_vs_MuPt2_0j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
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
	if (DEBUG) cout << __LINE__ << " "<< ih<< " " << zMuPt1_vs_MuPt2_1j.size() << " " << EvtHandle->weights().size() << endl;
        double dPhi_1J = mu1->phi()-mu2->phi();
	double mMt_1J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_1J));
        double mMt_1J = std::sqrt(mMt_1J2);
	zMuPt1_vs_MuPt2_1j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
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
        double dPhi_2J = mu1->phi()-mu2->phi();
        double mMt_2J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_2J));
        double mMt_2J = std::sqrt(mMt_2J2);
        zMuPt1_vs_MuPt2_2j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
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
        double dPhi_3J = mu1->phi()-mu2->phi();
	double mMt_3J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_3J));
        double mMt_3J = std::sqrt(mMt_3J2);
	zMuPt1_vs_MuPt2_3j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
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
        double dPhi_4J = mu1->phi()-mu2->phi();
        double mMt_4J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_4J));
        double mMt_4J = std::sqrt(mMt_4J2);
        zMuPt1_vs_MuPt2_4j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
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
        double dPhi_5J = mu1->phi()-mu2->phi();
        double mMt_5J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_5J));
        double mMt_5J = std::sqrt(mMt_5J2);
        zMuPt1_vs_MuPt2_5j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
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
        double dPhi_6J = mu1->phi()-mu2->phi();
        double mMt_6J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_6J));
        double mMt_6J = std::sqrt(mMt_6J2);
        zMuPt1_vs_MuPt2_6j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
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
        double dPhi_7J = mu1->phi()-mu2->phi();
        double mMt_7J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_7J));
        double mMt_7J = std::sqrt(mMt_7J2);
        zMuPt1_vs_MuPt2_7j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
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
        double dPhi_8J = mu1->phi()-mu2->phi();
        double mMt_8J2 = 2*mu1->pt()*mu2->pt()*(1-cos(dPhi_8J));
        double mMt_8J = std::sqrt(mMt_8J2);
        zMuPt1_vs_MuPt2_8j[ih] ->Fill (mu1->pt(),mu2->pt(), weight);
      }
    }
			 
    
    if (DEBUG) cout << __LINE__ << " " << njet << endl;
    if (mu1->pt() < 20) return;
    if (mu2->pt() < 20) return;
    //if (tau->mt() > 50 ) return; 
    zMuPt1->Fill(mu1->pt());
    zMu1Mt->Fill(mu1->mt());
    zMuPt2->Fill(mu2->pt());
    zMu2Mt->Fill(mu2->mt());
    zMass->Fill((mu1->p4()+mu2->p4()).M());
    for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
      std::string whichWeightId = std::to_string(1001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	//if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
      }
      if (DEBUG) cout << __LINE__ << endl;
      zMass_scale[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
    }
    for (unsigned int ih = 0 ; ih < 56; ih++){
      std::string whichWeightId = std::to_string(4002+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      zMass_pdf[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
    }
    for (unsigned int ih = 0 ; ih < 55; ih++){
      std::string whichWeightId = std::to_string(3001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      zMass_pdf[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
    }
    for (unsigned int ih = 0 ; ih < 102; ih++){
      std::string whichWeightId = std::to_string(2001+ih);
      double weight = theWeight;
      for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
      }
      zMass_pdf[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
    }

    if ((doExclusive && njet==0) || (!doExclusive && njet>=0)){
      //if (mu2->mt() > 50 ) return;
      //if (reco::deltaPhi(mu2->phi() , mu1->phi())<2.7)return;
      zMuPt1_0j->Fill(mu1->pt());
      zMu1Mt_0j->Fill(mu1->mt());
      zMuPt2_0j->Fill(mu2->pt());
      zMu2Mt_0j->Fill(mu2->mt());
      zMass_0j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_0j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_0j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_0j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_0j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==0) || (!doExclusive && njet>=0)
    if ((doExclusive && njet==1) || (!doExclusive && njet>=1)){
      zMuPt1_1j->Fill(mu1->pt());
      zMu1Mt_1j->Fill(mu1->mt());
      zMuPt2_1j->Fill(mu2->pt());
      zMu2Mt_1j->Fill(mu2->mt());
      zMass_1j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	  //if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl; 
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_1j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_1j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_1j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_1j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==1) || (!doExclusive && njet>=1)
    if ((doExclusive && njet==2) || (!doExclusive && njet>=2)){
      zMuPt1_2j->Fill(mu1->pt());
      zMu1Mt_2j->Fill(mu1->mt());
      zMuPt2_2j->Fill(mu2->pt());
      zMu2Mt_2j->Fill(mu2->mt());
      zMass_2j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_2j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_2j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_2j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_2j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==2) || (!doExclusive && njet>=2)


    if ((doExclusive && njet==3) || (!doExclusive && njet>=3)){
      zMuPt1_3j->Fill(mu1->pt());
      zMu1Mt_3j->Fill(mu1->mt());
      zMuPt2_3j->Fill(mu2->pt());
      zMu2Mt_3j->Fill(mu2->mt());
      zMass_3j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_3j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_3j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_3j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_3j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==3) || (!doExclusive && njet>=3)

    if ((doExclusive && njet==4) || (!doExclusive && njet>=4)){
      zMuPt1_4j->Fill(mu1->pt());
      zMu1Mt_4j->Fill(mu1->mt());
      zMuPt2_4j->Fill(mu2->pt());
      zMu2Mt_4j->Fill(mu2->mt());
      zMass_4j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_4j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_4j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_4j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_4j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==4) || (!doExclusive && njet>=4)

    if ((doExclusive && njet==5) || (!doExclusive && njet>=5)){
      zMuPt1_5j->Fill(mu1->pt());
      zMu1Mt_5j->Fill(mu1->mt());
      zMuPt2_5j->Fill(mu2->pt());
      zMu2Mt_5j->Fill(mu2->mt());
      zMass_5j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_5j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_5j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_5j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_5j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==5) || (!doExclusive && njet>=5)

    if ((doExclusive && njet==6) || (!doExclusive && njet>=6)){
      zMuPt1_6j->Fill(mu1->pt());
      zMu1Mt_6j->Fill(mu1->mt());
      zMuPt2_6j->Fill(mu2->pt());
      zMu2Mt_6j->Fill(mu2->mt());
      zMass_6j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_6j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_6j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_6j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_6j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==6) || (!doExclusive && njet>=6)


    if ((doExclusive && njet==7) || (!doExclusive && njet>=7)){
      zMuPt1_7j->Fill(mu1->pt());
      zMu1Mt_7j->Fill(mu1->mt());
      zMuPt2_7j->Fill(mu2->pt());
      zMu2Mt_7j->Fill(mu2->mt());
      zMass_7j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_7j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_7j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_7j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_7j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==7) || (!doExclusive && njet>=7)
    

    if ((doExclusive && njet==8) || (!doExclusive && njet>=8)){
      zMuPt1_8j->Fill(mu1->pt());
      zMu1Mt_8j->Fill(mu1->mt());
      zMuPt2_8j->Fill(mu2->pt());
      zMu2Mt_8j->Fill(mu2->mt());
      zMass_8j->Fill((mu1->p4()+mu2->p4()).M());
      for (unsigned int ih = 0 ; ih < zMass_scale.size(); ih++){
	std::string whichWeightId = std::to_string(1001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	if (DEBUG) cout << __LINE__ << endl;
	zMass_scale_8j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 56; ih++){
	std::string whichWeightId = std::to_string(4001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_8j[ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 55; ih++){
	std::string whichWeightId = std::to_string(3001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_8j[56+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
      for (unsigned int ih = 0 ; ih < 102; ih++){
	std::string whichWeightId = std::to_string(2001+ih);
	double weight = theWeight;
	for (unsigned int i=0; i<EvtHandle->weights().size(); i++) {
	  if (EvtHandle->weights()[i].id == whichWeightId) weight *= EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP();
	}
	zMass_pdf_8j[111+ih]->Fill((mu1->p4()+mu2->p4()).M(), weight); 
      }
    }//(doExclusive && njet==8) || (!doExclusive && njet>=8)


  }//hadDecay
  if (DEBUG) cout << __LINE__ << endl;

  
  



   
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenSystInclusiveZ::beginJob()
{

  std::stringstream name;
  
  
  zMass=fs->make<TH1F> ("ZMass",  "ZMass", 180, 0, 180); 
  zMass_0j=fs->make<TH1F> ("ZMass_0j",  "ZMass_0j", 180, 0, 180); 
  zMass_1j=fs->make<TH1F> ("ZMass_1j",  "ZMass_1j", 180, 0, 180); 
  zMass_2j=fs->make<TH1F> ("ZMass_2j",  "ZMass_2j", 180, 0, 180); 
  zMass_3j=fs->make<TH1F> ("ZMass_3j",  "ZMass_3j", 180, 0, 180);
  zMass_4j=fs->make<TH1F> ("ZMass_4j",  "ZMass_4j", 180, 0, 180);
  zMass_5j=fs->make<TH1F> ("ZMass_5j",  "ZMass_5j", 180, 0, 180);
  zMass_6j=fs->make<TH1F> ("ZMass_6j",  "ZMass_6j", 180, 0, 180);
  zMass_7j=fs->make<TH1F> ("ZMass_7j",  "ZMass_7j", 180, 0, 180);
  zMass_8j=fs->make<TH1F> ("ZMass_8j",  "ZMass_8j", 180, 0, 180);

  for (unsigned int ih = 0; ih < 9 ; ih++) {
    name.str("");
    name << "ZMass_scale_" << 1001+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale.push_back(histo);
    name.str("");
    name << "ZMass_scale_0j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_0j.push_back(histo);
    name.str("");
    name << "ZMass_scale_1j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_1j.push_back(histo);
    name.str("");
    name << "ZMass_scale_2j" << 1001+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_2j.push_back(histo);
    name.str("");
    name << "ZMass_scale_3j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_3j.push_back(histo);
    name.str("");
    name << "ZMass_scale_4j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_4j.push_back(histo);
    name.str("");
    name << "ZMass_scale_5j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_5j.push_back(histo);
    name.str("");
    name << "ZMass_scale_6j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_6j.push_back(histo);
    name.str("");
    name << "ZMass_scale_7j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_7j.push_back(histo);
    name.str("");
    name << "ZMass_scale_8j" << 1001+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_scale_8j.push_back(histo);
    name.str("");


    name.str("");
    name<< "hWeights_"<<1001+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    name.str("");
    name << "zMuPt1_vs_MuPt2_0j" << 1001+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_0j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_1j" << 1001+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_1j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_2j" << 1001+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_2j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_3j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_3j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_4j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_4j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_5j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_5j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_6j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_6j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_7j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_7j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_8j" << 1001+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_8j.push_back(histo2D);
    
   
  }
  for (unsigned int ih = 1; ih < 57 ; ih++) {
    name.str("");
    name << "ZMass_pdf_" << 4000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf.push_back(histo);
    name.str("");
    name << "ZMass_pdf_0j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_0j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_1j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_1j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_2j" << 4000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_2j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_3j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_3j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_4j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_4j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_5j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_5j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_6j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_6j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_7j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_7j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_8j" << 4000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_8j.push_back(histo);
    name.str("");


    name.str("");
    name<< "hWeights_"<<4000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    name.str("");
    name << "zMuPt1_vs_MuPt2_0j" << 4000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_0j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_1j" << 4000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_1j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_2j" << 4000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_2j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_3j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_3j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_4j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_4j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_5j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_5j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_6j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_6j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_7j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_7j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_8j" << 4000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_8j.push_back(histo2D);
    name.str("");
   
  }
  for (unsigned int ih = 1; ih < 56 ; ih++) {
    name.str("");
    name << "ZMass_pdf_" <<3000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf.push_back(histo);
    name.str("");
    name << "ZMass_pdf_0j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_0j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_1j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_1j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_2j" << 3000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_2j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_3j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_3j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_4j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_4j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_5j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_5j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_6j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_6j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_7j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_7j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_8j" << 3000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_8j.push_back(histo);
    name.str("");

    name.str("");
    name<< "hWeights_"<<3000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    name.str("");
    name << "zMuPt1_vs_MuPt2_0j" << 3000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_0j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_1j" << 3000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_1j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_2j" << 3000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_2j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_3j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_3j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_4j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_4j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_5j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_5j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_6j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_6j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_7j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_7j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_8j" << 3000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_8j.push_back(histo2D);
    name.str("");
   
  }
  for (unsigned int ih = 1; ih < 103 ; ih++) {
    name.str("");
    name << "ZMass_pdf_" <<2000+ih ; 
    TH1F* histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf.push_back(histo);
    name.str("");
    name << "ZMass_pdf_0j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_0j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_1j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_1j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_2j" <<2000+ih ; 
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_2j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_3j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_3j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_4j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_4j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_5j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_5j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_6j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_6j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_7j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_7j.push_back(histo);
    name.str("");
    name << "ZMass_pdf_8j" <<2000+ih ;
    histo = fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 180, 0, 180);
    zMass_pdf_8j.push_back(histo);
    name.str("");

    name.str("");
    name<< "hWeights_"<<2000+ih;
    histo=fs->make<TH1F> (name.str().c_str(),  name.str().c_str(), 3, 0, 3);
    hWeights.push_back(histo);
    name.str("");

    name.str("");
    name << "zMuPt1_vs_MuPt2_0j" << 2000+ih ; 
    TH2F* histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_0j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_1j" << 2000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_1j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_2j" << 2000+ih ; 
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_2j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_3j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_3j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_4j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_4j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_5j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_5j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_6j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_6j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_7j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_7j.push_back(histo2D);
    name.str("");
    name << "zMuPt1_vs_MuPt2_8j" << 2000+ih ;
    histo2D = fs->make<TH2F> (name.str().c_str(),  name.str().c_str(), 100, 0, 100, 100, 0, 100);
    zMuPt1_vs_MuPt2_8j.push_back(histo2D);
    name.str("");
   
  }
  zMuPt1 = fs->make<TH1F> ("zMuPt1", "zMuPt1", 100, 0, 100 ) ; 
  zMu1Mt = fs->make<TH1F> ("zMu1Mt", "zMu1Mt", 100, 0, 100 ) ; 
  zMuPt2= fs->make<TH1F> ("zMuPt2", "zMuPt2", 100, 0, 100 ) ; 
  zMu2Mt= fs->make<TH1F> ("zMu2Mt", "zMu2Mt", 100, 0, 100 ) ; 
  zMuPt1_0j = fs->make<TH1F> ("zMuPt1_0j", "zMuPt1_0j", 100, 0, 100 ) ; 
  zMu1Mt_0j = fs->make<TH1F> ("zMu1Mt_0j", "zMu1Mt_0j", 100, 0, 100 ) ; 
  zMuPt2_0j= fs->make<TH1F> ("zMuPt2_0j", "zMuPt2_0j", 100, 0, 100 ) ; 
  zMu2Mt_0j= fs->make<TH1F> ("zMu2Mt_0j", "zMu2Mt_0j", 100, 0, 100 ) ; 
  zMuPt1_1j = fs->make<TH1F> ("zMuPt1_1j", "zMuPt1_1j", 100, 0, 100 ) ; 
  zMu1Mt_1j = fs->make<TH1F> ("zMu1Mt_1j", "zMu1Mt_1j", 100, 0, 100 ) ; 
  zMuPt2_1j= fs->make<TH1F> ("zMuPt2_1j", "zMuPt2_1j", 100, 0, 100 ) ; 
  zMu2Mt_1j= fs->make<TH1F> ("zMu2Mt_1j", "zMu2Mt_1j", 100, 0, 100 ) ; 
  zMuPt1_2j = fs->make<TH1F> ("zMuPt1_2j", "zMuPt1_2j", 100, 0, 100 ) ; 
  zMu1Mt_2j = fs->make<TH1F> ("zMu1Mt_2j", "zMu1Mt_2j", 100, 0, 100 ) ; 
  zMuPt2_2j= fs->make<TH1F> ("zMuPt2_2j", "zMuPt2_2j", 100, 0, 100 ) ; 
  zMu2Mt_2j= fs->make<TH1F> ("zMu2Mt_2j", "zMu2Mt_2j", 100, 0, 100 ) ; 
  zMuPt1_3j = fs->make<TH1F> ("zMuPt1_3j", "zMuPt1_3j", 100, 0, 100 ) ;
  zMu1Mt_3j = fs->make<TH1F> ("zMu1Mt_3j", "zMu1Mt_3j", 100, 0, 100 ) ;
  zMuPt2_3j= fs->make<TH1F> ("zMuPt2_3j", "zMuPt2_3j", 100, 0, 100 ) ;
  zMu2Mt_3j= fs->make<TH1F> ("zMu2Mt_3j", "zMu2Mt_3j", 100, 0, 100 ) ;
  zMuPt1_4j = fs->make<TH1F> ("zMuPt1_4j", "zMuPt1_4j", 100, 0, 100 ) ;
  zMu1Mt_4j = fs->make<TH1F> ("zMu1Mt_4j", "zMu1Mt_4j", 100, 0, 100 ) ;
  zMuPt2_4j= fs->make<TH1F> ("zMuPt2_4j", "zMuPt2_4j", 100, 0, 100 ) ;
  zMu2Mt_4j= fs->make<TH1F> ("zMu2Mt_4j", "zMu2Mt_4j", 100, 0, 100 ) ;
  zMuPt1_5j = fs->make<TH1F> ("zMuPt1_5j", "zMuPt1_5j", 100, 0, 100 ) ;
  zMu1Mt_5j = fs->make<TH1F> ("zMu1Mt_5j", "zMu1Mt_5j", 100, 0, 100 ) ;
  zMuPt2_5j= fs->make<TH1F> ("zMuPt2_5j", "zMuPt2_5j", 100, 0, 100 ) ;
  zMu2Mt_5j= fs->make<TH1F> ("zMu2Mt_5j", "zMu2Mt_5j", 100, 0, 100 ) ;
  zMuPt1_6j = fs->make<TH1F> ("zMuPt1_6j", "zMuPt1_6j", 100, 0, 100 ) ;
  zMu1Mt_6j = fs->make<TH1F> ("zMu1Mt_6j", "zMu1Mt_6j", 100, 0, 100 ) ;
  zMuPt2_6j= fs->make<TH1F> ("zMuPt2_6j", "zMuPt2_6j", 100, 0, 100 ) ;
  zMu2Mt_6j= fs->make<TH1F> ("zMu2Mt_6j", "zMu2Mt_6j", 100, 0, 100 ) ;
  zMuPt1_7j = fs->make<TH1F> ("zMuPt1_7j", "zMuPt1_7j", 100, 0, 100 ) ;
  zMu1Mt_7j = fs->make<TH1F> ("zMu1Mt_7j", "zMu1Mt_7j", 100, 0, 100 ) ;
  zMuPt2_7j= fs->make<TH1F> ("zMuPt2_7j", "zMuPt2_7j", 100, 0, 100 ) ;
  zMu2Mt_7j= fs->make<TH1F> ("zMu2Mt_7j", "zMu2Mt_7j", 100, 0, 100 ) ;
  zMuPt1_8j = fs->make<TH1F> ("zMuPt1_8j", "zMuPt1_8j", 100, 0, 100 ) ;
  zMu1Mt_8j = fs->make<TH1F> ("zMu1Mt_8j", "zMu1Mt_8j", 100, 0, 100 ) ;
  zMuPt2_8j= fs->make<TH1F> ("zMuPt2_8j", "zMuPt2_8j", 100, 0, 100 ) ;
  zMu2Mt_8j= fs->make<TH1F> ("zMu2Mt_8j", "zMu2Mt_8j", 100, 0, 100 ) ;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenSystInclusiveZ::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenSystInclusiveZ::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
GenSystInclusiveZ::beginRun( edm::Run const & iRun,  edm::EventSetup const& iSetup) {
  nevent_run = 0;
}
void
GenSystInclusiveZ::endRun( const edm::Run& iRun, const edm::EventSetup& iSetup) {
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
 
double GenSystInclusiveZ::dPhi(math::XYZTLorentzVector& vec1, math::XYZTLorentzVector& vec2) {
  return reco::deltaPhi(vec1.phi(),vec2.phi());
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenSystInclusiveZ);
