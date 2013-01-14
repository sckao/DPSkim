// -*- C++ -*-
//
// Package:    DPSkim
// Class:      DPSkim
// 
/**\class DPSkim DPSkim.cc EXO/DPSkim/src/DPSkim.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri Jan 11 18:43:30 CST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//
// class declaration
//
using namespace std ;

class DPSkim : public edm::EDFilter {
   public:
      explicit DPSkim(const edm::ParameterSet&);
      ~DPSkim();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void TriggerTagging( edm::Handle<edm::TriggerResults> triggers, 
                           const edm::TriggerNames& trgNameList, 
			   int RunId, vector<int>& firedTrig ) ;
      bool TriggerSelection( edm::Handle<edm::TriggerResults> triggers, vector<int> firedTrigID ) ;
      bool PhotonSelection( edm::Handle<reco::PhotonCollection> photons ) ;
      bool JetSelection( edm::Handle<reco::PFJetCollection> jets ) ;

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag trigSource;
      edm::InputTag trigEvent;
      edm::InputTag photonSource;
      edm::InputTag metSource;
      edm::InputTag jetSource;

      std::vector<double> photonCuts ;
      std::vector<double> photonIso ;
      std::vector<double> metCuts ;
      std::vector<double> jetCuts ;
      std::vector<int>    firedTrig ;
      std::vector<string> triggerPatent ;
      int runID_ ;

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
DPSkim::DPSkim(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   trigSource           = iConfig.getParameter<edm::InputTag> ("trigSource");
   photonSource         = iConfig.getParameter<edm::InputTag> ("photonSource");
   metSource            = iConfig.getParameter<edm::InputTag> ("metSource");
   jetSource            = iConfig.getParameter<edm::InputTag> ("jetSource");

   jetCuts              = iConfig.getParameter<std::vector<double> >("jetCuts");
   metCuts              = iConfig.getParameter<std::vector<double> >("metCuts");
   photonCuts           = iConfig.getParameter<std::vector<double> >("photonCuts");
   photonIso            = iConfig.getParameter<std::vector<double> >("photonIso");
   triggerPatent        = iConfig.getParameter< std::vector<string> >("triggerName");

   firedTrig.clear() ;
   for ( size_t i=0; i< triggerPatent.size(); i++ ) firedTrig.push_back(-1) ;

   runID_ = 0 ;

}


DPSkim::~DPSkim() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool DPSkim::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   using namespace edm;

   bool pass = true ;

   // HLT Selection
   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( trigSource, triggers );
   const edm::TriggerNames& trgNameList = iEvent.triggerNames( *triggers ) ;

   int run_id    = iEvent.id().run()  ;
   TriggerTagging( triggers, trgNameList, run_id, firedTrig ) ;
   bool passHLT = TriggerSelection( triggers, firedTrig ) ;
   if ( !passHLT ) pass = false ;
   if ( !passHLT ) cout<<" fail HLT "<<endl ;

   // Object Selection
   Handle<reco::PhotonCollection>      photons;
   Handle<reco::PFJetCollection>       jets;
   Handle<reco::PFMETCollection>       met; 

   iEvent.getByLabel( trigSource,     triggers );
   iEvent.getByLabel( photonSource,   photons  );
   iEvent.getByLabel( jetSource,      jets  );
   iEvent.getByLabel( metSource,      met  );

   bool passPhoton = PhotonSelection( photons ) ;
   if ( !passPhoton ) pass = false ;
   if ( !passPhoton ) cout<<" fail Photon "<<endl ;

   bool passJet    = JetSelection( jets  );
   if ( !passJet ) pass = false ;
   if ( !passJet ) cout<<" fail Jet "<<endl ;

   const reco::PFMET pfMet = (*met)[0] ;
   bool passMET = ( pfMet.et() > metCuts[0]  ) ? true : false ; 
   if ( !passMET ) pass = false ;
   if ( !passMET ) cout<<" fail MET "<<endl ;

   return pass ;
}

void DPSkim::TriggerTagging( edm::Handle<edm::TriggerResults> triggers, const edm::TriggerNames& trgNameList, int RunId, vector<int>& firedTrig ) {

   if ( runID_ != RunId )  {
      for (size_t j=0; j< triggerPatent.size(); j++ ) firedTrig[j] = -1;

             // loop through trigger menu
	     for ( size_t i =0 ; i < trgNameList.size(); i++ ) {
		 string tName  = trgNameList.triggerName( i );
		 // loop through desired triggers
		 for ( size_t j=0; j< triggerPatent.size(); j++ )  {
			 if ( strncmp( tName.c_str(), triggerPatent[j].c_str(), triggerPatent[j].size() ) ==0 ) {
				 firedTrig[j] = i;
				 cout<<" Trigger Found ("<<j <<"):  "<<tName ;
				 cout<<" Idx: "<< i <<" triggers "<<endl;
			 }
		 }
	     }
	     runID_ = RunId ;
      }
}	

bool DPSkim::TriggerSelection( edm::Handle<edm::TriggerResults> triggers, vector<int> firedTrigID ) {

   bool pass =false ;
   uint32_t trgbits = 0 ;
   for ( size_t i=0; i< firedTrigID.size(); i++ ) {
	   if ( firedTrigID[i] == -1 ) continue ;
	   if ( triggers->accept( firedTrigID[i] ) == 1  ) trgbits |= ( 1 << i ) ;
   }

   if ( trgbits != 0 )   pass = true ;
   
   return pass ;
}

bool DPSkim::PhotonSelection( edm::Handle<reco::PhotonCollection> photons ) { 

   double maxPt = 0 ;
   int ng = 0 ; 
   for(reco::PhotonCollection::const_iterator it = photons->begin(); it != photons->end(); it++) {

       if ( it->pt()  < photonCuts[0] ) continue ;
       if ( fabs(it->eta()) > photonCuts[1] ) continue ;

       //if ( it->hasPixelSeed() ) continue ; 
       float ecalSumEt = it->ecalRecHitSumEtConeDR04();
       float hcalSumEt = it->hcalTowerSumEtConeDR04();
       float trkSumPt  = it->trkSumPtSolidConeDR04();

       bool trkIso  = ( ( trkSumPt / it->pt())     < photonIso[0] ) ;
       bool ecalIso = ( (ecalSumEt / it->energy()) < photonIso[2] && ecalSumEt < photonIso[1] ) ;
       bool hcalIso = ( (hcalSumEt / it->energy()) < photonIso[4] && hcalSumEt < photonIso[3] ) ;
       if ( !trkIso || !ecalIso || !hcalIso ) continue ;

       maxPt = ( it->pt() > maxPt ) ? it->pt() : maxPt ;
       ng++ ;
   }
   if ( ng < photonCuts[2] || maxPt < photonCuts[3] )  return false ;
   else                                                return true ;

}


bool DPSkim::JetSelection( edm::Handle<reco::PFJetCollection> jets ) {

   int nj = 0 ;
   for(reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); it++) {
       if ( it->pt() < jetCuts[0] || fabs( it->eta() ) > jetCuts[1] ) continue ;
       nj++ ;
   }

   if ( nj < jetCuts[2] )  return false ;
   else                    return true ;

}


// ------------ method called once each job just before starting event loop  ------------
void 
DPSkim::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DPSkim::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
DPSkim::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
DPSkim::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
DPSkim::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
DPSkim::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DPSkim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DPSkim);
