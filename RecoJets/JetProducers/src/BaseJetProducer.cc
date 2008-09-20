// File: BaseJetProducer.cc
// Author: F.Ratnikov UMd Aug 22, 2006
// $Id: BaseJetProducer.cc,v 1.37 2008/09/10 12:09:28 oehler Exp $
//--------------------------------------------
#include <memory>

#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "RecoJets/JetAlgorithms/interface/JetRecoTypes.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "RecoJets/JetAlgorithms/interface/JetMaker.h"
#include "RecoJets/JetAlgorithms/interface/ProtoJet.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "BaseJetProducer.h"
#include <iostream>

using namespace std;
using namespace reco;
using namespace JetReco;

namespace {
  const bool debug = false;

  bool makeCaloJet (const string& fTag) {
    return fTag == "CaloJet";
  }
  bool makePFJet (const string& fTag) {
    return fTag == "PFJet";
  }
  bool makeGenJet (const string& fTag) {
    return fTag == "GenJet";
  }
  bool makeBasicJet (const string& fTag) {
    return fTag == "BasicJet";
  }

  bool makeGenericJet (const string& fTag) {
    return !makeCaloJet (fTag) && !makePFJet (fTag) && !makeGenJet (fTag) && !makeBasicJet (fTag);
  }

  template <class T>  
  void dumpJets (const T& fJets) {
    for (unsigned i = 0; i < fJets.size(); ++i) {
      std::cout << "Jet # " << i << std::endl << fJets[i].print();
    }
  }

  void copyVariables (const ProtoJet& fProtojet, reco::Jet* fJet) {
    fJet->setJetArea (fProtojet.jetArea ());
    fJet->setPileup (fProtojet.pileup ());
    fJet->setNPasses (fProtojet.nPasses ());
  }

  void copyConstituents (const JetReco::InputCollection& fConstituents, const edm::View <Candidate>& fInput, reco::Jet* fJet) {
    // put constituents
    for (unsigned iConstituent = 0; iConstituent < fConstituents.size (); ++iConstituent) {
      fJet->addDaughter (fInput.ptrAt (fConstituents[iConstituent].index ()));
    }
  }

}

namespace cms
{

  // Constructor takes input parameters now: to be replaced with parameter set.
  BaseJetProducer::BaseJetProducer(const edm::ParameterSet& conf)
    : mSrc(conf.getParameter<edm::InputTag>( "src" )),
      mJetType (conf.getUntrackedParameter<string>( "jetType", "CaloJet")),
      mVerbose (conf.getUntrackedParameter<bool>("verbose", false)),
      mEtInputCut (conf.getParameter<double>("inputEtMin")),
      mEInputCut (conf.getParameter<double>("inputEMin")),
      mJetPtMin (conf.getParameter<double>("jetPtMin")),
      mVertexCorrectedInput(false)
  {
    std::string alias = conf.getUntrackedParameter<string>( "alias", conf.getParameter<std::string>("@module_label"));
    if (makeCaloJet (mJetType)) {
      produces<CaloJetCollection>().setBranchAlias (alias);
      mVertexCorrectedInput=conf.getParameter<bool>("correctInputToSignalVertex");
      if (mVertexCorrectedInput){
	mPVCollection=conf.getParameter<edm::InputTag>("pvCollection");
      }
    }
    else if (makePFJet (mJetType)) produces<PFJetCollection>().setBranchAlias (alias);
    else if (makeGenJet (mJetType)) produces<GenJetCollection>().setBranchAlias (alias);
    else if (makeBasicJet (mJetType)) produces<BasicJetCollection>().setBranchAlias (alias);
    using namespace std;
  }

  // Virtual destructor needed.
  BaseJetProducer::~BaseJetProducer() { } 

  // Functions that gets called by framework every event
  void BaseJetProducer::produce(edm::Event& e, const edm::EventSetup& fSetup)
  {
    //set default vertex for undefined cases:
    vertex=reco::Jet::Point(0,0,0);
    //getSignalVertex (when producing caloJets, and configuration wants it)
    if (makeCaloJet(mJetType)) {
      if (mVertexCorrectedInput){
	edm::Handle<reco::VertexCollection> thePrimaryVertexCollection;
	e.getByLabel(mPVCollection,thePrimaryVertexCollection);
	vertex = (*thePrimaryVertexCollection)[0].position();
      }
    }
    
    // get input
    edm::Handle<edm::View <Candidate> > inputHandle; 
    e.getByLabel( mSrc, inputHandle);
    // convert to input collection
    JetReco::InputCollection input;
    vector<Candidate*> garbageCollection;
    if (mVertexCorrectedInput) garbageCollection.reserve(inputHandle->size());
    input.reserve (inputHandle->size());
    for (unsigned i = 0; i < inputHandle->size(); ++i) {
      JetReco::InputItem tmpInput(&((*inputHandle)[i]),i);
      if (mVertexCorrectedInput){
	tmpInput.setOriginal(tmpInput.get());
	const CaloTower *tower=dynamic_cast<const CaloTower*>(&(*inputHandle)[i]);
	Candidate* tmpCandidate=new CaloTower(*tower);
	math::PtEtaPhiMLorentzVector newCaloTowerVector(tower->p4(vertex));
	reco::Particle::LorentzVector correctedP4(newCaloTowerVector.px(),newCaloTowerVector.py(),newCaloTowerVector.pz(),newCaloTowerVector.energy());
	garbageCollection.push_back(tmpCandidate);
	tmpCandidate->setP4(correctedP4);
	tmpInput.setBase(tmpCandidate);
      }
      if ((mEtInputCut <= 0 || tmpInput->et() > mEtInputCut) &&
        (mEInputCut <= 0 || tmpInput->energy() > mEInputCut)) {
	input.push_back (tmpInput);
      }
    }
    if (mVerbose) {
      std::cout << "BaseJetProducer::produce-> INPUT COLLECTION selected from" << mSrc 
		<< " with ET > " << mEtInputCut << " and/or E > " << mEInputCut << std::endl;
      std::cout << "correct input to vertex: "<<mVertexCorrectedInput<<std::endl;
      for (unsigned index = 0; index < input.size(); ++index) {
	std::cout << "  Input " << index << ", px/py/pz/pt/e: "
		  << input[index]->px() << '/' << input[index]->py() << '/' << input[index]->pz() << '/' 
		  << input[index]->pt() << '/' << input[index]->energy()
 		  << std::endl;
      }
    }
    
    // run algorithm
    vector <ProtoJet> output;
    if (input.empty ()) {
      edm::LogInfo ("Empty Event") << "empty input for jet algorithm: bypassing..." << std::endl;
    }
    else {
      runAlgorithm (input, &output);
    }

    // produce output collection
    if (mVerbose) {
      std::cout << "OUTPUT JET COLLECTION:" << std::endl;
    }
    // make sure protojets are sorted
    sortByPt (&output);
    if (makeCaloJet (mJetType)) {
      edm::ESHandle<CaloGeometry> geometry;
      fSetup.get<CaloGeometryRecord>().get(geometry);
      const CaloSubdetectorGeometry* towerGeometry = 
	geometry->getSubdetectorGeometry(DetId::Calo, CaloTowerDetId::SubdetId);
      auto_ptr<CaloJetCollection> jets (new CaloJetCollection);
      jets->reserve(output.size());
      for (unsigned iJet = 0; iJet < output.size (); ++iJet) {
	ProtoJet* protojet=0;
	protojet = &(output [iJet]);
	if (protojet->p4().pt()<mJetPtMin) continue;
	const JetReco::InputCollection& constituents = protojet->getTowerList();
	CaloJet::Specific specific;
	JetMaker::makeSpecific (constituents, *towerGeometry, &specific);
	jets->push_back (CaloJet (protojet->p4(), vertex, specific));
	Jet* newJet = &(jets->back());
        copyConstituents (constituents, *inputHandle, newJet);
	copyVariables (*protojet, newJet);
      }
      if (mVerbose) dumpJets (*jets);
      e.put(jets);
    }
    else if (makePFJet (mJetType)) {
      auto_ptr<PFJetCollection> jets (new PFJetCollection);
      jets->reserve(output.size());
      for (unsigned iJet = 0; iJet < output.size (); ++iJet) {
	ProtoJet* protojet = &(output [iJet]);
	if (protojet->p4().pt()<mJetPtMin) continue;
	const JetReco::InputCollection& constituents = protojet->getTowerList();
	PFJet::Specific specific;
	JetMaker::makeSpecific (constituents, &specific);
	jets->push_back (PFJet (protojet->p4(), vertex, specific));
	Jet* newJet = &(jets->back());
	copyConstituents (constituents, *inputHandle, newJet);
	copyVariables (*protojet, newJet);
      }
      if (mVerbose) dumpJets (*jets);
      e.put(jets);
    }
    else if (makeGenJet (mJetType)) {
      auto_ptr<GenJetCollection> jets (new GenJetCollection);
      jets->reserve(output.size());
      for (unsigned iJet = 0; iJet < output.size (); ++iJet) {
	ProtoJet* protojet = &(output [iJet]);
	if (protojet->p4().pt()<mJetPtMin) continue;
	const JetReco::InputCollection& constituents = protojet->getTowerList();
	GenJet::Specific specific;
	JetMaker::makeSpecific (constituents, &specific);
	jets->push_back (GenJet (protojet->p4(), vertex, specific));
	Jet* newJet = &(jets->back());
	copyConstituents (constituents, *inputHandle, newJet);
	copyVariables (*protojet, newJet);
      }
      if (mVerbose) dumpJets (*jets);
      e.put(jets);
    }
    else if (makeBasicJet (mJetType)) {
      auto_ptr<BasicJetCollection> jets (new BasicJetCollection);
      jets->reserve(output.size());
      for (unsigned iJet = 0; iJet < output.size (); ++iJet) {
	ProtoJet* protojet = &(output [iJet]);
	if (protojet->p4().pt()<mJetPtMin) continue;
	const JetReco::InputCollection& constituents = protojet->getTowerList();
	jets->push_back (BasicJet (protojet->p4(), vertex));
	Jet* newJet = &(jets->back());
	copyConstituents (constituents, *inputHandle, newJet);
	copyVariables (*protojet, newJet);
      }
      if (mVerbose) dumpJets (*jets);
      e.put(jets);
    }
      //clean up garbage from modified calojet input:
    if (mVertexCorrectedInput){
      for (vector<Candidate*>::iterator iter=garbageCollection.begin();iter!=garbageCollection.end();++iter){
	delete *iter;
      }
    }

//     else if (makeGenericJet (mJetType)) {
//       auto_ptr<GenericJetCollection> jets (new GenericJetCollection);
//       jets->reserve(output.size());
//       for (unsigned iJet = 0; iJet < output.size (); ++iJet) {
// 	ProtoJet* protojet = output [iJet];
// 	const JetReco::InputCollection& constituents = protojet->getTowerList();
// 	jets->push_back (GenericJet (protojet->p4()));
// 	Jet* newJet = &(jets->back());
// 	copyConstituents (constituents, *inputHandle, newJet);
// 	copyVariables (*protojet, newJet);
//       }
//       if (mVerbose) dumpJets (*jets);
//       e.put(jets);
//     }
  }
}

