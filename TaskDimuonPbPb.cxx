/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

// include root libraries
#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TChain.h"
#include "TMath.h"
#include "THnSparse.h"
// include AliRoot Libraries
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliVTrack.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "TaskDimuonPbPb.h"

class TaskDimuonPbPb;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(TaskDimuonPbPb) // classimp: necessary for root

TaskDimuonPbPb::TaskDimuonPbPb() : AliAnalysisTaskSE(), 
    fAODEvent(0),
    fVEvent(0),
    fArrayMCParticles(0),
    fMuonTrackCuts(0),
    fTriggerClass(0),
    fFirstRun(0),
    fLastRun(0),
    fListEventHistos(0), 
    fListGeneratedHistos(0),
    fListReconstructedHistos(0),
    fHistoTotalEventsPerRun(0),
    fHistoDiMuonOS(0),
    fHistoSingleMuon(0),
    fHistoJPsiGenerated(0),
    fHistoSingleMuonGenerated(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
TaskDimuonPbPb::TaskDimuonPbPb(const char* name,int firstRun, int lastRun, UInt_t triggerClass) : AliAnalysisTaskSE(name),
    fAODEvent(0),
    fVEvent(0),
    fArrayMCParticles(0),
    fMuonTrackCuts(0),
    fTriggerClass(triggerClass),
    fFirstRun(firstRun),
    fLastRun(lastRun),
    fListEventHistos(0), 
    fListGeneratedHistos(0),
    fListReconstructedHistos(0),
    fHistoTotalEventsPerRun(0),
    fHistoDiMuonOS(0),
    fHistoSingleMuon(0),
    fHistoJPsiGenerated(0),
    fHistoSingleMuonGenerated(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
    DefineOutput(2, TList::Class());    // you can add more output objects by calling DefineOutput(2, classname::Class())
    DefineOutput(3, TList::Class());    // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
TaskDimuonPbPb::~TaskDimuonPbPb()
{
    // destructor
    if(fListEventHistos && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fListEventHistos;     // at the end of your task, it is deleted from memory by calling this function
    }
    if(fListGeneratedHistos && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fListGeneratedHistos;
    }
    if(fListReconstructedHistos && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fListReconstructedHistos;
    }
}
//_____________________________________________________________________________
Bool_t TaskDimuonPbPb::IsMuonFromJPsi(int muonIndex){

  AliAODMCParticle* muonParticle = (AliAODMCParticle*)fArrayMCParticles->At(muonIndex);
  if(muonIndex <0) return kFALSE;
  if((muonParticle->GetPdgCode() == 13) || (muonParticle->GetPdgCode() == -13)){

    int iMother = muonParticle->GetMother();
    while ( iMother >= 0 ) {
      AliAODMCParticle *firstMother = (AliAODMCParticle*)fArrayMCParticles->At(iMother);
      int pdgCodeOfFirstMother = firstMother->GetPdgCode();
      if(pdgCodeOfFirstMother == 443)
      {
        return kTRUE;
      }
      iMother = firstMother->GetMother();
    }
  }
  return kFALSE;
}
//_____________________________________________________________________________
void TaskDimuonPbPb::NotifyRun()
{
  /// Set run number for cuts
  if ( fMuonTrackCuts ) fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void TaskDimuonPbPb::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file

     //Event histograms
    fListEventHistos = new TList();
    fListEventHistos->SetOwner(kTRUE);

    fHistoTotalEventsPerRun = new TH1I("fHistoTotalEventsPerRun","",fLastRun - fFirstRun,fFirstRun,fLastRun);
    fHistoTotalEventsPerRun->Sumw2();
    fHistoTotalEventsPerRun->GetXaxis()->SetTitle("Run Number");
    fHistoTotalEventsPerRun->GetYaxis()->SetTitle("# Total Events");
    fListEventHistos->Add(fHistoTotalEventsPerRun);
    
    fListGeneratedHistos = new TList();
    fListGeneratedHistos->SetOwner(kTRUE);
    fListReconstructedHistos = new TList();
    fListReconstructedHistos->SetOwner(kTRUE);

        //SingleMuon histograms
        Int_t nbinsSingleMuon[5]={1000,60,100, 20,fLastRun - fFirstRun}; //pT, Eta, Theta, Phi, run
        Double_t xminSingleMuon[5]={0,-5,0.75*TMath::Pi(),-TMath::Pi(), (Double_t)fFirstRun}, xmaxSingleMuon[5]={100,-2,1.25*TMath::Pi(),TMath::Pi(), (Double_t)fLastRun};
        fHistoSingleMuon = new THnSparseD("fHistoSingleMuon","",5, nbinsSingleMuon,xminSingleMuon,xmaxSingleMuon, 1024*16);
        fHistoSingleMuon->Sumw2();
        fHistoSingleMuon->GetAxis(0)->SetTitle("p_{T} GeV/c");
        fHistoSingleMuon->GetAxis(1)->SetTitle("#eta");
        fHistoSingleMuon->GetAxis(2)->SetTitle("#theta");
        fHistoSingleMuon->GetAxis(3)->SetTitle("#phi");
        fHistoSingleMuon->GetAxis(4)->SetTitle("Run Number");
        fListReconstructedHistos->Add(fHistoSingleMuon);

        //DiMuon histograms
        Int_t nbinsDiMuon[4]={1000,400,60, fLastRun - fFirstRun}; //Mmumu, pT, y, centrality
        Double_t xminDiMuon[4]={0,0,-5, (Double_t)fFirstRun}, xmaxDiMuon[4]={20,20,-2, (Double_t)fLastRun};
        fHistoDiMuonOS = new THnSparseD("fHistoDiMuonOS","",4,nbinsDiMuon,xminDiMuon,xmaxDiMuon, 1024*16);
        fHistoDiMuonOS->Sumw2();
        fHistoDiMuonOS->GetAxis(0)->SetTitle("M_{#mu#mu} GeV/c^{2}");
        fHistoDiMuonOS->GetAxis(1)->SetTitle("p_{T} GeV/c");
        fHistoDiMuonOS->GetAxis(2)->SetTitle("y");
        fHistoDiMuonOS->GetAxis(3)->SetTitle("Run Number");
        fListReconstructedHistos->Add(fHistoDiMuonOS);

        //Generated histograms
        fHistoJPsiGenerated = new THnSparseD("fHistoJPsiGenerated","",4,nbinsDiMuon,xminDiMuon,xmaxDiMuon, 1024*16);
        fHistoJPsiGenerated->Sumw2();
        fHistoJPsiGenerated->GetAxis(0)->SetTitle("M_{#mu#mu} GeV/c^{2}");
        fHistoJPsiGenerated->GetAxis(1)->SetTitle("p_{T} GeV/c");
        fHistoJPsiGenerated->GetAxis(2)->SetTitle("y");
        fHistoJPsiGenerated->GetAxis(3)->SetTitle("Run Number");
        fListGeneratedHistos->Add(fHistoJPsiGenerated);

        fHistoSingleMuonGenerated = new THnSparseD("fHistoSingleMuonGenerated","",5, nbinsSingleMuon,xminSingleMuon,xmaxSingleMuon, 1024*16);
        fHistoSingleMuonGenerated->Sumw2();
        fHistoSingleMuonGenerated->GetAxis(0)->SetTitle("p_{T} GeV/c");
        fHistoSingleMuonGenerated->GetAxis(1)->SetTitle("#eta");
        fHistoSingleMuonGenerated->GetAxis(2)->SetTitle("#theta");
        fHistoSingleMuonGenerated->GetAxis(3)->SetTitle("#phi");
        fHistoSingleMuonGenerated->GetAxis(4)->SetTitle("Run Number");
        fListGeneratedHistos->Add(fHistoSingleMuonGenerated);

        //The muon muonTrackCuts can be defined here. Hiwever it is better to defien it outside (in addTaskDimuonPPB.C). To be fixed
        fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts","StandardMuonTrackCuts");
        fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
        fMuonTrackCuts->SetFilterMask (AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs  | AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuPdca);//Set the cuts to be used for the muon selections. See all the available cuts in AliMuonTrackCuts.h
        fMuonTrackCuts->SetIsMC();
    

  //This is needed to save the outputs.
  PostData(1, fListEventHistos);
  PostData(2, fListGeneratedHistos);
  PostData(3, fListReconstructedHistos);

}
//_____________________________________________________________________________
void TaskDimuonPbPb::UserExec(Option_t *)
{
    // user exec this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());    
    if(!fAODEvent) {
        AliError("ERROR: Could not retrieve AOD event !!");
        return;
    }
    int runNumber;

    AliAODMCHeader* mcHeader = dynamic_cast<AliAODMCHeader*>(fAODEvent->FindListObject("mcHeader"));
    if(mcHeader) 
    {
        TString generatorName;
        generatorName.Form("%s",mcHeader->GetGeneratorName());
    }

    fVEvent = static_cast<AliVEvent *>(InputEvent());
    fArrayMCParticles = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
    runNumber = fAODEvent->GetRunNumber();
    Float_t muonMass2 = AliAnalysisMuonUtility::MuonMass2(); // the PDG rest mass of the muon (constante, used for getting the kinematics) en GeV

    //Verification event with two muons
    AliAODMCParticle* particle;
    int imuons = 0;
    if(fArrayMCParticles->GetEntries() > 2 )
    {
        cout << fArrayMCParticles->GetEntries() << endl;
        for (int ipart=0; ipart<fArrayMCParticles->GetEntries(); ipart++)
        {
            particle = (AliAODMCParticle*)fArrayMCParticles->At(ipart);
		    if(!particle) continue;

            Short_t PDGcode = particle->GetPdgCode();
            if (PDGcode == 13 || PDGcode == -13) imuons++;
        }
        if (imuons >2) return;
    }

    // Event Histos   
    fHistoTotalEventsPerRun->Fill(runNumber);

    //Generated
    AliAODMCParticle* particle1;
    AliAODMCParticle* particle2;
    TLorentzVector lvPart1, lvPart2, lvDiMuonGen;

	for(int iParticle=0;iParticle<fArrayMCParticles->GetEntries();iParticle++)
	{
		bool isMuonPlusFound=false, isMuonMinusFound=false;
		
	    particle1 = (AliAODMCParticle*)fArrayMCParticles->At(iParticle);
		if(!particle1) continue;

        Short_t muonPDG3 = particle1->GetPdgCode();
        if (!(muonPDG3 == 13 || muonPDG3 == -13)) continue;

        Float_t energy3 = TMath::Sqrt(particle1->P()*particle1->P() + muonMass2);
        lvPart1.SetPxPyPzE(particle1->Px(),particle1->Py(),particle1->Pz(),energy3); //def 4-vect muon1
        Short_t muonCharge3 = particle1->Charge();
        //singlemuon        
        Double_t propertiesSingleMuonGen[5]={lvPart1.Pt(),lvPart1.Eta(),lvPart1.Theta(),lvPart1.Phi(), (Double_t)runNumber};
        fHistoSingleMuonGenerated->Fill(propertiesSingleMuonGen,1);


		
        for(int iParticle2=iParticle+1;iParticle2<fArrayMCParticles->GetEntries();iParticle2++)
        {
            particle2 = (AliAODMCParticle*)fArrayMCParticles->At(iParticle2);
            if(!particle2) continue;

            Short_t muonPDG4 = particle2->GetPdgCode();
            if (!(muonPDG4 == 13 || muonPDG4 == -13)) continue;

            Float_t energy4 = TMath::Sqrt(particle2->P()*particle2->P() + muonMass2);
            lvPart2.SetPxPyPzE(particle2->Px(),particle2->Py(),particle2->Pz(),energy4); //def 4-vect muon1
            Short_t muonCharge4 = particle2->Charge();
            
            if(muonCharge3 == muonCharge4) continue;

            lvDiMuonGen = lvPart1 + lvPart2;
            Double_t propertiesDiMuonGen[4]={};
            propertiesDiMuonGen[0]=lvDiMuonGen.M();
            propertiesDiMuonGen[1]=lvDiMuonGen.Pt();
            propertiesDiMuonGen[2]=lvDiMuonGen.Rapidity();
            propertiesDiMuonGen[3]= (Double_t)runNumber;

            fHistoJPsiGenerated->Fill(propertiesDiMuonGen,1);
            
        }
    }    

    //Fill Single Muon and Dimuon properties histograms
    TLorentzVector lvMuon1, lvMuon2, lvDiMuonRec;
    AliVParticle* muonTrack1;
    AliVParticle* muonTrack2;
    int muonIndex1, muonIndex2;

    int numberOfTracks = AliAnalysisMuonUtility::GetNTracks(fVEvent); // get the number of muon tracks in the event
      
    for(Int_t iMuon1 = 0 ; iMuon1 < numberOfTracks ; iMuon1++) // loop ove rall these tracks
    {
        muonTrack1 = AliAnalysisMuonUtility::GetTrack(iMuon1,fVEvent);
        if( !muonTrack1 ) { AliError(Form("ERROR: Could not retrieve AOD or ESD track %d", iMuon1)); continue;}
        if( !fMuonTrackCuts->IsSelected(muonTrack1) ) continue;//include cuts on pDCA, Eta, Rabs
        muonIndex1 = muonTrack1->GetLabel();
        if(muonIndex1 <0) continue;
        AliAODMCParticle* muonParticle1 = (AliAODMCParticle*)fArrayMCParticles->At(muonIndex1);
        Short_t muonPDG1 = (muonParticle1->GetPdgCode());
        if (!(muonPDG1 == 13 || muonPDG1 ==-13)) continue;
        
        //single muon properties 
        Float_t energy1 = TMath::Sqrt(muonTrack1->P()*muonTrack1->P() + muonMass2);
        lvMuon1.SetPxPyPzE(muonTrack1->Px(),muonTrack1->Py(),muonTrack1->Pz(),energy1); //def 4-vect muon1
        Short_t muonCharge1 = muonTrack1->Charge();
                
        Double_t propertiesSingleMuonRec[5]={lvMuon1.Pt(),lvMuon1.Eta(),lvMuon1.Theta(),lvMuon1.Phi(), (Double_t)runNumber};
        fHistoSingleMuon->Fill(propertiesSingleMuonRec,1);
                
        for (Int_t iMuon2 = iMuon1+1; iMuon2 < numberOfTracks; iMuon2++)
        {
            muonTrack2 = AliAnalysisMuonUtility::GetTrack(iMuon2,fVEvent);
            if ( !muonTrack2 ) {AliError(Form("ERROR: Could not retrieve AOD or ESD track %d", iMuon2)); continue;}
            if ( ! fMuonTrackCuts->IsSelected(muonTrack2) ) continue;//include cuts on pDCA, Eta, Rabs
            muonIndex2 = muonTrack2->GetLabel();
            if(muonIndex2 <0) continue;
            AliAODMCParticle* muonParticle2 = (AliAODMCParticle*)fArrayMCParticles->At(muonIndex2);
            Short_t muonPDG2 = (muonParticle2->GetPdgCode());
            if (!(muonPDG2 == 13 || muonPDG2 ==-13)) continue;
            

            Float_t energy2 = TMath::Sqrt(muonTrack2->P()*muonTrack2->P() + muonMass2);
            lvMuon2.SetPxPyPzE(muonTrack2->Px(),muonTrack2->Py(),muonTrack2->Pz(),energy2); //def 4-vect muon1
            Short_t muonCharge2 = muonTrack2->Charge();

            if (muonCharge1 == muonCharge2) continue;
            //dimuon properties
            lvDiMuonRec = lvMuon1+lvMuon2;

            Double_t propertiesDiMuonRec[4]={};
            propertiesDiMuonRec[0]=lvDiMuonRec.M();
            propertiesDiMuonRec[1]=lvDiMuonRec.Pt();
            propertiesDiMuonRec[2]=lvDiMuonRec.Rapidity();
            propertiesDiMuonRec[3]= (Double_t)runNumber;

            fHistoDiMuonOS->Fill(propertiesDiMuonRec,1);
        }
    }//end for reconstructed muons

 	

    PostData(1, fListEventHistos);
    PostData(2, fListGeneratedHistos);
    PostData(3, fListReconstructedHistos);
}
//_____________________________________________________________________________
void TaskDimuonPbPb::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}

