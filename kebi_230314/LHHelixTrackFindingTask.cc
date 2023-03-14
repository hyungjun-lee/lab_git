#include "KBRun.hh"
#include "LHHelixTrackFindingTask.hh"

#include <iostream>

// #define DEBUG_STEP
#define FT

ClassImp(LHHelixTrackFindingTask)

    bool LHHelixTrackFindingTask::Init()
{
  auto run = KBRun::GetRun();
  fPar = run->GetParameterContainer();
  fTpc = (LHTpc *)(run->GetDetectorSystem()->GetTpc());
  // FTDetector = (KBFTDetector *)(run->GetDetectorSystem()->GetTpc());
  fPadPlane = (KBPadPlane *)fTpc->GetPadPlane();
  // FTPad = (KBPadPlane *)FTDetector->GetPadPlane();

  fHitArray = (TClonesArray *)run->GetBranch(fBranchNameHit);
#ifdef FT
  fHitArray_FT = (TClonesArray *)run->GetBranch(fBranchNameHit_FT);
#endif

  fTrackArray = new TClonesArray("KBHelixTrack");
  run->RegisterBranch(fBranchNameTracklet, fTrackArray, fPersistency);

  fTrackHits = new KBHitArray();
  fCandHits = new KBHitArray();
  fCandHits_FT = new KBHitArray();
  fGoodHits = new KBHitArray();
  fGoodHits_FT = new KBHitArray();
  fBadHits = new KBHitArray();
  fBadHits_FT = new KBHitArray();

  run->RegisterBranch("TrackHit", fTrackHits, false);
  run->RegisterBranch("CandHit", fCandHits, false);
  run->RegisterBranch("GoodHit", fGoodHits, false);
  run->RegisterBranch("BadHit", fBadHits, false);

  fDefaultScale = fPar->GetParDouble("LHTF_defaultScale");
  fTrackWCutLL = fPar->GetParDouble("LHTF_trackWCutLL");
  fTrackWCutHL = fPar->GetParDouble("LHTF_trackWCutHL");
  fTrackHCutLL = fPar->GetParDouble("LHTF_trackHCutLL");
  fTrackHCutHL = fPar->GetParDouble("LHTF_trackHCutHL");
  fReferenceAxis = fPar->GetParAxis("LHTF_refAxis");

  fNextStep = StepNo::kStepInitArray;

  return true;
}

void LHHelixTrackFindingTask::Exec(Option_t *)
{
  fNextStep = StepNo::kStepInitArray;
  while (ExecStep())
  {
  }
}

bool LHHelixTrackFindingTask::ExecStep()
{
  if (fNextStep == kStepEndOfEvent)
    return false;
  else if (fNextStep == kStepInitArray)
    fNextStep = StepInitArray();
  else if (fNextStep == kStepNewTrack)
    fNextStep = StepNewTrack();
  else if (fNextStep == kStepRemoveTrack)
    fNextStep = StepRemoveTrack();
  else if (fNextStep == kStepInitTrack)
    fNextStep = StepInitTrack();
  else if (fNextStep == kStepInitTrackAddHit)
    fNextStep = StepInitTrackAddHit();
  else if (fNextStep == kStepContinuum)
    fNextStep = StepContinuum();
  else if (fNextStep == kStepContinuumAddHit)
    fNextStep = StepContinuumAddHit();
  else if (fNextStep == kStepExtrapolation)
    fNextStep = StepExtrapolation();
  else if (fNextStep == kStepExtrapolationAddHit)
    fNextStep = StepExtrapolationAddHit();
  else if (fNextStep == kStepConfirmation)
    fNextStep = StepConfirmation();
  else if (fNextStep == kStepFinalizeTrack)
    fNextStep = StepFinalizeTrack();
  else if (fNextStep == kStepNextPhase)
    fNextStep = StepNextPhase();
  else if (fNextStep == kStepEndEvent)
    fNextStep = StepEndEvent();
  return true;
}

bool LHHelixTrackFindingTask::ExecStepUptoTrackNum(Int_t numTracks)
{
  if (fNextStep == kStepEndOfEvent)
    return false;

  while (ExecStep())
  {
    if (fNextStep == kStepNewTrack && fTrackArray->GetEntriesFast() >= numTracks)
      break;
  }

  return true;
}

int LHHelixTrackFindingTask::StepInitArray()
{
  fCurrentTrack = nullptr;

  fPhaseIndex = 0;

  fPadPlane->ResetHitMap();
  fPadPlane->SetHitArray(fHitArray);

  fTrackArray->Clear("C");
  fTrackHits->Clear();
  fCandHits->Clear();
  // fCandHits_FT -> Clear();
  fGoodHits->Clear();
  fGoodHits_FT->Clear();
  fBadHits->Clear();
  fBadHits_FT->Clear();

  if (firstinit)
  {
    cnt = 0;
    firstinit = false;
  }

  return kStepNewTrack;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LHHelixTrackFindingTask::StepNewTrack()
{
  fTrackHits->Clear();
  fCandHits->Clear();
  fGoodHits->Clear();
  fGoodHits_FT->Clear();
  ReturnBadHitsToPadPlane();

  KBTpcHit *hit = fPadPlane->PullOutNextFreeHit();
  if (hit == nullptr)
  {
    return kStepNextPhase;
  }

  Int_t idx = fTrackArray->GetEntries();
  fCurrentTrack = new ((*fTrackArray)[idx]) KBHelixTrack(idx);
  fCurrentTrack->SetReferenceAxis(fReferenceAxis);
  fCurrentTrack->AddHit(hit);
  fGoodHits->AddHit(hit);
#ifdef FT
  if (cnt == 0)
  {
    for (int i = 0; i < fHitArray_FT->GetEntries(); i++)
    {
      auto hit_FT = (KBHit *)fHitArray_FT->At(i);
      fGoodHits_FT->AddHit(hit_FT);
    }
    cnt++;
    kb_debug << "[hits in FT] :: " << fGoodHits_FT->GetNumHits() << endl;
  }

  for (int i = 0; i < fBadHits_FT->GetEntries(); i++)
  {
    auto usedhit_FT = (KBHit *)fBadHits_FT->At(i);
    // fCurrentTrack->AddHit(usedhit_FT); // quality check 필요;
    fGoodHits_FT->AddHit(usedhit_FT);
  }
  kb_debug << "[hits new track] :: " << fBadHits_FT->GetNumHits() << endl;
  fBadHits_FT->Clear();
#endif
  // kb_debug << "[NewTrack : ] " << hit->GetPadID() << "\t" << fGoodHits -> GetNumHits()<< endl;

  return kStepInitTrack;
}

int LHHelixTrackFindingTask::StepRemoveTrack()
{
  fGoodHits->Clear();
  fTrackHits->Clear();
  fCandHits->Clear();
  ReturnBadHitsToPadPlane();

  auto trackHits = fCurrentTrack->GetHitArray();
  Int_t numTrackHits = trackHits->GetNumHits();
  for (Int_t iTrackHit = 0; iTrackHit < numTrackHits; ++iTrackHit)
  {
    auto trackHit = (KBTpcHit *)trackHits->GetHit(iTrackHit);
    if (trackHit->GetHitID() >= 40000) // FT hit selection ..?
    {                                  // continue
      continue;
    }
    trackHit->AddTrackCand(-1);
    fPadPlane->AddHit(trackHit);
  }
  fTrackArray->Remove(fCurrentTrack);
  fCurrentTrack = nullptr;
  kb_debug << "[hits remove track ] :: " << fBadHits_FT->GetNumHits() << endl;
  firstinit = true;
  return kStepNewTrack;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LHHelixTrackFindingTask::StepInitTrack()
{
  // kb_debug << "[Init :: ]"<<fGoodHits -> GetNumHits() << endl;
  fCandHits->Clear();
  fPadPlane->PullOutNeighborHits(fGoodHits, fCandHits);
  // kb_debug << "[candhit :: ] " << fCandHits->GetNumHits() << endl;
  fGoodHits->MoveHitsTo(fTrackHits);
  fNumCandHits = fCandHits->GetEntriesFast();
  if (fNumCandHits == 0)
  {
    kb_debug << "no candhit remove!!" << endl;
    return kStepRemoveTrack;
  }
  fCandHits->SortByDistanceTo(fCurrentTrack->GetMean(), true);
  return kStepInitTrackAddHit;
}

int LHHelixTrackFindingTask::StepInitTrackAddHit()
{
  // kb_debug << "here " << endl;
  auto candHit = (KBTpcHit *)fCandHits->GetLastHit();
  fCandHits->RemoveLastHit();

  Double_t quality;

  if (fCurrentTrack->IsHelix())
    quality = CorrelateHitWithTrack(fCurrentTrack, candHit);
  else
    quality = CorrelateHitWithTrackCandidate(fCurrentTrack, candHit);

  if (quality > 0)
  {
    fGoodHits->AddHit(candHit);
    fCurrentTrack->AddHit(candHit);
    fCurrentTrack->FitPlane(); // XXX should comment out

    auto numHitsInTrack = fCurrentTrack->GetNumHits();

    if (numHitsInTrack > fCutMaxNumHitsInitTrack)
    {
      Int_t numCandHits2 = fCandHits->GetEntriesFast();
      for (Int_t iCand = 0; iCand < numCandHits2; ++iCand)
      {
        fPadPlane->AddHit((KBTpcHit *)fCandHits->GetHit(iCand));
      }

      fCandHits->Clear("C");

      kb_debug << "cutmaxnumhits remove!!" << endl;
      return kStepRemoveTrack;
    }

    if (numHitsInTrack >= fMinHitsToFitInitTrack)
    {
      fCurrentTrack->Fit();
      if (numHitsInTrack > fCutMinNumHitsInitTrack &&
          fCurrentTrack->GetHelixRadius() > fCutMinHelixRadius &&
          fCurrentTrack->TrackLength() > fTrackLengthCutScale * fCurrentTrack->GetRMST())
      {
        return kStepContinuum;
      }
      else
      {
        fCurrentTrack->FitPlane();
      }
    }
  }
  else
  {
    fBadHits->AddHit(candHit);
    // kb_debug << candHit-> GetPadID() << endl;
  }

  fNumCandHits = fCandHits->GetEntriesFast();
  if (fNumCandHits == 0)
    return kStepInitTrack;

  return kStepInitTrackAddHit;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LHHelixTrackFindingTask::StepContinuum()
{
  fPadPlane->PullOutNeighborHits(fGoodHits, fCandHits);
  fGoodHits->MoveHitsTo(fTrackHits);

#ifdef FT
  fGoodHits_FT->MoveHitsTo(fCandHits_FT);

  auto fNumCandHits_FT = fCandHits_FT->GetEntries();
#endif

  fNumCandHits = fCandHits->GetEntries();
  if (fNumCandHits == 0 && fNumCandHits_FT == 0) /// modified for FT.
  {
    return kStepExtrapolation;
  }

  fCandHits->SortByCharge(false);
  return kStepContinuumAddHit;
}

int LHHelixTrackFindingTask::StepContinuumAddHit()
{
  for (Int_t iHit = 0; iHit < fNumCandHits; iHit++)
  {
    KBTpcHit *candHit = (KBTpcHit *)fCandHits->GetLastHit();
    fCandHits->RemoveLastHit();

    Double_t quality = 0;
    if (CheckParentTrackID(candHit) == -2)
      quality = CorrelateHitWithTrack(fCurrentTrack, candHit);

    if (quality > 0)
    {
      fGoodHits->AddHit(candHit);
      fCurrentTrack->AddHit(candHit);
      fCurrentTrack->Fit();
    }
    else
    {
      fBadHits->AddHit(candHit);
    }
  }

#ifdef FT
  for (Int_t iHit = 0; iHit < fCandHits_FT->GetEntries(); iHit++)
  {
    KBTpcHit *candHit_FT = (KBTpcHit *)fCandHits_FT->GetLastHit();
    fCandHits_FT->RemoveLastHit(); // sort 해야하는지 추후에 고민 할 것

    Double_t quality = 0;
    if (CheckParentTrackID(candHit_FT) == -2)
      quality = CorrelateHitWithTrack(fCurrentTrack, candHit_FT);
    if (quality > 0)
    {
      fCurrentTrack->AddHit(candHit_FT); //
      fCurrentTrack->Fit();
    }
    else
    {
      fBadHits_FT->AddHit(candHit_FT); //
    }
  }

#endif
  return kStepContinuum;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LHHelixTrackFindingTask::StepExtrapolation()
{
  ReturnBadHitsToPadPlane();

  return kStepExtrapolationAddHit;
}

int LHHelixTrackFindingTask::StepExtrapolationAddHit()
{
  Int_t countIteration = 0;
  bool buildFromHeadOrTail = true;
  Double_t extrapolationLength = 0;
  while (AutoBuildByExtrapolation(fCurrentTrack, buildFromHeadOrTail, extrapolationLength))
  { // true
    if (++countIteration > 200)
      break;
  }

  countIteration = 0;
  buildFromHeadOrTail = !buildFromHeadOrTail;
  extrapolationLength = 0;
  while (AutoBuildByExtrapolation(fCurrentTrack, buildFromHeadOrTail, extrapolationLength))
  { // true
    if (++countIteration > 200)
      break;
  }

  ReturnBadHitsToPadPlane();

  bool isGood = CheckTrackQuality(fCurrentTrack);

  if (isGood)
  {
    return kStepConfirmation;
  }
  kb_debug << "extra polation remove!!" << endl;
  return kStepRemoveTrack;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int LHHelixTrackFindingTask::StepConfirmation()
{

  bool tailToHead = false;
  KBVector3 pTail(fCurrentTrack->PositionAtTail(), fReferenceAxis);
  KBVector3 pHead(fCurrentTrack->PositionAtHead(), fReferenceAxis);

  if (pHead.K() > pTail.K()) //?? (pHead.K() > pHead.K())->(pHead.K() > pTail.K())
    tailToHead = true;

  ReturnBadHitsToPadPlane();

  if (BuildAndConfirmTrack(fCurrentTrack, tailToHead) == false)
  {
    kb_debug << "confirmation remove 1!!" << endl;
    return kStepRemoveTrack;
  }

  tailToHead = !tailToHead;

  ReturnBadHitsToPadPlane();
  if (BuildAndConfirmTrack(fCurrentTrack, tailToHead) == false)
  {
    kb_debug << "confirmation remove 2!!" << endl;
    return kStepRemoveTrack;
  }
  ReturnBadHitsToPadPlane();

  return kStepFinalizeTrack;
}

int LHHelixTrackFindingTask::StepFinalizeTrack() // 정리
{
  auto trackHits = fCurrentTrack->GetHitArray();
  Int_t trackID = fCurrentTrack->GetTrackID();
  Int_t numTrackHits = trackHits->GetNumHits();
  for (Int_t iTrackHit = 0; iTrackHit < numTrackHits; ++iTrackHit)
  { // tpc Hit만 cand에 넣음
    auto trackHit = (KBTpcHit *)trackHits->GetHit(iTrackHit);
    if (trackHit->GetHitID() >= 40000) // FT hit selection ..?
    {                                  // continue
      continue;
    }
    trackHit->AddTrackCand(trackID);
    fPadPlane->AddHit(trackHit);
  } //
  fGoodHits->MoveHitsTo(fTrackHits);
  fGoodHits->Clear();
  fGoodHits_FT->Clear();
  firstinit = true;
  kb_debug << "[hits finalize track] :: " << fBadHits_FT->GetNumHits() << endl;
  return kStepNewTrack;
}

int LHHelixTrackFindingTask::StepNextPhase()
{
  if (fPhaseIndex == 0)
  {
    fPhaseIndex = 1;

    fPadPlane->ResetEvent();

    fTrackHits->Clear();
    fCandHits->Clear();
    fGoodHits->Clear();
    ReturnBadHitsToPadPlane();

    // fMinHitsToFitInitTrack = 7;
    // fCutMinNumHitsInitTrack = 10;
    fCutMaxNumHitsInitTrack = 25;
    // fCutMinNumHitsFinalTrack = 15;
    // fCutMinHelixRadius = 30.;
    // fTrackLengthCutScale = 2.5;
    // fCutdkInExpectedTrackPath = 4.;

    return kStepNewTrack;
  }
  else if (fPhaseIndex == 1)
    return kStepEndEvent;
}

int LHHelixTrackFindingTask::StepEndEvent()
{
  fTrackArray->Compress();

  Int_t numTracks = fTrackArray->GetEntriesFast();
  for (Int_t iTrack = 0; iTrack < numTracks; ++iTrack)
  {
    KBHelixTrack *track = (KBHelixTrack *)fTrackArray->At(iTrack);
    track->SetTrackID(iTrack);
    track->FinalizeHits();
  }

  kb_info << "Number of found tracks: " << fTrackArray->GetEntries() << endl;

  return kStepEndOfEvent;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void LHHelixTrackFindingTask::ReturnBadHitsToPadPlane()
{
  fNumBadHits = fBadHits->GetEntriesFast(); // badhit solved
  for (Int_t iBad = 0; iBad < fNumBadHits; ++iBad)
  {
    fPadPlane->AddHit((KBTpcHit *)fBadHits->GetHit(iBad));
  }
  fBadHits->Clear();
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

double LHHelixTrackFindingTask::CorrelateHitWithTrackCandidate(KBHelixTrack *track, KBTpcHit *hit)
{
  if (hit->GetNumTrackCands() != 0)
    return 0;

  Double_t quality = 0;

  auto row = hit->GetRow();
  auto layer = hit->GetLayer();

  KBVector3 qosHit(hit->GetPosition(), fReferenceAxis);

  auto trackHits = track->GetHitArray();
  bool passCutdk = false;
  Int_t numTrackHits = trackHits->GetNumHits();
  for (Int_t iTrackHit = 0; iTrackHit < numTrackHits; ++iTrackHit)
  {
    auto trackHit = (KBTpcHit *)trackHits->GetHit(iTrackHit);
    if (row == trackHit->GetRow() && layer == trackHit->GetLayer())
      continue; // XXX
    KBVector3 qosTrackHit(trackHit->GetPosition(), fReferenceAxis);
    auto di = qosHit.I() - qosTrackHit.I();
    auto dj = qosHit.J() - qosTrackHit.J();
    auto distPadCenter = sqrt(di * di + dj * dj);
    Double_t tangentDip = abs(qosTrackHit.K()) / sqrt(qosTrackHit.I() * qosTrackHit.I() + qosTrackHit.J() * qosTrackHit.J());
    Double_t dkInExpectedTrackPath = 1.2 * distPadCenter * tangentDip;
    Double_t dkBetweenTwoHits = abs(qosHit.K() - qosTrackHit.K());

    if (dkInExpectedTrackPath < fCutdkInExpectedTrackPath)
      dkInExpectedTrackPath = fCutdkInExpectedTrackPath;
    if (dkBetweenTwoHits < dkInExpectedTrackPath)
    {
      passCutdk = true;
    }
  }

  if (!passCutdk)
  {
    return 0;
  }

  if (track->IsBad())
  {
    quality = 1;
  }
  else if (track->IsLine())
  {
    KBVector3 perp = track->PerpLine(hit->GetPosition());

    Double_t rmsCut = track->GetRMST();
    if (rmsCut < fTrackHCutLL)
      rmsCut = fTrackHCutLL;
    if (rmsCut > fTrackHCutHL)
      rmsCut = fTrackHCutHL;
    rmsCut = 3 * rmsCut;

    if (perp.K() > rmsCut)
    {
      quality = 0;
    }
    else
    {
      perp.SetK(0);
      auto magcut = 15.;
      if (perp.Mag() < magcut)
      // if (perp.Mag() < 10*pos1.K()/sqrt(pos1.Mag()))
      {
        quality = 1;
      }
    }
  }
  else if (track->IsPlane())
  {
    Double_t dist = (track->PerpPlane(hit->GetPosition())).Mag();

    Double_t rmsCut = track->GetRMST();
    if (rmsCut < fTrackHCutLL)
      rmsCut = fTrackHCutLL;
    if (rmsCut > fTrackHCutHL)
      rmsCut = fTrackHCutHL;
    rmsCut = 3 * rmsCut;

    if (dist < rmsCut)
    {
      quality = 1;
    }
  }
  else
  {
  }

  return quality;
}

double LHHelixTrackFindingTask::CorrelateHitWithTrack(KBHelixTrack *track, KBTpcHit *hit, Double_t rScale)
{
  Double_t scale = rScale * fDefaultScale;
  Double_t trackLength = track->TrackLength();
  if (trackLength < 500.)
    scale = scale + (500. - trackLength) / 500.;

  auto trackHCutLL = fTrackHCutLL;
  auto trackHCutHL = fTrackHCutHL;
  auto trackWCutLL = fTrackWCutLL;
  auto trackWCutHL = fTrackWCutHL;

  Double_t rmsWCut = track->GetRMSR();

  if (rmsWCut < trackWCutLL)
    rmsWCut = trackWCutLL;
  if (rmsWCut > trackWCutHL)
    rmsWCut = trackWCutHL;
  rmsWCut = scale * rmsWCut;

  Double_t rmsHCut = track->GetRMST();

  if (rmsHCut < trackHCutLL)
    rmsHCut = trackHCutLL;
  if (rmsHCut > trackHCutHL)
    rmsHCut = trackHCutHL;
  rmsHCut = scale * rmsHCut;

  TVector3 qHead = track->Map(track->PositionAtHead());
  TVector3 qTail = track->Map(track->PositionAtTail());
  TVector3 q = track->Map(hit->GetPosition());

  if (qHead.Z() > qTail.Z())
  {
    if (CheckHitDistInAlphaIsLargerThanQuarterPi(track, q.Z() - qHead.Z()))
      return 0;
    if (CheckHitDistInAlphaIsLargerThanQuarterPi(track, qTail.Z() - q.Z()))
      return 0;
  }
  else
  {
    if (CheckHitDistInAlphaIsLargerThanQuarterPi(track, q.Z() - qTail.Z()))
      return 0;
    if (CheckHitDistInAlphaIsLargerThanQuarterPi(track, qHead.Z() - q.Z()))
      return 0;
  }

  Double_t dr = abs(q.X());

  Double_t quality = 0;
  if (dr < rmsWCut && abs(q.Y()) < rmsHCut)
    quality = sqrt((dr - rmsWCut) * (dr - rmsWCut)) / rmsWCut;

  return quality;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

bool LHHelixTrackFindingTask::BuildAndConfirmTrack(KBHelixTrack *track, bool &tailToHead)
{
  track->SortHits(!tailToHead);
  auto trackHits = track->GetHitArray();
  Int_t numHits = trackHits->GetNumHits();

  TVector3 q, m;

  Double_t extrapolationLength = 10.;
  for (Int_t iHit = 1; iHit < numHits; iHit++)
  {
    auto trackHit = (KBTpcHit *)trackHits->GetLastHit();

    Double_t quality = CorrelateHitWithTrack(track, trackHit);

    if (quality <= 0)
    {
      track->RemoveHit(trackHit);
      trackHit->RemoveTrackCand(trackHit->GetTrackID());
      Int_t helicity = track->Helicity();
      track->Fit();
      if (helicity != track->Helicity())
        tailToHead = !tailToHead;

      continue;
    }
  }

  extrapolationLength = 0;
  while (AutoBuildByExtrapolation(track, tailToHead, extrapolationLength))
  {
  }

  if (track->GetNumHits() < fCutMinNumHitsFinalTrack)
    return false;

  return true;
}

bool LHHelixTrackFindingTask::AutoBuildByExtrapolation(KBHelixTrack *track, bool &buildHead, Double_t &extrapolationLength)
{
  if (track->GetNumHits() < fCutMinNumHitsFinalTrack)
    return false;

  TVector3 p;
  if (buildHead)
    p = track->ExtrapolateHead(extrapolationLength);
  else
    p = track->ExtrapolateTail(extrapolationLength);

  return AutoBuildAtPosition(track, p, buildHead, extrapolationLength);
}

bool LHHelixTrackFindingTask::AutoBuildAtPosition(KBHelixTrack *track, TVector3 p, bool &tailToHead, Double_t &extrapolationLength, Double_t rScale)
{
  KBVector3 p2(p, fReferenceAxis);
  if (fPadPlane->IsInBoundary(p2.I(), p2.J()) == false)
    return false;

  Int_t helicity = track->Helicity();

  Double_t rms = 3 * track->GetRMSR();
  if (rms < 25)
    rms = 25;

  Int_t range = Int_t(rms / 8);
  fPadPlane->PullOutNeighborHits(p2.I(), p2.J(), range, fCandHits);
  fNumCandHits = fCandHits->GetEntriesFast();
  Bool_t foundHit = false;

  if (fNumCandHits != 0)
  {
    fCandHits->SortByCharge(false);

    for (Int_t iHit = 0; iHit < fNumCandHits; iHit++)
    {
      KBTpcHit *candHit = (KBTpcHit *)fCandHits->GetLastHit();
      fCandHits->RemoveLastHit();

      Double_t quality = 0;
      if (CheckParentTrackID(candHit) < 0)
        quality = CorrelateHitWithTrack(track, candHit, rScale);

      if (quality > 0)
      {
        track->AddHit(candHit);
        track->Fit();
        foundHit = true;
      }
      else
        fBadHits->AddHit(candHit);
    }
  }

  if (foundHit)
  {
    extrapolationLength = 10;
    if (helicity != track->Helicity())
      tailToHead = !tailToHead;
  }
  else
  {
    extrapolationLength += 10;
    if (extrapolationLength > 3 * track->TrackLength())
    {
      return false;
    }
  }

  return true;
}

int LHHelixTrackFindingTask::CheckParentTrackID(KBTpcHit *hit)
{
  vector<Int_t> *candTracks = hit->GetTrackCandArray();
  Int_t numCandTracks = candTracks->size();
  if (numCandTracks == 0)
    return -2;

  Int_t trackID = -1;
  for (Int_t iCand = 0; iCand < numCandTracks; ++iCand)
  {
    Int_t candTrackID = candTracks->at(iCand);
    if (candTrackID != -1)
    {
      trackID = candTrackID;
    }
  }

  return trackID;
}

bool LHHelixTrackFindingTask::CheckTrackQuality(KBHelixTrack *track)
{
  if (track->GetHelixRadius() < fCutMinHelixRadius)
    return false;

  return true;
}

double LHHelixTrackFindingTask::CheckTrackContinuity(KBHelixTrack *track)
{
  Int_t numHits = track->GetNumHits();
  if (numHits < 2)
    return -1;

  track->SortHits();

  Double_t total = 0;
  Double_t continuous = 0;

  TVector3 qPre, mPre, qCur, mCur; // I need q(position on helix)

  KBVector3 kqPre;
  KBVector3 kqCur;

  auto trackHits = track->GetHitArray();
  auto pPre = trackHits->GetHit(0)->GetPosition();
  track->ExtrapolateByMap(pPre, qPre, mPre);
  kqPre = KBVector3(qPre, fReferenceAxis);

  auto axis1 = fPadPlane->GetAxis1();
  auto axis2 = fPadPlane->GetAxis2();

  for (auto iHit = 1; iHit < numHits; iHit++)
  {
    auto pCur = trackHits->GetHit(iHit)->GetPosition();
    track->ExtrapolateByMap(pCur, qCur, mCur);

    kqCur = KBVector3(qCur, fReferenceAxis);

    auto val1 = kqCur.At(axis1) - kqPre.At(axis1);
    auto val2 = kqCur.At(axis2) - kqPre.At(axis2);

    auto length = sqrt(val1 * val1 + val2 * val2);

    total += length;
    if (length <= 1.2 * fPadPlane->PadDisplacement())
      continuous += length;

    kqPre = kqCur;
  }

  return continuous / total;
}

bool LHHelixTrackFindingTask::CheckHitDistInAlphaIsLargerThanQuarterPi(KBHelixTrack *track, Double_t dLength)
{
  if (dLength > 0)
  {
    if (dLength > .5 * track->TrackLength())
    {
      if (abs(track->AlphaAtTravelLength(dLength)) > .5 * TMath::Pi())
      {
        return true;
      }
    }
  }
  return false;
}
