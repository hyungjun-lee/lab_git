#include "LHFTHitTask.hh"
#include "FTHit.hh"
#include "KBRun.hh"

#include <vector>
#include <iostream>
#include <TH3D.h>
using namespace std;

LHFTHitTask::LHFTHitTask()
	: KBTask("LHFTHitTask", "")
{
}

LHFTHitTask::~LHFTHitTask()
{
}

bool LHFTHitTask::Init()
{
	KBRun *run = KBRun::GetRun();

	fMCStepArray = (TClonesArray *)run->GetBranch(Form("MCStep%d", fDetID));
	fStepArray = new TClonesArray("KBMCStep");
	run->RegisterBranch("FTStep", fStepArray, fPersistency);

	fFTHitArray = new TClonesArray("FTHit");
	run->RegisterBranch("FTHit", fFTHitArray, fPersistency);

	return true;
}

void LHFTHitTask::Exec(Option_t *)
{
	fStepArray -> Delete();
	fFTHitArray-> Delete();
	
	fStepArray -> Clear();
	fFTHitArray-> Clear();
	
	Long64_t nMCSteps = fMCStepArray->GetEntries();

	kb_info << " Number of found steps in FT : " << nMCSteps << endl;
	
	for (Long64_t iStep = 0; iStep < nMCSteps; iStep++)
	{
		KBMCStep *step = (KBMCStep *)fMCStepArray->At(iStep);
		double FT_step_x = step->GetX();
		double FT_step_y = step->GetY();
		double FT_step_z = step->GetZ();
		double FT_step_edep = step->GetEdep();
		double FT_step_time = step->GetTime();
		int FT_step_moduleID = step->GetModuleID();
		double FT_step_trackID = step->GetTrackID();
		double FT_step_copyNo = step->GetcopyNo();

		KBMCStep * SingleStep = new ((*fStepArray)[iStep]) KBMCStep();
		SingleStep->SetMCStep(FT_step_trackID, FT_step_moduleID, FT_step_x, FT_step_y, FT_step_z, FT_step_time, FT_step_edep);
		
		// if(FT_step_trackID != 1)
		// {
		// 	continue;
		// }	
		
		if (RAWSTEP)
		{
			FTHit *SingleHit = (FTHit *)fFTHitArray->ConstructedAt(iStep);
			SingleHit->SetFTHit(fDetID, FT_step_x, FT_step_y, FT_step_z);
			SingleHit->SetHitID(40000+iStep);
			// SingleHit->SetX(FT_step_x);
			// SingleHit->SetY(FT_step_y);
			// SingleHit->SetZ(FT_step_z);
			// SingleHit->SetHitID(iStep * 10);
			// SingleHit->SetCharge(1);
			// SingleHit->SetTrackID(FT_step_trackID);

		}
	}
}

void LHFTHitTask::SetDetID(int id)
{
	fDetID = id;
}