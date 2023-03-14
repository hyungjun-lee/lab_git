#include "FTHit.hh"

FTHit::FTHit()
{

}

FTHit::~FTHit()
{

}

void FTHit::SetFTHit(int _moduleID, double _reco_x, double _reco_y, double _reco_z)
{
	KBHit::SetDetID(_moduleID);
	KBHit::SetX(_reco_x);
	KBHit::SetY(_reco_y);
	KBHit::SetZ(_reco_z);
}
