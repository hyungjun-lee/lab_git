#include "KBTask.hh"
#include "TClonesArray.h"

#include <TH3D.h>

class LHFTHitTask : public KBTask
{
public:
	LHFTHitTask();
	virtual ~LHFTHitTask();

	bool Init();
	void Exec(Option_t *);

	void SetPersistency(bool val);
	void SetDetID(int id);

private:
	TClonesArray *fMCStepArray;
	TClonesArray *fStepArray;
	TClonesArray *fFTHitArray;
	TClonesArray *fFTHit;

	bool RECONSTRUCTION = false;
	bool RAWSTEP = true;

	bool fPersistency = true;
	int fDetID = 40;
	int layerN;

};
