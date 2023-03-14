// #include "KBHit.hh"
#include "KBTpcHit.hh"
#include "KBTask.hh"

class FTHit : public KBTpcHit // , public KBTpcHit // , public KBTask
{
	public :
		FTHit();
		virtual ~FTHit();
		void SetFTHit(int _moduleID, double _reco_x, double _reco_y, double _reco_z);

	private :
		int moduleID;
		double reco_x;
		double reco_y;
		double reco_z;
};
