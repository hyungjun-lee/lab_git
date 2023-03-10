int run_g4sim(const char *name);

void run(const char *name="g4event"){

	const bool bG4SIM = false;
	const bool bDIGI = true;
	const bool bRECO = true;

	if ( bG4SIM ){
		if ( !run_g4sim(name) ) return;
	}

	if ( bDIGI || bRECO ){

		auto run = KBRun::GetRun();
		run->SetIOFile(Form("out_%s_LH.mc", name), Form("out_%s_LH.conv", name));
		run->AddDetector(new LHTpc());

		auto fthit = new LHFTHitTask();
		fthit->SetDetID(40); // FT
		run->Add(fthit);

		if ( bDIGI ){
			auto drift = new LHDriftElectronTask();
			drift->SetPadPersistency(true);
			drift->SetDetID(10); //TPC
			run->Add(drift);

			auto electronics = new LHElectronicsTask(true);
			run->Add(electronics);
		}


		if ( bDIGI && bRECO ){
			auto psa = new KBPSATask();
			psa -> SetInputBranchName("TPCPad");
			psa -> SetOutputBranchName("TPCHit");
			psa -> SetPSA(new KBPSAFastFit());
			run->Add(psa);

			auto htfTask = new LHHelixTrackFindingTask();
			htfTask -> SetHitBranchName("TPCHit");
			htfTask -> SetHitBranchName_FT("FTHit");
			htfTask -> SetTrackletBranchName("Tracklet");
			run->Add(htfTask);

			// auto gfTask = new LHGenfitTask();
			// gfTask->SetDetID(10); //TPC
			// run->Add(gfTask);

			// run->Add(new LHNDFastDigiTask());

			/*
			//run->Add(new LHVertexFindingTask());
			*/
		}

		run->Init();
		run->Print();
		run->Run();
	}


}

int run_g4sim(const char *name = "LH")
{

	TString fname = Form("kbpar_%s.conf",name);
	gSystem->FindFile("./",fname);

	if ( fname.Length()==0 ){
		cout << "Can not find the macrofile" << endl;
		return 0;
	}else{
		fname = Form("kbpar_%s.conf",name);
		gSystem->Exec(Form("./run_mc.g4sim %s",fname.Data()));
	}

	return 1;
}

