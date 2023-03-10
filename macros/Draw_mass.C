#include "BSHelper.cxx"
#include "Filipad2.h"

enum {
  kECbegin = 0,
  kDATA = 1,
  kINEL,
  kECend
};
enum {
  kTrigbegin = 0,
  kMBAND = 1,
  kTrigend
};
enum {
  kSpeciesbegin = 0,
  kK0short = 1,
  kLambda,
  kAntilambda,
  kSpeciesend
};
enum {
  kSignbegin = 0,
  kPositive = 1,
  kNegative,
  kSignend
};
enum {
  kStepbegin = 0,
  kAll = 1,
  kBasiccut,
  kMasscut,
  kStepend
};

const char *MCFile = "~/cernbox/O2Mult/mc/AnalysisResults_mc_pbpb536_test_0201.root";
const char *DataFile = "~/cernbox/O2Mult/data/AnalysisResults_data_pbpb536_test_0201.root";
const char *MCdir = "multiplicity-counter";
const char *Datadir = "multiplicity-counter";
const char *MCDataset = "LHC22k3b2";
const char *DataDataset = "LHC22s_pass4";

Int_t n=0;
Int_t nn=0;
const char *mstring = "d#it{N}_{ch}/d#it{#eta}";
const char *avmstring = "LT d#it{N}_{ch}/d#it{#eta} #GT";
const char *avmstringincl = "#LT d#it{N}_{ch}/d#it{#eta} #GT / #LT d#it{N}_{ch}/d#it{#eta} #GT_{Inclusive}";
const char *mstringincl = "(d#it{N}_{ch}/d#it{#eta})/(d#it{N}_{ch}/d#it{#eta})_{Inclusive}";

void setpad(TVirtualPad *pad)
{
  pad->SetTopMargin(0.02);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.2);
  pad->SetBottomMargin(0.15);
  pad->SetName(Form("c%d", ++n));
}

void hset(TH1 &hid, TString xtit = "", TString ytit = "",
          double titoffx = 0.9, double titoffy = 1.2,
          double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05,
          int divx = 510, int divy = 510)
{
	hid.SetStats(0);

	hid.GetXaxis()->CenterTitle(1);
	hid.GetYaxis()->CenterTitle(1);

	hid.GetXaxis()->SetTitleOffset(titoffx);
	hid.GetYaxis()->SetTitleOffset(titoffy);

	hid.GetXaxis()->SetTitleSize(titsizex);
	hid.GetYaxis()->SetTitleSize(titsizey);

	hid.GetXaxis()->SetLabelOffset(labeloffx);
	hid.GetYaxis()->SetLabelOffset(labeloffy);

	hid.GetXaxis()->SetLabelSize(labelsizex);
	hid.GetYaxis()->SetLabelSize(labelsizey);

	hid.GetXaxis()->SetNdivisions(divx);
	hid.GetYaxis()->SetNdivisions(divy);

	hid.GetXaxis()->SetTitle(xtit);
	hid.GetYaxis()->SetTitle(ytit);
}

void Draw_mass() {

// MC File open and find objects
    auto fmc = TFile::Open(MCFile, "oepn");
    auto dir = (TDirectory *) fmc->Get(MCdir);
    dir->cd();

    // reconstructed z-vertex (MC)
    auto HRZ = (THnSparse *) gROOT->FindObject("hreczvtx");
    auto Hreczvtx = BSTHnSparseHelper(HRZ);
    auto hreczvtx = Hreczvtx.GetTH1("hrecz", 3, {kINEL, kMBAND, -1, -1});
    Double_t nrecevent = hreczvtx->Integral(hreczvtx->GetXaxis()->FindBin(-10), hreczvtx->GetXaxis()->FindBin(10));

    // V0Mass (MC)
    auto HV0MASS = (THnSparse *) gROOT->FindObject("hV0Mass");
    auto HV0Mass = BSTHnSparseHelper(HV0MASS);
    auto hK0shortMass = HV0Mass.GetTH1("hK0short", 2, {kINEL, kK0short, 2});
    auto hLambdaMass = HV0Mass.GetTH1("hLambda", 2, {kINEL, kLambda, -1});
    auto hAntiLambdaMass = HV0Mass.GetTH1("hAntiLambda", 2, {kINEL, kAntilambda, -1});

// Data file open and find objects
    auto fdata = TFile::Open(DataFile, "open");
    auto Ddir = (TDirectory *) fdata->Get(Datadir);
    Ddir->cd();

    // reconstructed zvtx and dndeta (Data)
    auto DHRZ = (THnSparse *) gROOT->FindObject("hreczvtx");
    auto DHreczvtx = BSTHnSparseHelper(DHRZ);
    auto Dhreczvtx = DHreczvtx.GetTH1("data_hrecz", 3, {kDATA, kMBAND, -1, -1});
    Double_t Dnrecevent = Dhreczvtx->Integral(Dhreczvtx->GetXaxis()->FindBin(-10), Dhreczvtx->GetXaxis()->FindBin(10));

    // V0Mass (Data)
    auto DHV0MASS = (THnSparse *) gROOT->FindObject("hV0Mass");
    auto DHV0Mass = BSTHnSparseHelper(DHV0MASS);
    auto DhK0shortMass = DHV0Mass.GetTH1("data_hK0short", 2, {kDATA, kK0short, 2});
    auto DhLambdaMass = DHV0Mass.GetTH1("data_hLambda", 2, {kDATA, kLambda, -1});
    auto DhAntiLambdaMass = DHV0Mass.GetTH1("data_hAntiLambda", 2, {kDATA, kAntilambda, -1});
//--------------------------------------------------------------------------------------------------------------------//

// draw K0shortMass

    auto canvas1 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas1->Draw();
    TPad *K0shortPad = canvas1->GetPad(1);
    K0shortPad->SetTickx();
    K0shortPad->SetGridx();
    K0shortPad->SetGridy();
    K0shortPad->SetTopMargin(0.06);
    K0shortPad->cd();

    // K0shortMass (MC)
    hK0shortMass->Scale(1. / nrecevent, "width");
    hK0shortMass->SetLineColor(kBlack);
    hK0shortMass->SetMarkerColor(kBlack);
    hK0shortMass->SetMarkerStyle(20);
    hK0shortMass->SetMarkerSize(0.4);
    hK0shortMass->GetYaxis()->SetMaxDigits(3);
    hK0shortMass->GetYaxis()->SetRangeUser(-2e3, 13e3);
    hset(*hK0shortMass, "K^{0}_{s} Mass (GeV/c^2)", "Normalized tracks", 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hK0shortMass->Draw();

    // K0shortMass (Data)
    DhK0shortMass->Scale(1. / Dnrecevent, "width");
    DhK0shortMass->SetLineColor(kRed);
    DhK0shortMass->SetMarkerColor(kRed);
    DhK0shortMass->SetMarkerStyle(24);
    DhK0shortMass->SetMarkerSize(0.4);
    DhK0shortMass->GetYaxis()->SetMaxDigits(3);
    DhK0shortMass->Draw("same");

    // legend of K0shortMass histogram
    TLegend *legend_K0short = new TLegend(0.51, 0.74, 0.8, 0.92, NULL, "brNDC");
    legend_K0short->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_K0short->AddEntry(hK0shortMass, "MC", "lep");
    legend_K0short->AddEntry(DhK0shortMass, "Data", "lep");
    legend_K0short->SetTextSize(0.04);
    legend_K0short->SetBorderSize(0);
    legend_K0short->SetTextAlign(22);
    legend_K0short->Draw();

    // Count the # of the V0 particles in the peack with bg
    float DataK0Count = DhK0shortMass->Integral(0.482, 0.509);
    float binwidthDATAK0Count = DataK0Count->GetXaxis()->GetBinWidth(1);
    float DataK0shortMassPeakCount = DataK0Count / binwidthDATAK0Count;
    std::cout << DataK0shortMassPeakCount << std::endl;
    TLegend *lk0 = new TLegend(0.51, 0.74, 0.8, 0.92, NULL, "brNDC");
    lk0->AddEntry(DhK0shortMass, Form("# of ptls: %.1f",DataK0shortMassPeakCount),"");
    lk0->SetBorderSize(0);
    lk0->SetTextAlign(12);
    lk0->SetTextSize(0.07);
    lk0->SetMargin(-0.04);
    lk0->Draw("same");

    // draw ratio of K0shortMass (MC/Data)
    TPad *K0short_ratioPad = canvas1->GetPad(2);
    K0short_ratioPad->SetTickx();
    K0short_ratioPad->SetGridx();
    K0short_ratioPad->SetGridy();
    K0short_ratioPad->cd();

    auto K0shortRatio = (TH1D *) hK0shortMass->Clone();
    K0shortRatio->Divide(DhK0shortMass);
    K0shortRatio->SetLineColor(kRed);
    K0shortRatio->SetMarkerColor(kRed);
    K0shortRatio->SetMarkerStyle(24);
    K0shortRatio->SetMarkerSize(0.4);
    K0shortRatio->GetYaxis()->SetMaxDigits(10);
    K0shortRatio->GetYaxis()->SetRangeUser(-5, 50);
    hset(*K0shortRatio, "K^{0}_{s} Mass (GeV/c^{2})", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    K0shortRatio->Draw();

// draw LambdaMass

    auto canvas2 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas2->Draw();
    TPad *LambdaMassPad = canvas2->GetPad(1);
    LambdaMassPad->SetTickx();
    LambdaMassPad->SetGridx();
    LambdaMassPad->SetGridy();
    LambdaMassPad->SetTopMargin(0.06);
    LambdaMassPad->cd();

    // LambdaMass (MC)
    hLambdaMass->Scale(1. / nrecevent, "width");
    hLambdaMass->SetLineColor(kBlack);
    hLambdaMass->SetMarkerColor(kBlack);
    hLambdaMass->SetMarkerStyle(20);
    hLambdaMass->SetMarkerSize(0.4);
    hLambdaMass->GetYaxis()->SetMaxDigits(3);
    hLambdaMass->GetXaxis()->SetRangeUser(1.05,1.3);
    hLambdaMass->GetYaxis()->SetRangeUser(-2e3,12e3);
    hset(*hLambdaMass, "#Lambda Mass (GeV/c^{2})", "Normalized tracks", 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hLambdaMass->Draw();

    // LambdaMass (Data)
    DhLambdaMass->Scale(1. / Dnrecevent, "width");
    DhLambdaMass->SetLineColor(kRed);
    DhLambdaMass->SetMarkerColor(kRed);
    DhLambdaMass->SetMarkerStyle(24);
    DhLambdaMass->SetMarkerSize(0.4);
    DhLambdaMass->GetYaxis()->SetMaxDigits(3);
    DhLambdaMass->Draw("same");

    // legend of LambdaMass histogram
    TLegend *legend_Lambda = new TLegend(0.51, 0.74, 0.8, 0.92,NULL,"brNDC");
    legend_Lambda->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_Lambda->AddEntry(hLambdaMass, "MC", "lep");
    legend_Lambda->AddEntry(DhLambdaMass, "Data", "lep");
    legend_Lambda->SetTextSize(0.04);
    legend_Lambda->SetBorderSize(0);
    legend_Lambda->SetTextAlign(22);
    legend_Lambda->Draw();

    // draw ratio of LambdaMass (MC/Data)
    TPad *Lambda_ratioPad = canvas2->GetPad(2);
    Lambda_ratioPad->SetTickx();
    Lambda_ratioPad->SetGridx();
    Lambda_ratioPad->SetGridy();
    Lambda_ratioPad->cd();

    auto LambdaRatio = (TH1D *) hLambdaMass->Clone();
    LambdaRatio->Divide(DhLambdaMass);
    LambdaRatio->SetLineColor(kRed);
    LambdaRatio->SetMarkerColor(kRed);
    LambdaRatio->SetMarkerStyle(24);
    LambdaRatio->SetMarkerSize(0.4);
    LambdaRatio->GetYaxis()->SetMaxDigits(10);
    LambdaRatio->GetYaxis()->SetRangeUser(0,20);
    hset(*LambdaRatio, "#Lambda Mass (GeV/c^{2})", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    LambdaRatio->Draw();

  // draw AntiLambdaMass

    auto canvas3 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas3->Draw();
    TPad *AntiLambdaMassPad = canvas3->GetPad(1);
    AntiLambdaMassPad->SetTickx();
    AntiLambdaMassPad->SetGridx();
    AntiLambdaMassPad->SetGridy();
    AntiLambdaMassPad->SetTopMargin(0.06);
    AntiLambdaMassPad->cd();

    // AntiLambdaMass (MC)
    hAntiLambdaMass->Scale(1. / nrecevent, "width");
    hAntiLambdaMass->SetLineColor(kBlack);
    hAntiLambdaMass->SetMarkerColor(kBlack);
    hAntiLambdaMass->SetMarkerStyle(20);
    hAntiLambdaMass->SetMarkerSize(0.4);
    hAntiLambdaMass->GetYaxis()->SetMaxDigits(3);
    hAntiLambdaMass->GetXaxis()->SetRangeUser(1.05, 1.3);
    hAntiLambdaMass->GetYaxis()->SetRangeUser(-2e3,12e3);
    hset(*hAntiLambdaMass, "#Lambda Mass (GeV/c^{2})", "Normalized tracks", 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hAntiLambdaMass->Draw();

    // AntiLambdaMass (Data)
    DhAntiLambdaMass->Scale(1. / Dnrecevent, "width");
    DhAntiLambdaMass->SetLineColor(kRed);
    DhAntiLambdaMass->SetMarkerColor(kRed);
    DhAntiLambdaMass->SetMarkerStyle(24);
    DhAntiLambdaMass->SetMarkerSize(0.4);
    DhAntiLambdaMass->GetYaxis()->SetMaxDigits(3);
    DhAntiLambdaMass->Draw("same");

    // legend of AntiLambdaMass histogram
    TLegend *legend_AntiLambda = new TLegend(0.51, 0.74, 0.8, 0.92, NULL, "brNDC");
    legend_AntiLambda->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_AntiLambda->AddEntry(hAntiLambdaMass, "MC", "lep");
    legend_AntiLambda->AddEntry(DhAntiLambdaMass, "Data", "lep");
    legend_AntiLambda->SetTextSize(0.04);
    legend_AntiLambda->SetBorderSize(0);
    legend_AntiLambda->SetTextAlign(22);
    legend_AntiLambda->Draw();

    // draw ratio of LambdaMass (MC/Data)
    TPad *AntiLambda_ratioPad = canvas3->GetPad(2);
    AntiLambda_ratioPad->SetTickx();
    AntiLambda_ratioPad->SetGridx();
    AntiLambda_ratioPad->SetGridy();
    AntiLambda_ratioPad->cd();

    auto AntiLambdaRatio = (TH1D *) hAntiLambdaMass->Clone();
    AntiLambdaRatio->Divide(DhAntiLambdaMass);
    AntiLambdaRatio->SetLineColor(kRed);
    AntiLambdaRatio->SetMarkerColor(kRed);
    AntiLambdaRatio->SetMarkerStyle(24);
    AntiLambdaRatio->SetMarkerSize(0.4);
    AntiLambdaRatio->GetYaxis()->SetMaxDigits(10);
    AntiLambdaRatio->GetYaxis()->SetRangeUser(0,20);
    hset(*AntiLambdaRatio, "#bar{#Lambda} Mass (GeV/c^{2})", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    AntiLambdaRatio->Draw();

}