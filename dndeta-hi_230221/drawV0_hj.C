#include "BSHelper.cxx"
#include "Filipad2.h"

// #define DrawMass
// #define DrawFitting
#define DrawMotherV0

enum
{
    kECbegin = 0,
    kDATA = 1,
    kINEL,
    kECend
};
enum
{
    kTrigbegin = 0,
    kMBAND = 1,
    kTrigend
};
enum
{
    kSpeciesbegin = 0,
    kK0short = 1,
    kLambda,
    kAntilambda,
    kSpeciesend
};
enum
{
    kSignbegin = 0,
    kPositive = 1,
    kNegative,
    kSignend
};
double Par[6] = {
    0,
};

/// @fit function
Double_t FitGaussian(Double_t *x, Double_t *par)
{
    return TMath::Gaus(x[0], par[0], par[1]) * par[2];
}

Double_t FitPol1(Double_t *x, Double_t *par)
{
    return -1 * x[0] * par[0] + par[1];
}

Double_t FitPol2(Double_t *x, Double_t *par)
{
    return -1 * x[0] * x[0] * par[0] + par[1] * x[1] + par[2];
}

Double_t FitBreitWigenr(Double_t *x, Double_t *par)
{
    Double_t dm = x[0] - par[1];
    return par[0] / (dm * dm + par[2] * par[2]);
}

Double_t FitCrystalBall(Double_t *x, Double_t *par)
{
    double *vx = x;
    double *vpar = par;
    return ROOT::Math::crystalball_function(*vx * (-1), par[0], par[1], par[2], par[3]) * par[4];
}

Double_t FitFuncGausPol1(Double_t *x, Double_t *par)
{
    return FitGaussian(x, par) + FitPol1(x, &par[3]);
}

Double_t FitFuncGausPol2(Double_t *x, Double_t *par)
{
    return FitGaussian(x, par) + FitPol1(x, &par[3]);
}

Double_t FitGausBreitWigenr(Double_t *x, Double_t *par)
{
    return FitGaussian(x, par) + FitBreitWigenr(x, &par[3]);
}

Double_t FitDoubleGaus(Double_t *x, Double_t *par)
{
    return FitGaussian(x, par) + FitGaussian(x, &par[3]);
}

void drawV0_hj()
{
    /// @root file load - - - - - - - - - - - - - - - -
    auto mc_file = TFile::Open("../mc/LHC22k3b/hyper/AnalysisResults_mc_hyper_PbPb_LHC22k3b2_with_V0Cuts.root", "open");
    auto data_file = TFile::Open("../data/LHC22s/hyper/AnalysisResults_data_hyper_PbPb_LHC22s_apass4_with_v0Cuts.root", "open");

    /// - - - - - - - - - - - - - - - - - - - - - - - -

    /// @mc object load - - - - - - - - - - - - - - - -
    auto mc_dir = (TDirectory *)mc_file->Get("multiplicity-counter/Tracks/ProcessMCCounting");
    mc_dir->cd();

    auto mc_RECZVTX = (THnSparse *)gROOT->FindObject("hreczvtx"); // mc reconstructed z-vertex
    auto mc_reczvtx = BSTHnSparseHelper(mc_RECZVTX);

    auto mc_V0MASS = (THnSparse *)gROOT->FindObject("hV0Mass"); // mc V0 particle mass
    auto mc_v0mass = BSTHnSparseHelper(mc_V0MASS);

    auto mc_V0DAUETA = (THnSparse *)gROOT->FindObject("hV0DauEta"); // mc V0 daughter eta
    auto mc_v0daueta = BSTHnSparseHelper(mc_V0DAUETA);

    auto mc_V0COUNT = (THnSparse *)gROOT->FindObject("hV0Count"); // mc V0 particle count
    auto mc_v0count = BSTHnSparseHelper(mc_V0COUNT);

    auto mc_MV0COUNT = (THnSparse *)gROOT->FindObject("hMotherV0Count"); // mc V0 particle count
    auto mc_mv0count = BSTHnSparseHelper(mc_MV0COUNT);
//     /// - - - - - - - - - - - - - - - - - - - - - - - - -

//     /// @data object load - - - - - - - - - - - - - - - -
//     // auto data_dir = (TDirectory *)data_file->Get("multiplicity-counter_id2661");
//     // data_dir->cd();

//     // auto data_RECZVTX = (THnSparse *)gROOT->FindObject("hreczvtx"); // data reconstructed z-vertex
//     // auto data_reczvtx = BSTHnSparseHelper(data_RECZVTX);

//     // auto data_V0MASS = (THnSparse *)gROOT->FindObject("hV0Mass"); // data V0 particle mass
//     // auto data_v0mass = BSTHnSparseHelper(data_V0MASS);

//     // auto data_V0DAUETA = (THnSparse *)gROOT->FindObject("hV0DauEta"); // data V0 daughter eta
//     // auto data_v0daueta = BSTHnSparseHelper(data_V0DAUETA);

//     // auto data_V0COUNT = (THnSparse *)gROOT->FindObject("hV0Count"); // data V0 particle count
//     // auto data_v0count = BSTHnSparseHelper(data_V0COUNT);
//     /// - - - - - - - - - - - - - - - - - - - - - - - -

//     /// @make histogram - - - - - - - - - - - - - - - -
    auto HistMCreczvtx =mc_reczvtx.GetTH1("HistMCreczvtx", 2, {kINEL, kMBAND, -1});

    auto HistMCK0ShortDauEta = mc_v0daueta.GetTH1("HistMCK0ShortDauEta", 2, {kINEL, -1, 1,-1});

//     auto HistMCK0Shortmass = mc_v0mass.GetTH1("HistMCK0Shortmass", 2, {kINEL, kK0short, -1});
//     auto HistMCLambdamass = mc_v0mass.GetTH1("HistMCLambdamass", 2, {kINEL, kLambda, -1});
//     auto HistMCAntilambdamass = mc_v0mass.GetTH1("HistMCAntilambdamass", 2, {kINEL, kAntilambda, -1});

//     auto HistMCK0ShortDauEta = mc_v0daueta.GetTH1("HistMCK0ShortDauEta", 3, {kINEL, -1, kK0short, -1});
//     auto HistMCLambdaDauEta = mc_v0daueta.GetTH1("HistMCLambdaDauEta", 3, {kINEL, -1, kLambda, -1});
//     auto HistMCAntilambdaDauEta = mc_v0daueta.GetTH1("HistMCAntilambdaDauEta", 3, {kINEL, -1, kAntilambda, -1});

//     auto HistMCK0ShortDauEtaPos = mc_v0daueta.GetTH1("HistMCK0ShortDauEtaPos", 3, {kINEL, kPositive, kK0short, -1});
//     auto HistMCK0ShortDauEtaNeg = mc_v0daueta.GetTH1("HistMCK0ShortDauEtaNeg", 3, {kINEL, kNegative, kK0short, -1});
//     auto HistMCLambdaDauEtaPos = mc_v0daueta.GetTH1("HistMCLambdaDauEtaPos", 3, {kINEL, kPositive, kLambda, -1});
//     auto HistMCLambdaDauEtaNeg = mc_v0daueta.GetTH1("HistMCLambdaDauEtaNeg", 3, {kINEL, kNegative, kLambda, -1});
//     auto HistMCAntilambdaDauEtaPos = mc_v0daueta.GetTH1("HistMCAntilambdaDauEtaPos", 3, {kINEL, kPositive, kAntilambda, -1});
//     auto HistMCAntilambdaDauEtaNeg = mc_v0daueta.GetTH1("HistMCAntilambdaDauEtaNeg", 3, {kINEL, kNegative, kAntilambda, -1});

//     // auto HistDATAreczvtx = data_reczvtx.GetTH1("HistDATAreczvtx", 3, {kDATA, kMBAND, -1, -1});

//     // auto HistDATAK0Shortmass = data_v0mass.GetTH1("HistDATAK0Shortmass", 2, {kDATA, kK0short, -1});
//     // auto HistDATALambdamass = data_v0mass.GetTH1("HistDATALambdamass", 2, {kDATA, kLambda, -1});
//     // auto HistDATAAntilambdamass = data_v0mass.GetTH1("HistDATAAntilambdamass", 2, {kDATA, kAntilambda, -1});

//     // auto HistDATAK0ShortDauEta = data_v0daueta.GetTH1("HistDATAK0ShortDauEta", 3, {kDATA, -1, kK0short, -1});
//     // auto HistDATALambdaDauEta = data_v0daueta.GetTH1("HistDATALambdaDauEta", 3, {kDATA, -1, kLambda, -1});
//     // auto HistDATAAntilambdaDauEta = data_v0daueta.GetTH1("HistDATAAntilambdaDauEta", 3, {kDATA, -1, kAntilambda, -1});

//     // auto HistDATAK0ShortDauEtaPos = data_v0daueta.GetTH1("HistDATAK0ShortDauEtaPos", 3, {kDATA, kPositive, kK0short, -1});
//     // auto HistDATAK0ShortDauEtaNeg = data_v0daueta.GetTH1("HistDATAK0ShortDauEtaNeg", 3, {kDATA, kNegative, kK0short, -1});
//     // auto HistDATALambdaDauEtaPos = data_v0daueta.GetTH1("HistDATALambdaDauEtaPos", 3, {kDATA, kPositive, kLambda, -1});
//     // auto HistDATALambdaDauEtaNeg = data_v0daueta.GetTH1("HistDATALambdaDauEtaNeg", 3, {kDATA, kNegative, kLambda, -1});
//     // auto HistDATAAntilambdaDauEtaPos = data_v0daueta.GetTH1("HistDATAAntilambdaDauEtaPos", 3, {kDATA, kPositive, kAntilambda, -1});
//     // auto HistDATAAntilambdaDauEtaNeg = data_v0daueta.GetTH1("HistDATAAntilambdaDauEtaNeg", 3, {kDATA, kNegative, kAntilambda, -1});

    auto HistMCK0ShortV0Count = mc_v0count.GetTH1("HistMCV0Count", 1, {kINEL, -1 , 3});
    auto HistMCLambdaV0Count = mc_v0count.GetTH1("HistMCV0Count", 1, {kINEL, kLambda, 2});
    auto HistMCAntilambdaV0Count = mc_v0count.GetTH1("HistMCV0Count", 1, {kINEL, kAntilambda, 2});

//     // auto HistDATAK0ShortV0Count = data_v0count.GetTH1("HistDATAV0Count", 2, {kDATA, kK0short, -1});
//     // auto HistDATALambdaV0Count = data_v0count.GetTH1("HistDATAV0Count", 2, {kDATA, kLambda, -1});
//     // auto HistDATAAntilambdaV0Count = data_v0count.GetTH1("HistDATAV0Count", 2, {kDATA, kAntilambda, -1});
    
    // auto HistMCMotherV0Count = mc_mv0count.GetTH1("HistMCMotherV0Count",1,{kINEL,-1});

//     /// - - - - - - - - - - - - - - - - - - - - - - -

//     /// @make histogram - - - - - - - - - - - - - - - -
    auto NormalizationValueMC = HistMCreczvtx->Integral(HistMCreczvtx->FindBin(-10), HistMCreczvtx->FindBin(10));
//     // auto NormalizationValueDATA = HistDATAreczvtx->Integral(HistDATAreczvtx->FindBin(-10), HistDATAreczvtx->FindBin(10));

//     HistMCK0Shortmass->Scale(1. / NormalizationValueMC, "width");
//     HistMCLambdamass->Scale(1. / NormalizationValueMC, "width");
//     HistMCAntilambdamass->Scale(1. / NormalizationValueMC, "width");

    HistMCK0ShortDauEta->Scale(1. / NormalizationValueMC, "width");

     HistMCK0ShortV0Count->Scale(1. / NormalizationValueMC, "width");
    HistMCLambdaV0Count->Scale(1. / NormalizationValueMC, "width");
    HistMCAntilambdaV0Count->Scale(1. / NormalizationValueMC, "width");

//     HistMCK0ShortDauEtaPos->Scale(1. / NormalizationValueMC, "width");
//     HistMCK0ShortDauEtaNeg->Scale(1. / NormalizationValueMC, "width");

//     HistMCLambdaDauEtaPos->Scale(1. / NormalizationValueMC, "width");
//     HistMCLambdaDauEtaNeg->Scale(1. / NormalizationValueMC, "width");

//     HistMCAntilambdaDauEtaPos->Scale(1. / NormalizationValueMC, "width");
//     HistMCAntilambdaDauEtaNeg->Scale(1. / NormalizationValueMC, "width");

    // HistMCMotherV0Count->Scale(1. / NormalizationValueMC,"width");
//     // HistDATAK0Shortmass->Scale(1. / NormalizationValueDATA, "width");
//     // HistDATALambdamass->Scale(1. / NormalizationValueDATA, "width");
//     // HistDATAAntilambdamass->Scale(1. / NormalizationValueDATA, "width");

//     // HistDATAK0ShortDauEtaPos->Scale(1. / NormalizationValueDATA, "width");
//     // HistDATAK0ShortDauEtaNeg->Scale(1. / NormalizationValueDATA, "width");

//     // HistDATALambdaDauEtaPos->Scale(1. / NormalizationValueDATA, "width");
//     // HistDATALambdaDauEtaNeg->Scale(1. / NormalizationValueDATA, "width");

//     // HistDATAAntilambdaDauEtaPos->Scale(1. / NormalizationValueDATA, "width");
//     // HistDATAAntilambdaDauEtaNeg->Scale(1. / NormalizationValueDATA, "width");
//     /// - - - - - - - - - - - - - - - - - - - - - - -

// /// @draw histogram - - - - - - - - - - - - - - - -
// TCanvas c99 = new TCanvas("c99","c99",800,600);
// HistMCK0ShortDauEta->Draw("hist][text0");
// auto c1 = new TCanvas("c1","c1",800,600);
gStyle->SetOptStat(0);
HistMCK0ShortV0Count->Draw("hist][text0");
HistMCK0ShortV0Count->SetTitle("V0 count after cut");
HistMCK0ShortV0Count->GetXaxis()->SetTitle("Species");
HistMCK0ShortV0Count->GetYaxis()->SetTitle("Count [#]");
HistMCK0ShortV0Count->GetXaxis()->ChangeLabel(1,-1,-1,-1,-1,-1," ");
HistMCK0ShortV0Count->GetXaxis()->ChangeLabel(3,-1,-1,-1,-1,-1," ");
HistMCK0ShortV0Count->GetXaxis()->ChangeLabel(5,-1,-1,-1,-1,-1," ");
HistMCK0ShortV0Count->GetXaxis()->ChangeLabel(7,-1,-1,-1,-1,-1," ");
HistMCK0ShortV0Count->GetXaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"K^{0}_{s}");
HistMCK0ShortV0Count->GetXaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"#Lambda");
HistMCK0ShortV0Count->GetXaxis()->ChangeLabel(6,-1,-1,-1,-1,-1,"#bar{#Lambda}");
// auto c2 = new TCanvas("c2","c2",800,600);
// HistMCLambdaV0Count->Draw("hist][text0");
// auto c3 = new TCanvas("c3","c3",800,600);
// HistMCAntilambdaV0Count->Draw("hist][text0");

/// mass
#ifdef DrawMass
    auto Filicanvas1 = new Filipad2(1, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas1->Draw();
    TPad *p1_1 = Filicanvas1->GetPad(1);
    p1_1->cd(1);
    p1_1->SetGrid();
    auto HistK0ShortRatio = (TH1D *)HistMCK0Shortmass->Clone();

    HistMCK0Shortmass->Draw();
    HistDATAK0Shortmass->Draw("same");
    HistMCK0Shortmass->GetYaxis()->SetTitle("K^{0}_{s} count[#]");
    HistMCK0Shortmass->GetYaxis()->SetTitleSize(0.05);
    HistMCK0Shortmass->GetYaxis()->CenterTitle();
    HistMCK0Shortmass->SetMarkerStyle(4);
    HistMCK0Shortmass->SetMarkerSize(0.45);
    HistMCK0Shortmass->SetMarkerColor(kRed);
    HistMCK0Shortmass->SetLineColor(kRed);

    HistDATAK0Shortmass->SetMarkerStyle(4);
    HistDATAK0Shortmass->SetMarkerSize(0.45);
    HistDATAK0Shortmass->SetMarkerColor(kBlue);
    HistDATAK0Shortmass->SetLineColor(kBlue);

    TLegend *l1 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l1->AddEntry(HistMCK0Shortmass, "MC : LHC22k3b2", "lp");
    l1->AddEntry(HistDATAK0Shortmass, "DATA : LHC22s pass4", "lp");
    l1->SetBorderSize(0);
    l1->SetTextAlign(12);
    l1->SetTextSize(0.04);
    l1->Draw("same");

    TPad *p1_2 = Filicanvas1->GetPad(2);
    p1_2->cd(2);
    p1_2->SetGrid();
    HistK0ShortRatio->Divide(HistDATAK0Shortmass);
    HistK0ShortRatio->Draw("same");
    HistK0ShortRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistK0ShortRatio->GetYaxis()->SetTitleSize(0.07);
    HistK0ShortRatio->GetYaxis()->CenterTitle();
    HistK0ShortRatio->GetXaxis()->SetTitleSize(0.07);
    HistK0ShortRatio->GetYaxis()->SetLabelSize(0.05);
    HistK0ShortRatio->GetXaxis()->SetLabelSize(0.06);
    HistK0ShortRatio->SetMarkerStyle(4);
    HistK0ShortRatio->SetMarkerSize(0.45);
    HistK0ShortRatio->SetLineColor(kBlack);
    HistK0ShortRatio->SetMarkerColor(kBlack);

    auto Filicanvas2 = new Filipad2(2, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas2->Draw();
    TPad *p2_1 = Filicanvas2->GetPad(1);
    p2_1->cd(1);
    p2_1->SetGrid();
    HistMCLambdamass->GetXaxis()->SetRangeUser(1.05, 1.3);
    auto HistLambdaRatio = (TH1D *)HistMCLambdamass->Clone();

    HistMCLambdamass->Draw();
    HistDATALambdamass->Draw("same");
    HistMCLambdamass->GetYaxis()->SetTitle("#Lambda count[#]");
    HistMCLambdamass->GetYaxis()->SetTitleSize(0.05);
    HistMCLambdamass->GetYaxis()->CenterTitle();
    HistMCLambdamass->GetYaxis()->SetRangeUser(0., 12000);

    HistMCLambdamass->SetMarkerStyle(4);
    HistMCLambdamass->SetMarkerSize(0.45);
    HistMCLambdamass->SetMarkerColor(kRed);
    HistMCLambdamass->SetLineColor(kRed);

    HistDATALambdamass->SetMarkerStyle(4);
    HistDATALambdamass->SetMarkerSize(0.45);
    HistDATALambdamass->SetMarkerColor(kBlue);
    HistDATALambdamass->SetLineColor(kBlue);

    TLegend *l2 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l2->AddEntry(HistMCLambdamass, "MC : LHC22k3b2", "lp");
    l2->AddEntry(HistDATALambdamass, "DATA : LHC22s pass4", "lp");
    l2->SetBorderSize(0);
    l2->SetTextAlign(12);
    l2->SetTextSize(0.04);
    l2->Draw("same");

    TPad *p2_2 = Filicanvas2->GetPad(2);
    p2_2->cd(2);
    p2_2->SetGrid();
    HistLambdaRatio->Divide(HistDATALambdamass);
    HistLambdaRatio->Draw("same");
    HistLambdaRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistLambdaRatio->GetYaxis()->SetTitleSize(0.07);
    HistLambdaRatio->GetYaxis()->CenterTitle();
    HistLambdaRatio->GetXaxis()->SetTitleSize(0.07);
    HistLambdaRatio->GetYaxis()->SetLabelSize(0.05);
    HistLambdaRatio->GetXaxis()->SetLabelSize(0.06);
    HistLambdaRatio->SetMarkerStyle(4);
    HistLambdaRatio->SetMarkerSize(0.45);
    HistLambdaRatio->SetLineColor(kBlack);
    HistLambdaRatio->SetMarkerColor(kBlack);

    auto Filicanvas3 = new Filipad2(3, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas3->Draw();
    TPad *p3_1 = Filicanvas3->GetPad(1);
    p3_1->cd(1);
    p3_1->SetGrid();
    HistMCAntilambdamass->GetXaxis()->SetRangeUser(1.05, 1.3);
    auto HistAntilambdaRatio = (TH1D *)HistMCAntilambdamass->Clone();

    HistMCAntilambdamass->Draw();
    HistDATAAntilambdamass->Draw("same");
    HistMCAntilambdamass->GetYaxis()->SetTitle("#bar{#Lambda} count[#]");
    HistMCAntilambdamass->GetYaxis()->SetTitleSize(0.05);
    HistMCAntilambdamass->GetYaxis()->CenterTitle();
    HistMCAntilambdamass->GetYaxis()->SetRangeUser(0., 12000);
    HistMCAntilambdamass->SetMarkerStyle(4);
    HistMCAntilambdamass->SetMarkerSize(0.45);
    HistMCAntilambdamass->SetMarkerColor(kRed);
    HistMCAntilambdamass->SetLineColor(kRed);

    HistDATAAntilambdamass->SetMarkerStyle(4);
    HistDATAAntilambdamass->SetMarkerSize(0.45);
    HistDATAAntilambdamass->SetMarkerColor(kBlue);
    HistDATAAntilambdamass->SetLineColor(kBlue);

    TLegend *l3 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l3->AddEntry(HistMCAntilambdamass, "MC : LHC22k3b2", "lp");
    l3->AddEntry(HistDATAAntilambdamass, "DATA : LHC22s pass4", "lp");
    l3->SetBorderSize(0);
    l3->SetTextAlign(12);
    l3->SetTextSize(0.04);
    l3->Draw("same");

    TPad *p3_2 = Filicanvas3->GetPad(2);
    p3_2->cd(2);
    p3_2->SetGrid();
    HistAntilambdaRatio->Divide(HistDATAAntilambdamass);
    HistAntilambdaRatio->Draw("same");
    HistAntilambdaRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistAntilambdaRatio->GetYaxis()->SetTitleSize(0.07);
    HistAntilambdaRatio->GetYaxis()->CenterTitle();
    HistAntilambdaRatio->GetXaxis()->SetTitleSize(0.07);
    HistAntilambdaRatio->GetYaxis()->SetLabelSize(0.05);
    HistAntilambdaRatio->GetXaxis()->SetLabelSize(0.06);
    HistAntilambdaRatio->SetMarkerStyle(4);
    HistAntilambdaRatio->SetMarkerSize(0.45);
    HistAntilambdaRatio->SetLineColor(kBlack);
    HistAntilambdaRatio->SetMarkerColor(kBlack);

    ///@Eta
    auto Filicanvas4 = new Filipad2(4, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas4->Draw();
    TPad *p4_1 = Filicanvas4->GetPad(1);
    p4_1->cd(1);
    p4_1->SetGrid();
    HistMCK0ShortDauEtaPos->GetXaxis()->SetRangeUser(-2, 2);
    auto HistK0ShortDauEtaPosRatio = (TH1D *)HistMCK0ShortDauEtaPos->Clone();

    HistMCK0ShortDauEtaPos->Draw();
    HistDATAK0ShortDauEtaPos->Draw("same");
    HistMCK0ShortDauEtaPos->GetYaxis()->SetTitle("K^{0}_{s} daughter positive track count[#]");
    HistMCK0ShortDauEtaPos->GetYaxis()->SetTitleSize(0.05);
    HistMCK0ShortDauEtaPos->GetYaxis()->CenterTitle();
    HistMCK0ShortDauEtaPos->SetMarkerStyle(4);
    HistMCK0ShortDauEtaPos->SetMarkerSize(0.45);
    HistMCK0ShortDauEtaPos->SetMarkerColor(kRed);
    HistMCK0ShortDauEtaPos->SetLineColor(kRed);

    HistDATAK0ShortDauEtaPos->SetMarkerStyle(4);
    HistDATAK0ShortDauEtaPos->SetMarkerSize(0.45);
    HistDATAK0ShortDauEtaPos->SetMarkerColor(kBlue);
    HistDATAK0ShortDauEtaPos->SetLineColor(kBlue);

    TLegend *l4 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l4->AddEntry(HistMCK0ShortDauEtaPos, "MC : LHC22k3b2", "lp");
    l4->AddEntry(HistDATAK0ShortDauEtaPos, "DATA : LHC22s pass4", "lp");
    l4->SetBorderSize(0);
    l4->SetTextAlign(12);
    l4->SetTextSize(0.04);
    l4->Draw("same");

    TPad *p4_2 = Filicanvas4->GetPad(2);
    p4_2->cd(2);
    p4_2->SetGrid();
    HistK0ShortDauEtaPosRatio->Divide(HistDATAK0ShortDauEtaPos);
    HistK0ShortDauEtaPosRatio->Draw("same");
    HistK0ShortDauEtaPosRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistK0ShortDauEtaPosRatio->GetYaxis()->SetTitleSize(0.07);
    HistK0ShortDauEtaPosRatio->GetYaxis()->CenterTitle();
    HistK0ShortDauEtaPosRatio->GetXaxis()->SetTitleSize(0.07);
    HistK0ShortDauEtaPosRatio->GetYaxis()->SetLabelSize(0.05);
    HistK0ShortDauEtaPosRatio->GetXaxis()->SetLabelSize(0.06);
    HistK0ShortDauEtaPosRatio->SetMarkerStyle(4);
    HistK0ShortDauEtaPosRatio->SetMarkerSize(0.45);
    HistK0ShortDauEtaPosRatio->SetLineColor(kBlack);
    HistK0ShortDauEtaPosRatio->SetMarkerColor(kBlack);

    auto Filicanvas5 = new Filipad2(5, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas5->Draw();
    TPad *p5_1 = Filicanvas5->GetPad(1);
    p5_1->cd(1);
    p5_1->SetGrid();
    HistMCK0ShortDauEtaNeg->GetXaxis()->SetRangeUser(-2, 2);
    auto HistK0ShortDauEtaNegRatio = (TH1D *)HistMCK0ShortDauEtaNeg->Clone();

    HistMCK0ShortDauEtaNeg->Draw();
    HistDATAK0ShortDauEtaNeg->Draw("same");
    HistMCK0ShortDauEtaNeg->GetYaxis()->SetTitle("K^{0}_{s} daughter negative track count[#]");
    HistMCK0ShortDauEtaNeg->GetYaxis()->SetTitleSize(0.05);
    HistMCK0ShortDauEtaNeg->GetYaxis()->CenterTitle();
    HistMCK0ShortDauEtaNeg->SetMarkerStyle(4);
    HistMCK0ShortDauEtaNeg->SetMarkerSize(0.45);
    HistMCK0ShortDauEtaNeg->SetMarkerColor(kRed);
    HistMCK0ShortDauEtaNeg->SetLineColor(kRed);

    HistDATAK0ShortDauEtaNeg->SetMarkerStyle(4);
    HistDATAK0ShortDauEtaNeg->SetMarkerSize(0.45);
    HistDATAK0ShortDauEtaNeg->SetMarkerColor(kBlue);
    HistDATAK0ShortDauEtaNeg->SetLineColor(kBlue);

    TLegend *l5 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l5->AddEntry(HistMCK0ShortDauEtaNeg, "MC : LHC22k3b2", "lp");
    l5->AddEntry(HistDATAK0ShortDauEtaNeg, "DATA : LHC22s pass4", "lp");
    l5->SetBorderSize(0);
    l5->SetTextAlign(12);
    l5->SetTextSize(0.04);
    l5->Draw("same");

    TPad *p5_2 = Filicanvas5->GetPad(2);
    p5_2->cd(2);
    p5_2->SetGrid();
    HistK0ShortDauEtaNegRatio->Divide(HistDATAK0ShortDauEtaNeg);
    HistK0ShortDauEtaNegRatio->Draw("same");
    HistK0ShortDauEtaNegRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistK0ShortDauEtaNegRatio->GetYaxis()->SetTitleSize(0.07);
    HistK0ShortDauEtaNegRatio->GetYaxis()->CenterTitle();
    HistK0ShortDauEtaNegRatio->GetXaxis()->SetTitleSize(0.07);
    HistK0ShortDauEtaNegRatio->GetYaxis()->SetLabelSize(0.05);
    HistK0ShortDauEtaNegRatio->GetXaxis()->SetLabelSize(0.06);
    HistK0ShortDauEtaNegRatio->SetMarkerStyle(4);
    HistK0ShortDauEtaNegRatio->SetMarkerSize(0.45);
    HistK0ShortDauEtaNegRatio->SetLineColor(kBlack);
    HistK0ShortDauEtaNegRatio->SetMarkerColor(kBlack);

    auto Filicanvas6 = new Filipad2(6, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas6->Draw();
    TPad *p6_1 = Filicanvas6->GetPad(1);
    p6_1->cd(1);
    p6_1->SetGrid();
    HistMCLambdaDauEtaPos->GetXaxis()->SetRangeUser(-2, 2);
    auto HistLambdaDauEtaPosRatio = (TH1D *)HistMCLambdaDauEtaPos->Clone();

    HistMCLambdaDauEtaPos->Draw();
    HistDATALambdaDauEtaPos->Draw("same");
    HistMCLambdaDauEtaPos->GetYaxis()->SetTitle("#Lambda daughter positive track  count[#]");
    HistMCLambdaDauEtaPos->GetYaxis()->SetTitleSize(0.05);
    HistMCLambdaDauEtaPos->GetYaxis()->CenterTitle();
    HistMCLambdaDauEtaPos->SetMarkerStyle(4);
    HistMCLambdaDauEtaPos->SetMarkerSize(0.45);
    HistMCLambdaDauEtaPos->SetMarkerColor(kRed);
    HistMCLambdaDauEtaPos->SetLineColor(kRed);

    HistDATALambdaDauEtaPos->SetMarkerStyle(4);
    HistDATALambdaDauEtaPos->SetMarkerSize(0.45);
    HistDATALambdaDauEtaPos->SetMarkerColor(kBlue);
    HistDATALambdaDauEtaPos->SetLineColor(kBlue);

    TLegend *l6 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l6->AddEntry(HistMCLambdaDauEtaPos, "MC : LHC22k3b2", "lp");
    l6->AddEntry(HistDATALambdaDauEtaPos, "DATA : LHC22s pass4", "lp");
    l6->SetBorderSize(0);
    l6->SetTextAlign(12);
    l6->SetTextSize(0.04);
    l6->Draw("same");

    TPad *p6_2 = Filicanvas6->GetPad(2);
    p6_2->cd(2);
    p6_2->SetGrid();
    HistLambdaDauEtaPosRatio->Divide(HistDATALambdaDauEtaPos);
    HistLambdaDauEtaPosRatio->Draw("same");
    HistLambdaDauEtaPosRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistLambdaDauEtaPosRatio->GetYaxis()->SetTitleSize(0.07);
    HistLambdaDauEtaPosRatio->GetYaxis()->CenterTitle();
    HistLambdaDauEtaPosRatio->GetXaxis()->SetTitleSize(0.07);
    HistLambdaRatio->GetYaxis()->SetLabelSize(0.05);
    HistK0ShortRatio->GetXaxis()->SetLabelSize(0.06);
    HistLambdaDauEtaPosRatio->SetMarkerStyle(4);
    HistLambdaDauEtaPosRatio->SetMarkerSize(0.45);
    HistLambdaDauEtaPosRatio->SetLineColor(kBlack);
    HistLambdaDauEtaPosRatio->SetMarkerColor(kBlack);

    auto Filicanvas7 = new Filipad2(7, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas7->Draw();
    TPad *p7_1 = Filicanvas7->GetPad(1);
    p7_1->cd(1);
    p7_1->SetGrid();
    HistMCLambdaDauEtaNeg->GetXaxis()->SetRangeUser(-2, 2);
    auto HistLambdaDauEtaNegRatio = (TH1D *)HistMCLambdaDauEtaNeg->Clone();

    HistMCLambdaDauEtaNeg->Draw();
    HistDATALambdaDauEtaNeg->Draw("same");
    HistMCLambdaDauEtaNeg->GetYaxis()->SetTitle("#Lambda daughter negative track count[#]");
    HistMCLambdaDauEtaNeg->GetYaxis()->SetTitleSize(0.05);
    HistMCLambdaDauEtaNeg->GetYaxis()->CenterTitle();
    HistMCLambdaDauEtaNeg->SetMarkerStyle(4);
    HistMCLambdaDauEtaNeg->SetMarkerSize(0.45);
    HistMCLambdaDauEtaNeg->SetMarkerColor(kRed);
    HistMCLambdaDauEtaNeg->SetLineColor(kRed);

    HistDATALambdaDauEtaNeg->SetMarkerStyle(4);
    HistDATALambdaDauEtaNeg->SetMarkerSize(0.45);
    HistDATALambdaDauEtaNeg->SetMarkerColor(kBlue);
    HistDATALambdaDauEtaNeg->SetLineColor(kBlue);

    TLegend *l7 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l7->AddEntry(HistMCLambdaDauEtaNeg, "MC : LHC22k3b2", "lp");
    l7->AddEntry(HistDATALambdaDauEtaNeg, "DATA : LHC22s pass4", "lp");
    l7->SetBorderSize(0);
    l7->SetTextAlign(12);
    l7->SetTextSize(0.04);
    l7->Draw("same");

    TPad *p7_2 = Filicanvas7->GetPad(2);
    p7_2->cd(2);
    p7_2->SetGrid();
    HistLambdaDauEtaNegRatio->Divide(HistDATALambdaDauEtaNeg);
    HistLambdaDauEtaNegRatio->Draw("same");
    HistLambdaDauEtaNegRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistLambdaDauEtaNegRatio->GetYaxis()->SetTitleSize(0.07);
    HistLambdaDauEtaNegRatio->GetYaxis()->CenterTitle();
    HistLambdaDauEtaNegRatio->GetXaxis()->SetTitleSize(0.07);
    HistLambdaDauEtaNegRatio->GetYaxis()->SetLabelSize(0.05);
    HistLambdaDauEtaNegRatio->GetXaxis()->SetLabelSize(0.06);
    HistLambdaDauEtaNegRatio->SetMarkerStyle(4);
    HistLambdaDauEtaNegRatio->SetMarkerSize(0.45);
    HistLambdaDauEtaNegRatio->SetLineColor(kBlack);
    HistLambdaDauEtaNegRatio->SetMarkerColor(kBlack);

    auto Filicanvas8 = new Filipad2(8, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas8->Draw();
    TPad *p8_1 = Filicanvas8->GetPad(1);
    p8_1->cd(1);
    p8_1->SetGrid();
    HistMCAntilambdaDauEtaPos->GetXaxis()->SetRangeUser(-2, 2);

    auto HistAntilambdaDauEtaPosRatio = (TH1D *)HistMCAntilambdaDauEtaPos->Clone();

    HistMCAntilambdaDauEtaPos->Draw();
    HistDATAAntilambdaDauEtaPos->Draw("same");
    HistMCAntilambdaDauEtaPos->GetYaxis()->SetTitle("#bar{#Lambda} daughter positive track count[#]");
    HistMCAntilambdaDauEtaPos->GetYaxis()->SetTitleSize(0.05);
    HistMCAntilambdaDauEtaPos->GetYaxis()->CenterTitle();
    HistMCAntilambdaDauEtaPos->SetMarkerStyle(4);
    HistMCAntilambdaDauEtaPos->SetMarkerSize(0.45);
    HistMCAntilambdaDauEtaPos->SetMarkerColor(kRed);
    HistMCAntilambdaDauEtaPos->SetLineColor(kRed);

    HistDATAAntilambdaDauEtaPos->SetMarkerStyle(4);
    HistDATAAntilambdaDauEtaPos->SetMarkerSize(0.45);
    HistDATAAntilambdaDauEtaPos->SetMarkerColor(kBlue);
    HistDATAAntilambdaDauEtaPos->SetLineColor(kBlue);

    TLegend *l8 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l8->AddEntry(HistMCAntilambdaDauEtaPos, "MC : LHC22k3b2", "lp");
    l8->AddEntry(HistDATAAntilambdaDauEtaPos, "DATA : LHC22s pass4", "lp");
    l8->SetBorderSize(0);
    l8->SetTextAlign(12);
    l8->SetTextSize(0.04);
    l8->Draw("same");

    TPad *p8_2 = Filicanvas8->GetPad(2);
    p8_2->cd(2);
    p8_2->SetGrid();
    HistAntilambdaDauEtaPosRatio->Divide(HistDATAAntilambdaDauEtaPos);
    HistAntilambdaDauEtaPosRatio->Draw("same");
    HistAntilambdaDauEtaPosRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistAntilambdaDauEtaPosRatio->GetYaxis()->SetTitleSize(0.07);
    HistAntilambdaDauEtaPosRatio->GetYaxis()->CenterTitle();
    HistAntilambdaDauEtaPosRatio->GetXaxis()->SetTitleSize(0.07);
    HistAntilambdaDauEtaPosRatio->GetYaxis()->SetLabelSize(0.05);
    HistAntilambdaDauEtaPosRatio->GetXaxis()->SetLabelSize(0.06);
    HistAntilambdaDauEtaPosRatio->SetMarkerStyle(4);
    HistAntilambdaDauEtaPosRatio->SetMarkerSize(0.45);
    HistAntilambdaDauEtaPosRatio->SetLineColor(kBlack);
    HistAntilambdaDauEtaPosRatio->SetMarkerColor(kBlack);

    auto Filicanvas9 = new Filipad2(9, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas9->Draw();
    TPad *p9_1 = Filicanvas9->GetPad(1);
    p9_1->cd(1);
    p9_1->SetGrid();
    HistMCAntilambdaDauEtaNeg->GetXaxis()->SetRangeUser(-2, 2);
    auto HistAntilambdaDauEtaNegRatio = (TH1D *)HistMCAntilambdaDauEtaNeg->Clone();

    HistMCAntilambdaDauEtaNeg->Draw();
    HistDATAAntilambdaDauEtaNeg->Draw("same");
    HistMCAntilambdaDauEtaNeg->GetYaxis()->SetTitle("#bar{#Lambda} daughter negative track count[#]");
    HistMCAntilambdaDauEtaNeg->GetYaxis()->SetTitleSize(0.05);
    HistMCAntilambdaDauEtaNeg->GetYaxis()->CenterTitle();
    HistMCAntilambdaDauEtaNeg->SetMarkerStyle(4);
    HistMCAntilambdaDauEtaNeg->SetMarkerSize(0.45);
    HistMCAntilambdaDauEtaNeg->SetMarkerColor(kRed);
    HistMCAntilambdaDauEtaNeg->SetLineColor(kRed);

    HistDATAAntilambdaDauEtaNeg->SetMarkerStyle(4);
    HistDATAAntilambdaDauEtaNeg->SetMarkerSize(0.45);
    HistDATAAntilambdaDauEtaNeg->SetMarkerColor(kBlue);
    HistDATAAntilambdaDauEtaNeg->SetLineColor(kBlue);

    TLegend *l9 = new TLegend(0.730769, 0.724638, 0.948161, 0.921739, NULL, "brNDC");
    l9->AddEntry(HistMCAntilambdaDauEtaNeg, "MC : LHC22k3b2", "lp");
    l9->AddEntry(HistDATAAntilambdaDauEtaNeg, "DATA : LHC22s pass4", "lp");
    l9->SetBorderSize(0);
    l9->SetTextAlign(12);
    l9->SetTextSize(0.04);
    l9->Draw("same");

    TPad *p9_2 = Filicanvas9->GetPad(2);
    p9_2->cd(2);
    p9_2->SetGrid();
    HistAntilambdaDauEtaNegRatio->Divide(HistDATAAntilambdaDauEtaNeg);
    HistAntilambdaDauEtaNegRatio->Draw("same");
    HistAntilambdaDauEtaNegRatio->GetYaxis()->SetTitle("Ratio (MC/DATA)");
    HistAntilambdaDauEtaNegRatio->GetYaxis()->SetTitleSize(0.07);
    HistAntilambdaDauEtaNegRatio->GetYaxis()->CenterTitle();
    HistAntilambdaDauEtaNegRatio->GetXaxis()->SetTitleSize(0.07);
    HistAntilambdaDauEtaNegRatio->GetYaxis()->SetLabelSize(0.05);
    HistAntilambdaDauEtaNegRatio->GetXaxis()->SetLabelSize(0.06);
    HistAntilambdaDauEtaNegRatio->SetMarkerStyle(4);
    HistAntilambdaDauEtaNegRatio->SetMarkerSize(0.45);
    HistAntilambdaDauEtaNegRatio->SetLineColor(kBlack);
    HistAntilambdaDauEtaNegRatio->SetMarkerColor(kBlack);
#endif

#ifdef DrawFitting
    /// @fitting
    auto Filicanvas10 = new Filipad2(10, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas10->Draw();
    TPad *p10_1 = Filicanvas10->GetPad(1);
    p10_1->cd(1);
    p10_1->SetGrid();
    auto FitDATAK0ShortMass = new TF1("FitDATAK0ShortMass", FitFuncGausPol1, 0.45, 0.55, 5);
    FitDATAK0ShortMass->SetParameters(0.49, 0.01, 1, 1, 1);
    FitDATAK0ShortMass->SetParNames("Gaus Mean", "Gaus Std.");
    FitDATAK0ShortMass->SetParLimits(0, 0.485, 0.50);
    FitDATAK0ShortMass->SetParLimits(1, 0, 0.01);
    FitDATAK0ShortMass->SetParLimits(2, 0, 150);

    auto HistFitDATAK0ShortMass = (TH1D *)HistDATAK0Shortmass->Clone();
    auto HistFitDATAK0ShortMassSubBg = (TH1D *)HistDATAK0Shortmass->Clone();

    gPad->SetGrid();
    HistFitDATAK0ShortMass->GetYaxis()->SetTitle("K^{0}_{s} count[#]");
    HistFitDATAK0ShortMass->GetXaxis()->SetRangeUser(0.45, 0.55);
    HistFitDATAK0ShortMass->Fit("FitDATAK0ShortMass", "R");
    HistFitDATAK0ShortMass->Draw("");
    FitDATAK0ShortMass->Draw("same");

    auto FitDATAK0ShortMassPeak = new TF1("FitDATAK0ShortMassGaus", FitGaussian, 0.45, 0.55, 3);

    FitDATAK0ShortMass->GetParameters(Par);
    FitDATAK0ShortMassPeak->SetParameters(Par);

    auto binwidthHistFitDATAK0ShortMass = HistFitDATAK0ShortMass->GetXaxis()->GetBinWidth(1);
    auto IntegralFitDATAK0ShortMassPeak = FitDATAK0ShortMassPeak->Integral(0.46, 0.54);
    IntegralFitDATAK0ShortMassPeak = IntegralFitDATAK0ShortMassPeak / binwidthHistFitDATAK0ShortMass;

    auto FitDATAK0ShortMassBg = new TF1("FitDATAK0ShortMassGaus", FitPol1, 0.45, 0.55, 2);
    FitDATAK0ShortMassBg->SetParameters(&Par[3]);
    FitDATAK0ShortMassBg->Draw("same");

    FitDATAK0ShortMass->SetLineColor(kGreen);
    FitDATAK0ShortMass->SetLineWidth(4);
    FitDATAK0ShortMassPeak->SetLineColor(kGreen);
    FitDATAK0ShortMassPeak->SetLineWidth(4);
    FitDATAK0ShortMassBg->SetLineColor(kBlack);
    FitDATAK0ShortMassBg->SetLineWidth(4);

    HistFitDATAK0ShortMass->GetYaxis()->SetTitleSize(0.05);
    HistFitDATAK0ShortMass->GetYaxis()->CenterTitle();

    TLegend *l10 = new TLegend(0.643813, 0.675362, 0.964883, 0.869565, NULL, "brNDC");
    l10->AddEntry(HistDATAK0Shortmass, "K^{0}_{s} DATA : LHC22s pass4", "lp");
    l10->AddEntry(FitDATAK0ShortMass, "Signal with background", "lp");
    l10->AddEntry(FitDATAK0ShortMassBg, "Background", "lp");
    l10->AddEntry(FitDATAK0ShortMassPeak, "Signal sub. background", "lp");
    l10->SetBorderSize(0);
    l10->SetTextAlign(12);
    l10->SetTextSize(0.04);
    l10->Draw("same");

    TPad *p10_2 = Filicanvas10->GetPad(2);
    p10_2->cd(2);
    p10_2->SetGrid();
    HistFitDATAK0ShortMassSubBg->Add(FitDATAK0ShortMassBg, -1);
    HistFitDATAK0ShortMassSubBg->GetXaxis()->SetRangeUser(0.45, 0.55);
    HistFitDATAK0ShortMassSubBg->GetYaxis()->SetRangeUser(-10, 100);
    HistFitDATAK0ShortMassSubBg->Draw();
    FitDATAK0ShortMassPeak->SetLineColor(kViolet);
    FitDATAK0ShortMassPeak->SetFillColorAlpha(kViolet, 0.25);
    FitDATAK0ShortMassPeak->SetFillStyle(3144);
    FitDATAK0ShortMassPeak->Draw("same");

    HistFitDATAK0ShortMassSubBg->GetYaxis()->SetTitle("Subtracintion Bg.");
    HistFitDATAK0ShortMassSubBg->GetYaxis()->SetTitleSize(0.07);
    HistFitDATAK0ShortMassSubBg->GetYaxis()->CenterTitle();

    HistFitDATAK0ShortMassSubBg->GetYaxis()->SetLabelSize(0.05);
    HistFitDATAK0ShortMassSubBg->GetXaxis()->SetTitleSize(0.08);
    HistFitDATAK0ShortMassSubBg->GetXaxis()->SetLabelSize(0.09);

    TLegend *l10b = new TLegend(0.605351, 0.665217, 0.938127, 0.852174, NULL, "brNDC");
    l10b->AddEntry(FitDATAK0ShortMassPeak, Form("# of particles : %.1f", IntegralFitDATAK0ShortMassPeak), "");
    l10b->SetBorderSize(0);
    l10b->SetTextAlign(12);
    l10b->SetTextSize(0.07);
    l10b->SetMargin(-0.04);
    l10b->Draw("same");

    auto Filicanvas11 = new Filipad2(11, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas11->Draw();
    TPad *p11_1 = Filicanvas11->GetPad(1);
    p11_1->cd(1);
    p11_1->SetGrid();
    auto FitMCK0ShortMass = new TF1("FitMCK0ShortMass", FitFuncGausPol1, 0.45, 0.55, 5);
    FitMCK0ShortMass->SetParameters(0.49, 0.01, 1, 1, 1);
    FitMCK0ShortMass->SetParNames("Gaus Mean", "Gaus Std.");
    FitMCK0ShortMass->SetParLimits(0, 0.485, 0.50);
    FitMCK0ShortMass->SetParLimits(1, 0, 0.01);
    FitMCK0ShortMass->SetParLimits(2, 0, 1e8);

    auto HistFitMCK0ShortMass = (TH1D *)HistMCK0Shortmass->Clone();
    auto HistFitMCK0ShortMassSubBg = (TH1D *)HistMCK0Shortmass->Clone();

    gPad->SetGrid();
    HistFitMCK0ShortMass->GetYaxis()->SetTitle("K^{0}_{s} count[#]");
    HistFitMCK0ShortMass->GetXaxis()->SetRangeUser(0.45, 0.55);
    HistFitMCK0ShortMass->Fit("FitMCK0ShortMass", "R");
    HistFitMCK0ShortMass->Draw("");
    FitMCK0ShortMass->Draw("same");

    auto FitMCK0ShortMassPeak = new TF1("FitMCK0ShortMassGaus", FitGaussian, 0.45, 0.55, 3);

    FitMCK0ShortMass->GetParameters(Par);
    FitMCK0ShortMassPeak->SetParameters(Par);

    auto binwidthHistFitMCK0ShortMass = HistFitMCK0ShortMass->GetXaxis()->GetBinWidth(1);
    auto IntegralFitMCK0ShortMassPeak = FitMCK0ShortMassPeak->Integral(0.46, 0.54);
    IntegralFitMCK0ShortMassPeak = IntegralFitMCK0ShortMassPeak / binwidthHistFitMCK0ShortMass;

    auto FitMCK0ShortMassBg = new TF1("FitMCK0ShortMassGaus", FitPol1, 0.45, 0.55, 2);
    FitMCK0ShortMassBg->SetParameters(&Par[3]);
    FitMCK0ShortMassBg->Draw("same");

    FitMCK0ShortMass->SetLineColor(kGreen);
    FitMCK0ShortMass->SetLineWidth(4);
    FitMCK0ShortMassPeak->SetLineColor(kGreen);
    FitMCK0ShortMassPeak->SetLineWidth(4);
    FitMCK0ShortMassBg->SetLineColor(kBlack);
    FitMCK0ShortMassBg->SetLineWidth(4);

    HistFitMCK0ShortMass->GetYaxis()->SetTitleSize(0.05);
    HistFitMCK0ShortMass->GetYaxis()->CenterTitle();

    TLegend *l11 = new TLegend(0.643813, 0.675362, 0.964883, 0.869565, NULL, "brNDC");
    l11->AddEntry(HistMCK0Shortmass, "K^{0}_{s} MC : LHC22k3b2", "lp");
    l11->AddEntry(FitMCK0ShortMass, "Signal with background", "lp");
    l11->AddEntry(FitMCK0ShortMassBg, "Background", "lp");
    l11->AddEntry(FitMCK0ShortMassPeak, "Signal sub. background", "lp");
    l11->SetBorderSize(0);
    l11->SetTextAlign(12);
    l11->SetTextSize(0.04);
    l11->Draw("same");

    TPad *p11_2 = Filicanvas11->GetPad(2);
    p11_2->cd(2);
    p11_2->SetGrid();
    HistFitMCK0ShortMassSubBg->Add(FitMCK0ShortMassBg, -1);
    HistFitMCK0ShortMassSubBg->GetXaxis()->SetRangeUser(0.45, 0.55);
    HistFitMCK0ShortMassSubBg->GetYaxis()->SetRangeUser(-500, 1.2e4);
    HistFitMCK0ShortMassSubBg->Draw("l");
    FitMCK0ShortMassPeak->SetLineColor(kViolet);
    FitMCK0ShortMassPeak->SetFillColorAlpha(kViolet, 0.25);
    FitMCK0ShortMassPeak->SetFillStyle(3144);
    FitMCK0ShortMassPeak->Draw("same");

    HistFitMCK0ShortMassSubBg->GetYaxis()->SetTitle("Subtracintion Bg.");
    HistFitMCK0ShortMassSubBg->GetYaxis()->SetTitleSize(0.07);
    HistFitMCK0ShortMassSubBg->GetYaxis()->CenterTitle();

    HistFitMCK0ShortMassSubBg->GetYaxis()->SetLabelSize(0.05);
    HistFitMCK0ShortMassSubBg->GetXaxis()->SetTitleSize(0.08);
    HistFitMCK0ShortMassSubBg->GetXaxis()->SetLabelSize(0.08);

    TLegend *l11b = new TLegend(0.605351, 0.665217, 0.938127, 0.852174, NULL, "brNDC");
    l11b->AddEntry(FitMCK0ShortMassPeak, Form("# of particles : %.1f", IntegralFitMCK0ShortMassPeak), "");
    l11b->SetBorderSize(0);
    l11b->SetTextAlign(12);
    l11b->SetTextSize(0.07);
    l11b->SetMargin(-0.04);
    l11b->Draw("same");

    auto Filicanvas12 = new Filipad2(12, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas12->Draw();
    TPad *p12_1 = Filicanvas12->GetPad(1);
    p12_1->cd(1);
    p12_1->SetGrid();
    auto FitDATALambdaMass = new TF1("FitDATALambdaMass", FitDoubleGaus, 1.08, 1.15, 6);
    FitDATALambdaMass->SetParameters(0.49, 0.01, 1, 1, 1, 1);
    FitDATALambdaMass->SetParNames("Gaus Mean 1", "Gaus Std. 1", "Amp. 1", "Gaus Mean 2", "Gaus Std. 2", "Amp. 2");
    FitDATALambdaMass->SetParLimits(0, 1.11, 1.115);
    FitDATALambdaMass->SetParLimits(1, 0, 0.002);
    FitDATALambdaMass->SetParLimits(2, 0, 1e4);
    FitDATALambdaMass->SetParLimits(3, 1.0, 1.15);
    FitDATALambdaMass->SetParLimits(4, 0, 0.03);
    FitDATALambdaMass->SetParLimits(5, 0, 1e4);

    auto HistFitDATALambdaMass = (TH1D *)HistDATALambdamass->Clone();
    auto HistFitDATALambdaMassSubBg = (TH1D *)HistDATALambdamass->Clone();

    gPad->SetGrid();
    HistFitDATALambdaMass->GetYaxis()->SetTitle("#Lambda count[#]");
    HistFitDATALambdaMass->GetXaxis()->SetRangeUser(1.10, 1.125);
    HistFitDATALambdaMass->GetYaxis()->SetRangeUser(7.5e2, 1e3);
    HistFitDATALambdaMass->Fit("FitDATALambdaMass", "R");
    HistFitDATALambdaMass->Draw("");
    FitDATALambdaMass->Draw("same");

    auto FitDATALambdaMassPeak = new TF1("FitDATALambdaMassGaus", FitGaussian, 1.08, 1.15, 3);

    FitDATALambdaMass->GetParameters(Par);
    FitDATALambdaMassPeak->SetParameters(Par);

    auto binwidthHistFitDATALambdaMass = HistFitDATALambdaMass->GetXaxis()->GetBinWidth(1);
    auto IntegralFitDATALambdaMassPeak = FitDATALambdaMassPeak->Integral(1.08, 1.15);
    IntegralFitDATALambdaMassPeak = IntegralFitDATALambdaMassPeak / binwidthHistFitDATALambdaMass;

    auto FitDATALambdaMassBg = new TF1("FitDATALambdaMassGaus", FitGaussian, 1.08, 1.15, 3);
    FitDATALambdaMassBg->SetParameters(&Par[3]);
    FitDATALambdaMassBg->Draw("same");

    FitDATALambdaMass->SetLineColor(kGreen);
    FitDATALambdaMass->SetLineWidth(4);
    FitDATALambdaMassPeak->SetLineColor(kGreen);
    FitDATALambdaMassPeak->SetLineWidth(4);
    FitDATALambdaMassBg->SetLineColor(kBlack);
    FitDATALambdaMassBg->SetLineWidth(4);

    HistFitDATALambdaMass->GetYaxis()->SetTitleSize(0.05);
    HistFitDATALambdaMass->GetYaxis()->CenterTitle();

    TLegend *l12 = new TLegend(0.643813, 0.675362, 0.964883, 0.869565, NULL, "brNDC");
    l12->AddEntry(HistDATALambdamass, "#Lambda DATA : LHC22s pass4", "lp");
    l12->AddEntry(FitDATALambdaMass, "Signal with background", "lp");
    l12->AddEntry(FitDATALambdaMassBg, "Background", "lp");
    l12->AddEntry(FitDATALambdaMassPeak, "Signal sub. background", "lp");
    l12->SetBorderSize(0);
    l12->SetTextAlign(12);
    l12->SetTextSize(0.04);
    l12->Draw("same");

    TPad *p12_2 = Filicanvas12->GetPad(2);
    p12_2->cd(2);
    p12_2->SetGrid();
    HistFitDATALambdaMassSubBg->Add(FitDATALambdaMassBg, -1);
    HistFitDATALambdaMassSubBg->GetXaxis()->SetRangeUser(1.10, 1.125);
    HistFitDATALambdaMassSubBg->GetYaxis()->SetRangeUser(-10, 100);
    HistFitDATALambdaMassSubBg->Draw();
    FitDATALambdaMassPeak->SetLineColor(kViolet);
    FitDATALambdaMassPeak->SetFillColorAlpha(kViolet, 0.25);
    FitDATALambdaMassPeak->SetFillStyle(3144);
    FitDATALambdaMassPeak->Draw("same");

    HistFitDATALambdaMassSubBg->GetYaxis()->SetTitle("Subtracintion Bg.");
    HistFitDATALambdaMassSubBg->GetYaxis()->SetTitleSize(0.07);
    HistFitDATALambdaMassSubBg->GetYaxis()->CenterTitle();

    HistFitDATALambdaMassSubBg->GetYaxis()->SetLabelSize(0.05);
    HistFitDATALambdaMassSubBg->GetXaxis()->SetTitleSize(0.08);
    HistFitDATALambdaMassSubBg->GetXaxis()->SetLabelSize(0.08);

    TLegend *l12b = new TLegend(0.695652,0.665217,0.936455,0.852174,NULL,"brNDC");
    l12b->AddEntry(FitDATALambdaMassPeak, Form("# of particles : %.1f", IntegralFitDATALambdaMassPeak), "");
    l12b->SetBorderSize(0);
    l12b->SetTextAlign(12);
    l12b->SetTextSize(0.07);
    l12b->SetMargin(-0.04);
    l12b->Draw("same");

    auto Filicanvas13 = new Filipad2(13, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas13->Draw();
    TPad *p13_1 = Filicanvas13->GetPad(1);
    p13_1->cd(1);
    p13_1->SetGrid();
    auto FitMCLambdaMass = new TF1("FitMCLambdaMass", FitDoubleGaus, 1.08, 1.15, 6);
    FitMCLambdaMass->SetParameters(0.49, 0.01, 1, 1, 1, 1);
    FitMCLambdaMass->SetParNames("Gaus Mean 1", "Gaus Std. 1", "Amp. 1", "Gaus Mean 2", "Gaus Std. 2", "Amp. 2");
    FitMCLambdaMass->SetParLimits(0, 1.11, 1.118);
    FitMCLambdaMass->SetParLimits(1, 0, 0.002);
    FitMCLambdaMass->SetParLimits(2, 0, 1e4);
    FitMCLambdaMass->SetParLimits(3, 1.0, 1.15);
    FitMCLambdaMass->SetParLimits(4, 0, 0.03);
    FitMCLambdaMass->SetParLimits(5, 0, 1e4);

    auto HistFitMCLambdaMass = (TH1D *)HistMCLambdamass->Clone();
    auto HistFitMCLambdaMassSubBg = (TH1D *)HistMCLambdamass->Clone();

    gPad->SetGrid();
    HistFitMCLambdaMass->GetYaxis()->SetTitle("#Lambda count[#]");
    HistFitMCLambdaMass->GetXaxis()->SetRangeUser(1.10, 1.125);
    HistFitMCLambdaMass->GetYaxis()->SetRangeUser(5e3, 1.1e4);
    HistFitMCLambdaMass->Fit("FitMCLambdaMass", "R");
    HistFitMCLambdaMass->Draw("");
    FitMCLambdaMass->Draw("same");

    auto FitMCLambdaMassPeak = new TF1("FitMCLambdaMassGaus", FitGaussian, 1.08, 1.15, 3);

    FitMCLambdaMass->GetParameters(Par);
    FitMCLambdaMassPeak->SetParameters(Par);

    auto binwidthHistFitMCLambdaMass = HistFitMCLambdaMass->GetXaxis()->GetBinWidth(1);
    auto IntegralFitMCLambdaMassPeak = FitMCLambdaMassPeak->Integral(1.08, 1.15);
    IntegralFitMCLambdaMassPeak = IntegralFitMCLambdaMassPeak / binwidthHistFitMCLambdaMass;

    auto FitMCLambdaMassBg = new TF1("FitMCLambdaMassGaus", FitGaussian, 1.08, 1.15, 3);
    FitMCLambdaMassBg->SetParameters(&Par[3]);
    FitMCLambdaMassBg->Draw("same");

    FitMCLambdaMass->SetLineColor(kGreen);
    FitMCLambdaMass->SetLineWidth(4);
    FitMCLambdaMassPeak->SetLineColor(kGreen);
    FitMCLambdaMassPeak->SetLineWidth(4);
    FitMCLambdaMassBg->SetLineColor(kBlack);
    FitMCLambdaMassBg->SetLineWidth(4);

    HistFitMCLambdaMass->GetYaxis()->SetTitleSize(0.05);
    HistFitMCLambdaMass->GetYaxis()->CenterTitle();

    TLegend *l13 = new TLegend(0.643813, 0.675362, 0.964883, 0.869565, NULL, "brNDC");
    l13->AddEntry(HistMCLambdamass, "#Lambda MC : LHC22k3b2", "lp");
    l13->AddEntry(FitMCLambdaMass, "Signal with background", "lp");
    l13->AddEntry(FitMCLambdaMassBg, "Background", "lp");
    l13->AddEntry(FitMCLambdaMassPeak, "Signal sub. background", "lp");
    l13->SetBorderSize(0);
    l13->SetTextAlign(12);
    l13->SetTextSize(0.04);
    l13->Draw("same");

    TPad *p13_2 = Filicanvas13->GetPad(2);
    p13_2->cd(2);
    p13_2->SetGrid();
    HistFitMCLambdaMassSubBg->Add(FitMCLambdaMassBg, -1);
    HistFitMCLambdaMassSubBg->GetXaxis()->SetRangeUser(1.10, 1.125);
    HistFitMCLambdaMassSubBg->GetYaxis()->SetRangeUser(-200, 4e3);
    HistFitMCLambdaMassSubBg->Draw();
    FitMCLambdaMassPeak->SetLineColor(kViolet);
    FitMCLambdaMassPeak->SetFillColorAlpha(kViolet, 0.25);
    FitMCLambdaMassPeak->SetFillStyle(3144);
    FitMCLambdaMassPeak->Draw("same");

    HistFitMCLambdaMassSubBg->GetYaxis()->SetTitle("Subtracintion Bg.");
    HistFitMCLambdaMassSubBg->GetYaxis()->SetTitleSize(0.07);
    HistFitMCLambdaMassSubBg->GetYaxis()->CenterTitle();

    HistFitMCLambdaMassSubBg->GetYaxis()->SetLabelSize(0.05);
    HistFitMCLambdaMassSubBg->GetXaxis()->SetTitleSize(0.08);
    HistFitMCLambdaMassSubBg->GetXaxis()->SetLabelSize(0.08);

    TLegend *l13b = new TLegend(0.695652,0.665217,0.936455,0.852174,NULL,"brNDC");
    l13b->AddEntry(FitMCLambdaMassPeak, Form("# of particles : %.1f", IntegralFitMCLambdaMassPeak), "");
    l13b->SetBorderSize(0);
    l13b->SetTextAlign(12);
    l13b->SetTextSize(0.07);
    l13b->SetMargin(-0.04);
    l13b->Draw("same");

    auto Filicanvas14 = new Filipad2(14, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas14->Draw();
    TPad *p14_1 = Filicanvas14->GetPad(1);
    p14_1->cd(1);
    p14_1->SetGrid();
    auto FitDATAAntilambdaMass = new TF1("FitDATAAntilambdaMass", FitDoubleGaus, 1.08, 1.15, 6);
    FitDATAAntilambdaMass->SetParameters(0.49, 0.01, 1, 1, 1, 1);
    FitDATAAntilambdaMass->SetParNames("Gaus Mean 1", "Gaus Std. 1", "Amp. 1", "Gaus Mean 2", "Gaus Std. 2", "Amp. 2");
    FitDATAAntilambdaMass->SetParLimits(0, 1.11, 1.118);
    FitDATAAntilambdaMass->SetParLimits(1, 0, 0.005);
    FitDATAAntilambdaMass->SetParLimits(2, 0, 1e3);
    FitDATAAntilambdaMass->SetParLimits(3, 1.0, 1.15);
    FitDATAAntilambdaMass->SetParLimits(4, 0, 0.03);
    FitDATAAntilambdaMass->SetParLimits(5, 0, 1e3);

    auto HistFitDATAAntilambdaMass = (TH1D *)HistDATAAntilambdamass->Clone();
    auto HistFitDATAAntilambdaMassSubBg = (TH1D *)HistDATAAntilambdamass->Clone();

    gPad->SetGrid();
    HistFitDATAAntilambdaMass->GetYaxis()->SetTitle("#bar{#Lambda} count[#]");
    HistFitDATAAntilambdaMass->GetXaxis()->SetRangeUser(1.10, 1.125);
    HistFitDATAAntilambdaMass->GetYaxis()->SetRangeUser(7.5e2, 1e3);
    HistFitDATAAntilambdaMass->Fit("FitDATAAntilambdaMass", "R");
    HistFitDATAAntilambdaMass->Draw("");
    FitDATAAntilambdaMass->Draw("same");

    auto FitDATAAntilambdaMassPeak = new TF1("FitDATAAntilambdaMassGaus", FitGaussian, 1.08, 1.15, 3);

    FitDATAAntilambdaMass->GetParameters(Par);
    FitDATAAntilambdaMassPeak->SetParameters(Par);

    auto binwidthHistFitDATAAntilambdaMass = HistFitDATAAntilambdaMass->GetXaxis()->GetBinWidth(1);
    auto IntegralFitDATAAntilambdaMassPeak = FitDATAAntilambdaMassPeak->Integral(1.08, 1.15);
    IntegralFitDATAAntilambdaMassPeak = IntegralFitDATAAntilambdaMassPeak / binwidthHistFitDATAAntilambdaMass;
    // auto IntegralHistFitMCAntilambdaMassSubBg = HistFitMCAntilambdaMassSubBg->Integral(HistFitMCAntilambdaMassSubBg->FindBin(1.108),HistFitMCAntilambdaMassSubBg->FindBin(1.122));

    auto FitDATAAntilambdaMassBg = new TF1("FitDATAAntilambdaMassGaus", FitGaussian, 1.08, 1.15, 3);
    FitDATAAntilambdaMassBg->SetParameters(&Par[3]);
    FitDATAAntilambdaMassBg->Draw("same");

    FitDATAAntilambdaMass->SetLineColor(kGreen);
    FitDATAAntilambdaMass->SetLineWidth(4);
    FitDATAAntilambdaMassPeak->SetLineColor(kGreen);
    FitDATAAntilambdaMassPeak->SetLineWidth(4);
    FitDATAAntilambdaMassBg->SetLineColor(kBlack);
    FitDATAAntilambdaMassBg->SetLineWidth(4);

    HistFitDATAAntilambdaMass->GetYaxis()->SetTitleSize(0.05);
    HistFitDATAAntilambdaMass->GetYaxis()->CenterTitle();

    TLegend *l14 = new TLegend(0.643813, 0.675362, 0.964883, 0.869565, NULL, "brNDC");
    l14->AddEntry(HistDATAAntilambdamass, "#bar{#Lambda} DATA : LHC22s pass4", "lp");
    l14->AddEntry(FitDATAAntilambdaMass, "Signal with background", "lp");
    l14->AddEntry(FitDATAAntilambdaMassBg, "Background", "lp");
    l14->AddEntry(FitDATAAntilambdaMassPeak, "Signal sub. background", "lp");
    l14->SetBorderSize(0);
    l14->SetTextAlign(12);
    l14->SetTextSize(0.04);
    l14->Draw("same");

    TPad *p14_2 = Filicanvas14->GetPad(2);
    p14_2->cd(2);
    p14_2->SetGrid();
    HistFitDATAAntilambdaMassSubBg->Add(FitDATAAntilambdaMassBg, -1);
    HistFitDATAAntilambdaMassSubBg->GetXaxis()->SetRangeUser(1.10, 1.125);
    HistFitDATAAntilambdaMassSubBg->GetYaxis()->SetRangeUser(-10, 100);
    HistFitDATAAntilambdaMassSubBg->Draw();
    FitDATAAntilambdaMassPeak->SetLineColor(kViolet);
    FitDATAAntilambdaMassPeak->SetFillColorAlpha(kViolet, 0.25);
    FitDATAAntilambdaMassPeak->SetFillStyle(3144);
    FitDATAAntilambdaMassPeak->Draw("same");

    HistFitDATAAntilambdaMassSubBg->GetYaxis()->SetTitle("Subtracintion Bg.");
    HistFitDATAAntilambdaMassSubBg->GetYaxis()->SetTitleSize(0.07);
    HistFitDATAAntilambdaMassSubBg->GetYaxis()->CenterTitle();

    HistFitDATAAntilambdaMassSubBg->GetYaxis()->SetLabelSize(0.05);
    HistFitDATAAntilambdaMassSubBg->GetXaxis()->SetTitleSize(0.08);
    HistFitDATAAntilambdaMassSubBg->GetXaxis()->SetLabelSize(0.08);

    TLegend *l14b = new TLegend(0.695652,0.665217,0.936455,0.852174,NULL,"brNDC");
    l14b->AddEntry(FitDATAAntilambdaMassPeak, Form("# of particles : %.1f", IntegralFitDATAAntilambdaMassPeak), "");
    l14b->SetBorderSize(0);
    l14b->SetTextAlign(12);
    l14b->SetTextSize(0.07);
    l14b->SetMargin(-0.04);
    l14b->Draw("same");

    auto Filicanvas15 = new Filipad2(15, 2, 0.4, 200, 500, 1, 1, 1);
    Filicanvas15->Draw();
    TPad *p15_1 = Filicanvas15->GetPad(1);
    p15_1->cd(1);
    p15_1->SetGrid();
    auto FitMCAntilambdaMass = new TF1("FitMCAntilambdaMass", FitDoubleGaus, 1.08, 1.15, 6);
    FitMCAntilambdaMass->SetParameters(0.49, 0.01, 1, 1, 1, 1);
    FitMCAntilambdaMass->SetParNames("Gaus Mean 1", "Gaus Std. 1", "Amp. 1", "Gaus Mean 2", "Gaus Std. 2", "Amp. 2");
    FitMCAntilambdaMass->SetParLimits(0, 1.11, 1.115);
    FitMCAntilambdaMass->SetParLimits(1, 0, 0.002);
    FitMCAntilambdaMass->SetParLimits(2, 0, 1e4);
    FitMCAntilambdaMass->SetParLimits(3, 1.0, 1.15);
    FitMCAntilambdaMass->SetParLimits(4, 0, 0.03);
    FitMCAntilambdaMass->SetParLimits(5, 0, 1e4);

    auto HistFitMCAntilambdaMass = (TH1D *)HistMCAntilambdamass->Clone();
    auto HistFitMCAntilambdaMassSubBg = (TH1D *)HistMCAntilambdamass->Clone();

    gPad->SetGrid();
    HistFitMCAntilambdaMass->GetYaxis()->SetTitle("#bar{#Lambda} count[#]");
    HistFitMCAntilambdaMass->GetXaxis()->SetRangeUser(1.10, 1.125);
    HistFitMCAntilambdaMass->GetYaxis()->SetRangeUser(5e3, 1.1e4);
    HistFitMCAntilambdaMass->Fit("FitMCAntilambdaMass", "R");
    HistFitMCAntilambdaMass->Draw("");
    FitMCAntilambdaMass->Draw("same");

    auto FitMCAntilambdaMassPeak = new TF1("FitMCAntilambdaMassGaus", FitGaussian, 1.08, 1.15, 3);

    FitMCAntilambdaMass->GetParameters(Par);
    FitMCAntilambdaMassPeak->SetParameters(Par);

    auto binwidthHistFitMCAntilambdaMass = HistFitMCAntilambdaMass->GetXaxis()->GetBinWidth(1);
    double IntegralFitMCAntilambdaMassPeak = FitMCAntilambdaMassPeak->Integral(1.08, 1.15);
    IntegralFitMCAntilambdaMassPeak = IntegralFitMCAntilambdaMassPeak / binwidthHistFitMCAntilambdaMass;
    // auto binwidthHistFitMCAntilambdaMass = HistFitMCAntilambdaMass->GetXaxis()->GetBinWidth(1);
    // auto IntegralHistFitMCAntilambdaMassSubBg = HistFitMCAntilambdaMassSubBg->Integral(HistFitMCAntilambdaMassSubBg->FindBin(1.108),HistFitMCAntilambdaMassSubBg->FindBin(1.122));

    auto FitMCAntilambdaMassBg = new TF1("FitMCAntilambdaMassGaus", FitGaussian, 1.08, 1.15, 3);
    FitMCAntilambdaMassBg->SetParameters(&Par[3]);
    FitMCAntilambdaMassBg->Draw("same");

    FitMCAntilambdaMass->SetLineColor(kGreen);
    FitMCAntilambdaMass->SetLineWidth(4);
    FitMCAntilambdaMassPeak->SetLineColor(kGreen);
    FitMCAntilambdaMassPeak->SetLineWidth(4);
    FitMCAntilambdaMassBg->SetLineColor(kBlack);
    FitMCAntilambdaMassBg->SetLineWidth(4);

    HistFitMCAntilambdaMass->GetYaxis()->SetTitleSize(0.05);
    HistFitMCAntilambdaMass->GetYaxis()->CenterTitle();

    TLegend *l15 = new TLegend(0.643813, 0.675362, 0.964883, 0.869565, NULL, "brNDC");
    l15->AddEntry(HistMCAntilambdamass, "#bar{#Lambda} MC : LHC22k3b2", "lp");
    l15->AddEntry(FitMCAntilambdaMass, "Signal with background", "lp");
    l15->AddEntry(FitMCAntilambdaMassBg, "Background", "lp");
    l15->AddEntry(FitMCAntilambdaMassPeak, "Signal sub. background", "lp");
    l15->SetBorderSize(0);
    l15->SetTextAlign(12);
    l15->SetTextSize(0.04);
    l15->Draw("same");

    TPad *p15_2 = Filicanvas15->GetPad(2);
    p15_2->cd(2);
    p15_2->SetGrid();
    HistFitMCAntilambdaMassSubBg->Add(FitMCAntilambdaMassBg, -1);
    HistFitMCAntilambdaMassSubBg->GetXaxis()->SetRangeUser(1.10, 1.125);
    HistFitMCAntilambdaMassSubBg->GetYaxis()->SetRangeUser(-200, 4e3);
    HistFitMCAntilambdaMassSubBg->Draw();
    FitMCAntilambdaMassPeak->SetLineColor(kViolet);
    FitMCAntilambdaMassPeak->SetFillColorAlpha(kViolet, 0.25);
    FitMCAntilambdaMassPeak->SetFillStyle(3144);
    FitMCAntilambdaMassPeak->Draw("same");

    HistFitMCAntilambdaMassSubBg->GetYaxis()->SetTitle("Subtracintion Bg.");
    HistFitMCAntilambdaMassSubBg->GetYaxis()->SetTitleSize(0.07);
    HistFitMCAntilambdaMassSubBg->GetYaxis()->CenterTitle();

    HistFitMCAntilambdaMassSubBg->GetYaxis()->SetLabelSize(0.05);
    HistFitMCAntilambdaMassSubBg->GetXaxis()->SetTitleSize(0.08);
    HistFitMCAntilambdaMassSubBg->GetXaxis()->SetLabelSize(0.08);

    TLegend *l15b = new TLegend(0.695652,0.665217,0.936455,0.852174,NULL,"brNDC");
    l15b->AddEntry(FitMCAntilambdaMassPeak, Form("# of particles : %.1f", IntegralFitMCAntilambdaMassPeak), "");
    l15b->SetBorderSize(0);
    l15b->SetTextAlign(12);
    l15b->SetTextSize(0.07);
    l15b->SetMargin(-0.04);
    l15b->Draw("same");

#endif

// #ifdef DrawMotherV0
// TCanvas c99 = new TCanvas("c99","c99",800,600);
// HistMCMotherV0Count->Draw("][");

// #endif
}
