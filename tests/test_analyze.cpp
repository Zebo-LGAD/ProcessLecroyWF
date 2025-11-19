#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TProfile.h>

#include <string>
#include <iostream>

#include "LinearRct.h"
#include "WFDataConverter.h"

const double gElectronicGain = 20.58; // charge (in fC) = signal / gElectronicGain;

std::map<int, int> gmTCTCh2Pad{{1, 0}, {2, 4}, {3, 3}, {4, 2}}; // Ch1 is used for trigger, so pad index is 0

// std::string gXReconstructionFile = "Linear-Fit-Data/20_38_47_width_86.5_190_Results";
// std::string gScopeDataFile = "../Beta-Sr90/ProcessedROOTFiles/AC-W1-190V-20C.root";
std::string gLinearPRFitFolder = "/mnt/e/Project-LGAD2/ACLGAD-Test/Position-Reconstruct-Algorithm/Linear-Fit-Data/";
// std::string gLinearPRFitFile = gLinearPRFitFolder + "10_12_18_width_86.5_180_Results";
std::string gLinearPRFitFile = gLinearPRFitFolder + "20_38_47_width_86.5_190_Results";
// std::string gScopeDataFile = "../Beta-Sr90/ProcessedROOTFiles/AC-W1-180V-20C.root";

std::string gExtractedFile = "Processed-Data/AC-W1-190V-20C-Extracted.root";

// Draw TOA difference histograms to determine time offsets between channels
std::map<int, double> gChTimeOffsets;
void NormalizeTOA(WFDataProcessor::WFDataTreeReader &reader)
{
    gSystem->mkdir("plots/test_analyze/NormalizeTOA/", true);
    TCanvas cTOA("cTOA", "cTOA", 800, 600);
    std::map<int, TH1D *> mhTOA;
    auto waveInfoMap = reader.GetChannelDataMap();
    for (auto &iter : waveInfoMap)
    {
        int ch = iter.first;
        std::string hname = "hTOA_ch" + std::to_string(ch);
        mhTOA[ch] = new TH1D(hname.c_str(), hname.c_str(), 1000, -10, 10);
    }

    for (int entry = 0; entry < reader.GetTree()->GetEntries(); entry++)
    {
        reader.GetTree()->GetEntry(entry);
        double toaTrigger = waveInfoMap[1]->toa; // Channel 1 is trigger

        for (auto &iter : waveInfoMap)
        {
            int ch = iter.first;
            double toa = iter.second->toa;
            double charge = iter.second->charge;
            if (charge >= 50.0) // Skip small signals
                mhTOA[ch]->Fill(toa - toaTrigger);
        }
    }
    // Determine time offsets
    TF1 fgaus("fgaus", "gaus", -5, 5);
    for (auto &iter : mhTOA)
    {
        int ch = iter.first;
        if (ch == 1)
        {
            gChTimeOffsets[ch] = 0.0;
            continue;
        }

        auto hTOA = iter.second;
        hTOA->Fit(&fgaus, "RQ");
        hTOA->Fit(&fgaus, "Q", "", fgaus.GetParameter(1) - 5 * fgaus.GetParameter(2), fgaus.GetParameter(1) + 5 * fgaus.GetParameter(2));
        double mean = fgaus.GetParameter(1);
        gChTimeOffsets[ch] = -mean;
        std::cout << "Channel " << ch << " TOA mean: " << mean << " ns" << std::endl;
        hTOA->GetXaxis()->SetRangeUser(mean - fgaus.GetParameter(2) * 10, mean + fgaus.GetParameter(2) * 10);
        cTOA.SaveAs(("plots/test_analyze/NormalizeTOA/Ch" + std::to_string(ch) + "_TOA.png").c_str());
    }
}

double processPosition(std::map<int, double> chCharges)
{
    int max_ch = 2;
    // Find maximum charge channel
    for (int ch = 2; ch <= 4; ch++)
        if (chCharges[ch] > chCharges[max_ch])
            max_ch = ch;
    // Calculate x position using chCharges ratio method
    double xMeasure = 0;

    if (max_ch == 2)
        xMeasure = (chCharges[2] * 2 + chCharges[3] * 3) / (chCharges[2] + chCharges[3]);
    if (max_ch == 3)
    {
        int neighbor_ch = (chCharges[2] > chCharges[4]) ? 2 : 4;
        xMeasure = (chCharges[3] * 3 + chCharges[neighbor_ch] * (neighbor_ch)) / (chCharges[3] + chCharges[neighbor_ch]);
    }
    if (max_ch == 4)
        xMeasure = (chCharges[4] * 4 + chCharges[3] * 3) / (chCharges[4] + chCharges[3]);
    return xMeasure;
}

int main()
{
    gSystem->mkdir("plots/test_analyze/", true);

    PRAlgorithms::LinearRct rct;
    // rct.AddCPMap(1, 0);
    // Add map in order of pad number
    rct.GetCaliDataMap().AddPCMap(1, 4);
    rct.GetCaliDataMap().AddPCMap(2, 3);
    rct.GetCaliDataMap().AddPCMap(3, 2);
    rct.GetCaliDataMap().AddPCMap(4, 1);
    rct.GetCaliDataMap().InitializeMap();
    rct.ReadLinearPRFits(gLinearPRFitFile);
    // rct.GetCaliDataMap().Print();

    rct.GetPRDataMap().AddPCMap(2, 4);
    rct.GetPRDataMap().AddPCMap(3, 3);
    rct.GetPRDataMap().AddPCMap(4, 2);
    rct.GetPRDataMap().InitializeMap();
    rct.SetPadNormFactorsFromCali();
    // rct.GetPRDataMap().Print();
    rct.SetUniformPadSize(125.0); // 125 um pad width
    rct.Print();

    auto valid = rct.JudgeCalibrationValid();
    std::cout << "Calibration valid: " << valid << std::endl;
    // return 1;

    WFDataProcessor::WFDataTreeReader reader({1, 2, 3, 4});
    auto rtn = reader.OpenFile(gExtractedFile);
    if (!rtn)
    {
        std::cerr << "Failed to open scope data file: " << gExtractedFile << std::endl;
        return -1;
    }
    NormalizeTOA(reader);

    std::map<int, double> toa;
    auto section2 = rct.GetPadRightSection(2);
    auto section3 = rct.GetPadRightSection(3);
    double sectionStart = (section2._edgeX.first + section2._edgeX.second) / 2.0;
    double sectionEnd = (section3._edgeX.first + section3._edgeX.second) / 2.0;
    // double sectionStart = (section2._edgeX.first + section2._edgeX.first) / 2.0;
    // double sectionEnd = (section3._edgeX.second + section3._edgeX.second) / 2.0;

    auto hitXCor = new TH2D("hitXCor", "hitXCor", 100, 0, 1000, 100, 0, 5);
    // auto hdeltaTx = new TH2D("hdeltaTx", "hdeltaTx", 40, 0, 1000, 100, -1, 1);
    // auto hdeltaTx2 = new TH2D("hdeltaTx2", "hdeltaTx2", 40, 0, 1000, 100, -1, 1);
    auto hdeltaTx = new TH2D("hdeltaTx", "hdeltaTx", 20, sectionStart, sectionEnd, 80, -0.5, 0.5);
    auto hdeltaTx2 = new TH2D("hdeltaTx2", "hdeltaTx2", 20, sectionStart, sectionEnd, 80, -0.5, 0.5);
    // auto hdeltaTx = new TH2D("hdeltaTx", "hdeltaTx", 40, 0, 1000, 80, -0.5, 0.5);
    // auto hdeltaTx2 = new TH2D("hdeltaTx2", "hdeltaTx2", 40, 0, 1000, 80, -0.5, 0.5);

    auto hxMeasure = new TH1D("hxMeasure", "hxMeasure", 100, 0, 5);
    auto hxPos = new TH1D("hxPos", "hxPos", 100, 0, 1000);

    auto mdata = reader.GetChannelDataMap();
    for (int entry = 0; entry < reader.GetTree()->GetEntries(); entry++)
    {

        reader.GetTree()->GetEntry(entry);
        auto waveInfoMap = reader.GetChannelDataMap();
        double doubleX = -9999.0;
        double doubleY = -9999.0;

        // Reconstruct X
        PRAlgorithms::LinearRct::PRinfo info;
        auto rtn = rct.PRfromWaveInfo(waveInfoMap, info);
        if (!rtn)
            continue;
        double reconX = info.reconX;
        int padLeft = info.pads.first;
        int padRight = info.pads.second;
        int chLeft = info.chs.first;
        int chRight = info.chs.second;
        double signalLeft = info.signals.first;
        double signalRight = info.signals.second;
        double normFactorLeft = rct.GetPadNormFactors().at(padLeft);
        double normFactorRight = rct.GetPadNormFactors().at(padRight);
        int maxCh = (signalLeft * normFactorLeft > signalRight * normFactorRight) ? chLeft : chRight;

        // Reconstrcut measured time
        for (auto &iter : waveInfoMap)
        {
            int ch = iter.first;
            toa[ch] = iter.second->toa;
        }

        // Calculate new toa
        double toa2 = 0;
        if (toa[chLeft] < -1e9 || toa[chRight] < -1e9 || toa[1] < -1e9)
            continue;

        double toaLeft = toa[chLeft] - toa[1] + gChTimeOffsets[chLeft];
        double toaRight = toa[chRight] - toa[1] + gChTimeOffsets[chRight];
        double normSignalLeft = signalLeft * normFactorLeft;
        double normSignalRight = signalRight * normFactorRight;
        toa2 = (toaLeft * normSignalLeft * normSignalLeft + toaRight * normSignalRight * normSignalRight) / (normSignalLeft * normSignalLeft + normSignalRight * normSignalRight);
        // toa2 = (toaLeft * leftSignal + toaRight * rightSignal) / (leftSignal + rightSignal);

        if (TMath::Abs(toa[maxCh] - toa[1] + gChTimeOffsets[maxCh]) < 0.4)
            hdeltaTx->Fill(reconX, toa[maxCh] - toa[1] + gChTimeOffsets[maxCh]);
        if (TMath::Abs(toa2) < 0.4)
            hdeltaTx2->Fill(reconX, toa2);

        hxPos->Fill(reconX);

        // Generate another x position using charge ratio method for comparison
        std::map<int, double> chCharges;
        for (auto &iter : waveInfoMap)
        {
            int ch = iter.first;
            double charge = iter.second->charge;
            chCharges[ch] = charge;
        }
        double xMeasure = processPosition(chCharges);
        hxMeasure->Fill(xMeasure);
        hitXCor->Fill(reconX, xMeasure);
    }

    TLatex latex;
    auto c = new TCanvas("cHitXCor", "cHitXCor", 800, 600);
    hitXCor->Draw("colz");
    c->SaveAs("plots/test_analyze/hitXCor.png");

    hxPos->Draw();
    c->SetLogy();
    c->SaveAs("plots/test_analyze/hxPos.png");
    c->SetLogy(0);
    hxMeasure->Draw();
    c->SaveAs("plots/test_analyze/hxMeasure.png");

    hdeltaTx->Draw("colz");
    c->SaveAs("plots/test_analyze/hdeltaTx.png");
    auto hpfx = hdeltaTx->ProfileX();
    hpfx->SetLineColor(kRed);
    hpfx->Draw("same");
    c->SaveAs("plots/test_analyze/hdeltaTxProfileX.png");
    hdeltaTx->FitSlicesY();
    auto hfpy = (TH1D *)gDirectory->Get("hdeltaTx_2");
    hfpy->SetLineColor(kBlue);
    hfpy->Draw("hist");
    // hfpy->GetYaxis()->SetRangeUser(0, 0.1);
    c->SaveAs("plots/test_analyze/hdeltaTxFitSlicesY.png");
    auto fGaus = new TF1("fGaus", "gaus", -1, 1);
    auto hdeltaTxSigmat = new TH1D(*hfpy);
    hdeltaTxSigmat->SetTitle("Time Resolution per X Position Bin;X (um);#sigma_{t} (ps)");
    for (int bin = 1; bin <= hfpy->GetNbinsX(); bin++)
    {
        // get projection for each bin
        auto proj = hdeltaTx->ProjectionY("proj", bin, bin);
        proj->Fit(fGaus, "Q", "", -0.2, 0.2);
        double mean = fGaus->GetParameter(1);
        double sigma = fGaus->GetParameter(2);
        sigma = proj->GetStdDev();
        // std::cout << "Bin " << bin << ": mean = " << mean << ", sigma = " << sigma << std::endl;
        proj->Draw();
        latex.DrawLatexNDC(0.6, 0.8, Form("Mean = %.3f ns", mean));
        latex.DrawLatexNDC(0.6, 0.7, Form("RMS = %.3f ps", sigma * 1000));
        hdeltaTxSigmat->SetBinContent(bin, sigma * 1000); // convert to ps
        c->SaveAs(Form("plots/test_analyze/hdeltaTx_ProjBin%d.png", bin));
    }
    hdeltaTxSigmat->Draw("hist");
    hdeltaTxSigmat->GetYaxis()->SetRangeUser(0, 100);
    c->SaveAs("plots/test_analyze/hdeltaTx_TimeResolution.png");
    auto hdeltaTxprjy = hdeltaTx->ProjectionY();
    hdeltaTxprjy->Draw();
    hdeltaTxprjy->Fit(fGaus, "Q", "", -0.2, 0.2);
    double mean1 = fGaus->GetParameter(1);
    double sigma1 = fGaus->GetParameter(2);
    latex.DrawLatexNDC(0.6, 0.8, Form("Mean = %.3f ns", mean1));
    latex.DrawLatexNDC(0.7, 0.7, Form("RMS = %.3f ps", sigma1 * 1000));
    c->SaveAs("plots/test_analyze/hdeltaTx_TimeResolution_ProjY.png");

    hdeltaTx2->Draw("colz");
    c->SaveAs("plots/test_analyze/hdeltaTx2.png");
    auto hpfx2 = hdeltaTx2->ProfileX();
    hpfx2->SetLineColor(kRed);
    hpfx2->Draw("same");
    c->SaveAs("plots/test_analyze/hdeltaTxProfileX2.png");
    hdeltaTx2->FitSlicesY();
    auto hfpy2 = (TH1D *)gDirectory->Get("hdeltaTx2_2");
    hfpy2->SetLineColor(kBlue);
    hfpy2->Draw("hist");
    // hfpy2->GetYaxis()->SetRangeUser(0, 0.1);
    c->SaveAs("plots/test_analyze/hdeltaTxFitSlicesY2.png");

    auto hdeltaTxSigmat2 = new TH1D(*hfpy2);
    hdeltaTxSigmat2->SetTitle("Time Resolution per X Position Bin;X (um);#sigma_{t} (ps)");
    for (int bin = 1; bin <= hfpy->GetNbinsX(); bin++)
    {
        // get projection for each bin
        auto proj = hdeltaTx2->ProjectionY("proj", bin, bin);
        proj->Fit(fGaus, "Q", "", -0.2, 0.2);
        double mean = fGaus->GetParameter(1);
        double sigma = fGaus->GetParameter(2);
        sigma = proj->GetStdDev();
        // std::cout << "Bin " << bin << ": mean = " << mean << ", sigma = " << sigma << std::endl;
        proj->Draw();
        latex.DrawLatexNDC(0.6, 0.8, Form("Mean = %.3f ns", mean));
        latex.DrawLatexNDC(0.6, 0.7, Form("RMS = %.3f ps", sigma * 1000));
        hdeltaTxSigmat2->SetBinContent(bin, sigma * 1000); // convert to ps
        c->SaveAs(Form("plots/test_analyze/hdeltaTx2_ProjBin%d.png", bin));
    }
    hdeltaTxSigmat2->Draw("hist");
    hdeltaTxSigmat2->GetYaxis()->SetRangeUser(0, 100);
    hdeltaTxSigmat->SetLineColor(kRed);
    hdeltaTxSigmat->Draw("hist same");
    c->SaveAs("plots/test_analyze/hdeltaTx2_TimeResolution.png");

    auto hdeltaTx2prjy = hdeltaTx2->ProjectionY();
    hdeltaTx2prjy->Draw();
    hdeltaTx2prjy->Fit(fGaus, "Q", "", -0.2, 0.2);
    double mean2 = fGaus->GetParameter(1);
    double sigma2 = fGaus->GetParameter(2);
    latex.DrawLatexNDC(0.6, 0.8, Form("Mean = %.3f ns", mean2));
    latex.DrawLatexNDC(0.7, 0.7, Form("RMS = %.3f ps", sigma2 * 1000));
    c->SaveAs("plots/test_analyze/hdeltaTx2_TimeResolution_ProjY.png");

    // rct.ReadCalibrationFile("Linear-Fit-Data/23_26_21_width_86.5_170_Results/Channel2and1.root");
    // rct.ReadCalibrationFile("Linear-Fit-Data/23_26_21_width_86.5_170_Results/Channel3and2.root");
    // rct.ReadCalibrationFile("Linear-Fit-Data/23_26_21_width_86.5_170_Results/Channel4and3.root");

    // rct.ReadCalibrationFiles("Linear-Fit-Data/23_26_21_width_86.5_170_Results");
    // rct.ReadCalibrationFiles("Linear-Fit-Data/20_38_47_width_86.5_190_Results");
    // rct.ReadCalibrationFiles("Linear-Fit-Data/10_12_18_width_86.5_180_Results");
    return 1;
}