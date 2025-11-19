#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TGraph.h>

#include <string>
#include <iostream>

#include "LinearRct.h"
#include "WFDataConverter.h"

const double gElectronicGain = 20.58; // charge (in fC) = signal / gElectronicGain;

// std::string gXReconstructionFile = "Linear-Fit-Data/20_38_47_width_86.5_190_Results";
// std::string gScopeDataFile = "../Beta-Sr90/ProcessedROOTFiles/AC-W1-190V-20C.root";
std::string gLinearPRFitFolder = "/mnt/e/Project-LGAD2/ACLGAD-Test/Position-Reconstruct-Algorithm/Linear-Fit-Data/";
// std::string gLinearPRFitFile = gLinearPRFitFolder + "10_12_18_width_86.5_180_Results";
std::string gLinearPRFitFile = gLinearPRFitFolder + "20_38_47_width_86.5_190_Results";
// std::string gScopeDataFile = "../Beta-Sr90/ProcessedROOTFiles/AC-W1-180V-20C.root";

// Draw TOA difference histograms to determine time offsets between channels
std::string gScopeDataFolder = "/mnt/e/Project-LGAD2/ACLGAD-Test/Laser/20_38_47_width_86.5_190";
std::string gScopeDataFile = "/mnt/e/Project-LGAD2/ACLGAD-Test/Laser/20_38_47_width_86.5_190";
std::string gOutputROOTFile = gScopeDataFolder + "AC-W1-190V-20C.root";

typedef struct _beam_xyz
{
    double x;
    double y;
    double z;
} BeamStep, BeamPosition;

typedef struct _scan_info
{
    double vBias = 0.0;
    double width = 86.5;
    double xStepLength = 2.5; // um
    int xMinStep = 0;
    double yStepLength = 2.5; // um
    int yMinStep = 0;
    double zStepLength = 2.5; // um
    int zMinStep = 0;
    BeamPosition GetBeamPosition(const BeamStep &step)
    {
        BeamPosition pos;
        pos.x = ((int)step.x - xMinStep) * xStepLength;
        pos.y = ((int)step.y - yMinStep) * yStepLength;
        pos.z = ((int)step.z - zMinStep) * zStepLength;
        return pos;
    }
} ScanInfo;

typedef struct _trc_file_info
{
    std::string mid_name;
    int index;
} TRCInfo;

BeamStep ExtractBeamXYZFromString(const std::string &folderName)
{
    BeamStep beamXYZ = {0.0, 0.0, 0.0};
    // folder name format: .../XX_YY_ZZ/, extract XX, YY, ZZ as x, y, z
    std::string folderNameLocal = folderName;
    if (folderName.back() == '/')
        folderNameLocal = folderName.substr(0, folderName.length() - 1);

    size_t pos1 = folderNameLocal.find_last_of('/');
    if (pos1 == std::string::npos)
        pos1 = -1;
    size_t pos2 = folderNameLocal.find_first_of('_', pos1 + 1);
    size_t pos3 = folderNameLocal.find_first_of('_', pos2 + 1);

    if (pos2 == std::string::npos || pos3 == std::string::npos)
    {
        std::cerr << "Error: Invalid folder name format: " << folderNameLocal << std::endl;
        return beamXYZ;
    }
    std::string xStr = folderNameLocal.substr(pos1 + 1, pos2 - pos1 - 1);
    std::string yStr = folderNameLocal.substr(pos2 + 1, pos3 - pos2 - 1);
    std::string zStr = folderNameLocal.substr(pos3 + 1, folderNameLocal.length() - pos3 - 2); // -2 to remove trailing '/'

    beamXYZ.x = std::stod(xStr);
    beamXYZ.y = std::stod(yStr);
    beamXYZ.z = std::stod(zStr);
    // std::cout << "Extracted beam position: (" << beamXYZ.x << ", " << beamXYZ.y << ", " << beamXYZ.z << ")" << std::endl;

    return beamXYZ;
}

ScanInfo ExtractScanInfoPart1(const std::string &fileName)
{
    ScanInfo scanInfo;
    // file name format: .../****_width_ww_vv/, extract ww as width, vv as vBias
    std::string fileNameLocal = fileName;
    if (fileName.back() == '/')
        fileNameLocal = fileName.substr(0, fileName.length() - 1);
    size_t pos1 = fileNameLocal.find_last_of('/');
    fileNameLocal = fileNameLocal.substr(pos1 + 1); // get the last part
    size_t posWidth = fileNameLocal.find("_width_");
    size_t posVBias = fileNameLocal.find("_", posWidth + 7);
    if (posWidth == std::string::npos || posVBias == std::string::npos)
    {
        std::cerr << "Error: Invalid file name format: " << fileNameLocal << std::endl;
        return scanInfo;
    }
    std::string widthStr = fileNameLocal.substr(posWidth + 7, posVBias - posWidth - 7);
    std::string vBiasStr = fileNameLocal.substr(posVBias + 1);
    scanInfo.width = std::stod(widthStr);
    scanInfo.vBias = std::stod(vBiasStr);
    std::cout << "Extracted scan info: width = " << scanInfo.width << ", vBias = " << scanInfo.vBias << std::endl;

    return scanInfo;
}

ScanInfo &ExtractScanInfoPart2(const std::map<std::string, BeamStep> &beamSteps, ScanInfo &scanInfo)
{
    if (beamSteps.empty())
        return scanInfo;

    // Determine min step and step length for x, y, z
    double xMin = beamSteps.begin()->second.x;
    double yMin = beamSteps.begin()->second.y;
    double zMin = beamSteps.begin()->second.z;
    double xMax = beamSteps.begin()->second.x;
    double yMax = beamSteps.begin()->second.y;
    double zMax = beamSteps.begin()->second.z;

    for (const auto &step : beamSteps)
    {
        if (step.second.x < xMin)
            xMin = step.second.x;
        if (step.second.y < yMin)
            yMin = step.second.y;
        if (step.second.z < zMin)
            zMin = step.second.z;
    }

    scanInfo.xMinStep = static_cast<int>(xMin);
    scanInfo.yMinStep = static_cast<int>(yMin);
    scanInfo.zMinStep = static_cast<int>(zMin);

    std::cout << "Determined scan info: xMinStep = " << scanInfo.xMinStep
              << ", yMinStep = " << scanInfo.yMinStep
              << ", zMinStep = " << scanInfo.zMinStep
              << ", xStepLength = " << scanInfo.xStepLength
              << ", yStepLength = " << scanInfo.yStepLength
              << ", zStepLength = " << scanInfo.zStepLength << std::endl;

    return scanInfo;
}

#include <fstream>
bool CheckTRCFileExists(const std::string &prefix, const std::string &suffix, const std::vector<int> &channels)
{
    std::string sFileName;
    std::ifstream infile;
    for (const auto &ch : channels)
    {
        sFileName = prefix + "C" + std::to_string(ch) + suffix;
        infile.open(sFileName);
        if (!infile.is_open())
            return false;
        infile.close();
    }
    return true;
}

std::vector<TRCInfo> ScanFolder(const std::string &folderName, const std::vector<int> &channels)
{
    std::vector<TRCInfo> trcFiles;

    std::string prefix = folderName;
    std::string suffix_pre = ".trc";
    std::string midName = "Trace";

    std::string suffix = Form("%s%05d", midName.c_str(), 0) + suffix_pre;
    if (CheckTRCFileExists(prefix, suffix, channels))
        trcFiles.push_back({_trc_file_info{midName, 0}});

    for (int index = 0; index < 20; index++)
    {
        midName = "--Trace--";
        suffix = Form("%s%05d", midName.c_str(), index) + suffix_pre;
        if (CheckTRCFileExists(prefix, suffix, channels))
            trcFiles.push_back({_trc_file_info{midName, index}});
        else if (index > 5)
            break;
    }
    return trcFiles;
}

void ConvertToROOT()
{
    WFDataProcessor::WFtrc2ROOT converter({1, 2, 3, 4});
    converter.SetForceMatch(true);

    // Search all sub directories in gScopeDataFolder
    auto dir = gSystem->OpenDirectory(gScopeDataFolder.c_str());
    std::cout << "Get DirEntry:" << gSystem->GetDirEntry(dir) << std::endl;
    const char *ptr_entry;

    auto scanInfo = ExtractScanInfoPart1(gScopeDataFolder);
    std::map<std::string, BeamStep> beamSteps;
    std::cout << "Generating beam steps from folder: " << gScopeDataFolder << std::endl;
    while (ptr_entry = gSystem->GetDirEntry(dir))
    {
        std::string entry = ptr_entry;
        if (entry == "." || entry == "..")
            continue;
        std::string fullPath = gScopeDataFolder + "/" + entry;
        if (gSystem->AccessPathName(fullPath.c_str()))
            continue;
        // std::cout << "Processing folder: " << fullPath << std::endl;
        auto step = ExtractBeamXYZFromString(fullPath);
        beamSteps[ptr_entry] = step;
        // WFDataProcessor::WFDataExtractor extractor({1, 2, 3, 4});
        // extractor.ExtractFromTRCFiles(fullPath, 0);
        static int nFolders = 0;
        nFolders++;
        // if (nFolders >= 5)
        // break;
    }

    scanInfo = ExtractScanInfoPart2(beamSteps, scanInfo);

    gSystem->mkdir("PRData/", true);

    int counter = 0, total = beamSteps.size();
    for (const auto &iter : beamSteps)
    {
        counter++;

        std::string sScanFolder = gScopeDataFolder + "/" + iter.first + "/";
        auto vector = ScanFolder(sScanFolder, {1, 2, 3, 4});

        gSystem->mkdir(("plots/PRData/" + iter.first + "/").c_str(), true);
        std::cout << Form("\rProcessing PRData/%s.root ", iter.first.c_str()) << counter << "/" << total << std::flush;
        if (counter == total)
            std::cout << std::endl;

        converter.OpenFile(Form("PRData/%s.root", iter.first.c_str()));

        for (int idx_file = 0; idx_file < vector.size(); idx_file++)
        {
            auto trcInfo = vector[idx_file];

            // converter.TurnOnPlots("plots/PRData/" + iter.first + "/");
            converter.ReadAllAndFill(sScanFolder, trcInfo.index, trcInfo.mid_name);
        }
        converter.CloseFile();
        // return 1;
    }
}

int main_forConvertROOT()
{
    ConvertToROOT();
    return 1;
}

int main_QACor()
{
    gSystem->Load("rootlogon.C");
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

    rct.GetPRDataMap().AddPCMap(1, 4);
    rct.GetPRDataMap().AddPCMap(2, 3);
    rct.GetPRDataMap().AddPCMap(3, 2);
    rct.GetPRDataMap().AddPCMap(4, 1);
    rct.GetPRDataMap().InitializeMap();
    rct.SetPadNormFactorsFromCali();
    // rct.GetPRDataMap().Print();
    rct.Print();
    auto valid = rct.JudgeCalibrationValid();
    std::cout << "Calibration valid: " << valid << std::endl;

    WFDataProcessor::WFROOTReader reader({1, 2, 3, 4});
    auto dir = gSystem->OpenDirectory("PRData/");
    const char *ptr_entry;

    WFDataProcessor::WFDataExtractor extractor({1, 2, 3, 4});
    WFDataProcessor::_extract_config rangeConfig{{0, 20.0}, 20.0, false, ""};
    for (const auto &ch : {1, 2, 3, 4})
        extractor.SetExtractConfig(ch, rangeConfig);

    // Generate output tgraph for charge - amplitude correlation
    gSystem->mkdir("plots/PRData/", true);
    gSystem->mkdir("plots/ChargeAmpCorrelation/", true);
    std::map<int, TH2D *> mhq10;
    std::map<int, TH2D *> mhq50;
    std::map<int, TH2D *> mhq90;
    std::map<int, TH2D *> mhqpm2ns;
    std::map<int, TH2D *> mhq;
    auto filew = new TFile("plots/ChargeAmpCorrelation/ChargeAmpCorrelation.root", "RECREATE");
    filew->cd();
    auto c = new TCanvas("c", "c", 800, 600);
    for (int ch = 1; ch <= 4; ch++)
    {

        mhq10[ch] = new TH2D(Form("ch%d_q10", ch), ";amplitude (mV);Charge 10% (mV*ns)", 300, 0, 300, 400, 0, 300);
        mhq10[ch]->SetName(Form("ch%d_q10", ch));
        mhq10[ch]->SetTitle(";amplitude (mV);Charge 10% (mV*ns)");
        mhq10[ch]->SetMarkerStyle(20);
        mhq10[ch]->SetMarkerSize(1.0);

        mhq50[ch] = new TH2D(Form("ch%d_q10", ch), ";amplitude (mV);Charge 10% (mV*ns)", 300, 0, 300, 400, 0, 300);
        mhq50[ch]->SetName(Form("ch%d_q50", ch));
        mhq50[ch]->SetTitle(";amplitude (mV);Charge 50% (mV*ns)");
        mhq50[ch]->SetMarkerStyle(20);
        mhq50[ch]->SetMarkerSize(1.0);

        mhq90[ch] = new TH2D(Form("ch%d_q10", ch), ";amplitude (mV);Charge 10% (mV*ns)", 300, 0, 300, 400, 0, 300);
        mhq90[ch]->SetName(Form("ch%d_q90", ch));
        mhq90[ch]->SetTitle(";amplitude (mV);Charge 90% (mV*ns)");
        mhq90[ch]->SetMarkerStyle(20);
        mhq90[ch]->SetMarkerSize(1.0);

        mhq[ch] = new TH2D(Form("ch%d_q10", ch), ";amplitude (mV);Charge 10% (mV*ns)", 300, 0, 300, 400, 0, 300);
        mhq[ch]->SetName(Form("ch%d_q", ch));
        mhq[ch]->SetTitle(";amplitude (mV);Charge (mV*ns)");
        mhq[ch]->SetMarkerStyle(20);
        mhq[ch]->SetMarkerSize(1.0);

        mhqpm2ns[ch] = new TH2D(Form("ch%d_q10", ch), ";amplitude (mV);Charge 10% (mV*ns)", 300, 0, 300, 400, 0, 300);
        mhqpm2ns[ch]->SetName(Form("ch%d_qpm2ns", ch));
        mhqpm2ns[ch]->SetTitle(";amplitude (mV);Charge +/-2ns (mV*ns)");
        mhqpm2ns[ch]->SetMarkerStyle(20);
        mhqpm2ns[ch]->SetMarkerSize(1.0);
    }

    int counter = 0;
    while (ptr_entry = gSystem->GetDirEntry(dir))
    {
        std::string entry = ptr_entry;
        if (entry == "." || entry == "..")
            continue;
        if (entry == "ConvertedROOT")
            continue;
        std::string fullPath = "PRData/" + entry;

        // Get base name without .root
        std::string scanName = fullPath.substr(0, fullPath.length() - 5);
        scanName = scanName.substr(7); // remove "PRData/"
        _beam_xyz beamPos = ExtractBeamXYZFromString(scanName);
        double beamX = (beamPos.x + 1300) * 2.5;
        double beamY = (beamPos.y - 1200) * 2.5;

        if (gSystem->AccessPathName(fullPath.c_str()))
            continue;
        std::cout << "***********Processing file: *************" << fullPath << std::endl;
        reader.OpenFile(fullPath);

        gSystem->mkdir("plots/PRData/", true);
        gSystem->mkdir(("plots/PRData/" + scanName + "/").c_str(), true);

        bool fAlreadyOutput = false;

        for (int evt = 0; evt < reader.GetTree()->GetEntries(); evt++)
        {
            reader.GetTree()->GetEntry(evt);

            auto dataMap = extractor.GetChannelDataMap();
            extractor.TurnOffPlots();
            extractor.ExtractFromWFROOTReader(reader);
            for (auto iter : extractor.GetChannelDataMap())
            {
                int ch = iter.first;
                auto data = iter.second;
                // mhq10[ch]->SetPoint(mhq10[ch]->GetN(), data->amp, data->charge_10);
                // mhq50[ch]->SetPoint(mhq50[ch]->GetN(), data->amp, data->charge_50);
                // mhq90[ch]->SetPoint(mhq90[ch]->GetN(), data->amp, data->charge_90);
                // mhq[ch]->SetPoint(mhq[ch]->GetN(), data->amp, data->charge);
                mhq10[ch]->Fill(data->amp, data->charge_10);
                mhq50[ch]->Fill(data->amp, data->charge_50);
                mhq90[ch]->Fill(data->amp, data->charge_90);
                mhq[ch]->Fill(data->amp, data->charge);
                mhqpm2ns[ch]->Fill(data->amp, data->charge_pm2ns);
            }
            continue;

            double outX = 0.0, outY = 0.0;
            PRAlgorithms::LinearRct::PRinfo info;
            rct.PRfromWaveInfo(dataMap, info);
            outX = info.reconX;

            if (TMath::Abs(outX - beamX) < 100.0)
            {
                if (fAlreadyOutput)
                    continue;
                else
                {
                    std::cout << "Normal event: " << evt << std::endl;
                    std::cout << "recX: " << outX << ", beamX: " << beamX << ", pads:<" << info.pads.first << ", " << info.pads.second << "> signals:<" << info.signals.first << ", " << info.signals.second << ">" << "chs:<" << info.chs.first << ", " << info.chs.second << ">" << std::endl;
                    fAlreadyOutput = true;
                    continue;
                }
            }
            std::cout << "Large deviation detected! Event: " << evt << std::endl;
            std::cout << "recX: " << outX << ", beamX: " << beamX << ", pads:<" << info.pads.first << ", " << info.pads.second << "> signals:<" << info.signals.first << ", " << info.signals.second << ">" << "chs:<" << info.chs.first << ", " << info.chs.second << ">" << std::endl;
            extractor.TurnOnPlots("plots/PRData/" + scanName + "/Event" + std::to_string(evt));
            extractor.ExtractFromWFROOTReader(reader);

            for (auto iter : dataMap)
                std::cout << " ch " << iter.first << ": a = " << iter.second->amp << ", q = " << iter.second->charge << ", q_10 = " << iter.second->charge_10 << ", q_50 = " << iter.second->charge_50 << ", q_90 = " << iter.second->charge_90 << ", q_pm2ns = " << iter.second->charge_pm2ns << std::endl;
            std::cout << std::endl;

            std::cout << std::endl;
        }
        reader.CloseFile();

        if (counter % 10 == 10 - 1)
        {
            // Draw and save charge - amplitude correlation graphs
            for (int ch = 1; ch <= 4; ch++)
            {
                c->cd();
                mhq10[ch]->Draw("COLZ");
                c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_q10.png", ch));
                c->cd();
                mhq50[ch]->Draw("COLZ");
                c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_q50.png", ch));
                c->cd();
                mhq90[ch]->Draw("COLZ");
                c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_q90.png", ch));
                c->cd();
                mhq[ch]->Draw("COLZ");
                c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_q.png", ch));
                mhqpm2ns[ch]->Draw("COLZ");
                c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_qpm2ns.png", ch));
            }
        }
        counter++;
    }

    for (int ch = 1; ch <= 4; ch++)
    {
        c->cd();
        mhq10[ch]->Draw("COLZ");
        c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_q10.png", ch));
        mhq10[ch]->Write();
        c->cd();
        mhq50[ch]->Draw("COLZ");
        c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_q50.png", ch));
        mhq50[ch]->Write();
        c->cd();
        mhq90[ch]->Draw("COLZ");
        c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_q90.png", ch));
        mhq90[ch]->Write();
        c->cd();
        mhq[ch]->Draw("COLZ");
        c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_q.png", ch));
        mhq[ch]->Write();
        mhqpm2ns[ch]->Draw("COLZ");
        c->SaveAs(Form("plots/ChargeAmpCorrelation/ch%d_qpm2ns.png", ch));
        mhqpm2ns[ch]->Write();
    }
    filew->Close();

    return 1;
}

// main for position reconstruction
int main()
{
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

    rct.GetPRDataMap().AddPCMap(1, 4);
    rct.GetPRDataMap().AddPCMap(2, 3);
    rct.GetPRDataMap().AddPCMap(3, 2);
    rct.GetPRDataMap().AddPCMap(4, 1);
    rct.GetPRDataMap().InitializeMap();
    rct.SetPadNormFactorsFromCali();
    // rct.GetPRDataMap().Print();
    rct.Print();
    auto valid = rct.JudgeCalibrationValid();
    std::cout << "Calibration valid: " << valid << std::endl;

    WFDataProcessor::WFDataExtractor extractor({1, 2, 3, 4});
    WFDataProcessor::_extract_config rangeConfig{{0, 20.0}, 20.0, false, ""};
    for (const auto &ch : {1, 2, 3, 4})
        extractor.SetExtractConfig(ch, rangeConfig);

    WFDataProcessor::WFROOTReader reader({1, 2, 3, 4});

    gSystem->mkdir("CheckPRAlgo/", true);
    auto filew = new TFile("CheckPRAlgo/test_prdata_output.root", "RECREATE");
    filew->cd();
    auto treew = new TTree("prtree", "PR Data Tree");
    double beamX, beamY;
    double recX;
    treew->Branch("beamX", &beamX, "beamX/D");
    treew->Branch("beamY", &beamY, "beamY/D");
    treew->Branch("recX", &recX, "recX/D");

    // Scan file in PRData/
    auto dir = gSystem->OpenDirectory("PRData/");
    const char *ptr_entry;
    while (ptr_entry = gSystem->GetDirEntry(dir))
    {
        std::string entry = ptr_entry;
        if (entry == "." || entry == "..")
            continue;
        std::string fullPath = "PRData/" + entry;
        if (gSystem->AccessPathName(fullPath.c_str()))
            continue;
        std::cout << "\rProcessing file: " << fullPath << std::flush;
        reader.OpenFile(fullPath);
        // full path subtract ".root"
        std::string scanName = fullPath.substr(0, fullPath.length() - 5);
        _beam_xyz beamPos = ExtractBeamXYZFromString(scanName);
        for (int evt = 0; evt < reader.GetTree()->GetEntries(); evt++)
        {
            reader.GetTree()->GetEntry(evt);
            extractor.ExtractFromWFROOTReader(reader);

            auto dataMap = extractor.GetChannelDataMap();

            double outX = 0.0, outY = 0.0;
            PRAlgorithms::LinearRct::PRinfo info;
            rct.PRfromWaveInfo(dataMap, info);
            outX = info.reconX;
            outY = info.reconY;
            // std::cout << "Reconstructed position: " << outX << ", beam X: " << (beamPos.x + 1300) * 2.5 << std::endl;
            beamX = (beamPos.x + 1300) * 2.5;
            beamY = (beamPos.y - 1200) * 2.5;
            recX = outX;
            treew->Fill();
            // std::cout << " Pads: " << info.pads.first << ", " << info.pads.second << "; Signals: " << info.signals.first << ", " << info.signals.second << "; channels: " << info.chs.first << ", " << info.chs.second << std::endl;
        }

        reader.CloseFile();
    }

    filew->cd();
    treew->Write();
    filew->Close();

    return 1;
}
