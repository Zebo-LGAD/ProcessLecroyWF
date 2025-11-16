#include <TTree.h>
#include <TFile.h>

#include <string>
#include <iostream>

#include "LinearRct.h"
#include "WFDataConverter.h"
std::map<int, int> gmTCTCh2Pad{{1, 0}, {2, 4}, {3, 3}, {4, 2}}; // Ch1 is used for trigger, so pad index is 0

// std::string gXReconstructionFile = "Linear-Fit-Data/20_38_47_width_86.5_190_Results";
// std::string gScopeDataFile = "../Beta-Sr90/ProcessedROOTFiles/AC-W1-190V-20C.root";
std::string gLinearPRFitFolder = "/mnt/e/Project-LGAD2/ACLGAD-Test/Position-Reconstruct-Algorithm/Linear-Fit-Data/";
std::string gLinearPRFitFile = gLinearPRFitFolder + "10_12_18_width_86.5_180_Results";
// std::string gScopeDataFile = "../Beta-Sr90/ProcessedROOTFiles/AC-W1-180V-20C.root";
std::string gScopeDataFolder = gLinearPRFitFolder + "../../Beta-Sr90/ProcessedROOTFiles2/";
std::string gOutputROOTFile = gScopeDataFolder + "AC-W1-190V-20C.root";

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

    rct.GetPRDataMap().AddPCMap(2, 4);
    rct.GetPRDataMap().AddPCMap(3, 3);
    rct.GetPRDataMap().AddPCMap(4, 2);
    rct.GetPRDataMap().InitializeMap();
    rct.SetPadNormFactorsFromCali();
    // rct.GetPRDataMap().Print();
    rct.Print();

    auto valid = rct.JudgeCalibrationValid();
    std::cout << "Calibration valid: " << valid << std::endl;
    return 1;

    WFDataProcessor::WFDataTreeReader reader({1, 2, 3, 4});
    auto rtn = reader.OpenFile(gOutputROOTFile);
    if (!rtn)
    {
        std::cerr << "Failed to open scope data file: " << gOutputROOTFile << std::endl;
        return -1;
    }

    auto mdata = gWFDataTreeReader->GetChannelDataMap();
    for (int entry = 0; entry < gWFDataTreeReader->GetTree()->GetEntries(); entry++)
    {
        gWFDataTreeReader->GetTree()->GetEntry(entry);

        std::vector<double> signals;
        std::cout << "Event " << entry << ":" << std::endl;
        for (const auto &pair : gmTCTCh2Pad)
        {
            int ch = pair.first;
            int pad = pair.second;
            auto winfo = mdata.at(ch);
            std::cout << "  " << "Ch " << ch << " valid: " << winfo->valid << std::endl;
            std::cout << "  " << "  " << "t1: " << winfo->t1 << ", amp: " << winfo->amp << std::endl;
        }
    }

    // rct.ReadCalibrationFile("Linear-Fit-Data/23_26_21_width_86.5_170_Results/Channel2and1.root");
    // rct.ReadCalibrationFile("Linear-Fit-Data/23_26_21_width_86.5_170_Results/Channel3and2.root");
    // rct.ReadCalibrationFile("Linear-Fit-Data/23_26_21_width_86.5_170_Results/Channel4and3.root");

    // rct.ReadCalibrationFiles("Linear-Fit-Data/23_26_21_width_86.5_170_Results");
    // rct.ReadCalibrationFiles("Linear-Fit-Data/20_38_47_width_86.5_190_Results");
    // rct.ReadCalibrationFiles("Linear-Fit-Data/10_12_18_width_86.5_180_Results");
    return 1;
}