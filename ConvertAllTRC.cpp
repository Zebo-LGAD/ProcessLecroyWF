// plot_waveform.cpp
#include "lcparser/lcparser.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <vector>

// #define DEBUG_DRAW

std::vector<std::string> gDirs{"../W2Die7/7-11_110V", "../W2Die7/7-11_110V_Multitrigger", "../W2Die7/7-11_110V_Reftrigger"};
std::vector<int> gMaxFiles{2020, 41492, 99999};

int ConvertAllTRC(std::string sWriteFile, std::string sDataFolder, int maxFiles)
{
    gSystem->mkdir("../plots", true);

    auto fileWriter = new TFile(Form("%s", sWriteFile.c_str()), "RECREATE");

    auto tree = new TTree("tree", "Waveforms from LeCroy .trc files");
    tree->SetDirectory(fileWriter);

    ScopeData ch2sd, ch3sd;
    tree->Branch("ch2", &ch2sd);
    tree->Branch("ch3", &ch3sd);

    std::ifstream testFile;

    int entryCounter = 0;
    for (int idx_file = 0; idx_file < maxFiles; idx_file++)
    {
        // Print progress
        if (idx_file % 100 == 0)
            std::cout << "\rProcessing file " << idx_file << "/" << maxFiles << std::flush;
        if (idx_file == maxFiles - 1)
            std::cout << std::endl;

        // Generate file names
        char cnum[6];
        sprintf(cnum, "%05d", idx_file);
        std::string filename1 = "C2--Trace--" + std::string(cnum) + ".trc";
        std::string filename2 = "C3--Trace--" + std::string(cnum) + ".trc";
        std::string filePath1 = sDataFolder + "/" + filename1;
        std::string filePath2 = sDataFolder + "/" + filename2;

        // Check if files exist
        testFile.open(filePath1);
        if (!testFile)
        {
            std::cerr << "File not found: " << filePath1 << ", skipping." << std::endl;
            testFile.close();
            continue;
        }
        testFile.close();
        testFile.open(filePath2);
        if (!testFile)
        {
            std::cerr << "File not found: " << filePath2 << ", skipping." << std::endl;
            testFile.close();
            continue;
        }
        testFile.close();

        // Start Converting
        ScopeData::errorCodes err = ch2sd.InitData(filePath1);
        if (err != ScopeData::kSUCCESS)
        {
            std::cerr << "Error parsing file: " << filePath1 << ", error code: " << err << std::endl;
            continue;
        }

        err = ch3sd.InitData(filePath2);
        if (err != ScopeData::kSUCCESS)
        {
            std::cerr << "Error parsing file: " << filePath2 << ", error code: " << err << std::endl;
            continue;
        }

        tree->Fill();
        entryCounter++;

#ifdef DEBUG_DRAW
        // Debug
        auto tg2 = ch2sd.createGraph();
        auto tg3 = ch3sd.createGraph();

        static TCanvas *c = new TCanvas("c", "Waveforms", 800, 600);
        c->cd();
        tg2->Draw("AL");
        c->SaveAs(Form("../plots/ch2_%05d.png", idx_file));
        tg3->Draw("AL");
        c->SaveAs(Form("../plots/ch3_%05d.png", idx_file));

        delete tg2;
        delete tg3;
#endif
    }

    fileWriter->Write();
    fileWriter->Close();
    delete fileWriter;

    std::cout << "Totally filled " << entryCounter << "waveforms into the tree." << std::endl;
    std::cout << "Data saved to " << sWriteFile << std::endl;
    std::cout << std::endl;

    return 0;
}

std::string ExtractLastDir(const std::string &path)
{
    // Extract last part of the path as the folder name
    std::string sDataFolder = path.substr(path.find_last_of("/") + 1);
    if (sDataFolder.empty())
        sDataFolder = "data";
    return sDataFolder;
}

void ConvertAllTRC()
{
    return;
    
    std::string storageDir = "../Processed-Data/";
    gSystem->mkdir(Form("%s/", storageDir.c_str()), true);

    // Example: only convert the first directory for testing
    // std::string sWriteFile = storageDir + "/" + ExtractLastDir(gDirs[0]) + ".root";
    // std::cout << sWriteFile << std::endl;
    // ConvertAllTRC(sWriteFile, gDirs[0], gMaxFiles[0] + 10);

    for (size_t i = 0; i < gDirs.size(); i++)
    // for (size_t i = 0; i < 1; i++)
    {
        std::string sWriteFile = storageDir + "/" + ExtractLastDir(gDirs[i]) + ".root";
        std::cout << sWriteFile << std::endl;
        ConvertAllTRC(sWriteFile, gDirs[i], gMaxFiles[i] + 1);
    }
}
