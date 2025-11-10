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

#include <iostream>

int test()
{
    gSystem->mkdir("../plots", true);
    std::string storageDir = "../Processed-Data/";
    gSystem->mkdir(Form("%s/", storageDir.c_str()), true);
    auto fileWriter = new TFile(Form("%s/waveforms.root", storageDir.c_str()), "RECREATE");

    auto c = new TCanvas("c", "Waveform", 800, 600);

    ScopeData scopeData;

    for (int i = 0; i < 10; i++)
    {
        std::string filePath = "../C2--Trace--0000" + std::to_string(i) + ".trc";
        ScopeData::errorCodes err = scopeData.InitData(filePath);
        if (err != ScopeData::kSUCCESS)
        {
            std::cerr << "Error parsing file: " << filePath << ", error code: " << err << std::endl;
            continue;
        }
        std::cout << scopeData << std::endl;

        // Draw plot
        c->cd();
        auto graph = scopeData.createGraph();
        graph->SetTitle(("Waveform " + std::to_string(i)).c_str());
        graph->SetLineColor(i + 1);
        graph->SetMarkerColor(i + 1);
        graph->SetMarkerStyle(20);
        graph->GetXaxis()->SetTitle("Time (s)");
        graph->GetYaxis()->SetTitle("Voltage (V)");
        graph->Draw("ALP");
        c->SaveAs(("../plots/waveform_" + std::to_string(i) + ".png").c_str());
        delete graph;
        graph = nullptr;
    }

    return 0;
}
