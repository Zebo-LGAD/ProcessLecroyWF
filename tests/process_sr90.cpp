#include "lcparser.h"
#include "WFDataConverter.h"

#include <iostream>
#include <TSystem.h>

std::string gOutputFolder = "/mnt/e/Project-LGAD2/ACLGAD-Test/Beta-Sr90/ProcessedROOTFiles2/";
std::string gWFFolder = gOutputFolder + "plots/waveforms/";
std::string gOutputROOTFile = gOutputFolder + "AC-W1-190V-20C.root";
std::string gInputTRCFolder = "/mnt/e/Project-LGAD2/ACLGAD-Test/Beta-Sr90/Source/AC-W1-190V-20C/";

WFDataProcessor::_extract_config ch1_config{{-3, 0.0}, 30.0, false, ""};
WFDataProcessor::_extract_config ch2_config{{-1, 3}, 30.0, false, ""};
WFDataProcessor::_extract_config ch3_config{{-1, 3}, 30.0, false, ""};
WFDataProcessor::_extract_config ch4_config{{-0.5, 3.5}, 30.0, false, ""};

void plot()
{
    using namespace WFDataProcessor;

    // Plots:
    gSystem->mkdir(gOutputFolder.c_str(), true);
    gSystem->mkdir(gWFFolder.c_str(), true);

    WFDataProcessor::WFDataExtractor extractor({1, 2, 3, 4});

    extractor.SetExtractConfig(1, ch1_config);
    extractor.SetExtractConfig(2, ch2_config);
    extractor.SetExtractConfig(3, ch3_config);
    extractor.SetExtractConfig(4, ch4_config);

    extractor.SetForceMatch(true);
    extractor.OpenFile(gOutputROOTFile);
    extractor.TurnOnPlots(gWFFolder.c_str(), true);

    int maxFiles = 20;
    for (int idx_file = 0; idx_file < maxFiles; idx_file++)
        extractor.ExtractFromTRCFiles(gInputTRCFolder, idx_file);
}

int main()
{
    WFDataProcessor::WFDataExtractor extractor({1, 2, 3, 4});
    extractor.SetExtractConfig(1, ch1_config);
    extractor.SetExtractConfig(2, ch2_config);
    extractor.SetExtractConfig(3, ch3_config);
    extractor.SetExtractConfig(4, ch4_config);

    extractor.SetForceMatch(true);
    extractor.OpenFile(gOutputROOTFile);

    int maxFiles = 3739;
    for (int idx_file = 0; idx_file < maxFiles; idx_file++)
    {
        if (idx_file % 100 == 0)
            std::cout << "\rProcessing file " << idx_file << "/" << maxFiles << std::flush;
        if (idx_file == maxFiles - 1)
            std::cout <<"\r"<< std::endl;
        extractor.ExtractFromTRCFiles(gInputTRCFolder, idx_file);
    }

    return 0;
}