#include "lcparser.h"
#include "WFDataConverter.h"

#include <iostream>
#include <TSystem.h>

int main()
{
    using namespace WFDataProcessor;
    gSystem->mkdir("../Processed-Data/plots/", true);

    WFtrc2ROOT converter({1, 2, 3, 4});
    converter.SetForceMatch(true);
    std::string wfSourceDir = "/mnt/e/Project-LGAD2/ACLGAD-Test/Beta-Sr90/Source/AC-W1-190V-20C/";
    int maxFiles = 50;

    // Test WFtrc2ROOT
    if (!converter.OpenFile("../Processed-Data/AC-W1-190V-20C.root"))
    {
        std::cerr << "Failed to initialize write file." << std::endl;
        return -1;
    }

    for (int idx_file = 0; idx_file < maxFiles; idx_file++)
    {
        auto rtn = converter.ReadAllAndFill(wfSourceDir, idx_file);
        if (!rtn)
            break;
    }
    converter.CloseFile();
    std::cout << "Data saved to ../Processed-Data/AC-W1-190V-20C.root" << std::endl;

    // Test WFROOTReader
    WFROOTReader reader({1, 2, 3, 4});
    if (!reader.OpenFile("../Processed-Data/AC-W1-190V-20C.root"))
    {
        std::cerr << "Failed to open read file." << std::endl;
        return -1;
    }
    TTree *tree = reader.GetTree();
    auto chDataMap = reader.GetChannelDataMap();
    for (int entry = 0; entry < tree->GetEntries(); entry++)
    {
        tree->GetEntry(entry);
        std::cout << "Entry " << entry << ":" << std::endl;
        for (const auto &pair : chDataMap)
        {
            int channel = pair.first;
            auto data = *pair.second;
            std::cout << "  Channel " << channel << ": "
                      << "Waveform Source = " << data->getWaveSource() << std::endl;
        }
    }
    reader.CloseFile();

    // Test WFDataExtractor
    WFDataExtractor extractor({1, 2, 3, 4});
    if (!extractor.OpenFile("../Processed-Data/AC-W1-190V-20C-Extracted.root"))
    {
        std::cerr << "Failed to open extract file." << std::endl;
        return -1;
    }
    // Test WFDataExtractor::ExtractFromWFtrc2ROOT
    converter.OpenFile("../Processed-Data/AC-W1-190V-20C.root");
    for (int idx_file = 0; idx_file < maxFiles; idx_file++)
    {
        auto rtn = converter.ReadAllAndFill(wfSourceDir, idx_file);
        if (!rtn)
            break;
        rtn = extractor.ExtractFromWFtrc2ROOT(converter);
        if (!rtn)
            std::cerr << "Extraction failed at file index " << idx_file << std::endl;
        auto chDataMap = extractor.GetChannelDataMap();
        for (const auto &pair : chDataMap)
        {
            int channel = pair.first;
            auto data = *pair.second;
            std::cout << "Extracted from file index " << idx_file << ", Channel " << channel << ": "
                      << "Amplitude = " << data.amp << " mV, "
                      << "Time of Amplitude = " << data.t_amp << " ns" << std::endl;
        }
    }
    extractor.CloseFile();
    converter.CloseFile();

    // Test WFDataTreeReader
    WFDataTreeReader treeReader({1, 2, 3, 4});
    if (!treeReader.OpenFile("../Processed-Data/AC-W1-190V-20C-Extracted.root"))
    {
        std::cerr << "Failed to open read file." << std::endl;
        return -1;
    }
    TTree *readTree = treeReader.GetTree();
    auto readTreeDataMap = treeReader.GetChannelDataMap();
    for (int entry = 0; entry < readTree->GetEntries(); entry++)
    {
        readTree->GetEntry(entry);
        std::cout << "Entry " << entry << ":" << std::endl;
        for (const auto &pair : readTreeDataMap)
        {
            int channel = pair.first;
            auto data = *pair.second;
            std::cout << "  Channel " << channel << ": "
                      << "Amplitude = " << data.amp << " mV, "
                      << "Time of Amplitude = " << data.t_amp << " ns" << std::endl;
        }
    }
    treeReader.CloseFile();

    return 0;
}