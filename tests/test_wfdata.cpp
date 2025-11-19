#include "lcparser.h"
#include "WFDataConverter.h"

#include <iostream>
#include <TSystem.h>

// Test all functionalities in WFDataConverter
int main_testall()
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

#include <filesystem>
WFDataProcessor::_extract_config ch1_config{{-3, 0.0}, 30.0, false, ""};
WFDataProcessor::_extract_config ch2_config{{-1, 3}, 30.0, false, ""};
WFDataProcessor::_extract_config ch3_config{{-1, 3}, 30.0, false, ""};
WFDataProcessor::_extract_config ch4_config{{-0.5, 3.5}, 30.0, false, ""};

int main()
{
    return 1;
    using namespace WFDataProcessor;
    gSystem->mkdir("Processed-Data/plots/", true);

    WFtrc2ROOT converter({1, 2, 3, 4});
    converter.SetForceMatch(true);
    
    WFDataExtractor extractor({1, 2, 3, 4});
    extractor.SetExtractConfig(1, ch1_config);
    extractor.SetExtractConfig(2, ch2_config);
    extractor.SetExtractConfig(3, ch3_config);
    extractor.SetExtractConfig(4, ch4_config);
    extractor.SetForceMatch(true);
    extractor.OpenFile("Processed-Data/AC-W1-190V-20C-Extracted.root");

    std::string wfSourceDir = "/mnt/e/Project-LGAD2/ACLGAD-Test/Beta-Sr90/Source/AC-W1-190V-20C/";
    // Get file list via <filesystem>
    int find_channel = converter.GetChannelDataMap().begin()->first;
    std::vector<std::string> wfFileList;
    for (const auto &entry : std::filesystem::directory_iterator(wfSourceDir))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".trc" && entry.path().string().find("C" + std::to_string(find_channel)) != std::string::npos)
        {
            std::cout << "\rFound file: " << entry.path().string() << std::flush;
            wfFileList.push_back(entry.path().string());
        }
    }
    std::cout << "Totally " << wfFileList.size() << " .trc files found." << std::endl;
    std::sort(wfFileList.begin(), wfFileList.end());
    int maxFiles = wfFileList.size();

    // Test WFtrc2ROOT
    if (!converter.OpenFile("Processed-Data/AC-W1-190V-20C.root"))
    {
        std::cerr << "Failed to initialize write file." << std::endl;
        return -1;
    }

    for (int idx_file = 0; idx_file < maxFiles; idx_file++)
    {
        if (idx_file < 20)
            extractor.TurnOnPlots("Processed-Data/plots/file" + std::to_string(idx_file));
        else
            extractor.TurnOffPlots();
        std::cout << "\rProcessing file " << idx_file + 1 << " / " << maxFiles << ": " << wfFileList[idx_file] << std::flush;
        auto rtn = converter.ReadAllAndFill(wfSourceDir, idx_file);
        if (!rtn)
            break;
        rtn = extractor.ExtractFromWFtrc2ROOT(converter);
        if (!rtn)
            std::cerr << "Extraction failed at file index " << idx_file << std::endl;
    }
    converter.CloseFile();
    std::cout << "Data saved to Processed-Data/AC-W1-190V-20C.root" << std::endl;
    extractor.CloseFile();
    std::cout << "Extracted data saved to Processed-Data/AC-W1-190V-20C-Extracted.root" << std::endl;
}