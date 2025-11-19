// plot_waveform.cpp
#include "WFDataConverter.h"
#include "lcparser.h"

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
using namespace WFDataProcessor;

std::vector<std::string> gDirs{"../W2Die7/7-11_110V", "../W2Die7/7-11_110V_Multitrigger", "../W2Die7/7-11_110V_Reftrigger"};
std::vector<int> gMaxFiles{2020, 41492, 99999};

int WFDataProcessor::ConvertAllTRC(std::string sWriteFile, std::string sDataFolder, int maxFiles)
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

std::string WFDataProcessor::ExtractLastDir(const std::string &path)
{
    // Extract last part of the path as the folder name
    std::string sDataFolder = path.substr(path.find_last_of("/") + 1);
    if (sDataFolder.empty())
        sDataFolder = "data";
    return sDataFolder;
}

void WFDataProcessor::ConvertAllTRC()
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

bool JudgeAllChannelsHaveData(const std::map<int, bool> &chDataHasData)
{
    bool allHaveData = true;
    for (const auto &pair : chDataHasData)
        allHaveData &= pair.second;
    return allHaveData;
}

#include <TMath.h>
#include <TLine.h>
void WFDataProcessor::processWave(const double *t, const double *a, int Nsamples, _waveinfo &winfo, _extract_config config)
{
    // Extract configuration
    double _threshold = config.threshold;
    _signal_range search_range = config.search_range;
    auto need_draw = config.need_draw;
    auto savePrefix = config.savePrefix;
    // double _threshold = 5 * _ped_std_dev; // 5 sigma above pedestal

    winfo.valid = VALID; // Assume valid until proven otherwise
    if (Nsamples <= 0 || t == nullptr || a == nullptr)
    {
        winfo.valid = NO_WAVEFORM;
        return;
    }

    double _t_start = t[0] * 1.0e9;          // convert to ns
    double _t_end = t[Nsamples - 1] * 1.0e9; // convert to ns

    // int _idx_ch = ch - 1; // index for g_low_range_t0 and g_high_range_t0 arrays
    // Process the waveform to find max amplitude, charge, and time of threshold crossing
    double _max_a = -1.0e9;
    double _max_t = 0;
    int _sample_max = 0;
    double _dt = (t[1] - t[0]) * 1.0e9; // assuming uniform sampling, units in ns
    double _charge_10 = 0;
    double _charge_50 = 0;
    double _charge_90 = 0;
    double _charge_pm2ns = 0;
    double _charge_full = 0;

    double _toa = -100e9;
    double _xMeasure = 0;

    // Generate new array inside
    auto _t = new double[Nsamples];
    auto _a = new double[Nsamples];

    // Find sample point for up and down range
    int _sample_up = 0;
    int _sample_down = Nsamples - 1;

    for (int _sample_point = 0; _sample_point < Nsamples; _sample_point++)
    {
        _t[_sample_point] = t[_sample_point] * 1.0e9; // convert to ns
        _a[_sample_point] = a[_sample_point] * 1.0e3; // convert to mV
        if (_sample_point > 0)
        {
            double _last_point_time = _t[_sample_point - 1]; // convert to ns
            double _this_point_time = _t[_sample_point];     // convert to ns
            if (_last_point_time < search_range.first && _this_point_time >= search_range.first)
                _sample_up = _sample_point;
            if (_last_point_time < search_range.second && _this_point_time >= search_range.second)
                _sample_down = _sample_point;
        }
    }

    // Calculate the pedestal, mean and standard deviation of the first part in the waveform
    double _ped_sum = 0;
    double _ped_square_sum = 0;
    int _ped_start_count = 0;

    for (int _sample_point = 0; _sample_point < _sample_up; _sample_point++)
    {
        _ped_sum += _a[_sample_point];
        _ped_square_sum += _a[_sample_point] * _a[_sample_point];
        _ped_start_count++;
    }
    double _ped_start = _ped_start_count > 0 ? _ped_sum / _ped_start_count : 0;
    double _ped_start_std_dev = _ped_start_count > 0 ? TMath::Sqrt((_ped_square_sum / _ped_start_count) - (_ped_start * _ped_start)) : 0;

    // _ped = 0, _ped_std_dev = 0; // Reset pedestal and std dev for debugging

    // Find max amplitude and calculate charge between _sample_up and _sample_down
    for (int _sample_point = _sample_up; _sample_point <= _sample_down; _sample_point++)
    {
        double _amp_temp = _a[_sample_point] - _ped_start; // units in mV
        double temp_time = _t[_sample_point];              // convert to ns

        if (temp_time >= search_range.first && temp_time <= search_range.second)
        {
            _charge_full += _amp_temp * _dt;
            if (_amp_temp > _max_a)
            {
                _max_a = _amp_temp;
                _max_t = temp_time;
                _sample_max = _sample_point;
            }
        }
    }

    // Calculate charge from t_amp-2ns to t_amp+2ns
    int _t_amp_minus_2ns_sample = _sample_max - static_cast<int>(2.0 / _dt);
    int _t_amp_plus_2ns_sample = _sample_max + static_cast<int>(2.0 / _dt);
    if (_t_amp_minus_2ns_sample < 0)
        _t_amp_minus_2ns_sample = 0;
    if (_t_amp_plus_2ns_sample >= Nsamples)
        _t_amp_plus_2ns_sample = Nsamples - 1;
    for (int _sample_point = _t_amp_minus_2ns_sample; _sample_point <= _t_amp_plus_2ns_sample; _sample_point++)
    {
        double _amp_temp = _a[_sample_point] - _ped_start; // units in mV
        _charge_pm2ns += _amp_temp * _dt;
    }

    // Scan from _t_max_ down to find rise time from 10% to 90%
    double _t1_10 = -100e9;
    double _t1_50 = -100e9;
    double _t1_90 = -100e9;
    double _amp_10 = 0.1 * _max_a;
    double _amp_90 = 0.9 * _max_a;
    double _amp_50 = 0.5 * _max_a;

    bool _found_10 = false;
    bool _found_50 = false;
    bool _found_90 = false;

    bool _found_t1 = false;
    double _t1 = -100e9;
    for (int _sample_point = _sample_max; _sample_point >= _sample_up; _sample_point--)
    {
        double _amp_temp = _a[_sample_point] - _ped_start; // units in mV
        double temp_time = _t[_sample_point];              // convert to ns

        if (!_found_90 && _amp_temp <= _amp_90)
        {
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _amp_90, _ped_start, _t1_90);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _amp_90, _ped_start, _t1_90);

            _found_90 = true;
        }

        if (!_found_50 && _amp_temp <= _amp_50)
        {
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _amp_50, _ped_start, _t1_50);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _amp_50, _ped_start, _t1_50);
            _found_50 = true;
        }

        if (!_found_10 && _amp_temp <= _amp_10)
        {
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _amp_10, _ped_start, _t1_10);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _amp_10, _ped_start, _t1_10);
            _found_10 = true;
        }

        if (!_found_t1 && _amp_temp < _threshold && (_threshold < _max_a)) // only search for t1 if threshold is less than max amplitude
        {
            _found_t1 = true;
            _t1 = temp_time;
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _threshold, _ped_start, _t1);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _threshold, _ped_start, _t1);
        }

        if (!_found_10)
            _charge_10 += _amp_temp * _dt;
        if (!_found_50)
            _charge_50 += _amp_temp * _dt;
        if (!_found_90)
            _charge_90 += _amp_temp * _dt;

        if (_found_10 && _found_50 && _found_90 && _found_t1)
            break; // all found, exit loop
    }
    if (!_found_10)
        winfo.valid |= NO_T1_10_FOUND;
    if (!_found_50)
        winfo.valid |= NO_T1_50_FOUND;
    if (!_found_90)
        winfo.valid |= NO_T1_90_FOUND;
    if (!_found_t1)
        winfo.valid |= NO_T1_FOUND;

    // Scan from _t_max_ up to find decreasing time from 10% to 90%
    double _t2_10 = -100e9;
    double _t2_50 = -100e9;
    double _t2_90 = -100e9;

    _found_10 = false;
    _found_50 = false;
    _found_90 = false;
    bool _found_t2 = false;
    double _t2 = -100e9;

    for (int _sample_point = _sample_max; _sample_point <= _sample_down; _sample_point++)
    {
        double _amp_temp = _a[_sample_point] - _ped_start; // units in mV
        double temp_time = _t[_sample_point];              // convert to ns

        if (!_found_90 && _amp_temp <= _amp_90)
        {
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _amp_90, _ped_start, _t1_90);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _amp_90, _ped_start, _t2_90);

            _found_90 = true;
        }

        if (!_found_50 && _amp_temp <= _amp_50)
        {
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _amp_50, _ped_start, _t1_50);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _amp_50, _ped_start, _t2_50);
            _found_50 = true;
        }

        if (!_found_10 && _amp_temp <= _amp_10)
        {
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _amp_10, _ped_start, _t1_10);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _amp_10, _ped_start, _t2_10);
            _found_10 = true;
        }

        if (!_found_t2 && _amp_temp < _threshold && (_threshold < _max_a)) // only search for t2 if threshold is less than max amplitude
        {
            _found_t2 = true;
            _t2 = temp_time;
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _threshold, _ped_start, _t2);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _threshold, _ped_start, _t2);
        }

        if (!_found_10)
            _charge_10 += _amp_temp * _dt;
        if (!_found_50)
            _charge_50 += _amp_temp * _dt;
        if (!_found_90)
            _charge_90 += _amp_temp * _dt;

        if (_found_10 && _found_50 && _found_90 && _found_t2)
            break; // all found, exit loop
    }
    if (!_found_10)
        winfo.valid |= NO_T2_10_FOUND;
    if (!_found_50)
        winfo.valid |= NO_T2_50_FOUND;
    if (!_found_90)
        winfo.valid |= NO_T2_90_FOUND;
    if (!_found_t2)
        winfo.valid |= NO_T2_FOUND;

    // // Scan from _sample_up to _sample_down to find time of arrival at 50% of max amplitude
    // double _toa_threshold = _amp_50; // 50% of max amplitude
    // for (int _sample_point = _sample_up; _sample_point < _sample_down; _sample_point++)
    // {
    //     double temp_time = _t[_sample_point]; // convert to ns
    //     if (temp_time < search_range.first)
    //         continue; // only search for toa in the region after the up_range_t0
    //     if (temp_time > search_range.second)
    //         break; // only search for toa in the region before the down_range_t0

    //     double _amp_temp = _a[_sample_point] - _ped_start;          // units in mV
    //     double _amp_temp_last = _a[_sample_point - 1] - _ped_start; // units in mV
    //     if (_amp_temp >= _toa_threshold && _amp_temp_last < _toa_threshold)
    //     {
    //         // auto flag = interpolateTOA(t, a, _sample_point, _toa_threshold, _ped_start, _toa);
    //         auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _toa_threshold, _ped_start, _toa);
    //         break;
    //     }
    // }
    _toa = _t1_50; // set toa to t1_50 for consistency

    // Calculate the pedestal, mean and standard deviation of the last 10% of the waveform
    double _ped_sum_end = 0;
    double _ped_square_sum_end = 0;
    int _ped_count_end = 0;

    for (int _sample_point = _sample_down; _sample_point < Nsamples; _sample_point++)
    {
        if (_t[_sample_point] > search_range.first)
        {
            _ped_sum_end += _a[_sample_point];
            _ped_square_sum_end += _a[_sample_point] * _a[_sample_point];
            _ped_count_end++;
        }
    }
    double _ped_end = _ped_count_end > 0 ? _ped_sum_end / _ped_count_end : 0;
    double _ped_end_std_dev = _ped_count_end > 0 ? TMath::Sqrt((_ped_square_sum_end / _ped_count_end) - (_ped_end * _ped_end)) : 0;
    // std::cout << "Channel " << ch << " pedestal: " << _ped_end << " mV, std dev: " << _ped_end_std_dev * 1.0e3 << " mV" << std::endl;

    winfo.nsamples = Nsamples;
    winfo.ped_start = _ped_start;
    winfo.ped_start_std_dev = _ped_start_std_dev;
    winfo.ped_end = _ped_end;
    winfo.ped_end_std_dev = _ped_end_std_dev;
    winfo.amp = _max_a;
    winfo.t_amp = _max_t;
    winfo.t1 = _t1;
    winfo.t1_10 = _t1_10;
    winfo.t1_50 = _t1_50;
    winfo.t1_90 = _t1_90;
    winfo.toa = _toa;
    winfo.t2 = _t2;
    winfo.t2_10 = _t2_10;
    winfo.t2_50 = _t2_50;
    winfo.t2_90 = _t2_90;
    winfo.charge_10 = _charge_10;
    winfo.charge_50 = _charge_50;
    winfo.charge_90 = _charge_90;
    winfo.charge_pm2ns = _charge_pm2ns;
    winfo.charge_full = _charge_full;
    winfo.charge = winfo.charge_50; // for backward compatibility

    if (!need_draw)
    {
        delete[] _t;
        delete[] _a;
        return;
    }
    // if (_toa == -100e9)
    // if (ch == 2 && amp > 50)
    {

        TGraph tg(Nsamples, _t, _a);
        tg.SetTitle(Form(";Time (ns);Amplitude (mV)"));

        double _tg_y_low = tg.GetYaxis()->GetXmin();
        double _tg_y_high = tg.GetYaxis()->GetXmax();
        // std::cout << "Y axis range: " << _tg_y_low << " to " << _tg_y_high << std::endl;

        static TCanvas *c1 = new TCanvas("c1", "Waveform", 800, 600);
        tg.Draw("AL");
        static int count = 0;
        // Gray for range & pedestal
        TLine liney1(search_range.first, _tg_y_low, search_range.first, _tg_y_high);
        liney1.SetLineColor(kGray);
        liney1.SetLineStyle(2);
        liney1.Draw("same");
        TLine liney2(search_range.second, _tg_y_low, search_range.second, _tg_y_high);
        liney2.SetLineColor(kGray);
        liney2.SetLineStyle(2);
        liney2.Draw("same");
        // Green for max
        TLine liney3(_max_t, _tg_y_low, _max_t, _tg_y_high);
        liney3.SetLineColor(kGreen + 2);
        liney3.SetLineStyle(2);
        liney3.Draw("same");
        // Red for t1 and t2 & threshold
        TLine liney4(_t1, _tg_y_low, _t1, _tg_y_high);
        liney4.SetLineColor(kRed);
        liney4.SetLineStyle(2);
        liney4.Draw("same");
        TLine liney5(_t2, _tg_y_low, _t2, _tg_y_high);
        liney5.SetLineColor(kRed);
        liney5.SetLineStyle(2);
        liney5.Draw("same");
        // Blue for toa & 50% threshold
        TLine liney6(_toa, _tg_y_low, _toa, _tg_y_high);
        liney6.SetLineColor(kBlue);
        liney6.SetLineStyle(2);
        liney6.Draw("same");

        TLine linex1(search_range.first, _ped_start, search_range.second, _ped_start); // convert to V
        linex1.SetLineColor(kGray);
        linex1.SetLineStyle(2);
        linex1.Draw("same");
        TLine linex1_up(search_range.first, _ped_start + 3 * _ped_start_std_dev, search_range.second, _ped_start + 3 * _ped_start_std_dev); // convert to V
        linex1_up.SetLineColor(kGray);
        linex1_up.SetLineStyle(2);
        linex1_up.Draw("same");
        TLine linex1_down(search_range.first, _ped_start - 3 * _ped_start_std_dev, search_range.second, _ped_start - 3 * _ped_start_std_dev); // convert to V
        linex1_down.SetLineColor(kGray);
        linex1_down.SetLineStyle(2);
        linex1_down.Draw("same");

        TLine linex2(search_range.first, _ped_start + _threshold, search_range.second, _ped_start + _threshold); // convert to V
        linex2.SetLineColor(kRed);
        linex2.SetLineStyle(2);
        linex2.Draw("same");

        TLine linex3(search_range.first, _ped_start + _amp_50, search_range.second, _ped_start + _amp_50); // convert to V
        linex3.SetLineColor(kBlue);
        linex3.SetLineStyle(2);
        linex3.Draw("same");

        TLine linex4(search_range.first, _ped_start + _max_a, search_range.second, _ped_start + _max_a); // convert to V
        linex4.SetLineColor(kGreen + 2);
        linex4.SetLineStyle(2);
        linex4.Draw("same");

        TLine linex5(search_range.second, _ped_end, _t_end, _ped_end); // convert to V
        linex5.SetLineColor(kGray + 2);
        linex5.SetLineStyle(2);
        linex5.SetLineWidth(3);
        linex5.Draw("same");

        if (savePrefix.empty())
            savePrefix = Form("../plots/waveform_%05d", count++);
        c1->SaveAs(Form("%s.png", savePrefix.c_str()));
        // std::cout << count << '\t' << _toa << '\t' << _t1 << '\t' << _t2 << '\t' << '\t' << _max_t << '\t' << t[_sample_up] * 1.0e9 << '\t' << t[_sample_down] * 1.0e9 << std::endl;
    }
    delete[] _t;
    delete[] _a;
    return;
}

std::string WFDataProcessor::GenerateScopeFileName(const std::string &folder, int channel, int idx_lecroy_wf, std::string mid_name, std::string ext)
{
    std::string filename = Form("C%d%s%05d%s", channel, mid_name.c_str(), idx_lecroy_wf, ext.c_str());
    std::string filepath = folder + "/" + filename;
    return filepath;
}

void WFDataProcessor::GenerateBranchForWaveInfo(TTree *tree, const std::string &branchNamePrefix, WFDataProcessor::_waveinfo *winfo)
{
    tree->Branch(Form("%s_valid", branchNamePrefix.c_str()), &winfo->valid, Form("%s_valid/i", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_nsamples", branchNamePrefix.c_str()), &winfo->nsamples, Form("%s_nsamples/I", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_ped_start", branchNamePrefix.c_str()), &winfo->ped_start, Form("%s_ped_start/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_ped_start_std_dev", branchNamePrefix.c_str()), &winfo->ped_start_std_dev, Form("%s_ped_start_std_dev/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_ped_end", branchNamePrefix.c_str()), &winfo->ped_end, Form("%s_ped_end/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_ped_end_std_dev", branchNamePrefix.c_str()), &winfo->ped_end_std_dev, Form("%s_ped_end_std_dev/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_amp", branchNamePrefix.c_str()), &winfo->amp, Form("%s_amp/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t_amp", branchNamePrefix.c_str()), &winfo->t_amp, Form("%s_t_amp/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t1", branchNamePrefix.c_str()), &winfo->t1, Form("%s_t1/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t1_10", branchNamePrefix.c_str()), &winfo->t1_10, Form("%s_t1_10/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t1_50", branchNamePrefix.c_str()), &winfo->t1_50, Form("%s_t1_50/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t1_90", branchNamePrefix.c_str()), &winfo->t1_90, Form("%s_t1_90/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_toa", branchNamePrefix.c_str()), &winfo->toa, Form("%s_toa/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_charge", branchNamePrefix.c_str()), &winfo->charge, Form("%s_charge/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t2", branchNamePrefix.c_str()), &winfo->t2, Form("%s_t2/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t2_10", branchNamePrefix.c_str()), &winfo->t2_10, Form("%s_t2_10/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t2_50", branchNamePrefix.c_str()), &winfo->t2_50, Form("%s_t2_50/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_t2_90", branchNamePrefix.c_str()), &winfo->t2_90, Form("%s_t2_90/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_q_10", branchNamePrefix.c_str()), &winfo->charge_10, Form("%s_q_10/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_q_50", branchNamePrefix.c_str()), &winfo->charge_50, Form("%s_q_50/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_q_90", branchNamePrefix.c_str()), &winfo->charge_90, Form("%s_q_90/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_q_pm2ns", branchNamePrefix.c_str()), &winfo->charge_pm2ns, Form("%s_q_pm2ns/D", branchNamePrefix.c_str()));
    tree->Branch(Form("%s_q_full", branchNamePrefix.c_str()), &winfo->charge_full, Form("%s_q_full/D", branchNamePrefix.c_str()));
}

void WFDataProcessor::SetBranchAddressToWaveInfo(TTree *tree, const std::string &branchNamePrefix, WFDataProcessor::_waveinfo *winfo)
{
    tree->SetBranchAddress(Form("%s_valid", branchNamePrefix.c_str()), &winfo->valid);
    tree->SetBranchAddress(Form("%s_nsamples", branchNamePrefix.c_str()), &winfo->nsamples);
    tree->SetBranchAddress(Form("%s_ped_start", branchNamePrefix.c_str()), &winfo->ped_start);
    tree->SetBranchAddress(Form("%s_ped_start_std_dev", branchNamePrefix.c_str()), &winfo->ped_start_std_dev);
    tree->SetBranchAddress(Form("%s_ped_end", branchNamePrefix.c_str()), &winfo->ped_end);
    tree->SetBranchAddress(Form("%s_ped_end_std_dev", branchNamePrefix.c_str()), &winfo->ped_end_std_dev);
    tree->SetBranchAddress(Form("%s_amp", branchNamePrefix.c_str()), &winfo->amp);
    tree->SetBranchAddress(Form("%s_t_amp", branchNamePrefix.c_str()), &winfo->t_amp);
    tree->SetBranchAddress(Form("%s_t1", branchNamePrefix.c_str()), &winfo->t1);
    tree->SetBranchAddress(Form("%s_t1_10", branchNamePrefix.c_str()), &winfo->t1_10);
    tree->SetBranchAddress(Form("%s_t1_50", branchNamePrefix.c_str()), &winfo->t1_50);
    tree->SetBranchAddress(Form("%s_t1_90", branchNamePrefix.c_str()), &winfo->t1_90);
    tree->SetBranchAddress(Form("%s_toa", branchNamePrefix.c_str()), &winfo->toa);
    tree->SetBranchAddress(Form("%s_charge", branchNamePrefix.c_str()), &winfo->charge);
    tree->SetBranchAddress(Form("%s_t2", branchNamePrefix.c_str()), &winfo->t2);
    tree->SetBranchAddress(Form("%s_t2_10", branchNamePrefix.c_str()), &winfo->t2_10);
    tree->SetBranchAddress(Form("%s_t2_50", branchNamePrefix.c_str()), &winfo->t2_50);
    tree->SetBranchAddress(Form("%s_t2_90", branchNamePrefix.c_str()), &winfo->t2_90);
    tree->SetBranchAddress(Form("%s_q_10", branchNamePrefix.c_str()), &winfo->charge_10);
    tree->SetBranchAddress(Form("%s_q_50", branchNamePrefix.c_str()), &winfo->charge_50);
    tree->SetBranchAddress(Form("%s_q_90", branchNamePrefix.c_str()), &winfo->charge_90);
    tree->SetBranchAddress(Form("%s_q_pm2ns", branchNamePrefix.c_str()), &winfo->charge_pm2ns);
    tree->SetBranchAddress(Form("%s_q_full", branchNamePrefix.c_str()), &winfo->charge_full);
}

WFDataProcessor::WFDataExtractor::WFDataExtractor(std::vector<int> channelsToRead) : VMultiChannelWriter<_waveinfo>(channelsToRead)
{
    for (int ch : channelsToRead)
        fmChExtractConfig[ch] = _extract_config();
}

bool WFDataProcessor::WFDataExtractor::AddChannel(int channel)
{
    auto rtn = VMultiIO::AddChannel(channel);
    if (rtn)
        fmChExtractConfig[channel] = _extract_config();
    return rtn;
}

bool WFDataProcessor::WFDataExtractor::GetExtractConfig(int channel, _extract_config &config) const
{
    if (fmChData.find(channel) == fmChData.end())
    {
        std::cerr << "Channel " << channel << " not found in data extractor." << std::endl;
        return false;
    }
    config = fmChExtractConfig.at(channel);
    return true;
}

bool WFDataProcessor::WFDataExtractor::SetExtractConfig(int channel, _extract_config config)
{
    if (fmChData.find(channel) == fmChData.end())
    {
        std::cerr << "Channel " << channel << " not found in data extractor." << std::endl;
        return false;
    }
    fmChExtractConfig[channel] = config;
    return true;
}

int WFDataProcessor::WFDataExtractor::SetExtractConfig(const std::map<int, _extract_config> &chRangeMap)
{
    int setCount = 0;
    for (const auto &pair : chRangeMap)
    {
        int channel = pair.first;
        _extract_config config = pair.second;
        auto rtn = SetExtractConfig(channel, config);
        if (rtn)
            setCount++;
    }
    return setCount;
}

bool WFDataProcessor::WFDataExtractor::ExtractFromScopeData(const std::map<int, ScopeData *> &chDataMap, const std::map<int, bool> &chDataHasDataMap)
{
    // Check if keys in three maps are consistent
    for (const auto &pair : fmChData)
    {
        int channel = pair.first;
        if (chDataMap.find(channel) == chDataMap.end() || chDataHasDataMap.find(channel) == chDataHasDataMap.end())
        {
            std::cerr << "Channel " << channel << " not found in input data maps." << std::endl;
            return false;
        }
    }
    // Already checked whether all channels have data when reading the data here
    if (fForceMatch)
    {
        // Check if all channels have data
        bool allHaveData = JudgeAllChannelsHaveData(chDataHasDataMap);
        if (!allHaveData)
            return false;
    }

    // Extract waveform information for each channel
    for (const auto &pair : fmChData)
    {
        int channel = pair.first;
        ScopeData *chData = chDataMap.at(channel);
        bool hasData = chDataHasDataMap.at(channel);
        if (!hasData)
        {
            fmChDataHasData[channel] = false;
            continue;
        }
        _extract_config config = fmChExtractConfig.at(channel);
        if (config.need_draw)
            config.savePrefix += "_counter_" + std::to_string(fExtractedCounter);

        WFDataProcessor::processWave(chData->getX().data(), chData->getY().data(), chData->getX().size(), *pair.second, config);
        fmChDataHasData[channel] = true;
    }
    // Increment extracted counter
    fExtractedCounter++;

    if (fTree)
        fTree->Fill();

    return true;
}

bool WFDataProcessor::WFDataExtractor::ExtractFromTRCFiles(const std::string &folder, int idx_file, std::string mid_name, std::string ext)
{
    std::map<int, ScopeData *> chDataMap;
    std::map<int, bool> chDataHasDataMap;
    for (const auto &pair : fmChData)
    {
        int channel = pair.first;
        chDataMap[channel] = new ScopeData();
        chDataHasDataMap[channel] = false;
    }

    ReadAllWF(folder, idx_file, chDataMap, chDataHasDataMap, mid_name, ext);
    auto rtn = ExtractFromScopeData(chDataMap, chDataHasDataMap);

    for (const auto &pair : chDataMap)
        delete pair.second;
    return rtn;
}

bool WFDataProcessor::WFDataExtractor::ExtractFromWFtrc2ROOT(const WFDataProcessor::WFtrc2ROOT &convertor)
{
    return ExtractFromScopeData(convertor.GetChannelDataMap(), convertor.GetChannelDataHasDataMap());
}

bool WFDataProcessor::WFDataExtractor::ExtractFromWFROOTReader(const WFDataProcessor::WFROOTReader &reader)
{
    std::map<int, ScopeData *> readDataMap;
    std::map<int, bool> readDataHasDataMap;
    for (auto &pair : reader.GetChannelDataMap())
    {
        int channel = pair.first;
        ScopeData *ptr = *(pair.second);
        readDataMap[channel] = ptr;
        readDataHasDataMap[channel] = true;
    }

    return ExtractFromScopeData(readDataMap, readDataHasDataMap);
}

bool WFDataProcessor::WFDataExtractor::TurnOnPlot(int channel, const std::string &savePrefix, bool bswitch)
{

    if (fmChExtractConfig.find(channel) == fmChExtractConfig.end())
    {
        std::cerr << "Channel " << channel << " not found in data extractor." << std::endl;
        return false;
    }
    fmChExtractConfig[channel].need_draw = bswitch;
    fmChExtractConfig[channel].savePrefix = savePrefix;
    return true;
}

int WFDataProcessor::WFDataExtractor::TurnOnPlots(const std::string &savePrefix, bool bswitch)
{
    int setCount = 0;
    for (const auto &pair : fmChExtractConfig)
    {
        std::string sPrefixForChannel = Form("%s_ch%d", savePrefix.c_str(), pair.first);
        auto rtn = TurnOnPlot(pair.first, sPrefixForChannel, bswitch);
        if (rtn)
            setCount++;
    }
    return setCount;
}

bool WFDataProcessor::WFDataExtractor::InitTree()
{
    auto rtn = VMultiIO::InitTree();
    if (!rtn)
        return false;

    fTree = new TTree("waveinfo", "Extracted waveform information");
    fTree->SetDirectory(fFile);

    for (const auto &pair : fmChData)
    {
        int channel = pair.first;
        GenerateBranchForWaveInfo(fTree, Form("ch%d", channel), pair.second);
    }
    return true;
}

void WFDataProcessor::WFDataExtractor::ClearMap()
{
    VMultiIO::ClearMap();
    fmChExtractConfig.clear();
    fExtractedCounter = 0;
}

bool WFDataProcessor::WFtrc2ROOT::ReadAllAndFill(const std::string &folder, int idx_file, std::string mid_name, std::string ext)
{
    auto rtn = ReadAllWF(folder, idx_file, fmChData, fmChDataHasData, mid_name, ext);

    // If force match is enabled, check if all channels have data
    if (fForceMatch && !rtn)
        return false;

    // Fill the tree
    if (fTree)
        fTree->Fill();
    return true;
}

bool WFDataProcessor::ReadAllWF(const std::string &folder, int idx_lecroy_wf, std::map<int, ScopeData *> &chDataMap, std::map<int, bool> &chDataHasDataMap, std::string mid_name, std::string ext)
{
    bool rtn = true;
    for (auto &pair : chDataMap)
    {
        int channel = pair.first;
        ScopeData *chData = pair.second;
        std::string filepath = GenerateScopeFileName(folder, channel, idx_lecroy_wf, mid_name, ext);
        ScopeData::errorCodes err = chData->InitData(filepath);

        if (err == ScopeData::kSUCCESS)
        {
            chDataHasDataMap[channel] = true;
            continue;
        }

        rtn = false;
        chDataHasDataMap[channel] = false;
        std::cerr << "Error parsing file: " << filepath << ", error code: " << err << std::endl;
        continue;
    }
    return rtn;
}

bool WFDataProcessor::WFtrc2ROOT::InitTree()
{
    auto rtn = VMultiIO::InitTree();
    if (!rtn)
        return false;

    fTree = new TTree("waveform", "Waveforms from LeCroy .trc files");
    fTree->SetDirectory(fFile);

    for (const auto &pair : fmChData)
    {
        int channel = pair.first;
        fTree->Branch(Form("ch%d", channel), pair.second);
    }
    return true;
}

bool WFDataProcessor::WFROOTReader::InitTree()
{
    auto rtn = VMultiIO::InitTree();
    if (!rtn)
        return false;

    fTree = (TTree *)fFile->Get("waveform");
    if (!fTree)
        return false;
    for (const auto &pair : fmChData)
    {
        int channel = pair.first;
        *pair.second = nullptr;
        fTree->SetBranchAddress(Form("ch%d", channel), pair.second);
    }
    return true;
}

bool WFDataProcessor::WFDataTreeReader::InitTree()
{
    auto rtn = VMultiIO::InitTree();
    if (!rtn)
        return false;

    fTree = (TTree *)fFile->Get("waveinfo");
    if (!fTree)
        return false;
    for (const auto &pair : fmChData)
    {
        int channel = pair.first;
        SetBranchAddressToWaveInfo(fTree, Form("ch%d", channel), pair.second);
    }
    return true;
}
