#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TROOT.h"
#include <iostream>
#include <vector>

#include "lcparser/lcparser.h"
// #define DRAW_DEBUG
// #define VERBOSE_OUTPUT

double g_low_range_t0[4] = {0.0, -3.0, -10, 0.0}; // units in ns
double g_high_range_t0[4] = {0.0, 3.0, 10, 0.0};  // units in ns

enum InvalidCode
{
    VALID = 0,
    NO_WAVEFORM = 1,
    NO_T1_FOUND = 2,
    NO_T1_10_FOUND = 4,
    NO_T1_50_FOUND = 8,
    NO_T1_90_FOUND = 16,
    NO_T2_FOUND = 32,
    NO_T2_10_FOUND = 64,
    NO_T2_50_FOUND = 128,
    NO_T2_90_FOUND = 256,

    RANGE_ERROR = 512,
};

struct waveinfo
{
    int nsamples;             ///  Number of samples in the waveform
    double ped_start;         ///  Pedestal at the start of the waveform (in mV)
    double ped_start_std_dev; ///  Standard deviation of the pedestal at the start (in mV)
    double amp;               ///  Maximum amplitude of the waveform (in mV)
    double t_amp;             /// Time when the waveform reaches maximum amplitude (in ns)
    double charge;            ///  Charge of the waveform (in mV*ns)
    double t1;                ///  Time when the waveform first crosses the threshold (in ns)
    double t1_10;             /// Time when the waveform reaches 10% of maximum amplitude (in ns)
    double t1_50;             ///  Time of arrival at 50% of maximum amplitude (in ns)
    double t1_90;             /// Time when the waveform reaches 90% of maximum amplitude (in ns)
    double t2;                ///  Time when the waveform last crosses the threshold (in ns)
    double t2_10;             /// Time when the waveform reaches 10% of maximum amplitude (in ns)
    double t2_50;             ///  Time of arrival at 50% of maximum amplitude (in ns)
    double t2_90;             /// Time when the waveform reaches 90% of maximum amplitude (in ns)
    double toa;               ///  Interpolated time of arrival at 50% of maximum amplitude (in ns)
    double ped_end;           ///  Pedestal at the end of the waveform (in mV)
    double ped_end_std_dev;   ///  Standard deviation of the pedestal at the end (in mV)
    uint16_t valid;           ///  Flag indicating if the waveform is valid (1) or not (0)
};

// Get interpolated time of arrival at given threshold, using linear interpolation
template <typename T>
bool interpolateTOA(const T *tarray, const T *aarray, int sample_point, T threshold, T pedestal, T &result)
{
    if (sample_point == 0)
    {
        std::cerr << "Warning: First sample point is above threshold, cannot interpolate." << std::endl;
        result = tarray[sample_point]; // Use the current time if no previous point exists
        return false;
    }

    T t_before = tarray[sample_point - 1] * 1.0e9;            // convert to ns
    T t_after = tarray[sample_point] * 1.0e9;                 // convert to ns
    T a_before = aarray[sample_point - 1] * 1.0e3 - pedestal; // units in mV
    T a_after = aarray[sample_point] * 1.0e3 - pedestal;      // units in mV
    result = t_before + (t_after - t_before) * (threshold - a_before) / (a_after - a_before);
    return false;
}

template <typename T>
bool rangeCheck(T value, T range1, T range2)
{
    T low = std::min(range1, range2);
    T high = std::max(range1, range2);
    return (value >= low && value <= high);
}

template <typename T>
bool interpolateTOA2(const T *tarray, const T *aarray, int NSamples, int sample_point_around_threshold, T threshold, T pedestal, T &result)
{
    T t_before, t_after, a_before, a_after;
    bool in_range = 0;

    if (sample_point_around_threshold == 0)
    {
        t_before = tarray[sample_point_around_threshold];
        a_before = aarray[sample_point_around_threshold] - pedestal; // units in mV
        t_after = tarray[sample_point_around_threshold + 1];
        a_after = aarray[sample_point_around_threshold + 1] - pedestal; // units in mV
        in_range = rangeCheck(threshold, a_before, a_after);
    }
    else if (sample_point_around_threshold >= NSamples - 1)
    {
        t_before = tarray[sample_point_around_threshold - 1];
        a_before = aarray[sample_point_around_threshold - 1] - pedestal; // units in mV
        t_after = tarray[sample_point_around_threshold];
        a_after = aarray[sample_point_around_threshold] - pedestal; // units in mV
        in_range = rangeCheck(threshold, a_before, a_after);
    }
    else
    {
        t_before = tarray[sample_point_around_threshold - 1];
        a_before = aarray[sample_point_around_threshold - 1] - pedestal; // units in mV
        t_after = tarray[sample_point_around_threshold];
        a_after = aarray[sample_point_around_threshold] - pedestal; // units in mV
        in_range = rangeCheck(threshold, a_before, a_after);
        if (!in_range)
        {
            t_before = tarray[sample_point_around_threshold];
            a_before = aarray[sample_point_around_threshold] - pedestal; // units in mV
            t_after = tarray[sample_point_around_threshold + 1];
            a_after = aarray[sample_point_around_threshold + 1] - pedestal; // units in mV
            in_range = rangeCheck(threshold, a_before, a_after);
        }
    }

    if (!in_range)
    {
#ifdef VERBOSE_OUTPUT
        std::cerr << "Warning: Cannot find points around threshold for interpolation." << std::endl;
        std::cerr << "  Sample point: " << sample_point_around_threshold << " / " << NSamples - 1 << std::endl;
        std::cerr << "  t_before: " << t_before << " ns, t_after: " << t_after << " ns" << std::endl;
        std::cerr << "  Threshold: " << threshold << " mV, a_before: " << a_before << " mV, a_after: " << a_after << " mV" << std::endl;
        std::cerr << std::endl;
#endif
        result = tarray[sample_point_around_threshold]; // Use the current time if no previous point exists
        return false;
    }

    result = t_before + (t_after - t_before) * (threshold - a_before) / (a_after - a_before);
    return true;
}

/// @brief Process a single waveform to extract features such as pedestal, amplitude, charge, and time of arrival.
/// Charge in mV*ns, Amplitude in mV, Time in ns, signal was search between low_range_t0 and high_range_t0 for each channel.
/// @param ch Channel number (1-based index)
/// @param t Time array (in seconds)
/// @param a Amplitude array (in volts)
/// @param Nsamples Number of samples in the waveform
/// @param threshold [in] Threshold for t1 and t2 calculation (in mV), default is 20.0 mV
void processWave(int ch, const double *t, const double *a, int Nsamples, waveinfo &winfo, double threshold = 20.0)
{
    winfo.valid = VALID; // Assume valid until proven otherwise
    if (Nsamples <= 0 || t == nullptr || a == nullptr)
    {
        winfo.valid = NO_WAVEFORM;
        return;
    }

    double _t_start = t[0] * 1.0e9;          // convert to ns
    double _t_end = t[Nsamples - 1] * 1.0e9; // convert to ns

    int _idx_ch = ch - 1; // index for g_low_range_t0 and g_high_range_t0 arrays
    // Process the waveform to find max amplitude, charge, and time of threshold crossing
    double _max_a = -1.0e9;
    double _max_t = 0;
    int _sample_max = 0;
    double _charge = 0;
    double _dt = (t[1] - t[0]) * 1.0e9; // assuming uniform sampling, units in ns

    double _toa = -100e9;
    double _xMeasure = 0;

    double _threshold = threshold;
    // double _threshold = 5 * _ped_std_dev; // 5 sigma above pedestal

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
            if (_last_point_time < g_low_range_t0[_idx_ch] && _this_point_time >= g_low_range_t0[_idx_ch])
                _sample_up = _sample_point;
            if (_last_point_time < g_high_range_t0[_idx_ch] && _this_point_time >= g_high_range_t0[_idx_ch])
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

        if (temp_time >= g_low_range_t0[_idx_ch] && temp_time <= g_high_range_t0[_idx_ch])
        {
            _charge += _amp_temp * _dt;
            if (_amp_temp > _max_a)
            {
                _max_a = _amp_temp;
                _max_t = temp_time;
                _sample_max = _sample_point;
            }
        }
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

        if (!_found_t1 && _amp_temp > _threshold) // only search for t1 if threshold is less than max amplitude
        {
            _found_t1 = true;
            _t1 = temp_time;
            // Linear interpolation to find more accurate crossing time
            // auto flag = interpolateTOA(t, a, _sample_point, _threshold, _ped_start, _t1);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _threshold, _ped_start, _t1);
        }
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

    // Scan from _sample_up to _sample_down to find time of arrival at 50% of max amplitude
    double _toa_threshold = _amp_50; // 50% of max amplitude
    for (int _sample_point = _sample_up; _sample_point < _sample_down; _sample_point++)
    {
        double temp_time = _t[_sample_point]; // convert to ns
        if (temp_time < g_low_range_t0[_idx_ch])
            continue; // only search for toa in the region after the up_range_t0
        if (temp_time > g_high_range_t0[_idx_ch])
            break; // only search for toa in the region before the down_range_t0

        double _amp_temp = _a[_sample_point] - _ped_start;          // units in mV
        double _amp_temp_last = _a[_sample_point - 1] - _ped_start; // units in mV
        if (_amp_temp >= _toa_threshold && _amp_temp_last < _toa_threshold)
        {
            // auto flag = interpolateTOA(t, a, _sample_point, _toa_threshold, _ped_start, _toa);
            auto flag = interpolateTOA2(_t, _a, Nsamples, _sample_point, _toa_threshold, _ped_start, _toa);
            break;
        }
    }

    // Calculate the pedestal, mean and standard deviation of the last 10% of the waveform
    double _ped_sum_end = 0;
    double _ped_square_sum_end = 0;
    int _ped_count_end = 0;

    for (int _sample_point = _sample_down; _sample_point < Nsamples; _sample_point++)
    {
        if (_t[_sample_point] > g_low_range_t0[_idx_ch])
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
    winfo.charge = _charge;
    winfo.t2 = _t2;
    winfo.t2_10 = _t2_10;
    winfo.t2_50 = _t2_50;
    winfo.t2_90 = _t2_90;

#ifdef DRAW_DEBUG
    // if (_toa == -100e9)
    // if (ch == 2 && amp > 50)
    {

        TGraph tg(Nsamples, _t, _a);
        tg.SetTitle(Form(";Time (ns);Amplitude (mV)"));

        double _tg_y_low = tg.GetYaxis()->GetXmin();
        double _tg_y_high = tg.GetYaxis()->GetXmax();
        std::cout << "Y axis range: " << _tg_y_low << " to " << _tg_y_high << std::endl;

        static TCanvas *c1 = new TCanvas("c1", "Waveform", 800, 600);
        tg.Draw("AL");
        static int count = 0;
        // Gray for range & pedestal
        TLine liney1(g_low_range_t0[_idx_ch], _tg_y_low, g_low_range_t0[_idx_ch], _tg_y_high);
        liney1.SetLineColor(kGray);
        liney1.SetLineStyle(2);
        liney1.Draw("same");
        TLine liney2(g_high_range_t0[_idx_ch], _tg_y_low, g_high_range_t0[_idx_ch], _tg_y_high);
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

        TLine linex1(g_low_range_t0[_idx_ch], _ped_start, g_high_range_t0[_idx_ch], _ped_start); // convert to V
        linex1.SetLineColor(kGray);
        linex1.SetLineStyle(2);
        linex1.Draw("same");
        TLine linex1_up(g_low_range_t0[_idx_ch], _ped_start + 3 * _ped_start_std_dev, g_high_range_t0[_idx_ch], _ped_start + 3 * _ped_start_std_dev); // convert to V
        linex1_up.SetLineColor(kGray);
        linex1_up.SetLineStyle(2);
        linex1_up.Draw("same");
        TLine linex1_down(g_low_range_t0[_idx_ch], _ped_start - 3 * _ped_start_std_dev, g_high_range_t0[_idx_ch], _ped_start - 3 * _ped_start_std_dev); // convert to V
        linex1_down.SetLineColor(kGray);
        linex1_down.SetLineStyle(2);
        linex1_down.Draw("same");

        TLine linex2(g_low_range_t0[_idx_ch], _ped_start + _threshold, g_high_range_t0[_idx_ch], _ped_start + _threshold); // convert to V
        linex2.SetLineColor(kRed);
        linex2.SetLineStyle(2);
        linex2.Draw("same");

        TLine linex3(g_low_range_t0[_idx_ch], _ped_start + _toa_threshold, g_high_range_t0[_idx_ch], _ped_start + _toa_threshold); // convert to V
        linex3.SetLineColor(kBlue);
        linex3.SetLineStyle(2);
        linex3.Draw("same");

        TLine linex4(g_low_range_t0[_idx_ch], _ped_start + _max_a, g_high_range_t0[_idx_ch], _ped_start + _max_a); // convert to V
        linex4.SetLineColor(kGreen + 2);
        linex4.SetLineStyle(2);
        linex4.Draw("same");

        TLine linex5(g_high_range_t0[_idx_ch], _ped_end, _t_end, _ped_end); // convert to V
        linex5.SetLineColor(kGray + 2);
        linex5.SetLineStyle(2);
        linex5.SetLineWidth(3);
        linex5.Draw("same");

        c1->SaveAs(Form("../plots/waveform_ch%d_toa_error_count_%d.png", ch, count++));
        std::cout << count << '\t' << _toa << '\t' << _t1 << '\t' << _t2 << '\t' << '\t' << _max_t << '\t' << t[_sample_up] * 1.0e9 << '\t' << t[_sample_down] * 1.0e9 << std::endl;
    }
#endif
    delete[] _t;
    delete[] _a;
    return;
}

void ProcessTree(std::string sFileName = "../Processed-Data/7-11_110V_Reftrigger.root")
// void ProcessTree(std::string sFileName = "../Processed-Data/7-11_110V.root")
// void ProcessTree(std::string sFileName = "../Processed-Data/7-11_110V_Multitrigger.root")
{
    // Open the ROOT file
    auto fileRead = new TFile(sFileName.c_str(), "READ");
    if (!fileRead || fileRead->IsZombie())
    {
        std::cerr << "Error: Could not open file " << sFileName << std::endl;
        return;
    }

    auto treeRead = static_cast<TTree *>(fileRead->Get("tree"));
    if (!treeRead)
    {
        std::cerr << "Error: Tree not found in file." << std::endl;
        fileRead->Close();
        return;
    }

    // Set up READ branches
    ScopeData *ch2 = NULL, *ch3 = NULL;
    treeRead->SetBranchAddress("ch2", &ch2);
    treeRead->SetBranchAddress("ch3", &ch3);

    // Generate output file path
    std::string sWriteDir = "../Processed-Data";

    gSystem->mkdir(Form("%s", sWriteDir.c_str()), true);
    // Get Last part of the input file name after the last '/'
    std::string baseName = sFileName.substr(sFileName.find_last_of("/\\") + 1);
    // Remove the file extension
    baseName = baseName.substr(0, baseName.find_last_of('.'));

    gSystem->mkdir(Form("%s/%s", sWriteDir.c_str(), baseName.c_str()), true);
    std::string sWriteFile = baseName + "-Extract.root";
    std::string sFullWritePath = sWriteDir + "/" + baseName + "/" + sWriteFile;

    auto fileWrite = new TFile(sFullWritePath.c_str(), "RECREATE");
    auto treeWrite = new TTree("tree", "Extracted Waveform Data");

    waveinfo winfo2, winfo3;
    // Set up WRITE branches
    treeWrite->Branch("ch2_nsamples", &winfo2.nsamples, "ch2_nsamples/I");
    treeWrite->Branch("ch2_ped_start", &winfo2.ped_start, "ch2_ped_start/D");
    treeWrite->Branch("ch2_ped_start_sigma", &winfo2.ped_start_std_dev, "ch2_ped_start_sigma/D");
    treeWrite->Branch("ch2_amp", &winfo2.amp, "ch2_amp/D");
    treeWrite->Branch("ch2_t_amp", &winfo2.t_amp, "ch2_t_amp/D");
    treeWrite->Branch("ch2_t1", &winfo2.t1, "ch2_t1/D");
    treeWrite->Branch("ch2_charge", &winfo2.charge, "ch2_charge/D");
    treeWrite->Branch("ch2_t1_10", &winfo2.t1_10, "ch2_t1_10/D");
    treeWrite->Branch("ch2_t1_50", &winfo2.t1_50, "ch2_t1_50/D");
    treeWrite->Branch("ch2_t1_90", &winfo2.t1_90, "ch2_t1_90/D");
    treeWrite->Branch("ch2_t2", &winfo2.t2, "ch2_t2/D");
    treeWrite->Branch("ch2_t2_10", &winfo2.t2_10, "ch2_t2_10/D");
    treeWrite->Branch("ch2_t2_50", &winfo2.t2_50, "ch2_t2_50/D");
    treeWrite->Branch("ch2_t2_90", &winfo2.t2_90, "ch2_t2_90/D");
    treeWrite->Branch("ch2_toa", &winfo2.toa, "ch2_toa/D");
    treeWrite->Branch("ch2_ped_end", &winfo2.ped_end, "ch2_ped_end/D");
    treeWrite->Branch("ch2_ped_end_sigma", &winfo2.ped_end_std_dev, "ch2_ped_end_sigma/D");
    treeWrite->Branch("ch2_valid", &winfo2.valid, "ch2_valid/s");

    treeWrite->Branch("ch3_nsamples", &winfo3.nsamples, "ch3_nsamples/I");
    treeWrite->Branch("ch3_ped_start", &winfo3.ped_start, "ch3_ped_start/D");
    treeWrite->Branch("ch3_ped_start_sigma", &winfo3.ped_start_std_dev, "ch3_ped_start_sigma/D");
    treeWrite->Branch("ch3_amp", &winfo3.amp, "ch3_amp/D");
    treeWrite->Branch("ch3_t_amp", &winfo3.t_amp, "ch3_t_amp/D");
    treeWrite->Branch("ch3_charge", &winfo3.charge, "ch3_charge/D");
    treeWrite->Branch("ch3_t1", &winfo3.t1, "ch3_t1/D");
    treeWrite->Branch("ch3_t1_10", &winfo3.t1_10, "ch3_t1_10/D");
    treeWrite->Branch("ch3_t1_50", &winfo3.t1_50, "ch3_t1_50/D");
    treeWrite->Branch("ch3_t1_90", &winfo3.t1_90, "ch3_t1_90/D");
    treeWrite->Branch("ch3_t2", &winfo3.t2, "ch3_t2/D");
    treeWrite->Branch("ch3_t2_10", &winfo3.t2_10, "ch3_t2_10/D");
    treeWrite->Branch("ch3_t2_50", &winfo3.t2_50, "ch3_t2_50/D");
    treeWrite->Branch("ch3_t2_90", &winfo3.t2_90, "ch3_t2_90/D");
    treeWrite->Branch("ch3_toa", &winfo3.toa, "ch3_toa/D");
    treeWrite->Branch("ch3_ped_end", &winfo3.ped_end, "ch3_ped_end/D");
    treeWrite->Branch("ch3_ped_end_sigma", &winfo3.ped_end_std_dev, "ch3_ped_end_sigma/D");
    treeWrite->Branch("ch3_valid", &winfo3.valid, "ch3_valid/s");

    Long64_t nEntries = treeRead->GetEntries();
#ifdef DRAW_DEBUG
    gSystem->mkdir("../plots", true);
    // auto c = new TCanvas("c", "Waveform Canvas", 800, 600);
#endif

    // Record CPU time for only first 100 entries
    clock_t start_time = clock();

#ifdef DRAW_DEBUG
    // for (Long64_t entry = 0; entry < 100000; ++entry)
    for (Long64_t entry = 0; entry < 100; ++entry)
#else
    for (Long64_t entry = 0; entry < nEntries; ++entry)
#endif
    {
        if (entry % 1000 == 0)
            std::cout << "Processing entry " << entry << " / " << nEntries - 1 << "\r" << std::flush;
        if (entry == nEntries - 1)
            std::cout << "Processing entry " << entry << " / " << nEntries - 1 << "\r" << std::endl;

        treeRead->GetEntry(entry);
        if (ch2)
            processWave(2, ch2->getX().data(), ch2->getY().data(), ch2->getX().size(), winfo2, 50);
        else
            winfo2.valid = NO_WAVEFORM;
        if (ch3)
            processWave(3, ch3->getX().data(), ch3->getY().data(), ch3->getX().size(), winfo3, 100);
        else
            winfo3.valid = NO_WAVEFORM;
        treeWrite->Fill();
    }

    fileWrite->cd();
    treeWrite->Write();
    fileWrite->Close();
    fileRead->Close();

    clock_t end_time = clock();
    double cpu_time_used = double(end_time - start_time) / CLOCKS_PER_SEC;
    std::cout << "Processed " << nEntries << " entries in " << cpu_time_used << " seconds." << std::endl;

    delete fileRead;
    delete fileWrite;
}