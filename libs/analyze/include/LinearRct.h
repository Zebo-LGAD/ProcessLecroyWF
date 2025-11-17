// #include <TFile.h>
// #include <TParameter.h>
// #include <TH1D.h>
// #include <TF1.h>

#ifndef LINEARRCT_HH
#define LINEARRCT_HH

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include "VPRAlgorithm.h"
class TF1;
class TH1D;

namespace WFDataProcessor
{
    struct _waveinfo;
}

// Position Reconstruction Algorithms namespace
namespace PRAlgorithms
{

    // typedef std::pair<int, int> _neighbor_pads;
    bool operator==(const _neighbor_pads &lhs, const _neighbor_pads &rhs);

    typedef struct _linear_fit_data
    {
        _linear_fit_data() = default;
        _linear_fit_data(std::string config_file_path);
        int _padIndex = -1;
        TF1 *_pol1 = NULL;
        std::pair<double, double> _edgeX;
        _neighbor_pads _edgePads{0, 0};
        std::vector<double> _factors;
        TH1D *_hist = NULL;
    } SectionData;

    typedef struct _draw_config
    {
        _draw_config() = default;
        bool _need_draw = false;
        std::string _save_folder = "";
        std::string _save_prefix = "";
    } DrawConfig;

    bool ReadLinearPRFit(const std::string &filePath, SectionData &section, _draw_config draw_config = {});

    class LinearRct : public VPRAlgorithm
    {
    public:
        LinearRct() = default;
        ~LinearRct() override;

        void Clear();

        // Read results from calibration files
        bool ReadLinearPRFits(const std::string &folderPath);
        bool JudgeCalibrationValid() const;

        // Input signals from oscilloscope channels
        double ReconstructX(const std::vector<double> &signals) const override { return 0.0; };
        double ReconstructY(const std::vector<double> &signals) const override { return 0.0; };
        double ReconstructX(PRAlgorithms::_neighbor_pads pads, _neighbor_pads_signals signals) const;

        typedef struct _lr_pr_info
        {
            _neighbor_pads pads;
            _neighbor_chs chs;
            _neighbor_pads_signals signals;
            double reconX;
            double reconY;
        } PRinfo;
        
        // Input waveinfo from oscilloscope channels
        bool PRfromWaveInfo(const std::map<int, WFDataProcessor::_waveinfo *> &chWaveInfoMap, double &outX, double &outY) const;
        bool PRfromWaveInfo(const std::map<int, WFDataProcessor::_waveinfo *> &chWaveInfoMap, PRinfo & info) const;

        const SectionData& GetPadRightSection(int pad) const;
        const std::map<_neighbor_pads, SectionData> &GetSectionMap() const { return fmSectionDB; };
        const std::map<int, double> &GetPadNormFactors() const { return fmPadNormFactors; };
        bool SetPadNormFactors(const std::map<int, double> &padNormFactors);
        bool SetPadNormFactorsFromCali();

        ChannelPadMap &GetCaliDataMap() { return fmCalibrationData; };
        const ChannelPadMap &GetCaliDataMapConst() const { return fmCalibrationData; };

        void Print(std::ostream &os = std::cout) const { os << (*this); };
        friend std::ostream &operator<<(std::ostream &os, const LinearRct &rct);

    private:
        ChannelPadMap fmCalibrationData;                   // This map is used for Calibration data
        std::map<_neighbor_pads, SectionData> fmSectionDB; // _neighbo_pads means (left pad, right pad)
        std::map<int, double> fmPadNormFactors;         // Data channel normalization factors, use different methods to calibrate amplitude differences among channels
    };

    std::ostream &operator<<(std::ostream &os, const PRAlgorithms::LinearRct &rct);

}

#endif