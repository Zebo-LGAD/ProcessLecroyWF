#include "LinearRct.h"
#include <iostream>
#include <TFile.h>
#include <TParameter.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>

// #define DEBUG_DRAW

PRAlgorithms::LinearRct::~LinearRct()
{
    Clear();
}

void PRAlgorithms::LinearRct::Clear()
{
    for (auto &pair : fmSectionDB)
    {
        if (pair.second._pol1)
            delete pair.second._pol1;
        if (pair.second._hist)
            delete pair.second._hist;
    }
    fmSectionDB.clear();
    fmPadNormFactors.clear();
    fmCalibrationData.Clear();
    VPRAlgorithm::GetPRDataMap().Clear();
}

bool PRAlgorithms::LinearRct::ReadLinearPRFits(const std::string &folderPath)
{
    auto caliMap = GetCaliDataMap();
    // Get all pads with right neighbor
    for (auto iter : caliMap.GetValidPadX())
    {
        auto pad_xy = iter;
        auto pads = caliMap.GetPadXYToNeighborPadsX().at(pad_xy);
        int leftPad = pads.first;
        int rightPad = pads.second;
        int leftCh = caliMap.GetChannel(leftPad);
        int rightCh = caliMap.GetChannel(rightPad);

        // Try first time with "Channel${leftCh}and${rightCh}.root"
        std::string filePath = folderPath + "/Channel" + std::to_string(leftCh) + "and" + std::to_string(rightCh) + ".root";
        PRAlgorithms::SectionData temp_section;
        temp_section._edgePads = pads;

        auto rtn = ReadLinearPRFit(filePath, temp_section);
        if (rtn)
        {
            fmSectionDB[temp_section._edgePads] = (temp_section);
            continue;
        }

        // Try second time with "Channel${rightCh}and${leftCh}.root"
        filePath = folderPath + "/Channel" + std::to_string(rightCh) + "and" + std::to_string(leftCh) + ".root";
        rtn = ReadLinearPRFit(filePath, temp_section);
        if (rtn)
        {
            fmSectionDB[temp_section._edgePads] = (temp_section);
            continue;
        }
        std::cout << "Failed to read Linear PR Fit file for pads " << leftPad << " and " << rightPad << std::endl;
        return false;
    }

    if (fmSectionDB.empty())
    {
        std::cout << "No section data read from calibration files!" << std::endl;
        return false;
    }

    return true;
}

bool PRAlgorithms::LinearRct::JudgeCalibrationValid() const
{
    auto &datamap = GetPRDataMapConst();
    auto &calimap = GetCaliDataMapConst();
    auto rtn = (datamap <= calimap);
    if (!rtn)
    {
        std::cout << "Calibration data map is not compatible with PR data map!" << std::endl;
        Print();
    }

    // Judge whether all pad normalization factors are set
    auto &dataMap = GetPRDataMapConst();
    for (const auto &pair : dataMap.GetMapP2C())
    {
        int pad = pair.first;
        if (fmPadNormFactors.find(pad) == fmPadNormFactors.end())
        {
            int channel = dataMap.GetChannel(pad);
            std::cout << "Normalization factor for pad " << pad << " channel: " << channel << " not found!" << std::endl;
            return false;
        }
    }

    return rtn;
}

double PRAlgorithms::LinearRct::ReconstructX(PRAlgorithms::_neighbor_pads pads, _neighbor_pads_signals signals) const
{
    auto iter = fmSectionDB.find(pads);
    if (iter == fmSectionDB.end())
    {
        std::cout << "Section data for pads (" << pads.first << ", " << pads.second << ") not found!" << std::endl;
        return 0.0;
    }

    auto section = iter->second;

    // auto fmPadFactors = section._factors;
    auto leftPad = pads.first;
    auto rightPad = pads.second;

    double factorLeft = fmPadNormFactors.at(leftPad);
    double factorRight = fmPadNormFactors.at(rightPad);
    double signalLeft = signals.first;
    double signalRight = signals.second;
    double normedSignalLeft = signalLeft * factorLeft;
    double normedSignalRight = signalRight * factorRight;

    // Avoid too small signals
    if (normedSignalLeft < 10 && normedSignalRight < 10) // approximately 10/20=0.5 fC
        return 0.0;

    static int count_left_small = 0;
    static int count_right_small = 0;
    static int count_normal = 0;
    if (normedSignalLeft < 0.15 * (normedSignalLeft + normedSignalRight))
    {
        double posrightEdge = section._edgeX.second;
        double padwidth = fmCalibrationData.GetPadSize(rightPad).first;
        return posrightEdge + padwidth / 2.0;
    }

    if (normedSignalRight < 0.15 * (normedSignalLeft + normedSignalRight))
    {
        double posleftEdge = section._edgeX.first;
        double padwidth = fmCalibrationData.GetPadSize(leftPad).first;
        return posleftEdge - padwidth / 2.0;
    }

    double ratio = (normedSignalLeft) / (normedSignalLeft + normedSignalRight);
    double _a0 = section._pol1->GetParameter(0);
    double _a1 = section._pol1->GetParameter(1);
    double x = (ratio - _a0) / _a1;

    return x;
}

#include "WFDataConverter.h"
// TODO: complete implementation
bool PRAlgorithms::LinearRct::PRfromWaveInfo(const std::map<int, WFDataProcessor::_waveinfo *> &chWaveInfoMap, double &outX, double &outY) const
{
    PRinfo info;
    auto rtn = PRfromWaveInfo(chWaveInfoMap, info);
    if (!rtn)
        return false;
    outX = info.reconX;
    outY = info.reconY;
    return true;
}

bool PRAlgorithms::LinearRct::PRfromWaveInfo(const std::map<int, WFDataProcessor::_waveinfo *> &chWaveInfoMap, PRinfo &info) const
{
    // Generate neighbor pads signals map
    std::map<_neighbor_pads, _neighbor_pads_signals> neighborPadsSignalsMap;
    std::map<_neighbor_pads, _neighbor_pads_signals> neighborPadsNormSignalsMap;
    auto dataMap = GetPRDataMapConst();
    for (auto iter : dataMap.GetValidPadX())
    {
        auto pad_xy = iter;
        auto pads = dataMap.GetPadXYToNeighborPadsX().at(pad_xy);
        int leftPad = pads.first;
        int rightPad = pads.second;
        double leftNormFactor = fmPadNormFactors.at(leftPad);
        double rightNormFactor = fmPadNormFactors.at(rightPad);

        auto channels = dataMap.GetNeighborPadsToChsX().at(pads);
        int leftCh = channels.first;
        int rightCh = channels.second;
        double leftSignal = chWaveInfoMap.at(leftCh)->charge;
        double rightSignal = chWaveInfoMap.at(rightCh)->charge;

        neighborPadsSignalsMap[pads] = std::make_pair(leftSignal, rightSignal);
        neighborPadsNormSignalsMap[pads] = std::make_pair(leftSignal * leftNormFactor, rightSignal * rightNormFactor);
    }

    // Compare which neighbor has the maximum normalized signal sum
    auto iter = std::max_element(neighborPadsNormSignalsMap.begin(), neighborPadsNormSignalsMap.end(),
                                 [](const std::pair<_neighbor_pads, _neighbor_pads_signals> &a, const std::pair<_neighbor_pads, _neighbor_pads_signals> &b)
                                 {
                                     double sumA = a.second.first + a.second.second;
                                     double sumB = b.second.first + b.second.second;
                                     return sumA < sumB;
                                 });
    if (iter == neighborPadsNormSignalsMap.end())
    {
        std::cout << "No valid neighbor pads found for position reconstruction!" << std::endl;
        return false;
    }
    auto selectedPads = iter->first;
    auto selectedSignals = neighborPadsSignalsMap.at(selectedPads);
    auto channels = dataMap.GetNeighborPadsToChsX().at(selectedPads);

    info.pads = selectedPads;
    info.chs = channels;
    info.signals = selectedSignals;
    info.reconX = ReconstructX(selectedPads, selectedSignals);
    info.reconY = 0.0; // Y reconstruction not implemented yet

    return true;
}

const PRAlgorithms::SectionData &PRAlgorithms::LinearRct::GetPadRightSection(int pad) const
{
    auto caliMap = GetCaliDataMapConst();
    int padRight;
    auto rtn = caliMap.GetRightPad(pad, padRight);
    auto iter = fmSectionDB.find(std::make_pair(pad, padRight));
    if (!rtn || iter == fmSectionDB.end())
    {
        std::cout << "Section data for pad " << pad << " and its right pad not found!" << std::endl;
        throw std::runtime_error("Section data not found");
    }
    return iter->second;
}

bool PRAlgorithms::LinearRct::SetPadNormFactors(const std::map<int, double> &padNormFactors)
{
    // Search all pads in data map, make sure each channel has a normalization factor
    auto &dataMap = GetPRDataMapConst();
    for (const auto &pair : dataMap.GetMapP2C())
    {
        int pad = pair.first;
        if (padNormFactors.find(pad) == padNormFactors.end())
        {
            int channel = dataMap.GetChannel(pad);
            std::cout << "Normalization factor for pad " << pad << " channel: " << channel << " not found!" << std::endl;
            return false;
        }
        fmPadNormFactors[pad] = padNormFactors.at(pad);
    }

    return true;
}

bool PRAlgorithms::LinearRct::SetPadNormFactorsFromCali()
{
    auto caliMap = GetCaliDataMap();
    // Generate pad normalization factors
    std::map<int, double> padNormFactors;
    auto datavec = fmSectionDB.begin()->second._factors;
    for (size_t idx = 0; idx < datavec.size(); idx++)
    {
        int channel = idx + 1;
        int pad = caliMap.GetPad(channel);
        padNormFactors[pad] = datavec[idx];
    }
    auto rtn = SetPadNormFactors(padNormFactors);
    return rtn;
}

bool PRAlgorithms::LinearRct::SetUniformPadSize(double padcolumnwidth, double padrowwidth)
{
    auto rtn1 = fmCalibrationData.SetUniformPadSize(padcolumnwidth, padrowwidth);
    auto rtn2 = VPRAlgorithm::SetUniformPadSize(padcolumnwidth, padrowwidth);
    return rtn1 && rtn2;
}

#include <TSystem.h>
#include <TLatex.h>
bool PRAlgorithms::operator==(const _neighbor_pads &lhs, const _neighbor_pads &rhs)
{
    return lhs.first == rhs.first && lhs.second == rhs.second;
}
bool PRAlgorithms::ReadLinearPRFit(const std::string &filePath, SectionData &section, _draw_config draw_config)
{
    auto need_draw = draw_config._need_draw;
    auto save_folder = draw_config._save_folder;
    auto save_prefix = draw_config._save_prefix;

    auto file = new TFile(filePath.c_str(), "READ");
    if (!file->IsOpen())
    {
        delete file;
        std::cout << ("Cannot open calibration file: " + filePath) << std::endl;
        return false;
    }

    section._pol1 = (TF1 *)file->Get("linearFitLine;1");
    double ledge = ((TParameter<double> *)file->Get("LeftPadEdge;1"))->GetVal();
    double redge = ((TParameter<double> *)file->Get("RightPadEdge;1"))->GetVal();
    section._edgeX = std::make_pair(ledge, redge);

    section._factors = (std::vector<double>)(*(std::vector<double> *)file->Get("Factors;1"));
    section._hist = (TH1D *)file->Get("linearRatioHist;1");
    section._hist->SetDirectory(nullptr); // Detach histogram from file

    delete file;

    if (!need_draw)
        return true;
    // file->ls();
    gSystem->mkdir(Form("%s", save_folder.c_str()), true);
    TLatex latex;
    std::cout << "Read section: " << ledge << " to " << redge << std::endl;
    std::cout << "  with " << section._factors.size() << " factors." << std::endl;
    for (size_t i = 0; i < section._factors.size(); i++)
        std::cout << "    Factor " << i << ": " << section._factors[i] << std::endl;
    std::cout << " and pol1: y = " << section._pol1->GetParameter(0) << " + " << section._pol1->GetParameter(1) << " * x" << std::endl;

    auto c = new TCanvas("c", "c", 800, 600);
    section._hist->Draw();
    section._pol1->Draw("same");
    latex.DrawLatexNDC(0.15, 0.85, Form("y = %.4f + %.4f * x/1000", section._pol1->GetParameter(0), section._pol1->GetParameter(1) * 1000));
    latex.DrawLatexNDC(0.15, 0.75, Form("File: %s", gSystem->BaseName(filePath.c_str())));
    std::string sSaveName = save_folder + "/" + save_prefix + ".png";
    c->SaveAs(sSaveName.c_str());
    delete c;

    std::cout << std::endl;
    return true;
}

std::ostream &PRAlgorithms::operator<<(std::ostream &os, const PRAlgorithms::LinearRct &rct)
{
    os << "Base class VPRAlgorithm info:" << std::endl;
    os << static_cast<const PRAlgorithms::VPRAlgorithm &>(rct);
    os << std::endl;

    os << "Section calibration map info:" << std::endl;
    os << rct.GetCaliDataMapConst();
    os << "LinearRct calibrated with " << rct.fmSectionDB.size() << " sections." << std::endl;
    os << "Section details:" << std::endl;
    auto caliMap = rct.GetCaliDataMapConst();
    for (auto iter : caliMap.GetValidPadX())
    {
        auto pad_xy = iter;
        auto pads = caliMap.GetPadXYToNeighborPadsX().at(pad_xy);
        auto chs = caliMap.GetNeighborPadsToChsX().at(pads);

        if (rct.fmSectionDB.find(pads) == rct.fmSectionDB.end())
        {
            os << " Section data for pads (" << pads.first << ", " << pads.second << ") not read!" << std::endl;
            continue;
        }
        auto sec = rct.fmSectionDB.at(pads);

        os << " Pads (" << sec._edgePads.first << ", " << sec._edgePads.second << "), EdgeX (" << sec._edgeX.first << ", " << sec._edgeX.second << "), scope channels (";
        os << rct.GetCaliDataMapConst().GetChannel(sec._edgePads.first) << ", ";
        os << rct.GetCaliDataMapConst().GetChannel(sec._edgePads.second) << ")";
        os << std::endl;

        os << "  with " << sec._factors.size() << " factors." << std::endl;
        for (size_t j = 0; j < sec._factors.size(); j++)
        {
            // here j represents [channel index - 1], if channel index match current chs pair, print with bold and blue color
            if (j == (size_t)(chs.first - 1) || j == (size_t)(chs.second - 1))
                // os << "    \033[1mFactor " << j << ": " << sec._factors[j] << "\033[0m" << std::endl;
                os << "    \033[1;34mFactor for pad " << caliMap.GetPad(j + 1) << ": " << sec._factors[j] << "\033[0m" << std::endl;
            else
                os << "    Factor for pad " << caliMap.GetPad(j + 1) << ": " << sec._factors[j] << std::endl;
        }
        os << "  and pol1: y = " << sec._pol1->GetParameter(0) << " + " << sec._pol1->GetParameter(1) << " * x" << std::endl;
    }
    os << std::endl;

    os << "Pad normalization factors:" << std::endl;
    for (const auto &pair : rct.fmPadNormFactors)
        os << " Pad " << pair.first << ": " << pair.second << std::endl;

    return os;
}

PRAlgorithms::_linear_fit_data::_linear_fit_data(std::string config_file_path)
{
    ReadLinearPRFit(config_file_path, *this);
}
