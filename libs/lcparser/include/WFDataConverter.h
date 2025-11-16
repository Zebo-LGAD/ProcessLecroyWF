#ifndef WFDataConverter_H
#define WFDataConverter_H
#include <string>
#include <vector>
#include <map>

#include "lcparser.h"

#include <TFile.h>
#include <iostream>
#include <TTree.h>

class TFile;
class TTree;
class ScopeData;

namespace WFDataProcessor
{

    typedef std::pair<double, double> _signal_range;
    struct _extract_config
    {
        _signal_range search_range{-10.0, 10.0}; // ns
        double threshold{20.0};                  // mV
        bool need_draw{false};
        std::string savePrefix{""};
    };

    /// @brief Invalid code enumeration for waveform analysis
    enum InvalidCode
    {
        VALID = 0,
        NO_WAVEFORM = 1,
        NO_T1_FOUND = 1<<1,
        NO_T1_10_FOUND = 1<<2,
        NO_T1_50_FOUND = 1<<3,
        NO_T1_90_FOUND = 1<<4,
        NO_T2_FOUND = 1<<5,
        NO_T2_10_FOUND = 1<<6,
        NO_T2_50_FOUND = 1<<7,
        NO_T2_90_FOUND = 1<<8,
        RANGE_ERROR = 1<<9,
    };

    /// @brief Structure to hold extracted waveform information
    struct _waveinfo
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
        uint32_t valid;           ///  Flag indicating if the waveform is valid (1) or not (0)
    };

    /// @brief Get interpolated time of arrival at given threshold, using linear interpolation
    /// @tparam T template type (e.g., float, double)
    /// @param tarray array of time samples
    /// @param aarray array of amplitude samples
    /// @param sample_point left point index for interpolation, i.e., aarray[sample_point] >= threshold > aarray[sample_point - 1]
    /// @param threshold threshold for interpolation, to calculate time of arrival
    /// @param pedestal pedestal value to subtract from amplitude
    /// @param result interpolated time of arrival (output)
    /// @return whether interpolation was successful, if sample_point is 0, return false
    template <typename T>
    bool interpolateTOA(const T *tarray, const T *aarray, int sample_point, T threshold, T pedestal, T &result);

    /// @brief check if a value is within a specified range
    /// @tparam T
    /// @param value specified value to check
    /// @param range1 range boundary 1
    /// @param range2 range boundary 2
    /// @return whether the value is within the range [min(range1, range2), max(range1, range2)]
    template <typename T>
    bool rangeCheck(T value, T range1, T range2);

    /// @brief Anohter version of interpolateTOA with enhanced range checking, much more robust
    /// @tparam T template type (e.g., float, double)
    /// @param tarray time array
    /// @param aarray amplitude array
    /// @param NSamples number of samples in the waveform
    /// @param sample_point_around_threshold closet point index around threshold, will search both neighboring points for interpolation
    /// @param threshold threshold for interpolation, to calculate time of arrival
    /// @param pedestal pedestal value to subtract from amplitude
    /// @param result interpolated time of arrival (output)
    /// @return
    template <typename T>
    bool interpolateTOA2(const T *tarray, const T *aarray, int NSamples, int sample_point_around_threshold, T threshold, T pedestal, T &result);

    /// @brief Generate the standard filename for LeCroy Scope waveform files
    /// @param folder Folder path that contains the waveform files
    /// @param channel which channel number
    /// @param idx_lecroy_wf which waveform index
    /// @param mid_name what is the middle name part
    /// @param ext extension of the file, default is .trc
    /// @return generated full file path string
    std::string GenerateScopeFileName(const std::string &folder, int channel, int idx_lecroy_wf, std::string mid_name = "--Trace--", std::string ext = ".trc");

    /// @brief Read all waveform files for specified channels from a folder, if a file is missing or cannot be parsed, mark the corresponding channel as having no data, and return as false, but will continue reading other channels.
    /// @param folder folder path that contains the waveform files
    /// @param idx_lecroy_wf which waveform index
    /// @param chDataMap [out] map of channel number to ScopeData pointer
    /// @param chDataHasDataMap [out] map of channel number to boolean indicating whether data was successfully read
    /// @param mid_name what is the middle name part
    /// @param ext extension of the file, default is .trc
    /// @return return true if all specified channels have data, false otherwise
    bool ReadAllWF(const std::string &folder, int idx_lecroy_wf, std::map<int, ScopeData *> &chDataMap, std::map<int, bool> &chDataHasDataMap, std::string mid_name = "--Trace--", std::string ext = ".trc");

    /// @brief Process a single waveform to extract features such as pedestal, amplitude, charge, and time of arrival.
    /// Charge in mV*ns, Amplitude in mV, Time in ns, signal was search between low_range_t0 and high_range_t0 for each channel.
    /// @param ch Channel number (1-based index)
    /// @param t Time array (in seconds)
    /// @param a Amplitude array (in volts)
    /// @param Nsamples Number of samples in the waveform
    /// @param threshold [in] Threshold for t1 and t2 calculation (in mV), default is 20.0 mV
    void processWave(const double *t, const double *a, int Nsamples, _waveinfo &winfo, _extract_config config);

    void GenerateBranchForWaveInfo(TTree *tree, const std::string &branchNamePrefix, WFDataProcessor::_waveinfo *winfo);
    void SetBranchAddressToWaveInfo(TTree *tree, const std::string &branchNamePrefix, WFDataProcessor::_waveinfo *winfo);

}

namespace WFDataProcessor
{
    /// @brief Base class for multi-channel waveform data reader and writer
    template <typename T>
    class VMultiIO
    {
    public:
        virtual ~VMultiIO();

        virtual bool AddChannel(int channel);
        virtual bool OpenFile(const std::string &filename);
        virtual void CloseFile();

        virtual TTree *GetTree() { return fTree; }
        virtual TFile *GetFile() { return fFile; }

        const std::map<int, T *> &GetChannelDataMap() const { return fmChData; };
        const std::map<int, bool> &GetChannelDataHasDataMap() const { return fmChDataHasData; };
        static VMultiIO<T> *&CurrentInstance();

    protected:
        VMultiIO(bool isRead) : fIsRead(isRead) { CurrentInstance() = this; };
        VMultiIO(bool isRead, std::vector<int> channelsToRead);
        VMultiIO(const VMultiIO &) = delete;
        VMultiIO &operator=(const VMultiIO &) = delete;

        std::map<int, T *> fmChData;         // channel -> waveform data pointer
        std::map<int, bool> fmChDataHasData; // channel -> whether data has been read

        TFile *fFile = nullptr;
        TTree *fTree = nullptr;
        virtual void ClearMap();

        virtual bool InitTree() = 0;
        bool fIsRead = false;
    };

    template <typename T>
    class VMultiChannelReader : public VMultiIO<T>
    {
    public:
        VMultiChannelReader() : VMultiIO<T>(true) {}
        VMultiChannelReader(std::vector<int> channelsToRead) : VMultiIO<T>(true, channelsToRead) {}
    };

    template <typename T>
    class VMultiChannelWriter : public VMultiIO<T>
    {
    public:
        VMultiChannelWriter() : VMultiIO<T>(false) {}
        VMultiChannelWriter(std::vector<int> channelsToRead) : VMultiIO<T>(false, channelsToRead) {};

        bool GetForceMatch() const { return fForceMatch; }
        void SetForceMatch(bool forceMatch) { fForceMatch = forceMatch; }

    protected:
        bool fForceMatch = true; // Whether to force match all channels when reading data
    };

}

namespace WFDataProcessor
{

    class WFtrc2ROOT;
    /// @brief Extract information from a ScopeData object
    class WFDataExtractor : public VMultiChannelWriter<_waveinfo>
    {
    public:
        WFDataExtractor() = default;
        WFDataExtractor(std::vector<int> channelsToRead);

        bool AddChannel(int channel) override;

        const std::map<int, _extract_config> &GetExtractConfig() const { return fmChExtractConfig; };
        bool GetExtractConfig(int channel, _extract_config &config) const;
        bool SetExtractConfig(int channel, _extract_config range);
        int SetExtractConfig(const std::map<int, _extract_config> &chRangeMap);

        bool ExtractFromScopeData(const std::map<int, ScopeData *> &chDataMap, const std::map<int, bool> &chDataHasDataMap);
        bool ExtractFromTRCFiles(const std::string &folder, int idx_file, std::string mid_name = "--Trace--", std::string ext = ".trc");
        bool ExtractFromWFtrc2ROOT(const WFDataProcessor::WFtrc2ROOT &convertor);

        bool TurnOnPlot(int channel, const std::string &savePrefix, bool bswitch = true);
        bool TurnOffPlot(int channel) { return TurnOnPlot(channel, "", false); };
        int TurnOnPlots(const std::string &savePrefix = "", bool bswitch = true);
        int TurnOffPlots() { return TurnOnPlots("", false); };

    private:
        bool InitTree() override;

        virtual void ClearMap() override;
        int fExtractedCounter = 0;

        std::map<int, _extract_config> fmChExtractConfig; // channel -> search range
    };

    class WFDataTreeReader : public VMultiChannelReader<_waveinfo>
    {
    public:
        WFDataTreeReader(std::vector<int> channelsToRead) : VMultiChannelReader<_waveinfo>(channelsToRead) {};

    private:
        bool InitTree() override;
    };

    /// @brief Convert LeCroy .trc waveform files into ROOT TTree format, supporting multiple channels and optional channel matching.
    class WFtrc2ROOT : public VMultiChannelWriter<ScopeData>
    {
    public:
        WFtrc2ROOT() = default;
        WFtrc2ROOT(std::vector<int> channelsToRead) : VMultiChannelWriter<ScopeData>(channelsToRead) {};

        // Read scope data from a folder for all channels, and fill into tree, but check if all channels have data when fForceMatch is true
        bool ReadAllAndFill(const std::string &folder, int idx_file, std::string mid_name = "--Trace--", std::string ext = ".trc");

    private:
        bool InitTree() override;
    };

    class WFROOTReader : public VMultiChannelReader<ScopeData *>
    {
    public:
        WFROOTReader(std::vector<int> channelsToRead) : VMultiChannelReader<ScopeData *>(channelsToRead) {};

    private:
        bool InitTree() override;
    };

    int ConvertAllTRC(std::string sWriteFile, std::string sDataFolder, int maxFiles);
    void ConvertAllTRC();
    std::string ExtractLastDir(const std::string &path);

}
#define gWFDataExtractor ((WFDataProcessor::WFDataExtractor *&)WFDataProcessor::WFDataExtractor::CurrentInstance())
#define gWFtrc2ROOT ((WFDataProcessor::WFtrc2ROOT *&)WFDataProcessor::WFtrc2ROOT::CurrentInstance())
#define gWFROOTReader ((WFDataProcessor::WFROOTReader *&)WFDataProcessor::WFROOTReader::CurrentInstance())
#define gWFDataTreeReader ((WFDataProcessor::WFDataTreeReader *&)WFDataProcessor::WFDataTreeReader::CurrentInstance())

namespace WFDataProcessor
{

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
        return true;
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

    template <typename T>
    inline void VMultiIO<T>::ClearMap()
    {
        fmChData.clear();
        fmChDataHasData.clear();
    }

    template <typename T>
    inline bool VMultiIO<T>::InitTree()
    {
        if (!fFile)
            return false;
        if (fTree)
        {
            std::cerr << "Tree already initialized!" << std::endl;
            return false;
        }
        return true;
    }

    template <typename T>
    inline VMultiIO<T> *&VMultiIO<T>::CurrentInstance()
    {
        // TODO: insert return statement here
        static VMultiIO<T> *instance = nullptr;
        return instance;
    }

    template <typename T>
    inline VMultiIO<T>::VMultiIO(bool isRead, std::vector<int> channelsToRead) : VMultiIO<T>(isRead)
    {
        for (int ch : channelsToRead)
            AddChannel(ch);
        CurrentInstance() = this;
    }
    template <typename T>
    inline VMultiIO<T>::~VMultiIO()
    {
        CloseFile();
        for (auto &pair : fmChData)
        {
            delete pair.second;
            pair.second = nullptr;
        }
        if (CurrentInstance() == this)
            CurrentInstance() = nullptr;
    }
    template <typename T>
    inline bool VMultiIO<T>::AddChannel(int channel)
    {
        if (fTree)
        {
            std::cerr << "Cannot add channel after the write tree has been initialized." << std::endl;
            return false;
        }
        if (channel < 0)
            return false;
        if (fmChData.find(channel) != fmChData.end())
            return true; // already exists

        T *chData = new T;
        fmChData[channel] = chData;
        fmChDataHasData[channel] = false;
        return true;
    }

    template <typename T>
    inline bool VMultiIO<T>::OpenFile(const std::string &filename)
    {
        if (fFile)
        {
            std::cerr << "File already opened. Close it before opening a new file." << std::endl;
            return false;
        }

        if (fIsRead)
            fFile = TFile::Open(filename.c_str());
        else
            fFile = TFile::Open(filename.c_str(), "RECREATE");

        if (!fFile || fFile->IsZombie())
            return false;

        return InitTree();
    }

    template <typename T>
    inline void VMultiIO<T>::CloseFile()
    {
        if (fFile)
        {
            if (fTree && fFile->IsWritable())
            {
                fFile->cd();
                fTree->Write();
            }

            fFile->Close();
            delete fFile;
            fFile = nullptr;
            fTree = nullptr;
        }
    }
}
#endif // WFDataConverter_H