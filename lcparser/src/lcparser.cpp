// Read data from lecroy oscilloscope waveform .trc files (converted from lecroyparser)

#include "lcparser.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <sstream>
#include <glob.h>
#include <stdexcept>

#include "TGraph.h"

ClassImp(ScopeData);

namespace Endian
{
    // 检查系统是否为小端序
    bool isLittleEndian()
    {
        int num = 1;
        return (*(char *)&num == 1);
    }

    // 交换16位整数字节序
    int16_t swap16(int16_t value)
    {
        return (value >> 8) | (value << 8);
    }

    // 交换32位整数字节序
    int32_t swap32(int32_t value)
    {
        return ((value >> 24) & 0xff) |
               ((value << 8) & 0xff0000) |
               ((value >> 8) & 0xff00) |
               ((value << 24) & 0xff000000);
    }

    // 交换64位整数字节序
    int64_t swap64(int64_t value)
    {
        return ((value >> 56) & 0xff) |
               ((value >> 40) & 0xff00) |
               ((value >> 24) & 0xff0000) |
               ((value >> 8) & 0xff000000) |
               ((value << 8) & 0xff00000000) |
               ((value << 24) & 0xff0000000000) |
               ((value << 40) & 0xff000000000000) |
               ((value << 56) & 0xff00000000000000);
    }

    // 交换32位浮点数字节序
    float swapFloat(float value)
    {
        uint32_t temp = swap32(*(uint32_t *)&value);
        return *(float *)&temp;
    }

    // 交换64位浮点数字节序
    double swapDouble(double value)
    {
        uint64_t temp = swap64(*(uint64_t *)&value);
        return *(double *)&temp;
    }
}

struct TimeStamp
{
    double second;
    UChar_t minute;
    UChar_t hour;
    UChar_t day;
    UChar_t month;
    UChar_t year;
};

const std::vector<std::string> ScopeData::waveSourceList = {
    "Channel 1", "Channel 2", "Channel 3", "Channel 4", "Unknown"};

const std::vector<std::string> ScopeData::verticalCouplingList = {
    "DC50", "GND", "DC1M", "GND", "AC1M"};

const std::vector<std::string> ScopeData::bandwidthLimitList = {
    "off", "on"};

const std::vector<std::string> ScopeData::recordTypeList = {
    "single_sweep", "interleaved", "histogram", "graph",
    "filter_coefficient", "complex", "extrema", "sequence_obsolete",
    "centered_RIS", "peak_detect"};

const std::vector<std::string> ScopeData::processingList = {
    "No Processing", "FIR Filter", "interpolated", "sparsed",
    "autoscaled", "no_resulst", "rolling", "cumulative"};

ScopeData::ScopeData(const std::string &path, bool parseAll, int sparse)
{
    if (!path.empty())
        InitData(path, parseAll, sparse);
}

void ScopeData::Clean()
{
    fPath = "";
    x.clear();
    y.clear();
    yAll.clear();

    endianness = "";
    templateName = "";
    instrumentName = "";
    instrumentNumber = 0;
    waveSource = "";
    waveArrayCount = 0;
    verticalCoupling = "";
    bandwidthLimit = "";
    recordType = "";
    processingDone = "";
    timeBase = "";
    triggerTime = "";

    verticalGain = 0;
    verticalOffset = 0;
    horizInterval = 0;
    horizOffset = 0;
    nominalBits = 0;

    posWAVEDESC = 0;
    commOrder = 0;
    commType = 0;
    waveDescriptor = 0;
    userText = 0;
    trigTimeArray = 0;
    waveArray1 = 0;

    needSwap = 0;
}

ScopeData::errorCodes ScopeData::InitData(const std::string &path, bool parseAll, int sparse)
{
    Clean();

    if (!path.empty())
    {
        this->fPath = path;

        if (parseAll)
        {
            // 查找所有匹配的文件
            std::string basePath = path.substr(0, path.find_last_of("/\\"));
            std::string core_filename = path.substr(path.find_last_of("/\\") + 1);

            // 提取文件名中的数字部分
            size_t digit_pos = core_filename.find_first_of("0123456789");
            if (digit_pos != std::string::npos)
                core_filename = core_filename.substr(digit_pos);

            std::vector<std::string> files;
            globFiles(basePath + "/C*" + core_filename, files);
            std::sort(files.begin(), files.end());

            // 解析所有文件
            for (size_t i = 0; i < files.size(); i++)
            {
                std::vector<double> xTemp, yTemp;
                auto rtn = parseFile(files[i], xTemp, yTemp, sparse);
                if (rtn != kSUCCESS)
                    return rtn;

                if (i == 0)
                    x = xTemp;
                yAll.push_back(yTemp);
            }
        }
        else
            return parseFile(path, x, y, sparse);
    }

    return kSUCCESS;
}

ScopeData::errorCodes ScopeData::parseFile(const std::string &path, std::vector<double> &x, std::vector<double> &y, int sparse)
{
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file.is_open())
    {
        std::runtime_error("Cannot open file: " + path);
        return kFILE_NOT_FOUND;
    }

    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size))
    {
        std::runtime_error("Cannot read file: " + path);
        return kREAD_ERROR;
    }

    file.close();

    return parseData(buffer, x, y, sparse);
}

TGraph *ScopeData::createGraph() const
{
    if (!x.empty() && !y.empty() && x.size() == y.size())
    {
        // return new TGraph(x.size(), x.data(), y.data());

        std::vector<double> t_in_ns(x.size()), y_in_mV(x.size());
        for (size_t i = 0; i < x.size(); i++)
        {
            t_in_ns[i] = (x[i] * 1e9); // Convert to ns
            y_in_mV[i] = (y[i] * 1e3); // Convert to mV
        }

        return new TGraph(t_in_ns.size(), t_in_ns.data(), y_in_mV.data());
    }
    return nullptr;
}

ScopeData::errorCodes ScopeData::parseData(const std::vector<char> &data, std::vector<double> &x, std::vector<double> &y, int sparse)
{
    // 查找WAVEDESC位置
    std::string header(data.begin(), data.begin() + std::min(static_cast<size_t>(50), data.size()));
    posWAVEDESC = header.find("WAVEDESC");
    if (posWAVEDESC == std::string::npos)
    {
        std::runtime_error("WAVEDESC not found in file");
        return kINVALID_FORMAT;
    }

    // 解析通信顺序（字节序）
    needSwap = 0;
    commOrder = parseInt16(data, 34);
    endianness = (commOrder == 0) ? ">" : "<";
    // needSwap = (Endian::isLittleEndian() != (commOrder == 1));
    needSwap = !(Endian::isLittleEndian() == (endianness == "<")); // If both are little or big endian, no swap needed

    // 解析其他元数据
    templateName = parseString(data, 16, 16);
    commType = parseInt16(data, 32);
    waveDescriptor = parseInt32(data, 36);
    userText = parseInt32(data, 40);
    trigTimeArray = parseInt32(data, 48);
    waveArray1 = parseInt32(data, 60);
    instrumentName = parseString(data, 76, 16);
    instrumentNumber = parseInt32(data, 92);
    waveArrayCount = parseInt32(data, 116);
    verticalGain = parseFloat(data, 156);
    verticalOffset = parseFloat(data, 160);
    nominalBits = parseInt16(data, 172);
    horizInterval = parseFloat(data, 176);
    horizOffset = parseDouble(data, 180);

    // 解析时间戳
    triggerTime = parseTimeStamp(data, 296);

    // 解析枚举类型
    int recordTypeIdx = parseInt16(data, 316);
    recordType = (recordTypeIdx >= 0 && recordTypeIdx < static_cast<int>(recordTypeList.size()))
                     ? recordTypeList[recordTypeIdx]
                     : "Unknown";

    int processingDoneIdx = parseInt16(data, 318);
    processingDone = (processingDoneIdx >= 0 && processingDoneIdx < static_cast<int>(processingList.size()))
                         ? processingList[processingDoneIdx]
                         : "Unknown";

    timeBase = parseTimeBase(data, 324);

    int verticalCouplingIdx = parseInt16(data, 326);
    verticalCoupling = (verticalCouplingIdx >= 0 && verticalCouplingIdx < static_cast<int>(verticalCouplingList.size()))
                           ? verticalCouplingList[verticalCouplingIdx]
                           : "Unknown";

    int bandwidthLimitIdx = parseInt16(data, 334);
    bandwidthLimit = (bandwidthLimitIdx >= 0 && bandwidthLimitIdx < static_cast<int>(bandwidthLimitList.size()))
                         ? bandwidthLimitList[bandwidthLimitIdx]
                         : "Unknown";

    int waveSourceIdx = parseInt16(data, 344);
    waveSource = (waveSourceIdx >= 0 && waveSourceIdx < static_cast<int>(waveSourceList.size()))
                     ? waveSourceList[waveSourceIdx]
                     : "Unknown";

    // 解析波形数据
    size_t start = posWAVEDESC + waveDescriptor + userText + trigTimeArray;

    if (start + waveArray1 > data.size())
    {
        std::runtime_error("Wave data exceeds file size");
        return kINVALID_FORMAT;
    }

    if (waveArrayCount == 0)
    {
        std::runtime_error("WaveArrayCount is zero");
        return kZERO_ARRAY_COUNT;
    }

    if (commType == 0)
    { // 8位整数数据
        for (size_t i = start; i < start + waveArray1; i++)
        {
            int8_t value = static_cast<int8_t>(data[i]);
            y.push_back(verticalGain * value - verticalOffset);
        }
    }
    else
    { // 16位整数数据
        for (size_t i = start; i < start + waveArray1; i += 2)
        {
            int16_t value;
            memcpy(&value, &data[i], sizeof(int16_t));
            if (needSwap)
            {
                value = Endian::swap16(value);
            }
            y.push_back(verticalGain * value - verticalOffset);
        }
    }

    // 生成X数据
    for (int i = 0; i < static_cast<int>(y.size()); i++)
    {
        x.push_back(i * horizInterval + horizOffset);
    }

    // 稀疏采样
    if (sparse > 0 && static_cast<int>(x.size()) > sparse)
    {
        std::vector<double> xSparse, ySparse;
        int step = static_cast<int>(x.size()) / sparse;

        for (int i = 0; i < sparse; i++)
        {
            int idx = i * step;
            if (idx < static_cast<int>(x.size()))
            {
                xSparse.push_back(x[idx]);
                ySparse.push_back(y[idx]);
            }
        }

        x = xSparse;
        y = ySparse;
    }
    return kSUCCESS;
}

int16_t ScopeData::parseInt16(const std::vector<char> &data, size_t pos) const
{
    if (pos + 1 + posWAVEDESC >= data.size())
    {
        return 0;
    }

    int16_t value;
    memcpy(&value, &data[pos + posWAVEDESC], sizeof(int16_t));
    return needSwap ? Endian::swap16(value) : value;
}

int32_t ScopeData::parseInt32(const std::vector<char> &data, size_t pos) const
{
    if (pos + 3 + posWAVEDESC >= data.size())
    {
        return 0;
    }

    int32_t value;
    memcpy(&value, &data[pos + posWAVEDESC], sizeof(int32_t));
    return needSwap ? Endian::swap32(value) : value;
}

float ScopeData::parseFloat(const std::vector<char> &data, size_t pos) const
{
    if (pos + 3 + posWAVEDESC >= data.size())
    {
        return 0.0f;
    }

    float value;
    memcpy(&value, &data[pos + posWAVEDESC], sizeof(float));

    return needSwap ? Endian::swapFloat(value) : value;
}

double ScopeData::parseDouble(const std::vector<char> &data, size_t pos) const
{
    if (pos + 7 + posWAVEDESC >= data.size())
    {
        return 0.0;
    }

    double value;
    memcpy(&value, &data[pos + posWAVEDESC], sizeof(double));

    return needSwap ? Endian::swapDouble(value) : value;
}

std::string ScopeData::parseString(const std::vector<char> &data, size_t pos, size_t length) const
{
    if (pos + length + posWAVEDESC > data.size())
    {
        return "";
    }

    const char *strStart = &data[pos + posWAVEDESC];
    return std::string(strStart, strnlen(strStart, length));
}

uint8_t ScopeData::parseByte(const std::vector<char> &data, size_t pos) const
{
    if (pos + posWAVEDESC >= data.size())
    {
        return 0;
    }

    return static_cast<uint8_t>(data[pos + posWAVEDESC]);
}

int16_t ScopeData::parseWord(const std::vector<char> &data, size_t pos) const
{
    return parseInt16(data, pos);
}

std::string ScopeData::parseTimeStamp(const std::vector<char> &data, size_t pos) const
{
    double second = parseDouble(data, pos);
    uint8_t minute = parseByte(data, pos + 8);
    uint8_t hour = parseByte(data, pos + 9);
    uint8_t day = parseByte(data, pos + 10);
    uint8_t month = parseByte(data, pos + 11);
    int16_t year = parseWord(data, pos + 12);

    std::ostringstream oss;
    oss << year << "-"
        << std::setw(2) << std::setfill('0') << static_cast<int>(month) << "-"
        << std::setw(2) << std::setfill('0') << static_cast<int>(day) << " "
        << std::setw(2) << std::setfill('0') << static_cast<int>(hour) << ":"
        << std::setw(2) << std::setfill('0') << static_cast<int>(minute) << ":"
        << std::fixed << std::setprecision(2) << second;

    return oss.str();
}

std::string ScopeData::parseTimeBase(const std::vector<char> &data, size_t pos) const
{
    int timeBaseNumber = parseInt16(data, pos);

    if (timeBaseNumber < 48)
    {
        const std::vector<std::string> units = {"ps", "ps", "ps", "ps", "ps", "ps", "ps", "ps", "ps",
                                                "ns", "ns", "ns", "ns", "ns", "ns", "ns", "ns", "ns",
                                                "us", "us", "us", "us", "us", "us", "us", "us", "us",
                                                "ms", "ms", "ms", "ms", "ms", "ms", "ms", "ms", "ms",
                                                "s", "s", "s", "s", "s", "s", "s", "s", "s",
                                                "ks", "ks", "ks"};

        const std::vector<int> values = {1, 2, 5, 10, 20, 50, 100, 200, 500,
                                         1, 2, 5, 10, 20, 50, 100, 200, 500,
                                         1, 2, 5, 10, 20, 50, 100, 200, 500,
                                         1, 2, 5, 10, 20, 50, 100, 200, 500,
                                         1, 2, 5, 10, 20, 50, 100, 200, 500,
                                         1, 2, 5};

        if (timeBaseNumber < static_cast<int>(values.size()))
        {
            return std::to_string(values[timeBaseNumber]) + " " + units[timeBaseNumber] + "/div";
        }
    }
    else if (timeBaseNumber == 100)
    {
        return "EXTERNAL";
    }

    return "Unknown";
}

void ScopeData::globFiles(const std::string &pattern, std::vector<std::string> &files)
{
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value != 0)
    {
        globfree(&glob_result);
        return;
    }

    for (size_t i = 0; i < glob_result.gl_pathc; ++i)
    {
        files.push_back(std::string(glob_result.gl_pathv[i]));
    }

    globfree(&glob_result);
}

std::ostream &operator<<(std::ostream &os, const ScopeData &data)
{
    os << "Le Croy Scope Data" << std::endl;
    os << "Path: " << data.fPath << std::endl;
    os << "Endianness: " << data.endianness << std::endl;
    os << "Instrument: " << data.instrumentName << std::endl;
    os << "Instrument Number: " << data.instrumentNumber << std::endl;
    os << "Template Name: " << data.templateName << std::endl;
    os << "Channel: " << data.waveSource << std::endl;
    os << "WaveArrayCount: " << data.waveArrayCount << std::endl;
    os << "Vertical Coupling: " << data.verticalCoupling << std::endl;
    os << "Bandwidth Limit: " << data.bandwidthLimit << std::endl;
    os << "Record Type: " << data.recordType << std::endl;
    os << "Processing: " << data.processingDone << std::endl;
    os << "TimeBase: " << data.timeBase << std::endl;
    os << "TriggerTime: " << data.triggerTime << std::endl;
    return os;
}
