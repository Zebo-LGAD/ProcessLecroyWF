

#ifndef LCPARSER_H
#define LCPARSER_H

#include "TObject.h"

class TGraph;

// 用于字节序转换的工具函数
namespace Endian
{
    // 检查系统是否为小端序
    bool isLittleEndian();

    // 交换16位整数字节序
    int16_t swap16(int16_t value);

    // 交换32位整数字节序
    int32_t swap32(int32_t value);

    // 交换64位整数字节序
    int64_t swap64(int64_t value);

    // 交换32位浮点数字节序
    float swapFloat(float value);
    // 交换64位浮点数字节序
    double swapDouble(double value);
}

class ScopeData : public TObject
{
private:
    std::string fPath = "";
    std::vector<double> x{};
    std::vector<double> y{};
    std::vector<std::vector<double>> yAll{}; // 用于存储多个通道的数据

    // 波形文件元数据
    std::string endianness = "";
    std::string templateName = "";
    std::string instrumentName = "";
    int instrumentNumber = 0;
    std::string waveSource = "";
    int waveArrayCount = 0;
    std::string verticalCoupling = "";
    std::string bandwidthLimit = "";
    std::string recordType = "";
    std::string processingDone = "";
    std::string timeBase = "";
    std::string triggerTime = "";

    float verticalGain = 0;
    float verticalOffset = 0;
    double horizInterval = 0;
    double horizOffset = 0;
    int nominalBits = 0;

    size_t posWAVEDESC = 0;
    int commOrder = 0;
    int commType = 0;
    int waveDescriptor = 0;
    int userText = 0;
    int trigTimeArray = 0;
    int waveArray1 = 0;

    // 常量列表
    static const std::vector<std::string> waveSourceList;       //!
    static const std::vector<std::string> verticalCouplingList; //!
    static const std::vector<std::string> bandwidthLimitList;   //!
    static const std::vector<std::string> recordTypeList;       //!
    static const std::vector<std::string> processingList;       //!

    // 添加一个标志来表示是否需要交换字节序
    bool needSwap = 0; //!

public:
    enum errorCodes
    {
        kSUCCESS = 0,
        kFILE_NOT_FOUND = -1,
        kREAD_ERROR = -2,
        kINVALID_FORMAT = -3,
        kZERO_ARRAY_COUNT = -4,
    };

    ScopeData(const std::string &path = "", bool parseAll = false, int sparse = -1);
    void Clean();
    errorCodes InitData(const std::string &path = "", bool parseAll = false, int sparse = -1);

    // 获取X数据
    const std::vector<double> &getX() const { return x; }

    // 获取Y数据
    const std::vector<double> &getY() const { return y; }

    // 获取所有通道的Y数据
    const std::vector<std::vector<double>> &getAllY() const { return yAll; }

    // 获取元数据
    const std::string &getInstrumentName() const { return instrumentName; }
    const std::string &getWaveSource() const { return waveSource; }
    int getWaveArrayCount() const { return waveArrayCount; }
    const std::string &getVerticalCoupling() const { return verticalCoupling; }
    const std::string &getBandwidthLimit() const { return bandwidthLimit; }
    const std::string &getRecordType() const { return recordType; }
    const std::string &getProcessingDone() const { return processingDone; }
    const std::string &getTimeBase() const { return timeBase; }
    const std::string &getTriggerTime() const { return triggerTime; }

    // 重载输出运算符
    friend std::ostream &operator<<(std::ostream &os, const ScopeData &data);

    TGraph *createGraph() const;

private:
    // 解析单个文件
    errorCodes parseFile(const std::string &path, std::vector<double> &x, std::vector<double> &y, int sparse);

    // 解析数据
    errorCodes parseData(const std::vector<char> &data, std::vector<double> &x, std::vector<double> &y, int sparse);

    // 从数据中解析16位整数
    int16_t parseInt16(const std::vector<char> &data, size_t pos) const;

    // 从数据中解析32位整数
    int32_t parseInt32(const std::vector<char> &data, size_t pos) const;

    // 从数据中解析浮点数
    float parseFloat(const std::vector<char> &data, size_t pos) const;

    // 从数据中解析双精度浮点数
    double parseDouble(const std::vector<char> &data, size_t pos) const;

    // 从数据中解析字符串
    std::string parseString(const std::vector<char> &data, size_t pos, size_t length) const;

    // 从数据中解析字节
    uint8_t parseByte(const std::vector<char> &data, size_t pos) const;

    // 从数据中解析字(16位)
    int16_t parseWord(const std::vector<char> &data, size_t pos) const;

    // 解析时间戳
    std::string parseTimeStamp(const std::vector<char> &data, size_t pos) const;

    // 解析时间基
    std::string parseTimeBase(const std::vector<char> &data, size_t pos) const;

    // 模拟Python的glob功能
    void globFiles(const std::string &pattern, std::vector<std::string> &files);

    ClassDef(ScopeData, 1);
};

#endif
