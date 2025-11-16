#include <vector>
#include <map>
#include <iostream>

#include <TMatrix.h>

namespace PRAlgorithms
{
    typedef std::pair<int, int> _neighbor_pads;               // left, right pads
    typedef std::pair<double, double> _neighbor_pads_signals; // two signals corresponding to _neighbor_pads
    typedef std::pair<int, int> _neighbor_chs;                // two channels corresponding to _neighbor_pads
    typedef std::pair<int, int> _neighbor_chs_signals;        // two signals corresponding to _neighbor_chs
    typedef std::pair<int, int> _pad_xy;                      // pad column and row

    // Used for sensor oscilloscope test, mapping scope channels to sensor pad.
    // Pad is important here, but channel may change in different oscilloscope setups.
    // While judging two maps, only pad to channel map is considered.
    class ChannelPadMap
    {
    public:
        ChannelPadMap() = default;

        bool AddCPMap(int channel, int pad);
        bool AddCPMap(int channel, int pad, int padcolumn, int padrow = 0);
        bool AddPCMap(int pad, int channel);
        bool AddPCMap(int pad, int channel, int padcolumn, int padrow = 0);
        bool SetMapC2P(std::map<int, int> channelToPad);
        bool SetMapP2C(std::map<int, int> padToChannel);

        int GetPad(int channel) const;
        int GetChannel(int pad) const;
        int GetChannel(int pad, int &padcolumn, int &padrow) const;
        int GetPad(int channel, int &padcolumn, int &padrow) const;
        bool GetPadXY(int pad, int &column, int &row) const;
        int GetNPads() const { return fmPadToChannel.size(); };

        const std::map<int, int> &GetMapC2P() const { return fmChannelToPad; };
        const std::map<int, int> &GetMapP2C() const { return fmPadToChannel; };

        bool InsertPadInMatrix(int pad, int column, int row = 0);
        void Clear();
        /// @brief Close the map for further insertion, and generate pad column-row neighbor pads map, and neighbor pads to channels maps
        void InitializeMap();

        bool GetLeftPad(int pad, int &leftpad) const;
        bool GetRightPad(int pad, int &rightpad) const;
        bool GetUpPad(int pad, int &uppad) const;
        bool GetDownPad(int pad, int &downpad) const;

        friend std::ostream &operator<<(std::ostream &os, const ChannelPadMap &map);
        friend bool operator<=(const ChannelPadMap &lhs, const ChannelPadMap &rhs);

        void Print(std::ostream &os = std::cout) { os << (*this); };

        const TMatrix &GetPadMatrix() const { return fPadMatrix; };

        /// @brief Get the pad (column,row) to neighbor pads map in X and Y directions, X: (current, right), Y: (current, down)
        /// @return
        const std::vector<_pad_xy> &GetValidPadX() const { return fvValidPadX; };
        const std::vector<_pad_xy> &GetValidPadY() const { return fvValidPadY; };
        const std::map<_pad_xy, _neighbor_pads> &GetPadXYToNeighborPadsX() const { return fmPadXYToNeighborPadsX; };
        const std::map<_pad_xy, _neighbor_pads> &GetPadXYToNeighborPadsY() const { return fmPadXYToNeighborPadsY; };
        const std::map<_neighbor_pads, _neighbor_chs> &GetNeighborPadsToChsX() const { return fmNeighborPadsToChsX; };
        const std::map<_neighbor_pads, _neighbor_chs> &GetNeighborPadsToChsY() const { return fmNeighborPadsToChsY; };

    private:
        std::map<int, int> fmChannelToPad;
        std::map<int, int> fmPadToChannel;
        TMatrix fPadMatrix; // Pad arrangement matrix for 2D sensor (or 1D line)

        bool fMapClosed = false;                                      // Whether the map is closed for further insertion
        std::vector<_pad_xy> fvValidPadX;                             // Pads has right neighbor
        std::vector<_pad_xy> fvValidPadY;                             // Pads has down neighbor
        std::map<_pad_xy, _neighbor_pads> fmPadXYToNeighborPadsX;     // Map from pad (column,row) to its neighbor pads <current, right>
        std::map<_pad_xy, _neighbor_pads> fmPadXYToNeighborPadsY;     // Map from pad (column,row) to its neighbor pads <current, down>
        std::map<_neighbor_pads, _neighbor_chs> fmNeighborPadsToChsX; // Neighbor pads to channels map in X(row) direction
        std::map<_neighbor_pads, _neighbor_chs> fmNeighborPadsToChsY; // Neighbor pads to channels map in Y(column) direction
    };

    bool JudgeSubMatrix(const TMatrix &submatrix, const TMatrix &matrix);


    
    /// @brief pads in lhs is subset of pads in rhs, channels in lhs is not considered, and lhs's pad matrix is submatrix of rhs's pad matrix, so as to ensure the relative positions among pads in lhs is same as in rhs
    bool operator<=(const ChannelPadMap &lhs, const ChannelPadMap &rhs);
    bool operator>=(const ChannelPadMap &lhs, const ChannelPadMap &rhs);
    bool operator==(const ChannelPadMap &lhs, const ChannelPadMap &rhs);

    std::ostream &operator<<(std::ostream &os, const ChannelPadMap &map);

    class VPRAlgorithm
    {
    public:
        virtual ~VPRAlgorithm() {}

        // Input signals from oscilloscope channels
        virtual double ReconstructX(const std::vector<double> &signals) const = 0;
        virtual double ReconstructY(const std::vector<double> &signals) const = 0;

        /// @brief Get the position reconstruction data map, which maps channel numbers to pad indices
        /// @return
        ChannelPadMap &GetPRDataMap() { return fmCPReconstruct; };
        const ChannelPadMap &GetPRDataMapConst() const { return fmCPReconstruct; };

        friend std::ostream &operator<<(std::ostream &os, const VPRAlgorithm &algo);
        void Print(std::ostream &os = std::cout) { os << (*this); };

    private:
        ChannelPadMap fmCPReconstruct; // This map is used for reconstruction data
    };
    std::ostream &operator<<(std::ostream &os, const VPRAlgorithm &algo);

}
