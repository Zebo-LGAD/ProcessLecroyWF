
#include "VPRAlgorithm.h"

bool PRAlgorithms::ChannelPadMap::SetMapC2P(std::map<int, int> channelToPad)
{
    for (const auto &pair : channelToPad)
    {
        auto rtn = AddCPMap(pair.first, pair.second);
        if (!rtn)
        {
            fmChannelToPad.clear();
            fmPadToChannel.clear();
            fmPadToXYWidth.clear();
            return false;
        }
    }
    return true;
}

bool PRAlgorithms::ChannelPadMap::SetMapP2C(std::map<int, int> padToChannel)
{
    for (const auto &pair : padToChannel)
    {
        auto rtn = AddCPMap(pair.second, pair.first);
        if (!rtn)
        {
            fmChannelToPad.clear();
            fmPadToChannel.clear();
            fmPadToXYWidth.clear();
            return false;
        }
    }
    return true;
}

bool PRAlgorithms::ChannelPadMap::AddCPMap(int channel, int pad)
{
    int n_pads = GetNPads();
    return AddCPMap(channel, pad, n_pads, 0);
}

bool PRAlgorithms::ChannelPadMap::AddCPMap(int channel, int pad, int padcolumn, int padrow)
{
    if (fMapClosed)
    {
        std::cerr << "Map is closed for further insertion." << std::endl;
        return false;
    }

    if (channel <= 0 || pad <= 0)
        return false;

    // Search in both maps to avoid conflict
    auto itC2P = fmChannelToPad.find(channel);
    if (itC2P != fmChannelToPad.end())
    {
        std::cerr << "Channel " << channel << " already mapped to Pad " << itC2P->second << std::endl;
        return false;
    }
    auto itP2C = fmPadToChannel.find(pad);
    if (itP2C != fmPadToChannel.end())
    {
        std::cerr << "Pad " << pad << " already mapped to Channel " << itP2C->second << std::endl;
        return false;
    }

    fmChannelToPad[channel] = pad;
    fmPadToChannel[pad] = channel;
    InsertPadInMatrix(pad, padcolumn, padrow);

    return true;
}

bool PRAlgorithms::ChannelPadMap::AddPCMap(int pad, int channel)
{
    return AddCPMap(channel, pad);
}

bool PRAlgorithms::ChannelPadMap::AddPCMap(int pad, int channel, int padcolumn, int padrow)
{
    return AddCPMap(channel, pad, padcolumn, padrow);
}

int PRAlgorithms::ChannelPadMap::GetPad(int channel) const
{
    auto it = fmChannelToPad.find(channel);
    if (it != fmChannelToPad.end())
        return it->second;
    return -1; // Not found
}

int PRAlgorithms::ChannelPadMap::GetChannel(int pad) const
{
    auto it = fmPadToChannel.find(pad);
    if (it != fmPadToChannel.end())
        return it->second;
    return -1; // Not found
}

int PRAlgorithms::ChannelPadMap::GetChannel(int pad, int &padcolumn, int &padrow) const
{
    auto channel = GetChannel(pad);
    if (channel != -1)
        GetPadXY(pad, padcolumn, padrow);
    return channel;
}

int PRAlgorithms::ChannelPadMap::GetPad(int channel, int &padcolumn, int &padrow) const
{
    auto pad = GetPad(channel);
    if (pad != -1)
        GetPadXY(pad, padcolumn, padrow);
    return pad;
}

bool PRAlgorithms::ChannelPadMap::GetPadXY(int pad, int &column, int &row) const
{
    for (int r = 0; r < fPadMatrix.GetNrows(); r++)
        for (int c = 0; c < fPadMatrix.GetNcols(); c++)
            if (fPadMatrix(r, c) == (float)pad)
            {
                column = c;
                row = r;
                return true;
            }
    return false; // Not found
}

bool PRAlgorithms::ChannelPadMap::InsertPadInMatrix(int pad, int column, int row)
{
    if (fMapClosed)
    {
        std::cerr << "Map is closed for further insertion." << std::endl;
        return false;
    }
    // Pad must be positive
    if (pad <= 0)
        return false;
    // Search if pad already exists
    for (int r = 0; r < fPadMatrix.GetNrows(); r++)
        for (int c = 0; c < fPadMatrix.GetNcols(); c++)
            if (fPadMatrix(r, c) == (float)pad)
            {
                std::cout << "Pad " << pad << " already exists in pad matrix at (" << r << ", " << c << ")" << std::endl;
                return false;
            }

    if (fPadMatrix.GetNrows() < row + 1 || fPadMatrix.GetNcols() < column + 1)
    {
        int nrows = std::max(fPadMatrix.GetNrows(), row + 1);
        int ncols = std::max(fPadMatrix.GetNcols(), column + 1);
        fPadMatrix.ResizeTo(nrows, ncols);
    }
    fPadMatrix(row, column) = pad;

    return false;
}

void PRAlgorithms::ChannelPadMap::Clear()
{
    fmChannelToPad.clear();
    fmPadToChannel.clear();
    fmPadToXYWidth.clear();
    fvValidPadX.clear();
    fvValidPadY.clear();
    fmPadXYToNeighborPadsX.clear();
    fmPadXYToNeighborPadsY.clear();
    fmNeighborPadsToChsX.clear();
    fmNeighborPadsToChsY.clear();
    fPadMatrix.ResizeTo(0, 0);
    fMapClosed = false;
}

bool PRAlgorithms::ChannelPadMap::GetLeftPad(int pad, int &leftpad) const
{
    int column, row;
    if (!GetPadXY(pad, column, row))
        return false;
    if (column == 0)
        return false; // No left pad
    leftpad = (int)fPadMatrix(row, column - 1);
    if (leftpad <= 0)
        return false;
    return true;
}

bool PRAlgorithms::ChannelPadMap::GetRightPad(int pad, int &rightpad) const
{
    int column, row;
    if (!GetPadXY(pad, column, row))
        return false;
    if (column >= fPadMatrix.GetNcols() - 1)
        return false; // No right pad
    rightpad = (int)fPadMatrix(row, column + 1);
    if (rightpad <= 0)
        return false;
    return true;
}

bool PRAlgorithms::ChannelPadMap::GetUpPad(int pad, int &uppad) const
{
    int column, row;
    if (!GetPadXY(pad, column, row))
        return false;
    if (row == 0)
        return false; // No up pad
    uppad = (int)fPadMatrix(row - 1, column);
    if (uppad <= 0)
        return false;
    return true;
}

bool PRAlgorithms::ChannelPadMap::GetDownPad(int pad, int &downpad) const
{
    int column, row;
    if (!GetPadXY(pad, column, row))
        return false;
    if (row >= fPadMatrix.GetNrows() - 1)
        return false; // No down pad
    downpad = (int)fPadMatrix(row + 1, column);
    if (downpad <= 0)
        return false;
    return true;
}

bool JudgeWidthValid(double columnwidth, double rowwidth)
{
    if (columnwidth == 0 && rowwidth == 0)
        return false;
    if (columnwidth < 0 || rowwidth < 0)
        return false;

    return true;
}

bool PRAlgorithms::ChannelPadMap::SetUniformPadSize(double padcolumnwidth, double padrowwidth)
{
    if (!JudgeWidthValid(padcolumnwidth, padrowwidth))
        return false;

    for (const auto &pair : fmPadToChannel)
    {
        int pad = pair.first;
        fmPadToXYWidth[pad] = _pad_xy_width(padcolumnwidth, padrowwidth);
    }
    fIsUniformPadSize = true;
    fUniformPadXYWidth = _pad_xy_width(padcolumnwidth, padrowwidth);
    return true;
}

bool PRAlgorithms::ChannelPadMap::SetPadSize(int pad, double padcolumnwidth, double padrowwidth)
{
    if (!JudgeWidthValid(padcolumnwidth, padrowwidth))
        return false;

    auto it = fmPadToChannel.find(pad);
    if (it == fmPadToChannel.end())
        return false;

    if (fIsUniformPadSize)
    {
        if (padcolumnwidth == fUniformPadXYWidth.first && padrowwidth == fUniformPadXYWidth.second)
        {
            std::cout << "Warning: Pad size is the same as uniform pad size, ignoring resetting." << std::endl;
            return false; // Same as uniform size
        }
        fIsUniformPadSize = false;
        fUniformPadXYWidth = _pad_xy_width(0, 0);
    }

    fmPadToXYWidth[pad] = _pad_xy_width(padcolumnwidth, padrowwidth);
    return true;
}

PRAlgorithms::_pad_xy_width PRAlgorithms::ChannelPadMap::GetPadSize(int pad) const
{
    if (fIsUniformPadSize)
        return fUniformPadXYWidth;

    auto it = fmPadToXYWidth.find(pad);
    if (it != fmPadToXYWidth.end())
        return it->second;
    return _pad_xy_width(0, 0); // Not found
}

void PRAlgorithms::ChannelPadMap::InitializeMap()
{
    // First find all neighbor pads from pad matrix
    // Scan for right neigbors
    int nrows = fPadMatrix.GetNrows();
    int ncols = fPadMatrix.GetNcols();

    for (int r = 0; r < fPadMatrix.GetNrows(); r++)
    {
        for (int c = 0; c < fPadMatrix.GetNcols(); c++)
        {
            int pad = (int)fPadMatrix(r, c);
            if (pad <= 0)
                continue;

            // Right neighbor
            if (c < ncols - 1)
            {
                // Right neighbor
                int rightpad = (int)fPadMatrix(r, c + 1);
                if (rightpad > 0)
                {
                    fvValidPadX.push_back(_pad_xy(c, r));
                    _neighbor_pads npads(pad, rightpad);
                    fmPadXYToNeighborPadsX[_pad_xy(c, r)] = npads;
                    _neighbor_chs nchs(GetChannel(pad), GetChannel(rightpad));
                    fmNeighborPadsToChsX[npads] = nchs;
                }
            }

            // Down neighbor
            if (r < nrows - 1)
            {
                int downpad = (int)fPadMatrix(r + 1, c);
                if (downpad > 0)
                {
                    fvValidPadY.push_back(_pad_xy(c, r));
                    _neighbor_pads npads(pad, downpad);
                    fmPadXYToNeighborPadsY[_pad_xy(c, r)] = npads;
                    _neighbor_chs nchs(GetChannel(pad), GetChannel(downpad));
                    fmNeighborPadsToChsY[npads] = nchs;
                }
            }
        }
    }

    // Set map as closed
    fMapClosed = true;
}

bool PRAlgorithms::JudgeSubMatrix(const TMatrix &submatrix, const TMatrix &matrix)
{
    if (submatrix.GetNrows() > matrix.GetNrows() || submatrix.GetNcols() > matrix.GetNcols())
        return false;

    for (int startRow = 0; startRow <= matrix.GetNrows() - submatrix.GetNrows(); startRow++)
    {
        for (int startCol = 0; startCol <= matrix.GetNcols() - submatrix.GetNcols(); startCol++)
        {
            bool match = true;
            for (int r = 0; r < submatrix.GetNrows(); r++)
            {
                for (int c = 0; c < submatrix.GetNcols(); c++)
                {
                    if (submatrix(r, c) != matrix(startRow + r, startCol + c))
                    {
                        match = false;
                        break;
                    }
                }
                if (!match)
                    break;
            }
            if (match)
                return true;
        }
    }
    return false;
}

bool PRAlgorithms::operator<=(const ChannelPadMap &lhs, const ChannelPadMap &rhs)
{
    // Check whether rhs contains all pads in lhs
    for (const auto &pair : lhs.fmPadToChannel)
    {
        auto it = rhs.fmPadToChannel.find(pair.first);
        if (it == rhs.fmPadToChannel.end())
            return false;
    }

    // Check wheter lhs is submatrix of rhs's pad matrix
    auto isSub = JudgeSubMatrix(lhs.fPadMatrix, rhs.fPadMatrix);
    if (isSub)
        return true;

    return false;
}

bool PRAlgorithms::operator>=(const ChannelPadMap &lhs, const ChannelPadMap &rhs)
{
    if (rhs <= lhs)
        return true;

    return false;
}

bool PRAlgorithms::operator==(const ChannelPadMap &lhs, const ChannelPadMap &rhs)
{
    return (lhs <= rhs) && (lhs >= rhs);
}

std::ostream &PRAlgorithms::operator<<(std::ostream &os, const ChannelPadMap &map)
{
    if (map.fMapClosed)
        os << "ChannelPadMap is closed for insertion." << std::endl;
    else
        os << "ChannelPadMap is open for insertion." << std::endl;

    // Print Channel to Pad Map information
    os << "Channel to Pad Map:" << std::endl;
    for (const auto &pair : map.fmChannelToPad)
        os << " Channel " << pair.first << " -> Pad " << pair.second << std::endl;

    // Print Pad to Channel Map information
    os << "Pad to Channel Map:" << std::endl;
    for (const auto &pair : map.fmPadToChannel)
        os << " Pad " << pair.first << " -> Channel " << pair.second << std::endl;

    // Print pad matrix
    os << "Pad Matrix (" << map.fPadMatrix.GetNrows() << " x " << map.fPadMatrix.GetNcols() << "):" << std::endl;
    for (int r = 0; r < map.fPadMatrix.GetNrows(); r++)
    {
        for (int c = 0; c < map.fPadMatrix.GetNcols(); c++)
        {
            int pad = (int)map.fPadMatrix(r, c);
            int channel = map.GetChannel(pad);
            os << "[Pad " << pad << ", Ch " << channel << "] ";
        }
        os << std::endl;
    }

    // Print neighbor pads information
    if (map.fmPadXYToNeighborPadsX.empty() && map.fmPadXYToNeighborPadsY.empty())
    {
        os << "No neighbor pads information available." << std::endl;
        return os;
    }
    if (!map.fmPadXYToNeighborPadsX.empty())
    {
        os << "Neighbor Pads Map (X direction):" << std::endl;
        for (const auto &pair : map.fmPadXYToNeighborPadsX)
        {
            os << " PadXY (" << pair.first.first << ", " << pair.first.second << ") -> Neighbor Pads (" << pair.second.first << ", " << pair.second.second << ")" << std::endl;
            // Print corresponding channels
            auto it = map.fmNeighborPadsToChsX.find(pair.second);
            if (it != map.fmNeighborPadsToChsX.end())
                os << "   Corresponding Channels: (" << it->second.first << ", " << it->second.second << ")" << std::endl;
        }
    }

    if (!map.fmPadXYToNeighborPadsY.empty())
    {
        os << "Neighbor Pads Map (Y direction):" << std::endl;
        for (const auto &pair : map.fmPadXYToNeighborPadsY)
        {
            os << " PadXY (" << pair.first.first << ", " << pair.first.second << ") -> Neighbor Pads (" << pair.second.first << ", " << pair.second.second << ")" << std::endl;
            // Print corresponding channels
            auto it = map.fmNeighborPadsToChsY.find(pair.second);
            if (it != map.fmNeighborPadsToChsY.end())
                os << "   Corresponding Channels: (" << it->second.first << ", " << it->second.second << ")" << std::endl;
        }
    }

    return os;
}

std::ostream &PRAlgorithms::operator<<(std::ostream &os, const VPRAlgorithm &algo)
{
    os << "VPRAlgorithm with reconstruction map of " << algo.fmCPReconstruct.GetNPads() << " pads." << std::endl;
    os << algo.fmCPReconstruct;

    return os;
}

bool PRAlgorithms::VPRAlgorithm::SetUniformPadSize(double padcolumnwidth, double padrowwidth)
{
    return fmCPReconstruct.SetUniformPadSize(padcolumnwidth, padrowwidth);
}
