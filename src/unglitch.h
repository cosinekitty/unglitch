/*
    unglitch.h  -  Don Cross  -  13 December 2017.
*/
#ifndef __DDC_UNGLITCH_H
#define __DDC_UNGLITCH_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <string>
#include <deque>
#include <vector>

#include "tinyxml2.h"

namespace unglitch
{
    class Chunk;
    class GlitchRemover;

    typedef std::vector<float> FloatVector;

    struct GapInfo
    {        
        size_t offset;      // post-processed sample offset
        size_t length;      // number of silent samples
        size_t front;       // original, pre-filtered sample offset

        GapInfo(size_t _offset, size_t _length, size_t _front)
            : offset(_offset)
            , length(_length)
            , front(_front)
            {}

        size_t OriginalCenter() const { return front + (length / 2); }
        size_t FilteredCenter() const { return offset + (length / 2); }
    };

    struct PreGapInfo
    {
        size_t front;       // index of first silent sample
        size_t length;      // number of silent samples

        PreGapInfo(size_t _front, size_t _length)
            : front(_front)
            , length(_length)
            {}

        size_t Center() const { return front + (length / 2); }
    };

    typedef std::vector<GapInfo> GapList;    
    typedef std::vector<PreGapInfo> PreGapList;    

    class GapFinder
    {
    private:
        GapList gaplist;
        size_t silenceFront;    // pre-filtered data offset (before removing glitches)
        size_t silenceOffset;   // post-filtered data offset (glitches removed)
        size_t silentSamples;
        size_t totalSamples;

    public:
        GapFinder()
        {
            Reset();
        }

        void Reset()
        {
            gaplist.clear();
            silenceFront = 0;
            silenceOffset = 0;
            silentSamples = 0;
            totalSamples = 0;
        }

        void Process(const Chunk &chunk, long programStartPosition);
        size_t TotalSamples() const { return totalSamples; }
        const GapList& SilentGaps() const { return gaplist; }
    };

    class PreScanGapFinder
    {
    private:
        PreGapList gaplist;
        size_t totalSamples;
        const size_t sampleDelay;
        double leftBias;
        double rightBias;
        size_t silentSamples;
        size_t silenceFront;

    public:
        PreScanGapFinder(size_t _sampleDelay)
            : totalSamples(0)
            , sampleDelay(_sampleDelay)
            , leftBias(0.0)
            , rightBias(0.0)
            , silentSamples(0)
            , silenceFront(0)
            {}
        
        void Process(const FloatVector& leftBlock, const FloatVector &rightBlock);
        const PreGapList& SilentGaps() const { return gaplist; }
    };

    inline float PeakValue(const FloatVector & buffer)
    {
        using namespace std;
        float peak = 0.0f;
        for (float data : buffer)
            peak = max(peak, abs(data));
        return peak;
    }

    class Error
    {
    private:
        const std::string message;

    public:
        Error(std::string _message) 
            : message(_message)
            {}

        const std::string& Message() const { return message; }
    };

    class WaveBlock
    {
    private:
        std::string filename;       // .au filename for this block
        long start;     // sample offset of the start of the block
        long length;    // number of samples in the block
        float ymin;    // minimum sample value
        float ymax;    // maximum sample value

    public:
        WaveBlock(std::string _filename, long _start, long _length, float _ymin, float _ymax)
            : filename(_filename)
            , start(_start)
            , length(_length)
            , ymin(_ymin)
            , ymax(_ymax)
            {}

        const std::string& Filename() const { return filename; }
        long Start() const { return start; }
        long Length() const { return length; }
        float Peak() const { return std::max(fabs(ymin), fabs(ymax)); }
    };

    class WaveTrack
    {
    private:
        int rate;           // sampling rate
        long nsamples;      // total number of samples in the track
        std::vector<WaveBlock> blockList;

    public:
        void Parse(tinyxml2::XMLElement *trackElem);
        int NumBlocks() const { return static_cast<int>(blockList.size()); }
        long NumSamples() const { return nsamples; }
        const WaveBlock& Block(int index) const { return blockList.at(index); }
    };

    class Project
    {
    private:
        std::vector<WaveTrack> channelList;
        std::string dataPath;   // e.g. "/home/don/radio/edit/2017-12-10_data/"

        struct DcBias
        {
            float left;
            float right;
            float limit;

            DcBias(float _left, float _right, float _limit)
                : left(_left)
                , right(_right)
                , limit(_limit)
                {}
        };

        struct ScanInfo
        {
            DcBias  bias;
            PreGapList gaplist;

            ScanInfo(DcBias _bias, const PreGapList & _gaplist)
                : bias(_bias)
                , gaplist(_gaplist)
                {}
        };

    public:
        void Load(const char *inFileName);
        void Convert(std::string outFilePrefix);

    private:
        void InitDataPath(const char *inFileName);
        std::string BlockFileName(const std::string& rawBlockFileName) const;

        static int FoldTally(
            std::vector<int>& folded, 
            const std::vector<int>& histogram,
            double bias
        );

        static double FindLimit(
            const std::vector<int> & leftHistogram, 
            double leftBias, 
            const std::vector<int> & rightHistogram, 
            double rightBias
        );

        ScanInfo PreScan() const;
        static std::string OutProgramFileName(std::string prefix, int hour);

        bool IsStartingNextProgram(
            long &boundary,     // out: offset into block to split program (if function returns true)
            int hour, 
            long recordingPosition,
            long programPosition, 
            long blockLength, 
            const PreGapList &gaplist) const;

        static FloatVector SplitBuffer(FloatVector& buffer, long offset);
        static void Append(FloatVector &target, const FloatVector& source);
        static bool Overlap(double a, double b, double x, double y)
        {
            return 
                (a >= x && a <= y) ||
                (b >= x && b <= y) ||
                (x >= a && x <= b) ||
                (y >= a && y <= b)
            ;
        }

        static void PrintProgramSummary(const std::string& filename, const GlitchRemover &remover);
    };

    class AudioReader   // reads .au files
    {
    private:
        FILE *infile;
        std::string filename;
        long nsamples;
        uint32_t rate;
        uint32_t channels;
        uint32_t dataOffset;

    public:
        AudioReader(std::string inFileName);
        ~AudioReader();

        void AssertChannels(uint32_t requiredNumChannels) const;
        void Seek(long sampleIndex);
        void Read(float *buffer, int length);
        const std::string& Filename() const { return filename; }

        double DurationSeconds() const
        {
            return static_cast<double>(nsamples) / static_cast<double>(rate);
        }

    private:
        static uint32_t DecodeInt(const char *buffer, int offset);
    };

    struct Chunk
    {
        FloatVector left;
        FloatVector right;
        long position;

        Chunk()
            : position(-1L)
            {}

        Chunk(const FloatVector& _left, const FloatVector& _right, int offset, int length, long front)
            : left(&_left[offset], &_left[offset + length])
            , right(&_right[offset], &_right[offset + length])
            , position(front + offset)
            {}

        void Clear()
        {
            left.clear();
            right.clear();
            position = -1L;
        }

        int Length() const
        {
            return static_cast<int>(left.size());
        }

        float Peak() const
        {
            return std::max(PeakValue(left), PeakValue(right));
        }

        FloatVector& First(bool flip) 
        {
            return flip ? right : left;
        }

        const FloatVector& Second(bool flip) const
        {
            return flip ? left : right;
        }

        void Extend(const FloatVector& _left, const FloatVector& _right, int offset, int extra)
        {
            for (int i=0; i < extra; ++i)
            {
                left.push_back(_left[offset+i]);
                right.push_back(_right[offset+i]);
            }
        }
    };

    class AudioWriter   // writes .au files
    {
    private:
        FILE *outfile;
        std::string outFileName;
        FloatVector buffer;
        int rate;
        int channels;

    public:
        AudioWriter(std::string _outFileName, int _rate, int _channels);
        ~AudioWriter()
        {
            Close();
        }

        void WriteChunk(const Chunk &chunk)
        {
            WriteStereo(chunk.left, chunk.right);
        }

        void Close();
        void StartNewFile(std::string _outFileName);

        std::string OutFileName() const
        {
            return outFileName;
        }

        static void DeleteFrontSamples(std::string filename, long numSamples);

    private:
        void WriteStereo(const FloatVector& left, const FloatVector& right);
        void WriteData(const void *data, size_t nbytes);
    };

    class GlitchGraph       // tracks glitches per minute, so we can create a density graph
    {
    private:
        static const int HEIGHT_LIMIT;
        std::vector<int> tally;

    public:
        void Reset();
        void Increment(long programSampleOffset);
        std::string Format() const;

    private:
        int Height(size_t minute) const;
        static char MinuteTickMark(size_t minute);
    };

    class GlitchRemover
    {
    private:
        long position;
        long programStartPosition;
        long glitchStartSample;
        int glitchCount;
        AudioWriter& writer;
        Chunk lastGoodChunk;
        std::deque<Chunk> chunklist;
        Chunk partial;
        const float sampleLimit;
        float programPeak;
        GlitchGraph graph;
        GapFinder gapFinder;

        static const int ChunkSamples = 1000;
        static const size_t MaxGlitchChunks = 3;
        static const int CrossFadeSamples = 250;

    public:
        GlitchRemover(AudioWriter& _writer, float _sampleLimit)
            : position(0)
            , programStartPosition(0)
            , glitchStartSample(0)
            , glitchCount(0)
            , writer(_writer)
            , sampleLimit(_sampleLimit)
            , programPeak(0.0f)
            {}

        void Fix(FloatVector& left, FloatVector& right);
        void WriteChunk(const Chunk &chunk);
        void Flush();
        void AdjustProgramLength(const std::string& filename, size_t programLengthSamples);

        int GlitchCount() const
        {
            return glitchCount;
        }

        float ProgramPeak() const
        {
            return programPeak;
        }

        float HeadroomDecibels() const
        {
            if (programPeak <= 0.0)
                return 144.0;

            return -20.0 * log10(programPeak);
        }

        std::string FormatGlitchGraph() const
        {
            return graph.Format();
        }

        void ResetProgram()
        {
            programStartPosition = position;
            glitchCount = 0;
            programPeak = 0.0f;
            graph.Reset();
        }

    private:
        void ProcessChunk(Chunk chunk);
        int ChunkListSampleCount() const;
        void CrossFade(Chunk &chunk);
        static bool IsHeartsOfSpace(const std::string& filename);
        static int GetInteger(const std::string& text, int start, int length);
    };
}

#endif // __DDC_UNGLITCH_H
