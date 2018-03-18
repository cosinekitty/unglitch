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
#include <iostream>
#include <iomanip>
#include <string>
#include <deque>
#include <vector>

#include "tinyxml2.h"

namespace unglitch
{
    bool IsLittleEndian();
    void ToggleFloatEndian(float *buffer, int length);

    typedef std::vector<float> FloatVector;

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

        DcBias FindBias() const;
        static std::string OutProgramFileName(std::string prefix, int hour);
        bool IsStartingNextProgram(int hour, long programPosition, const FloatVector& leftBuffer, const FloatVector& rightBuffer, long &boundary) const;
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

    enum class ChunkStatus
    {
        Unknown,
        Keep,
        Discard,
        CancelGlitch,
    };

    struct Chunk
    {
        ChunkStatus status;
        FloatVector left;
        FloatVector right;
        long position;

        Chunk()
            : status(ChunkStatus::Unknown)
            , position(-1L)
            {}

        Chunk(const FloatVector& _left, const FloatVector& _right, int offset, int length, long front)
            : status(ChunkStatus::Unknown)
            , left(&_left[offset], &_left[offset + length])
            , right(&_right[offset], &_right[offset + length])
            , position(front + offset)
            {}

        void Clear()
        {
            status = ChunkStatus::Unknown;
            left.clear();
            right.clear();
            position = -1L;
        }

        int Length() const
        {
            return static_cast<int>(left.size());
        }

        bool IsResolved() const
        {
            return status != ChunkStatus::Unknown;
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

    private:
        void WriteStereo(const FloatVector& left, const FloatVector& right);
        void WriteData(const void *data, size_t nbytes);
    };

    struct GlitchChannelState
    {
        int runLength;
        float prevPeak;

        GlitchChannelState()
            : runLength(0)
            , prevPeak(0.0f)
            {}
    };

    class GlitchRemover
    {
    private:
        long position;
        long glitchStartSample;
        int glitchCount;
        AudioWriter& writer;
        GlitchChannelState leftState;
        GlitchChannelState rightState;
        Chunk lastGoodChunk;
        std::deque<Chunk> chunklist;
        Chunk partial;
        const float sampleLimit;

        static const int ChunkSamples = 1000;
        static const int MaxGlitchChunks = 3;
        static const int CrossFadeSamples = 250;

    public:
        GlitchRemover(AudioWriter& _writer, float _sampleLimit)
            : position(0)
            , glitchStartSample(0)
            , glitchCount(0)
            , writer(_writer)
            , sampleLimit(_sampleLimit)
            {}

        void Fix(FloatVector& left, FloatVector& right);
        void Flush();
        void WriteChunk(const Chunk &chunk);

        int GlitchCount() const
        {
            return glitchCount;
        }

        void ResetGlitchCount()
        {
            glitchCount = 0;
        }

    private:
        void ProcessChunk(Chunk chunk);

        void ProcessChunkChannel(
            long sample,
            GlitchChannelState &state, 
            const FloatVector &first, 
            const FloatVector &second,
            ChunkStatus &status);

        int ChunkListSampleCount() const;
        void CrossFade(Chunk &first, Chunk &last);
        void WarnExceedSampleLimit(const Chunk& chunk) const;
    };
}

#endif // __DDC_UNGLITCH_H
