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
#include <string>
#include <deque>
#include <vector>

#include "tinyxml2.h"

namespace unglitch
{
    bool IsLittleEndian();
    void ToggleFloatEndian(float *buffer, int length);

    typedef std::vector<float> FloatVector;

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
        WaveBlock& Block(int index) { return blockList.at(index); }
    };

    class Project
    {
    private:
        std::vector<WaveTrack> channelList;
        std::string dataPath;   // e.g. "/home/don/radio/edit/2017-12-10_data/"

    public:
        void Load(const char *inFileName);
        void Convert(const char *outFileName);

    private:
        void InitDataPath(const char *inFileName);
        std::string BlockFileName(const std::string& rawBlockFileName) const;
    };

    class AudioReader   // reads .au files
    {
    private:
        FILE *infile;
        std::string filename;

    public:
        AudioReader(std::string inFileName);
        ~AudioReader();

        void Read(float *buffer, int length);
        const std::string& Filename() const { return filename; }

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

        Chunk()
            : status(ChunkStatus::Unknown)
            {}

        Chunk(const FloatVector& _left, const FloatVector& _right, int offset, int length)
            : status(ChunkStatus::Unknown)
            , left(&_left[offset], &_left[offset + length])
            , right(&_right[offset], &_right[offset + length])
            {}

        void Clear()
        {
            status = ChunkStatus::Unknown;
            left.clear();
            right.clear();
        }

        int Length() const
        {
            return static_cast<int>(left.size());
        }

        bool IsResolved() const
        {
            return status != ChunkStatus::Unknown;
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

    public:
        AudioWriter(std::string _outFileName, int _rate, int _channels);
        ~AudioWriter();

        void WriteStereo(const FloatVector& left, const FloatVector& right);

        void WriteChunk(const Chunk &chunk)
        {
            WriteStereo(chunk.left, chunk.right);
        }

    private:
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
        AudioWriter& writer;
        GlitchChannelState leftState;
        GlitchChannelState rightState;
        std::deque<Chunk> chunklist;
        Chunk partial;

        static const int ChunkSamples = 500;
        static const int MaxGlitchChunks = 4;
        static const int WindowChunks = 1 + MaxGlitchChunks + 1;

    public:
        GlitchRemover(AudioWriter& _writer)
            : writer(_writer)
            {}

        void Fix(FloatVector& left, FloatVector& right);
        void Flush();

    private:
        static float PeakValue(const FloatVector& vect);
        void ProcessChunk(Chunk chunk);

        void ProcessChunkChannel(
            GlitchChannelState &state, 
            const FloatVector &first, 
            const FloatVector &second,
            ChunkStatus &status);
    };
}

#endif // __DDC_UNGLITCH_H