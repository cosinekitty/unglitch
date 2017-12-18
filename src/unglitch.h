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
        float Threshold() const;
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

    private:
        static uint32_t DecodeInt(const char *buffer, int offset);
    };

    void Consume(FloatVector &buffer, int nsamples);

    class AudioWriter   // writes .au files
    {
    private:
        FILE *outfile;
        std::string outFileName;
        FloatVector buffer;

    public:
        AudioWriter(std::string _outFileName, int _rate, int _channels);
        ~AudioWriter();

        void WriteStereo(FloatVector& left, FloatVector& right);

    private:
        void WriteData(const void *data, size_t nbytes);
    };

    class GlitchFilter
    {
    private:
        const int minGlitchSamples;   // shortest run of large sample values to fix
        const int maxGlitchSamples;   // maximum number of consecutive samples to fix
        const int gapSamples;         // number of quiet samples after a glitch to trigger ending the glitch
        float threshold;              // absolute value above which we consider a glitch
        bool inGlitch;
        int quietSampleCount;
        float peak;
        FloatVector glitch;           // holds raw samples for a glitch in progress
        long sampleOffset;
        long glitchOffset;

    public:
        GlitchFilter(int _minGlitchSamples, int _maxGlitchSamples, int _gapSamples, float _threshold);
        void FixGlitches(const FloatVector& inBuffer, FloatVector& outBuffer);
        void Flush(FloatVector& outBuffer);
    };
}

#endif // __DDC_UNGLITCH_H