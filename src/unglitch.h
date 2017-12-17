/*
    unglitch.h  -  Don Cross  -  13 December 2017.
*/
#ifndef __DDC_UNGLITCH_H
#define __DDC_UNGLITCH_H

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

    private:
        static uint32_t DecodeInt(const char *buffer, int offset);
    };

    class AudioWriter   // writes .au files
    {
    private:
        FILE *outfile;

    public:
        AudioWriter(std::string outFileName);        
        ~AudioWriter();
    };
}

#endif // __DDC_UNGLITCH_H