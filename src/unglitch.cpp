/*
    unglitch.cpp  -  Don Cross  -  13 December 2017.

    Removes certain glitches that appear when I record radio programs.
*/

#include "unglitch.h"

/*
	<wavetrack name="2017-12-10" channel="0" linked="1" mute="0" solo="0" height="160" minimized="0" isSelected="1" rate="44100" gain="1.0" pan="0.0">
		<waveclip offset="0.00000000">
			<sequence maxsamples="262144" sampleformat="262159" numsamples="474554741">
				<waveblock start="0">
					<simpleblockfile filename="e0031f11.au" len="194051" min="-0.28178" max="0.288315" rms="0.068274"/>
				</waveblock>
                <!-- ... more waveblocks ... -->
            </sequence>
        </waveclip>
    </wavetrack>
*/

namespace unglitch
{
    long NumericAttribute(tinyxml2::XMLElement *elem, const char *name)
    {
        using namespace std;

        const char *attr = elem->Attribute(name);
        if (!attr)
            throw Error(string("Cannot find attribute ") + name);

        return atol(attr);
    }

    double FloatAttribute(tinyxml2::XMLElement *elem, const char *name)
    {
        using namespace std;

        const char *attr = elem->Attribute(name);
        if (!attr)
            throw Error(string("Cannot find attribute ") + name);
        
        return atof(attr);
    }

    void WaveTrack::Parse(tinyxml2::XMLElement *trackElem)
    {
        using namespace std;
        using namespace tinyxml2;

        // extract sample rate attribute from <wavetrack ...>
        rate = NumericAttribute(trackElem, "rate");
        if (rate != 44100)
            throw Error(string("Unexpected sampling rate ") + to_string(rate));

        XMLElement *clipElem = trackElem->FirstChildElement("waveclip");
        if (!clipElem)
            throw Error("Cannot find <waveclip> inside <wavetrack>");

        XMLElement *seqElem = clipElem->FirstChildElement("sequence");
        if (!seqElem)
            throw Error("Cannot find <sequence> inside <waveclip>");

        nsamples = NumericAttribute(seqElem, "numsamples");

        for (XMLElement *blockElem = seqElem->FirstChildElement("waveblock"); blockElem; blockElem = blockElem->NextSiblingElement("waveblock"))
        {
            /*
                <waveblock start="0">
					<simpleblockfile 
                        filename="e0031f11.au" 
                        len="194051" 
                        min="-0.28178" 
                        max="0.288315" 
                        rms="0.068274"
                    />
				</waveblock>
            */

            long start = NumericAttribute(blockElem, "start");

            XMLElement *fileElem = blockElem->FirstChildElement("simpleblockfile");
            if (!fileElem)
                throw Error("Cannot find <simpleblockfile> inside <waveblock>");

            const char *filename = fileElem->Attribute("filename");
            if (!filename)
                throw Error("Missing filename attribute in <simpleblockfile>");

            long length = NumericAttribute(fileElem, "len");
            float ymin = FloatAttribute(fileElem, "min");
            float ymax = FloatAttribute(fileElem, "max");
            blockList.push_back(WaveBlock(filename, start, length, ymin, ymax));
        }
    }

    float WaveTrack::Threshold() const
    {
        // Find the smallest non-negative value that includes the vast majority of the block peaks.
        // For simplicity we assume peak values are always in the range [0.001 .. 1.000]
        // and we keep a histogram of bands rounded up to the next higher increment of 0.001.
        const int NUMBANDS = 1000;
        std::vector<int> histogram(1 + NUMBANDS);

        for (const WaveBlock& block : blockList)
        {
            float peak = block.Peak();
            float fband = ceil(peak * NUMBANDS);
            if (fband < 0 || fband > NUMBANDS)
                throw Error("Invalid fband = " + std::to_string(fband));

            int band = static_cast<int>(fband);
            ++histogram[band];
        }

        const float KeepRatio = 0.98;
        int sum = 0;
        for (int band=0; band <= NUMBANDS; ++band)
        {
            sum += histogram[band];
            float ratio = static_cast<float>(sum) / static_cast<float>(blockList.size());
            if (ratio >= KeepRatio)
                return static_cast<float>(band) / static_cast<float>(NUMBANDS);
        }

        throw Error("Could not find threshold");
    }

    void Project::InitDataPath(const char *inFileName)
    {
        using namespace std;

        // Given inFileName   = "/home/don/radio/edit/2017-12-10.aup",
        // calculate dataPath = "/home/don/radio/edit/2017-12-10_data/"
        const char * const find = ".aup";
        int length = strlen(inFileName);
        int flen = strlen(find);
        if ((length < flen) || strcmp(&inFileName[length-flen], find))
            throw Error("Invalid Audacity project filename. Must end with .aup");

        dataPath = string(inFileName).substr(0, length-flen) + "_data/";
    }

    void Project::Load(const char *inFileName)
    {
        using namespace std;
        using namespace tinyxml2;

        InitDataPath(inFileName);

        XMLDocument doc;
        XMLError status = doc.LoadFile(inFileName);
        if (status != XML_SUCCESS)
            throw Error(string("Cannot open Audacity project file ") + inFileName);

        cout << "Loaded xml: " << inFileName << endl;

        XMLElement *root = doc.RootElement();
        for (XMLElement *trackElem = root->FirstChildElement("wavetrack"); trackElem; trackElem = trackElem->NextSiblingElement("wavetrack"))  
        {
            cout << "Found wavetrack" << endl;
            channelList.push_back(WaveTrack());
            WaveTrack &track = channelList.back();
            track.Parse(trackElem);
        }
    }

    void Project::Convert(const char *outFileName)
    {
        using namespace std;

        // Convert multiple pairs of single-channel audio into a single stereo audio file.
        // Assume that left and right channes come in equal-size blocks.
        if (channelList.size() != 2)
            throw Error("Input must be stereo.");

        WaveTrack &leftTrack = channelList[0];
        WaveTrack &rightTrack = channelList[1];

        const int nblocks = leftTrack.NumBlocks();
        if (nblocks != rightTrack.NumBlocks())
            throw Error("Left and right tracks have different number of blocks.");

        float threshold = std::max(leftTrack.Threshold(), rightTrack.Threshold());
        cout << "Threshold = " << threshold << endl;
        AudioWriter writer(outFileName, 44100, 2);

        const int maxGlitchSamples = 2000;
        const int gapSamples = 1000;
        GlitchFilter leftFilter(maxGlitchSamples, gapSamples, threshold);
        GlitchFilter rightFilter(maxGlitchSamples, gapSamples, threshold);

        FloatVector leftBuffer;
        FloatVector rightBuffer;
        FloatVector leftCleaned;
        FloatVector rightCleaned;
        long position = 0;
        for (int b=0; b < nblocks; ++b)
        {
            WaveBlock& leftBlock = leftTrack.Block(b);
            WaveBlock& rightBlock = rightTrack.Block(b);            

            if (leftBlock.Start() != position)
                throw Error("Left block not at expected position.");

            if (rightBlock.Start() != position)
                throw Error("Right block not at expected position.");

            if (leftBlock.Length() != rightBlock.Length())
                throw Error("Left and right blocks have different lengths.");            

            AudioReader leftReader(BlockFileName(leftBlock.Filename()));
            AudioReader rightReader(BlockFileName(rightBlock.Filename()));

            if (static_cast<long>(leftBuffer.size()) < leftBlock.Length())
                leftBuffer.resize(leftBlock.Length());

            if (static_cast<long>(rightBuffer.size()) < rightBlock.Length())
                rightBuffer.resize(rightBlock.Length());

            leftReader.Read(leftBuffer.data(), leftBlock.Length());
            rightReader.Read(rightBuffer.data(), rightBlock.Length());

            leftFilter.FixGlitches(leftBuffer, leftCleaned);
            rightFilter.FixGlitches(rightBuffer, rightCleaned);

            writer.WriteStereo(leftCleaned, rightCleaned);

            position += leftBlock.Length();
        }

        leftFilter.Flush(leftCleaned);
        rightFilter.Flush(rightCleaned);
        if (leftCleaned.size() != rightCleaned.size())
            throw Error("Bug detected: post-flush left and right channel lengths do not match!");

        writer.WriteStereo(leftCleaned, rightCleaned);
    }

    std::string Project::BlockFileName(const std::string& fn) const
    {
        // Input:  "e0020ade.au"
        // Output: "/home/don/radio/edit/2017-12-10_data/e00/d20/e0020ade.au"
        return dataPath + fn.substr(0, 3) + "/d" + fn.substr(3, 2) + "/" + fn;
    }

    AudioReader::AudioReader(std::string inFileName)
        : infile(fopen(inFileName.c_str(), "rb"))
        , filename(inFileName)
    {
        if (!infile)
            throw Error("Cannot open input audio file " + inFileName);

        // Decode .au file header
        // https://en.wikipedia.org/wiki/Au_file_format
        // 00000000  64 6e 73 2e 5c 30 00 00  ff ff ff ff 06 00 00 00  |dns.\0..........|
        // 00000010  44 ac 00 00 01 00 00 00  41 75 64 61 63 69 74 79  |D.......Audacity|
        // 00000020  42 6c 6f 63 6b 46 69 6c  65 31 31 32 1b 00 83 be  |BlockFile112....|
        // 00000030  cf cb 8f 3e 78 5c af 3d  58 b6 8a be 38 71 8f 3e  |...>x\.=X...8q.>|
        // 00000040  8b b9 ae 3d ae dc 84 be  dc 0d 9c 3e d4 6a ab 3d  |...=.......>.j.=|
        // 00000050  62 0b 93 be 55 6f 8e 3e  38 e7 af 3d 25 3a 81 be  |b...Uo.>8..=%:..|
        // 00000060  97 13 11 3e ae 32 b6 3d  78 b3 50 be 5f 34 53 3e  |...>.2.=x.P._4S>|
        // 00000070  cf 4a c7 3d 2b ed 22 be  92 7b 15 3e 31 dc 96 3d  |.J.=+."..{.>1..=|
        // 00000080  96 f1 3d be 78 8f e7 3d  72 24 a0 3d a4 e8 0d be  |..=.x..=r$.=....|

        // Detect Big-Endian versus little endian
        char header[24];
        int nread = fread(header, 1, sizeof(header), infile);
        if (nread != sizeof(header))
            throw Error("Cannot read file header for " + inFileName);

        if (memcmp(header, "dns.", 4))
            throw Error("Incorrect au file header in " + inFileName);        

        uint32_t dataOffset = DecodeInt(header, 1*4);
        uint32_t encoding = DecodeInt(header, 3*4);
        if (encoding != 6)
            throw Error("Unsupported data encoding: expected 32-bit IEEE floating point in " + inFileName);

        uint32_t rate = DecodeInt(header, 4*4);
        if (rate != 44100)
            throw Error("Unsupported sampling rate in file " + inFileName);

        uint32_t channels = DecodeInt(header, 5*4);
        if (channels != 1)
            throw Error("Expected mono audio in file " + inFileName);

        // Get ready to start reading audio data.
        if (fseek(infile, dataOffset, SEEK_SET))
            throw Error("Unable to seek to data offset in " + inFileName);
    }

    AudioReader::~AudioReader()
    {
        if (infile)
        {
            fclose(infile);
            infile = nullptr;
        }
    }

    uint32_t AudioReader::DecodeInt(const char *buffer, int offset)
    {
        uint32_t x = 0xffu & static_cast<uint32_t>(buffer[offset+3]);
        x = (x << 8) | (0xffu & static_cast<uint32_t>(buffer[offset+2]));
        x = (x << 8) | (0xffu & static_cast<uint32_t>(buffer[offset+1]));
        x = (x << 8) | (0xffu & static_cast<uint32_t>(buffer[offset]));
        return x;
    }

    void AudioReader::Read(float *buffer, int length)
    {
        int bytesWanted = sizeof(float) * length;
        int bytesRead = fread(buffer, 1, bytesWanted, infile);
        if (bytesRead != bytesWanted)
            throw Error("Could not read desired number of bytes from input file " + filename);
    }

    AudioWriter::AudioWriter(std::string _outFileName, int _rate, int _channels)
        : outfile(fopen(_outFileName.c_str(), "wb"))
        , outFileName(_outFileName)
    {
        if (!outfile)
            throw Error("Cannot open output file " + outFileName);

        // Write file header
        // 00000000  64 6e 73 2e 5c 30 00 00  ff ff ff ff 06 00 00 00  |dns.\0..........|
        // 00000010  44 ac 00 00 01 00 00 00  41 75 64 61 63 69 74 79  |D.......Audacity|
        // 00000020  42 6c 6f 63 6b 46 69 6c  65 31 31 32 1b 00 83 be  |BlockFile112....|
        // 00000030  cf cb 8f 3e 78 5c af 3d  58 b6 8a be 38 71 8f 3e  |...>x\.=X...8q.>|
        // 00000040  8b b9 ae 3d ae dc 84 be  dc 0d 9c 3e d4 6a ab 3d  |...=.......>.j.=|
        // 00000050  62 0b 93 be 55 6f 8e 3e  38 e7 af 3d 25 3a 81 be  |b...Uo.>8..=%:..|
        // 00000060  97 13 11 3e ae 32 b6 3d  78 b3 50 be 5f 34 53 3e  |...>.2.=x.P._4S>|
        // 00000070  cf 4a c7 3d 2b ed 22 be  92 7b 15 3e 31 dc 96 3d  |.J.=+."..{.>1..=|
        // 00000080  96 f1 3d be 78 8f e7 3d  72 24 a0 3d a4 e8 0d be  |..=.x..=r$.=....|

        uint32_t header[6];

        memcpy(header, "dns.", 4);      // audio file signature
        header[1] = sizeof(header);     // data begins immediately after this header
        header[2] = 0xffffffffu;        // total data size unknown
        header[3] = 6;                  // 32-bit IEEE floating point format
        header[4] = _rate;              // sampling rate in Hz
        header[5] = _channels;          // number of channels

        WriteData(header, sizeof(header));
    }

    AudioWriter::~AudioWriter()
    {
        if (outfile)
        {
            fclose(outfile);
            outfile = nullptr;
        }
    }

    void AudioWriter::WriteStereo(FloatVector& left, FloatVector& right)
    {
        const int leftLength = left.size();
        const int rightLength = right.size();
        const int length = std::min(leftLength, rightLength);

        // Merge and interleave the sample data from both input channels.
        buffer.resize(2 * length);
        int k = 0;
        for (int i=0; i < length; ++i)
        {
            buffer[k++] = left[i];
            buffer[k++] = right[i];
        }

        WriteData(buffer.data(), 2 * sizeof(float) * length);

        // Remove 'length' samples from both 'left' and 'right'.
        Consume(left, length);
        Consume(right, length);
    }

    void AudioWriter::WriteData(const void *data, size_t nbytes)
    {
        size_t written = fwrite(data, 1, nbytes, outfile);
        if (written != nbytes)
            throw Error("Error writing to file " + outFileName);
    }

    void Consume(FloatVector &buffer, int nsamples)
    {
        const int length = static_cast<int>(buffer.size());

        if (nsamples > length)
            throw Error("Attempt to consume more data than buffer contains!");

        for (int i=0; i + nsamples < length; ++i)
            buffer[i] = buffer[i + nsamples];

        buffer.resize(static_cast<size_t>(length - nsamples));
    }

    GlitchFilter::GlitchFilter(int _maxGlitchSamples, int _gapSamples, float _threshold)
        : maxGlitchSamples(_maxGlitchSamples)
        , gapSamples(_gapSamples)
        , threshold(_threshold)
        , inGlitch(false)
        , quietSampleCount(0)
    {        
    }

    void GlitchFilter::FixGlitches(const FloatVector& inBuffer, FloatVector& outBuffer)
    {        
        // State machine, processing one sample at a time.
        for (float x : inBuffer)
        {
            float ax = std::abs(x);

            if (!inGlitch && (ax > threshold))
            {
                inGlitch = true;
                peak = ax;
                quietSampleCount = 0;
            }

            if (inGlitch)
            {
                glitch.push_back(x);
                if (ax > threshold)
                {
                    quietSampleCount = 0;
                    if (ax > peak)
                        peak = ax;
                }
                else if (++quietSampleCount >= gapSamples)
                    Flush(outBuffer);
            }
            else
            {
                outBuffer.push_back(x);
            }
        }
    }

    void GlitchFilter::Flush(FloatVector& outBuffer)
    {        
        if (inGlitch)
        {
            // Leaving the glitch state.
            inGlitch = false;

            // Calculate glitch attentuation factor.
            float atten = threshold / peak;

            // Make the glitchy part of the glitch buffer quieter.
            int glitchSampleCount = static_cast<int>(glitch.size()) - quietSampleCount;
            for (int i=0; i < glitchSampleCount; ++i)
                glitch[i] *= atten;

            // Append the corrected glitch buffer to the output buffer.
            for (float g : glitch)
                outBuffer.push_back(g);

            // Empty out the glitch buffer so it is ready for another glitch.
            glitch.clear();

            std::cout << "GlitchFilter: fixed " << glitchSampleCount << " samples." << std::endl;
        }
    }
}

int main(int argc, const char *argv[])
{
    using namespace std;
    using namespace unglitch;

    if (argc != 3)
    {
        cerr << "USAGE: unglitch infile.aup outfile.au" << endl;
        return 1;
    }

    const char *inAudacityProjectFileName = argv[1];
    const char *outAudioFileName = argv[2];
    try
    {
        Project project;
        project.Load(inAudacityProjectFileName);
        cout << "Loaded project." << endl;
        project.Convert(outAudioFileName);
        cout << "Finished converting audio." << endl;
    }
    catch (const Error &error)
    {
        cerr << "ERROR: " << error.Message() << endl;
        return 1;
    }

    return 0;
}
