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
    bool IsLittleEndian()
    {
        unsigned n = 0x30313233;
        return ((const char *)&n)[0] == '3';
    }

    void ToggleFloatEndian(float *buffer, int length)
    {
        // Reverse bytes for little-endian machines.
        if (IsLittleEndian())
        {
            uint32_t *raw = reinterpret_cast<uint32_t *>(buffer);
            for (int i=0; i < length; ++i)
            {
                raw[i] = 
                    ((raw[i] & 0x000000ffu) << 24) |
                    ((raw[i] & 0x0000ff00u) << 8) |
                    ((raw[i] & 0x00ff0000u) >> 8) |
                    ((raw[i] & 0xff000000u) >> 24);
            }
        }
    }

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

        AudioWriter writer(outFileName);

        vector<float> leftBuffer;
        vector<float> rightBuffer;
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

            leftReader.Read(&leftBuffer.front(), leftBlock.Length());
            rightReader.Read(&rightBuffer.front(), rightBlock.Length());

            position += leftBlock.Length();
        }
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

        ToggleFloatEndian(buffer, length);
    }

    AudioWriter::AudioWriter(std::string outFileName)
        : outfile(fopen(outFileName.c_str(), "wb"))
    {
        if (!outfile)
            throw Error("Cannot open output file " + outFileName);
    }

    AudioWriter::~AudioWriter()
    {
        if (outfile)
        {
            fclose(outfile);
            outfile = nullptr;
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
