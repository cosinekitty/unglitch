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

        AudioWriter writer(outFileName, 44100, 2);

        GlitchRemover remover(writer);

        FloatVector leftBuffer;
        FloatVector rightBuffer;
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

            leftBuffer.resize(leftBlock.Length());
            leftReader.Read(leftBuffer.data(), leftBlock.Length());

            rightBuffer.resize(rightBlock.Length());
            rightReader.Read(rightBuffer.data(), rightBlock.Length());

            remover.Fix(leftBuffer, rightBuffer);

            position += leftBlock.Length();
        }

        remover.Flush();
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
        int nread = fread(buffer, sizeof(float), length, infile);
        if (nread != length)
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

    void AudioWriter::WriteStereo(const FloatVector& left, const FloatVector& right)
    {
        if (left.size() != right.size())
            throw Error("WriteStereo: buffer size mismatch.");

        const int length = static_cast<int>(left.size());

        // Merge and interleave the sample data from both input channels.
        buffer.resize(2 * length);
        int k = 0;
        for (int i=0; i < length; ++i)
        {
            buffer[k++] = left[i];
            buffer[k++] = right[i];
        }

        WriteData(buffer.data(), 2 * sizeof(float) * length);
    }

    void AudioWriter::WriteData(const void *data, size_t nbytes)
    {
        size_t written = fwrite(data, 1, nbytes, outfile);
        if (written != nbytes)
            throw Error("Error writing to file " + outFileName);
    }

    std::string Format(long value, int width)
    {
        using namespace std;

        string text = to_string(value);
        while (static_cast<int>(text.length()) < width)
            text = "0" + text;

        return text;
    }

    std::string TimeStamp(long offset)
    {
        using namespace std;

        const long rate = 44100L;

        long seconds = offset / rate;
        offset %= rate;
        long millis = (offset * 1000L) / rate;
        long minutes = seconds / 60L;
        seconds %= 60L;
        long hours = minutes / 60L;
        minutes %= 60L;

        return        
            to_string(hours) + ":" + 
            Format(minutes, 2) + ":" + 
            Format(seconds, 2) + "." + 
            Format(millis, 3);
    }

    void GlitchRemover::Fix(FloatVector& left, FloatVector& right)
    {
        // Break both channels into chunks of ChunkSamples samples.
        // Characterize the chunk by its max absolute value.
        // Crudest algorithm: optionally exclude entire run of 1..MaxGlitchChunks.
        // Improvement: on outermost chunks in a run to be excluded,
        // find exact point in middle to start/stop chopping first/last.
        // Improvement: after removing a section, cross-fade to eliminate popping.

        // Convert flat (left, right) buffers into chunks.        
        // There may be a partial chunk left over from the last iteration.
        // If so, extend it to ChunkSamples in length.
        int length = static_cast<int>(left.size());

        int offset = 0;
        if (partial.Length() > 0)
        {
            int extra = std::min(ChunkSamples - partial.Length(), length);
            partial.Extend(left, right, 0, extra);
            
            if (partial.Length() < ChunkSamples)
                goto done;     // we ran out of data this time... nothing left to do

            ProcessChunk(position - extra, partial);
            partial.Clear();
            offset += extra;
        }

        while (offset + ChunkSamples <= length)
        {
            ProcessChunk(position + offset, Chunk(left, right, offset, ChunkSamples));
            offset += ChunkSamples;
        }

        // There may be a partial chunk left over from this iteration.
        if (offset < length)
            partial = Chunk(left, right, offset, length-offset);

    done:
        position += length;
    }

    void GlitchRemover::ProcessChunk(long sample, Chunk chunk)
    {
        using namespace std;

        // Called once for every new chunk.
        // If we are already in a glitch in either channel, check for end of glitch.
        // A glitch can end when we exceed maximum chunk count for a glitch,
        // or when it looks like the glitch is naturally over.
        // If we were not in a glitch and everything still looks fine,
        // just write this chunk to the output and keep going.

        chunk.status = ChunkStatus::Keep;

        ProcessChunkChannel(sample, leftState, chunk.left, chunk.right, chunk.status);
        ProcessChunkChannel(sample, rightState, chunk.right, chunk.left, chunk.status);

        switch (chunk.status)
        {
        case ChunkStatus::Discard:
            chunklist.push_back(chunk);
            break;

        case ChunkStatus::CancelGlitch:
            Flush();
            writer.WriteChunk(chunk);
            break;

        case ChunkStatus::Keep:
            if (!chunklist.empty())
            {
                cout << ++glitchCount << ". Discarding " << ChunkListSampleCount() << " samples at " << TimeStamp(glitchStartSample) << endl;
                CrossFade();
            }
            writer.WriteChunk(chunk);
            break;

        default:
            throw Error("Invalid chunk.status");
        }
    }

    void GlitchRemover::CrossFade()
    {
        if (chunklist.empty())
            return;

        // Keep the first few samples of the first bad chunk (fading out)
        // and the last few samples of the last bad chunk (fading in).
        // This eliminates an abrupt popping sound for the audio we are removing.
        // Note that the first and last bad chunks may be the same chunk.

        Chunk &first = chunklist.front();
        Chunk &last = chunklist.back();

        if (first.Length() != ChunkSamples)
            throw Error("CrossFade: first chunk is wrong length.");

        if (last.Length() != ChunkSamples)
            throw Error("CrossFade: second chunk is wrong length.");

        int k = ChunkSamples - CrossFadeSamples;
        for (int i=0; i < CrossFadeSamples; ++i)
        {
            float fadeIn = static_cast<float>(i) / static_cast<float>(CrossFadeSamples-1);
            float fadeOut = 1.0f - fadeIn;
            first.left[i] = fadeOut*first.left[i] + fadeIn*last.left[i+k];
            first.right[i] = fadeOut*first.right[i] + fadeIn*last.right[i+k];
        }

        first.left.resize(CrossFadeSamples);
        first.right.resize(CrossFadeSamples);

        writer.WriteStereo(first.left, first.right);
        chunklist.clear();
    }

    void GlitchRemover::ProcessChunkChannel(
        long sample,
        GlitchChannelState &state, 
        const FloatVector &first, 
        const FloatVector &second,
        ChunkStatus &status)
    {
        // There are 4 possible cases:
        // 1. (most common) start in good state, end in good state.
        // 2. go from good state to glitch state
        // 3. stay in glitch state
        // 4. end glitch state (no more glitch, or glitch too long)

        const float BadJump = 0.1f;
        float peak1 = PeakValue(first);
        float peak2 = PeakValue(second);
        bool inglitch = (peak1 - state.prevPeak > BadJump) && (peak1 - peak2 > BadJump);

        if (state.runLength == 0)
        {
            // Starting in good state.
            if (inglitch)
            {
                // Entering a glitch. Keep prevPeak value as-is.
                state.runLength = 1;
                status = ChunkStatus::Discard;
                glitchStartSample = sample;
            }
            else
            {
                // Staying in good state. Update prevPeak.
                state.prevPeak = peak1;
            }
        }
        else
        {
            // In bad state. Do we end the bad state now or keep going?
            if (inglitch)
            {
                if (++state.runLength <= MaxGlitchChunks)
                {
                    // Staying in glitch state.
                    status = ChunkStatus::Discard;
                }
                else
                {
                    // Glitch is too long. Cancel glitch.
                    state.runLength = 0;
                    state.prevPeak = peak1;
                    status = ChunkStatus::CancelGlitch;
                }
            }
            else
            {
                // Leaving glitch state.
                state.prevPeak = peak1;
                state.runLength = 0;
            }
        }
    }

    void GlitchRemover::Flush()
    {
        // Write any chunks remaining in the chunklist.
        while (!chunklist.empty())
        {
            writer.WriteChunk(chunklist.front());
            chunklist.pop_front();
        }

        // Write partial chunk if it contains any data.
        if (partial.Length() > 0)
        {
            writer.WriteChunk(partial);
            partial.Clear();
        }
    }

    float GlitchRemover::PeakValue(const FloatVector& vect)
    {
        float peak = 0.0f;
        for (float x : vect)
            peak = std::max(peak, std::abs(x));

        return peak;
    }

    int GlitchRemover::ChunkListSampleCount() const
    {
        int sum = 0;
        for (const Chunk &chunk : chunklist)
            sum += chunk.Length();
        return sum;
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
