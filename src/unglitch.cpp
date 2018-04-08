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
    const double HEADROOM_DEFAULT = 9.2;
    const int SamplingRate = 44100;
    bool Verbose = false;

    inline double MinutesFromSamples(size_t samples)
    {
        return samples / (60.0 * SamplingRate);
    }

    inline double SecondsFromSamples(size_t samples)
    {
        return samples / static_cast<double>(SamplingRate);
    }

    inline size_t SamplesFromSeconds(double seconds)
    {
        if (seconds < 0.0)
            throw Error("Attempt to convert negative time to samples.");

        return static_cast<size_t>(seconds * SamplingRate);
    }

    inline size_t SamplesFromMinutes(double minutes)
    {
        return SamplesFromSeconds(60.0 * minutes);
    }

    inline float PeakFromHeadroom(double headroom_dB)
    {
        return static_cast<float>(pow(10.0, -headroom_dB/20.0));
    }

    std::string TimeStamp(long offset);

    long NumericAttribute(tinyxml2::XMLElement *elem, const char *name)
    {
        using namespace std;

        const char *attr = elem->Attribute(name);
        if (!attr)
            throw Error(string("Cannot find attribute ") + name);

        return atol(attr);
    }

    float FloatAttribute(tinyxml2::XMLElement *elem, const char *name)
    {
        using namespace std;

        const char *attr = elem->Attribute(name);
        if (!attr)
            throw Error(string("Cannot find attribute ") + name);

        return static_cast<float>(atof(attr));
    }

    void WaveTrack::Parse(tinyxml2::XMLElement *trackElem)
    {
        using namespace std;
        using namespace tinyxml2;

        // extract sample rate attribute from <wavetrack ...>
        rate = NumericAttribute(trackElem, "rate");
        if (rate != SamplingRate)
            throw Error(string("Unexpected sampling rate ") + to_string(rate));

        XMLElement *clipElem = trackElem->FirstChildElement("waveclip");
        if (!clipElem)
            throw Error("Cannot find <waveclip> inside <wavetrack>");

        XMLElement *seqElem = clipElem->FirstChildElement("sequence");
        if (!seqElem)
            throw Error("Cannot find <sequence> inside <waveclip>");

        nsamples = NumericAttribute(seqElem, "numsamples");
        long checkSamples = 0;

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
            checkSamples += length;
        }

        if (nsamples != checkSamples)
            throw Error("WaveTrack::Parse - numsamples does not match sum(len)");
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

        //cout << "Loaded xml: " << inFileName << endl;

        XMLElement *root = doc.RootElement();
        for (XMLElement *trackElem = root->FirstChildElement("wavetrack"); trackElem; trackElem = trackElem->NextSiblingElement("wavetrack"))
        {
            //cout << "Found wavetrack" << endl;
            channelList.push_back(WaveTrack());
            WaveTrack &track = channelList.back();
            track.Parse(trackElem);
        }
    }

    size_t Length(const FloatVector& left, const FloatVector& right)
    {
        using namespace std;

        if (left.size() != right.size())
            throw Error("Left and right buffers are different lengths: " + to_string(left.size()) + ", " + to_string(right.size()));

        return left.size();
    }

    void PreScanGapFinder::Process(const FloatVector& left, const FloatVector& right)
    {
        using namespace std;

        static const double MinSilenceSeconds  = 0.5;    // least duration we consider a silent period
        static const float  AmplitudeThreshold = 0.006f; // how loud can we get before not considered silent
        const size_t MinSilenceSamples = SamplesFromSeconds(MinSilenceSeconds);

        size_t length = Length(left, right);
        size_t preSamples = totalSamples;       // total number of samples processed before this block
        totalSamples += length;                 // total number of samples processed, including this block

        if (preSamples < sampleDelay)
        {
            // We are not yet ready to find silent periods because we
            // are still in the initial phase of estimating DC bias.

            for (float data : left)
                leftBias += data;

            for (float data : right)
                rightBias += data;

            if (totalSamples < sampleDelay)
                return;     // not yet ready to calculate bias or find silent periods

            // Initialize bias values and fall through to start finding silent periods
            leftBias /= totalSamples;
            rightBias /= totalSamples;
            cout << "ESTBIAS: " << TimeStamp(totalSamples) << " ==> left=" << leftBias << ", right=" << rightBias << endl;
        }

        // Look for silent periods centered around estimated DC bias for each channel.
        for (size_t i=0; i < length; ++i)
        {
            double signal = max(abs(left[i] - leftBias), abs(right[i] - rightBias));
            if (signal < AmplitudeThreshold)
            {
                if (silentSamples == 0)
                    silenceFront = preSamples + i;
                ++silentSamples;
            }
            else
            {
                if (silentSamples >= MinSilenceSamples)
                {
                    gaplist.push_back(PreGapInfo(silenceFront, silentSamples));

                    cout << "PREGAP #" << setw(3) << gaplist.size()
                        << " : front=" << TimeStamp(silenceFront)
                        << ", length=" << TimeStamp(silentSamples)
                        << endl;
                }
                silentSamples = 0;
            }
        }
    }

    Project::Project(const std::vector<double>& manualSplitPointsInSeconds)
    {
        using namespace std;

        // Use -1 as a sentinel that the split point should be automatic.
        for (int i=0; i < MAX_SPLIT_POINTS; ++i)
            manualSplitSamples[i] = -1L;

        // We can split between hours (1,2), (2,3), or (3,4).
        // Classify the split points by which program boundary they correspond to.
        for (double seconds : manualSplitPointsInSeconds)
        {
            // Round to the nearest hour mark.
            int hour = static_cast<int>(floor(0.5 + (seconds / 3600.0)));
            if (hour < 1 || hour > MAX_SPLIT_POINTS)
                throw Error("Invalid split point at offset = " + to_string(seconds) + " seconds.");

            if (manualSplitSamples[hour-1] != -1L)
                throw Error("Redundant split point after hour " + to_string(hour));

            manualSplitSamples[hour-1] = static_cast<long>(SamplesFromSeconds(seconds));
        }
    }

    Project::ScanInfo Project::PreScan() const
    {
        // Measure the DC offset in both channels.
        if (channelList.size() != 2)
            throw Error("Input must be stereo.");

        const WaveTrack &leftTrack = channelList[0];
        const WaveTrack &rightTrack = channelList[1];

        const int nblocks = leftTrack.NumBlocks();
        if (nblocks != rightTrack.NumBlocks())
            throw Error("Left and right tracks have different number of blocks.");

        PreScanGapFinder gapFinder(SamplesFromMinutes(20.0));

        double leftSum = 0.0;
        double rightSum = 0.0;
        long position = 0;
        FloatVector leftBuffer;
        FloatVector rightBuffer;
        for (int b=0; b < nblocks; ++b)
        {
            const WaveBlock& leftBlock = leftTrack.Block(b);
            const WaveBlock& rightBlock = rightTrack.Block(b);

            if (leftBlock.Start() != position)
                throw Error("Left block not at expected position.");

            if (rightBlock.Start() != position)
                throw Error("Right block not at expected position.");

            const long length = leftBlock.Length();
            if (length != rightBlock.Length())
                throw Error("Left and right blocks have different lengths.");

            AudioReader leftReader(BlockFileName(leftBlock.Filename()));
            leftReader.AssertChannels(1);

            AudioReader rightReader(BlockFileName(rightBlock.Filename()));
            rightReader.AssertChannels(1);

            leftBuffer.resize(length);
            leftReader.Read(leftBuffer.data(), length);

            rightBuffer.resize(length);
            rightReader.Read(rightBuffer.data(), length);

            for (float data : leftBuffer)
                leftSum += data;

            for (float data : rightBuffer)
                rightSum += data;

            gapFinder.Process(leftBuffer, rightBuffer);

            position += length;
        }

        double leftBias = leftSum / position;
        double rightBias = rightSum / position;

        return ScanInfo(
            DcBias(static_cast<float>(leftBias), static_cast<float>(rightBias)),
            gapFinder.SilentGaps()
        );
    }

    std::string Project::OutProgramFileName(std::string prefix, int hour)
    {
        return prefix + "-" + std::to_string(hour) + ".au";
    }

    void Project::PrintProgramSummary(
        const std::string& filename,
        const GlitchRemover &remover,
        long programDurationSamples)
    {
        using namespace std;

        cout << filename << " : "
            << "duration=" << TimeStamp(programDurationSamples) << ", "
            << remover.GlitchCount() << " glitches, peak="
            << setprecision(5) << remover.ProgramPeak()
            << ", headroom=" << setprecision(1) << remover.HeadroomDecibels() << " dB."
            << endl;

        const size_t roundedUpMinutes = static_cast<size_t>(
            ceil(MinutesFromSamples(static_cast<size_t>(programDurationSamples))));

        cout << remover.FormatGlitchGraph(roundedUpMinutes) << endl;
    }

    void Project::Convert(std::string projname, double headroom_dB)
    {
        using namespace std;

        string outFilePrefix = string("cleaned-") + projname;

        float limit = PeakFromHeadroom(headroom_dB);

        cout << "Prescanning..." << endl;
        ScanInfo scan = PreScan();
        cout << fixed
            << setprecision(5) << "Bias: left=" << scan.bias.left << ", right=" << scan.bias.right
            << setprecision(1) << ", headroom=" << headroom_dB << " dB"
            << " (limit=" << setprecision(4) << limit << ")"
            << endl;

        // Convert multiple pairs of single-channel audio into a single stereo audio file.
        // Assume that left and right channes come in equal-size blocks.
        if (channelList.size() != 2)
            throw Error("Input must be stereo.");

        WaveTrack &leftTrack = channelList[0];
        WaveTrack &rightTrack = channelList[1];

        const int nblocks = leftTrack.NumBlocks();
        if (nblocks != rightTrack.NumBlocks())
            throw Error("Left and right tracks have different number of blocks.");

        if (leftTrack.NumSamples() != rightTrack.NumSamples())
            throw Error("Left and right tracks have different number of samples.");

        int hour = 1;
        AudioWriter writer(OutProgramFileName(outFilePrefix, hour), SamplingRate, 2);

        GlitchRemover remover(writer, limit);

        FloatVector leftBuffer;
        FloatVector rightBuffer;
        long position = 0;
        long boundary = 0;
        long programPosition = 0;
        for (int b=0; b < nblocks; ++b)
        {
            const WaveBlock& leftBlock = leftTrack.Block(b);
            const WaveBlock& rightBlock = rightTrack.Block(b);

            if (leftBlock.Start() != position)
                throw Error("Left block not at expected position.");

            if (rightBlock.Start() != position)
                throw Error("Right block not at expected position.");

            if (leftBlock.Length() != rightBlock.Length())
                throw Error("Left and right blocks have different lengths.");

            const long blockLength = leftBlock.Length();

            AudioReader leftReader(BlockFileName(leftBlock.Filename()));
            AudioReader rightReader(BlockFileName(rightBlock.Filename()));

            leftBuffer.resize(blockLength);
            leftReader.Read(leftBuffer.data(), blockLength);

            rightBuffer.resize(blockLength);
            rightReader.Read(rightBuffer.data(), blockLength);

            for (float &data : leftBuffer)
                data -= scan.bias.left;

            for (float &data : rightBuffer)
                data -= scan.bias.right;

            if (IsStartingNextProgram(boundary, hour, position, programPosition, blockLength, scan.gaplist))
            {
                ++hour;
                cout << "SPLIT: " << TimeStamp(position + boundary) << endl;

                // Split buffers into before/after the boundary.
                FloatVector leftBefore = SplitBuffer(leftBuffer, boundary);
                FloatVector rightBefore = SplitBuffer(rightBuffer, boundary);
                remover.Fix(leftBefore, rightBefore);

                string prevFileName = writer.OutFileName();
                remover.Flush();
                PrintProgramSummary(prevFileName, remover, programPosition + boundary);
                remover.ResetProgram();

                writer.StartNewFile(OutProgramFileName(outFilePrefix, hour));
                remover.AdjustProgramLength(prevFileName, programPosition + boundary);

                programPosition = 0;
            }

            remover.Fix(leftBuffer, rightBuffer);

            position += blockLength;
            programPosition += leftBuffer.size();   // works whether or not we just split the program
        }

        string finalFileName = writer.OutFileName();
        remover.Flush();
        PrintProgramSummary(finalFileName, remover, programPosition);
        remover.ResetProgram();
        writer.Close();
        remover.AdjustProgramLength(finalFileName, programPosition);
    }

    bool Project::IsStartingNextProgram(
        long &boundary,     // out: offset into block to split program (if function returns true)
        int hour,
        long recordingPosition,
        long programPosition,
        long blockLength,
        const PreGapList &gaplist) const
    {
        using namespace std;

        boundary = 0;   // always initialize output parameter

        if ((hour >= 1) && (hour <= MAX_SPLIT_POINTS))
        {
            long split = manualSplitSamples[hour-1];
            if (split > 0)
            {
                if (split >= recordingPosition && split < recordingPosition + blockLength)
                {
                    // Manual override for the split position.
                    boundary = split - recordingPosition;
                    return true;
                }
            }
            else
            {
                // Look for automatic split position.
                // We are starting a new program if we find a silent period
                // somewhere between 58 and 61 minutes into the current program.
                const double minProgramMinutes = (hour==1) ? 58.25 : 59.25;
                const double maxProgramMinutes = 2.0 + minProgramMinutes;
                double minutesBegin = MinutesFromSamples(programPosition);
                double minutesEnd = MinutesFromSamples(programPosition + blockLength - 1);
                if (Overlap(minutesBegin, minutesEnd, minProgramMinutes, maxProgramMinutes))
                {
                    // This block is inside the window of times where we expect a program boundary.
                    // If we find a silent period whose center is inside this block,
                    // assume that is the program boundary.
                    for (const PreGapInfo & gap : gaplist)
                    {
                        long center = gap.Center();
                        if (center >= recordingPosition && center <= recordingPosition + blockLength)
                        {
                            boundary = center - recordingPosition;
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    FloatVector Project::SplitBuffer(FloatVector& buffer, long offset)
    {
        // Extract the beginning 'offset' samples from 'buffer' and return them.
        // Remove those samples from the beginning of 'buffer' itself.
        const long length = static_cast<long>(buffer.size());
        if (offset<0 || offset>length)
            throw Error("SplitBuffer: invalid offset");

        FloatVector front(offset);
        for (long i=0; i < offset; ++i)
            front[i] = buffer[i];

        for (long i=offset; i < length; ++i)
            buffer[i-offset] = buffer[i];

        buffer.resize(length-offset);
        return front;
    }

    void Project::Append(FloatVector &target, const FloatVector& source)
    {
        const size_t offset = target.size();
        const size_t increase = source.size();
        target.resize(offset + increase);
        for (size_t index=0; index < increase; ++index)
            target[offset + index] = source[index];
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

        char header[24];
        int nread = fread(header, 1, sizeof(header), infile);
        if (nread != sizeof(header))
            throw Error("Cannot read file header for " + inFileName);

        if (memcmp(header, "dns.", 4))
            throw Error("Incorrect au file header in " + inFileName);

        dataOffset = DecodeInt(header, 1*4);
        uint32_t encoding = DecodeInt(header, 3*4);
        if (encoding != 6)
            throw Error("Unsupported data encoding: expected 32-bit IEEE floating point in " + inFileName);

        rate = DecodeInt(header, 4*4);
        if (rate != SamplingRate)
            throw Error("Unsupported sampling rate in file " + inFileName);

        channels = DecodeInt(header, 5*4);

        // Figure out total number of samples.
        if (fseek(infile, 0, SEEK_END))
            throw Error("Unable to seek to EOF in " + inFileName);

        long size = ftell(infile);
        if (size < static_cast<long>(dataOffset))
            throw Error("Unable to determine file position in " + inFileName);

        nsamples = (size - dataOffset) / (channels * sizeof(float));

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

    void AudioReader::AssertChannels(uint32_t requiredNumChannels) const
    {
        using namespace std;

        if (channels != requiredNumChannels)
            throw Error("Expected " + to_string(requiredNumChannels) + " in file '" + filename + "', but found " + to_string(channels));
    }

    uint32_t AudioReader::DecodeInt(const char *buffer, int offset)
    {
        uint32_t x = 0xffu & static_cast<uint32_t>(buffer[offset+3]);
        x = (x << 8) | (0xffu & static_cast<uint32_t>(buffer[offset+2]));
        x = (x << 8) | (0xffu & static_cast<uint32_t>(buffer[offset+1]));
        x = (x << 8) | (0xffu & static_cast<uint32_t>(buffer[offset]));
        return x;
    }

    void AudioReader::Seek(long sampleIndex)
    {
        long fileOffset = (channels * sizeof(float))*sampleIndex + dataOffset;
        if (fseek(infile, fileOffset, SEEK_SET))
            throw Error("Unable to seek to sample index " + std::to_string(sampleIndex));
    }

    void AudioReader::Read(float *buffer, int length)
    {
        int nread = fread(buffer, sizeof(float), length, infile);
        if (nread != length)
            throw Error("Could not read desired number of bytes from input file " + filename);
    }

    AudioWriter::AudioWriter(std::string _outFileName, int _rate, int _channels)
        : outfile(nullptr)
        , rate(_rate)
        , channels(_channels)
    {
        StartNewFile(_outFileName);
    }

    void AudioWriter::Close()
    {
        if (outfile)
        {
            fclose(outfile);
            outfile = nullptr;
        }
    }

    void AudioWriter::DeleteFrontSamples(std::string filename, long numSamples)
    {
        FILE *file = fopen(filename.c_str(), "r+b");
        if (!file)
            throw Error("AudioWriter::DeleteFrontSamples: cannot open file: " + filename);

        const int NumChannels = 2;
        uint32_t oldDataOffset;
        if (fseek(file, 4, SEEK_SET))
        {
            fclose(file);
            throw Error("AudioWriter::DeleteFrontSamples: cannot seek for old data offset");
        }

        if (1 != fread(&oldDataOffset, sizeof(oldDataOffset), 1, file))
        {
            fclose(file);
            throw Error("AudioWriter::DeleteFrontSamples: cannot read old data offset");
        }

        uint32_t newDataOffset = oldDataOffset + (4 * NumChannels * numSamples);

        if (fseek(file, 4, SEEK_SET))
        {
            fclose(file);
            throw Error("AudioWriter::DeleteFrontSamples: cannot seek for new data offset");
        }

        if (1 != fwrite(&newDataOffset, sizeof(newDataOffset), 1, file))
        {
            fclose(file);
            throw Error("AudioWriter::DeleteFrontSamples: cannot write new data offset");
        }

        fclose(file);
    }

    void AudioWriter::StartNewFile(std::string _outFileName)
    {
        if (outfile)
            fclose(outfile);

        outFileName = _outFileName;
        outfile = fopen(_outFileName.c_str(), "wb");
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
        header[4] = rate;               // sampling rate in Hz
        header[5] = channels;           // number of channels

        WriteData(header, sizeof(header));
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

        long seconds = offset / SamplingRate;
        offset %= SamplingRate;
        long millis = ((offset * 1000L) + (SamplingRate/2)) / SamplingRate;
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

        // Convert flat (left, right) buffers into chunks.
        // There may be a partial chunk left over from the last iteration.
        // If so, extend it to ChunkSamples in length.
        int length = static_cast<int>(left.size());

        int offset = 0;
        int plen = partial.Length();
        if (plen > 0)
        {
            if (plen >= ChunkSamples)
                throw Error("Bug in GlitchRemover::Fix: plen=" + std::to_string(plen));

            // Calculate the number of samples we need to fill up the partial chunk.
            int extra = ChunkSamples - plen;
            if (length < extra)
            {
                // We don't have enough new data to make a complete chunk.
                // So append what we have and bail out now.
                partial.Extend(left, right, 0, length);
                goto done;
            }

            partial.Extend(left, right, 0, extra);
            ProcessChunk(partial);
            partial.Clear();
            offset += extra;
        }

        while (offset + ChunkSamples <= length)
        {
            ProcessChunk(Chunk(left, right, offset, ChunkSamples, position));
            offset += ChunkSamples;
        }

        // We always process chunks in units of ChunkSamples.
        // There may be a partial chunk left over from this iteration.
        if (offset < length)
            partial = Chunk(left, right, offset, length-offset, position);

    done:
        position += length;
    }

    void GlitchRemover::ProcessChunk(Chunk chunk)
    {
        using namespace std;

        // Called once for every new chunk.
        // If we are already in a glitch in either channel, check for end of glitch.
        // A glitch can end when we exceed maximum chunk count for a glitch,
        // or when it looks like the glitch is naturally over.
        // If we were not in a glitch and everything still looks fine,
        // just write this chunk to the output and keep going.

        if (PeakValue(chunk.left) > sampleLimit || PeakValue(chunk.right) > sampleLimit)
        {
            // This chunk is glitchy.
            // Tally each glitchy chunk we find in the graph, whether or not we end up removing it.
            graph.Increment(chunk.position - programStartPosition);

            if (chunklist.size() < MaxGlitchChunks)
            {
                if (chunklist.empty())
                    glitchStartSample = chunk.position;     // this is the very first glitchy chunk in a run

                // We are in a glitch. Keep saving bad chunks until
                // we decide whether to commit or cancel the glitch removal.
                chunklist.push_back(chunk);
            }
            else
            {
                // Too many glitches. Cancel the glitch removal.
                cout << "WARNING: Canceling removal of glitch at " << TimeStamp(glitchStartSample) << endl;
                Flush();
                lastGoodChunk = chunk;
            }
        }
        else
        {
            // This chunk is not glitchy.
            if (chunklist.empty())
            {
                // We are continuing a good section. This is the most common case. (We hope.)
                WriteChunk(lastGoodChunk);
                lastGoodChunk = chunk;
            }
            else
            {
                // We just found the end of a glitch that is short enough to remove.
                ++glitchCount;

                if (Verbose)
                {
                    cout << "GLITCH # "
                        << setw(4) << glitchCount
                        << " @ " << TimeStamp(glitchStartSample)
                        << " ==> removing " << chunklist.size()
                        << endl;
                }

                chunklist.clear();
                CrossFade(chunk);
            }
        }
    }

    void GlitchRemover::CrossFade(Chunk &chunk)
    {
        static_assert(CrossFadeSamples > 1, "CrossFadeSamples must be greater than 1");
        static_assert(CrossFadeSamples < ChunkSamples, "CrossFadeSamples must be less than ChunkSamples");

        if (chunk.Length() != ChunkSamples)
            throw Error("CrossFade: chunk has wrong length.");

        if (lastGoodChunk.Length() == 0)
        {
            lastGoodChunk = chunk;
            return;
        }

        if (lastGoodChunk.Length() != ChunkSamples)
            throw Error("CrossFade: lastGoodChunk has wrong length.");

        // Crossfade the front of 'chunk' into the tail of 'lastGoodChunk'.
        int k = ChunkSamples - CrossFadeSamples;
        for (int i=0; i < CrossFadeSamples; ++i)
        {
            float fadeIn = static_cast<float>(i) / static_cast<float>(CrossFadeSamples-1);
            float fadeOut = 1.0f - fadeIn;
            lastGoodChunk.left[i+k] = fadeOut*lastGoodChunk.left[i+k] + fadeIn*chunk.left[i];
            lastGoodChunk.right[i+k] = fadeOut*lastGoodChunk.right[i+k] + fadeIn*chunk.right[i];
        }

        // Write first chunk with crossfade at tail.
        WriteChunk(lastGoodChunk);

        // Prepend the part of 'chunk' after the crossfade to 'partial'.
        Chunk remainder(chunk.left, chunk.right, CrossFadeSamples, k, chunk.position);
        remainder.Extend(partial.left, partial.right, 0, partial.Length());
        partial = remainder;

        if (partial.Length() >= ChunkSamples)
        {
            // There is enough to form a new chunk, so that front part becomes the new lastGoodChunk.
            lastGoodChunk = Chunk(partial.left, partial.right, 0, ChunkSamples, partial.position);

            // Chop off the full chunk from the front of 'partial'.
            partial = Chunk(partial.left, partial.right, ChunkSamples, partial.Length() - ChunkSamples, partial.position);
        }
        else
        {
            // We don't have enough to create a new lastGoodChunk yet.
            // We will have to fill up partial later.
            lastGoodChunk.Clear();
        }
    }

    void GlitchRemover::WriteChunk(const Chunk &chunk)
    {
        using namespace std;

        float peak = chunk.Peak();
        if (peak > sampleLimit)
        {
            cout << "WARNING: Chunk at " << TimeStamp(chunk.position) << " has peak="
                << setprecision(5) << peak
                << " above limit=" << sampleLimit
                << endl;
        }

        if (peak > programPeak)
            programPeak = peak;

        gapFinder.Process(chunk, programStartPosition);
        writer.WriteChunk(chunk);
    }

    void GlitchRemover::Flush()
    {
        if (lastGoodChunk.Length() > 0)
        {
            WriteChunk(lastGoodChunk);
            lastGoodChunk.Clear();
        }

        // Write any chunks remaining in the chunklist.
        while (!chunklist.empty())
        {
            WriteChunk(chunklist.front());
            chunklist.pop_front();
        }

        // Write partial chunk if it contains any data.
        if (partial.Length() > 0)
        {
            WriteChunk(partial);
            partial.Clear();
        }
    }

    int GlitchRemover::GetInteger(const std::string& text, int start, int length)
    {
        int value = 0;

        for (int i=start; i < start + length; ++i)
        {
            if (text[i] < '0' || text[i] > '9')
                throw Error("Expected digit at offset " + std::to_string(i) + " in string '" + text + "'");

            value = (10 * value) + (text[i] - '0');
        }

        return value;
    }

    bool GlitchRemover::IsHeartsOfSpace(const std::string& filename)
    {
        using namespace std;

        // "Hearts of Space" programs are the fourth hour on a Saturday night.
        // Such a filename looks like: "cleaned-2018-03-17-4.au"
        //                              012345678901234567890123

        if (filename.length() != 23)
            throw Error("Unexpected output file format: " + filename);

        int year  = GetInteger(filename,  8, 4);
        int month = GetInteger(filename, 13, 2);
        int day   = GetInteger(filename, 16, 2);
        int hour  = GetInteger(filename, 19, 1);

        if (hour == 4)
        {
            // Is this date on a Saturday?
            tm time_in;
            memset(&time_in, 0, sizeof(time_in));
            time_in.tm_hour = 12;
            time_in.tm_year = year - 1900;
            time_in.tm_mon = month - 1;
            time_in.tm_mday = day;

            time_t time_temp = mktime(&time_in);

            const tm * time_out = localtime(&time_temp);
            cout << "HOS: day of week=" << time_out->tm_wday << " for file " << filename << endl;
            return time_out->tm_wday == 6;
        }

        return false;
    }

    void GlitchRemover::AdjustProgramLength(const std::string& filename, size_t programLengthSamples)
    {
        using namespace std;

        const double ToleranceSeconds = 5.0;
        const double IdealProgramMinutes = IsHeartsOfSpace(filename) ? 59.0 : 58.5;
        const double MinProgramMinutes = IdealProgramMinutes - (ToleranceSeconds / 60.0);  // never chop audio shorter than this many minutes

        double bestError = -1.0;
        const GapInfo *bestGap = nullptr;

        // Find the center-of-silence point that most closely fits a program length of 0:58:30.
        for (const GapInfo& gap : gapFinder.SilentGaps())
        {
            // Calculate how long the resulting audio would be if we chopped off the front here.
            double candidateLengthMinutes = MinutesFromSamples(programLengthSamples - gap.OriginalCenter());
            if (candidateLengthMinutes < MinProgramMinutes)
                break;

            // Find error from ideal program length
            double errorSeconds = 60.0 * abs(IdealProgramMinutes - candidateLengthMinutes);
            if (!bestGap || (errorSeconds < bestError))
            {
                bestGap = &gap;
                bestError = errorSeconds;
            }
        }

        if (bestGap)
        {
            cout << "ADJUST: Best candidate program length="
                << TimeStamp(programLengthSamples - bestGap->OriginalCenter())
                << ", error=" << setprecision(3) << bestError << " seconds."
                << endl;

            cout << "BESTGAP: offset=" << TimeStamp(bestGap->offset)
                << ", front=" << TimeStamp(bestGap->front)
                << ", length=" << TimeStamp(bestGap->length)
                << endl;

            cout << "POSTFILTER: Actual candidate length = "
                << TimeStamp(gapFinder.TotalSamples())
                << " - "
                << TimeStamp(bestGap->FilteredCenter())
                << " = "
                << TimeStamp(gapFinder.TotalSamples() - bestGap->FilteredCenter())
                << endl;

            if (bestError < ToleranceSeconds)
            {
                cout << "ADJUST: Deleting " << TimeStamp(bestGap->FilteredCenter())
                    << " of filtered data from front of file: "
                    << filename
                    << endl;

                AudioWriter::DeleteFrontSamples(filename, bestGap->FilteredCenter());
            }
            else
                cout << "ADJUST: Not adjusting program length because error is outside tolerance limit." << endl;
        }
        else
            cout << "ADJUST: Did not find a suitable silence gap for removing filler." << endl;

        gapFinder.Reset();
    }

    int GlitchRemover::ChunkListSampleCount() const
    {
        int sum = 0;
        for (const Chunk &chunk : chunklist)
            sum += chunk.Length();
        return sum;
    }

    const int GlitchGraph::HEIGHT_LIMIT = 30;

    void GlitchGraph::Reset()
    {
        tally.clear();
    }

    void GlitchGraph::Increment(long programSampleOffset)
    {
        // Convert sample offset to minutes (rounding down).
        size_t minutes = static_cast<size_t>(programSampleOffset / (60 * SamplingRate));
        const size_t MaxMinutes = 90;
        if (minutes > MaxMinutes)
            throw Error("GlitchGraph;:Increment -- time offset too large: minutes=" + std::to_string(minutes));

        if (minutes >= tally.size())
            tally.resize(1 + minutes);

        ++tally.at(minutes);
    }

    std::string GlitchGraph::Format(size_t duration) const
    {
        using namespace std;

        const size_t horizontal = duration + 1;

        int maxHeight = 0;
        for (size_t minute=0; minute < horizontal; ++minute)
        {
            int height = Height(minute);
            if (height > maxHeight)
                maxHeight = height;
        }

        const int vertical = maxHeight + 1;       // always include horizontal time axis "0----+----1---"

        vector<char> image(vertical * horizontal, ' ');

        for (size_t minute=0; minute < horizontal; ++minute)
        {
            int height = Height(minute);
            char marker = (height == HEIGHT_LIMIT) ? '@' : '|';
            image.at(minute) = MinuteTickMark(minute);
            for (int y=1; y <= height; ++y)
                image.at(minute + y*horizontal) = marker;
        }

        string text;
        for (int y = vertical-1; y >= 0; --y)
        {
            for (size_t x = 0; x < horizontal; ++x)
                text.push_back(image.at(x + y*horizontal));

            // Trim trailing whitespace.
            while (!text.empty() && (text.back() == ' '))
                text.pop_back();

            text.push_back('\n');
        }

        return text;
    }

    char GlitchGraph::MinuteTickMark(size_t minute)
    {
        if (minute % 10 == 0)
            return static_cast<char>('0' + ((minute / 10) % 10));

        if (minute % 5 == 0)
            return '+';

        return '-';
    }

    int GlitchGraph::Height(size_t minute) const
    {
        if (minute < tally.size())
            return std::min(GlitchGraph::HEIGHT_LIMIT, tally.at(minute));

        return 0;
    }

    void GapFinder::Process(const Chunk &chunk, long programStartPosition)
    {
        using namespace std;

        static const double MinSilenceSeconds  = 0.6;    // least duration we consider a silent period
        static const float  AmplitudeThreshold = 0.01f;  // how loud can we get before not considered silent

        const size_t MinSilenceSamples = SamplesFromSeconds(MinSilenceSeconds);

        const size_t buflen = chunk.left.size();
        const size_t front = chunk.position - programStartPosition;

        for (size_t i=0; i < buflen; ++i)
        {
            float mono = max(abs(chunk.left[i]), abs(chunk.right[i]));
            if (mono < AmplitudeThreshold)
            {
                if (silentSamples == 0)
                {
                    silenceOffset = totalSamples + i;
                    silenceFront = front + i;
                }
                ++silentSamples;
            }
            else
            {
                if (silentSamples >= MinSilenceSamples)
                {
                    cout << "GAP: offset=" << TimeStamp(silenceOffset)
                        << ", front=" << TimeStamp(silenceFront)
                        << ", length=" << TimeStamp(silentSamples) << endl;

                    gaplist.push_back(GapInfo(silenceOffset, silentSamples, silenceFront));
                }
                silentSamples = 0;
            }
        }

        totalSamples += buflen;
    }

    double ParseTimeOffset(const char *text)
    {
        using namespace std;

        // Parse a string of the format "h:mm:ss.mmm"
        //                               012345678901
        int hours, minutes;
        double seconds;
        if (3 != sscanf(text, "%d:%d:%lf", &hours, &minutes, &seconds))
            throw Error("Split point after '-s' must be formatted as h:mm:ss.sss");

        if (hours<0 || hours>4)
            throw Error("Hours value after '-s' must be 0..4");

        if (minutes<0 || minutes>59)
            throw Error("Minutes value after '-s' must be 0..59");

        if (seconds<0.0 || seconds>=60.0)
            throw Error("Seconds value after '-s' must be 0..59.99999");

        return seconds + 60.0*(minutes + 60*hours);
    }
}

int main(int argc, const char *argv[])
{
    using namespace std;
    using namespace unglitch;

    if (argc < 2)
    {
        cerr <<
            "USAGE:\n"
            "\n"
            "    unglitch projname [options...]\n"
            "        Given an Audacity project projname.aup, creates\n"
            "        a series of approximately hour-long audio files\n"
            "        that are normalized and have glitches removed.\n"
            "\n"
            "OPTIONS:\n"
            "\n"
            "-h dB\n"
            "    Specify desired headroom in decibels, where 0.0 <= dB <= 40.0.\n"
            "    Any samples that extend into the headroom are considered glitches.\n"
            "    The default value is " << HEADROOM_DEFAULT << " dB.\n"
            "\n"
            "-s h:mm:ss.mmm\n"
            "    Force program split at specified time index.\n"
            "\n"
            "-v\n"
            "    Generate verbose debug output.\n"
            "\n"
            << endl;

        return 1;
    }

    try
    {
        vector<double> manualSplitPointsInSeconds;
        double headroom = HEADROOM_DEFAULT;

        for (int i=2; i < argc; ++i)
        {
            const char *opt = argv[i];
            if (!strcmp(opt, "-v"))
                Verbose = true;
            else if (!strcmp(opt, "-s"))
            {
                if (i+1 < argc)
                    manualSplitPointsInSeconds.push_back(ParseTimeOffset(argv[++i]));
                else
                {
                    cerr << "ERROR: Missing time offset after '-s'" << endl;
                    return 1;
                }
            }
            else if (!strcmp(opt, "-h"))
            {
                if (i+1 < argc)
                {
                    const char *text = argv[++i];
                    if (1 != sscanf(text, "%lf", &headroom))
                    {
                        cerr << "ERROR: Invalid headroom '" << text << "' after '-h'" << endl;
                        return 1;
                    }

                    if (headroom<0.0 || headroom>40.0)
                    {
                        cerr << "ERROR: Headroom value is out of range." << endl;
                        return 1;
                    }
                }
                else
                {
                    cerr << "ERROR: Missing headroom value after '-h'" << endl;
                    return 1;
                }
            }
            else
            {
                cerr << "ERROR: Unknown option '" << opt << "'" << endl;
                return 1;
            }
        }

        string projname(argv[1]);
        string inAudacityProjectFileName = projname + ".aup";
        Project project(manualSplitPointsInSeconds);
        project.Load(inAudacityProjectFileName.c_str());
        cout << "Loaded project: " << inAudacityProjectFileName << endl;
        project.Convert(projname, headroom);
        cout << "Finished converting audio." << endl;
    }
    catch (const Error &error)
    {
        cerr << "ERROR: " << error.Message() << endl;
        return 1;
    }

    return 0;
}
