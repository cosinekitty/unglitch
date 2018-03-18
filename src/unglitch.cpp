/*
    unglitch.cpp  -  Don Cross  -  13 December 2017.

    Removes certain glitches that appear when I record radio programs.
*/

#include "unglitch.h"

#define DEBUG_DUMP_HISTOGRAMS 0

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
    const int SamplingRate = 44100;

    std::string TimeStamp(long offset);

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

    const int NUM_FOLDED_BINS = 1000;       // 3 decimal places of resolution in the range [0.0, +1.0].
    const int NUM_HISTOGRAM_BINS = 2 * NUM_FOLDED_BINS;

    void HistogramTally(std::vector<int> &histogram, const FloatVector &buffer)
    {
        for (float data : buffer)
        {
            if (data < -1.0 || data > +1.0)
                throw Error("HistogramTally: data out of range.");

            int bin = std::floor(((data + 1.0) / 2.0) * (NUM_HISTOGRAM_BINS - 1));
            if (bin < 0 || bin >= NUM_HISTOGRAM_BINS)
                throw Error("HistogramTally: bin out of range.");

            ++histogram.at(bin);
        }
    }

#if DEBUG_DUMP_HISTOGRAMS
    void WriteHistogram(
        const char *outCsvFileName,
        const std::vector<int> &histogram)
    {
        FILE *outfile = fopen(outCsvFileName, "wt");
        if (!outfile)
            throw Error("Cannot open histogram output csv file.");

        fprintf(outfile,"\"floor\",\"count\",\"fraction\"\n");

        int population = 0;
        int front = -1;
        int back = -1;
        for (int i=0; i < NUM_HISTOGRAM_BINS; ++i)
        {
            int count = histogram.at(i);
            if (count != 0)
            {
                population += count;
                back = i;
                if (front < 0)
                    front = i;
            }
        }

        if (back >= 0)
        {
            double delta = 2.0 / NUM_HISTOGRAM_BINS;
            for (int i = front; i <= back; ++i)
            {
                double low = -1.0 + (delta * i);
                int count = histogram.at(i);
                double fraction = static_cast<double>(count) / static_cast<double>(population);
                fprintf(outfile, "%0.3lf,%d,%0.15lf\n", low, count, fraction);
            }
        }

        fclose(outfile);
    }
#endif

    int Project::FoldTally(
        std::vector<int>& folded, 
        const std::vector<int>& histogram,
        double bias)
    {
        int population = 0;
        double delta = 2.0 / NUM_HISTOGRAM_BINS;
        for (int i=0; i < NUM_HISTOGRAM_BINS; ++i)
        {
            int count = histogram.at(i);
            if (count > 0)
            {
                // Calculate the median value of the raw bin.
                double median = -1.0 + (delta * (i + 0.5));

                // Subtract out the DC bias to get the corrected median value
                // and take absolute value.
                double fvalue = std::abs(median - bias);

                // convert to the folded bin index
                int bin = static_cast<int>(std::floor(0.5 + fvalue*(NUM_FOLDED_BINS - 1)));

                if (bin < 0 || bin >= NUM_FOLDED_BINS)
                    throw Error("FoldTally: bin is out of bounds: bin=" + std::to_string(bin) + ", i=" + std::to_string(i));

                population += count;
                folded.at(bin) += count;
            }
        }

        return population;
    }

    double Project::FindLimit(
        const std::vector<int> & leftHistogram, 
        double leftBias, 
        const std::vector<int> & rightHistogram, 
        double rightBias)
    {
        // The input histograms for left and right channels span the range -1.0 to +1.0,
        // and refer to raw values that still include DC bias.
        // Subtract the respective biases from both histograms and fold them to
        // absolute value range 0.0 to +1.0. This makes it easier to determine the
        // absolute glitch cutoff that keeps the vast majority of sample values.
        std::vector<int> folded(NUM_FOLDED_BINS);
        int population = FoldTally(folded, leftHistogram, leftBias);
        population += FoldTally(folded, rightHistogram, rightBias);

        const int DENOMINATOR = 800000;     // denominator of fraction of extreme samples to discard
        int keep = population - (population / DENOMINATOR);
        int total = 0;
        for (int i=0; i < NUM_FOLDED_BINS; ++i)
        {
            total += folded.at(i);
            if (total >= keep)
                return static_cast<double>(i) / static_cast<double>(NUM_FOLDED_BINS);
        }

        throw Error("FindLimit: Could not find keep limit.");
    }

    Project::DcBias Project::FindBias() const
    {
        // Measure the DC offset in both channels.
        if (channelList.size() != 2)
            throw Error("Input must be stereo.");

        const WaveTrack &leftTrack = channelList[0];
        const WaveTrack &rightTrack = channelList[1];

        const int nblocks = leftTrack.NumBlocks();
        if (nblocks != rightTrack.NumBlocks())
            throw Error("Left and right tracks have different number of blocks.");

        std::vector<int> leftHistogram(NUM_HISTOGRAM_BINS);
        std::vector<int> rightHistogram(NUM_HISTOGRAM_BINS);

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

            HistogramTally(leftHistogram, leftBuffer);
            HistogramTally(rightHistogram, rightBuffer);

            position += length;
        }

        double leftBias = leftSum / position;
        double rightBias = rightSum / position;

#if DEBUG_DUMP_HISTOGRAMS
        WriteHistogram("left.csv", leftHistogram);
        WriteHistogram("right.csv", rightHistogram);
#endif

        double limit = FindLimit(leftHistogram, leftBias, rightHistogram, rightBias);

        return DcBias(leftBias, rightBias, limit);
    }

    std::string Project::OutProgramFileName(std::string prefix, int hour)
    {
        return prefix + "-" + std::to_string(hour) + ".au";
    }

    void Project::Convert(std::string outFilePrefix)
    {
        using namespace std;

        cout << "Finding DC bias..." << endl;
        DcBias bias = FindBias();
        double headroom_dB = -20.0 * log10(bias.limit);
        cout << fixed 
            << setprecision(5) << "Bias: left=" << bias.left << ", right=" << bias.right 
            << ", limit=" << bias.limit 
            << setprecision(1) << ", headroom=" << headroom_dB << " dB"
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

        const long nsamples = leftTrack.NumSamples();
        if (nsamples != rightTrack.NumSamples())
            throw Error("Left and right tracks have different number of samples.");

        const long MinSplitSamples = 60L * SamplingRate;  // do not split within final minute of audio

        int hour = 1;
        AudioWriter writer(OutProgramFileName(outFilePrefix, hour), SamplingRate, 2);

        GlitchRemover remover(writer, bias.limit);

        FloatVector leftBuffer;
        FloatVector rightBuffer;
        bool skippingFiller = false;    // are we skipping commercials and station ID stuff between programs?
        FloatVector leftFiller;
        FloatVector rightFiller;
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
                data -= bias.left;

            for (float &data : rightBuffer)
                data -= bias.right;

            if ((nsamples-position > MinSplitSamples) && IsStartingNextProgram(hour, programPosition, leftBuffer, rightBuffer, boundary))
            {
                programPosition = 0;    // program position includes any skipped filler at front; needed by IsStartingNextProgram()
                ++hour;
                skippingFiller = true;

                // Split buffers into before/after the boundary.
                FloatVector leftBefore = SplitBuffer(leftBuffer, boundary);
                FloatVector rightBefore = SplitBuffer(rightBuffer, boundary);
                remover.Fix(leftBefore, rightBefore);
                remover.Flush();
                cout << "Removed " << remover.GlitchCount() << " glitches from " << writer.OutFileName() << endl;
                remover.ResetGlitchCount();

                cout << "Splitting program at " << TimeStamp(position + boundary) << endl;
                writer.StartNewFile(OutProgramFileName(outFilePrefix, hour));
            }
            else
            {
                programPosition += blockLength;     // FIXFIXFIX - should always happen at bottom of loop, but need to fix IsStartingNewProgram time offset thresholds
            }

            if (skippingFiller)
            {
                // Detect end of filler. There should be at between 1.0 and 2.5 minutes of filler.
                // Filler ends with WMFE announcer speaking "... and Daytona Beach" or "... for all devices".
                // Detect the generic pattern: [speech, silence, music] after at least 1.0 minutes.
                // Buffer up all suspected filler so that if we can't find the transition between filler
                // and program, we give up and dump all the filler to the output and keep going.
                const double fillerMinutes = leftFiller.size() / (SamplingRate * 60.0);
                if (fillerMinutes > 2.5)
                {
                    cout << "Could not determine end of filler. Keeping entire beginning of program." << endl;
                    remover.Fix(leftFiller, rightFiller);
                    leftFiller.clear();
                    rightFiller.clear();
                    skippingFiller = false;
                }
                else
                {
                    Append(leftFiller, leftBuffer);
                    Append(rightFiller, rightBuffer);
                }
            }

            if (!skippingFiller)
                remover.Fix(leftBuffer, rightBuffer);

            position += blockLength;
        }

        remover.Flush();
        cout << "Removed " << remover.GlitchCount() << " glitches from " << writer.OutFileName() << endl;
        remover.ResetGlitchCount();
    }

    bool Project::IsStartingNextProgram(
        int hour,
        long programPosition, 
        const FloatVector& leftBuffer, 
        const FloatVector& rightBuffer, 
        long &boundary) const
    {
        using namespace std;

        boundary = 0;

        if (leftBuffer.size() != rightBuffer.size())
            throw Error("Buffers different lengths in IsStartingNextProgram");

        const long length = static_cast<long>(leftBuffer.size());

        // We are starting a new program if we find a silent period
        // somewhere between 58 and 61 minutes into the current program.
        const double minProgramMinutes = (hour==1) ? 58.25 : 59.25;
        const double maxProgramMinutes = 2.0 + minProgramMinutes;
        const double samplesPerSecond = static_cast<double>(SamplingRate);
        const double samplesPerMinute = samplesPerSecond * 60.0;
        double minutesBegin = programPosition / samplesPerMinute;
        double minutesEnd = minutesBegin + (length / samplesPerMinute);
        if (Overlap(minutesBegin, minutesEnd, minProgramMinutes, maxProgramMinutes))
        {
            const double minSilenceSeconds = 0.4;       // least duration we consider a silent period
            const float amplitudeThreshold = 0.01;      // how loud can we get before not considered silent
            bool inSilence = false;
            for (long i=0; i < length; ++i)
            {
                float height = max(abs(leftBuffer[i]), abs(rightBuffer[i]));
                if (height < amplitudeThreshold)
                {
                    if (!inSilence)
                    {
                        boundary = i;
                        inSilence = true;
                    }
                }
                else
                {
                    if (inSilence)
                    {
                        inSilence = false;
                        long gap = i - boundary;
                        double silentSeconds = gap / samplesPerSecond;
                        if (silentSeconds >= minSilenceSeconds)
                        {
                            boundary += gap/2;   // split right in middle of silence
                            return true;
                        }
                    }
                }
            }
            if (inSilence)
            {
                long gap = length - boundary;
                double silentSeconds = gap / samplesPerSecond;
                if (silentSeconds >= minSilenceSeconds)
                {
                    boundary += gap/2;   // split right in middle of silence
                    return true;
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
        if (offset<0 || offset>=length)
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
        if (size < dataOffset)
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
        long millis = (offset * 1000L) / SamplingRate;
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
                ++glitchCount;
                //cout << glitchCount << ". Discarding " << ChunkListSampleCount() << " samples at " << TimeStamp(glitchStartSample) << endl;
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

        const float Threshold = 0.27f;
        const float BadJump = 0.08f;
        float peak1 = PeakValue(first);
        float peak2 = PeakValue(second);
        bool inglitch = (peak1 > Threshold) && (peak1 - state.prevPeak > BadJump) && (peak1 - peak2 > BadJump);

        inglitch = inglitch || (std::max(peak1, peak2) > sampleLimit);

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
                    std::cout << "Glitch too long at " << TimeStamp(glitchStartSample) << std::endl;
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

    if (argc != 2)
    {
        cerr << 
            "USAGE:\n"
            "\n"
            "    unglitch projname\n"
            "        Given an Audacity project projname.aup, creates\n"
            "        a series of approximately hour-long audio files\n"
            "        that are normalized and have glitches removed.\n"
            << endl;

        return 1;
    }

    string projname(argv[1]);
    string inAudacityProjectFileName = projname + ".aup";
    string outAudioPrefix = string("cleaned-") + projname;

    try
    {
        Project project;
        project.Load(inAudacityProjectFileName.c_str());
        cout << "Loaded project: " << inAudacityProjectFileName << endl;
        project.Convert(outAudioPrefix);
        cout << "Finished converting audio." << endl;
    }
    catch (const Error &error)
    {
        cerr << "ERROR: " << error.Message() << endl;
        return 1;
    }

    return 0;
}
