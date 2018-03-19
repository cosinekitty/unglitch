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

    void Project::PrintProgramSummary(const std::string& filename, const GlitchRemover &remover)
    {
        using namespace std;

        cout << filename << " : " << remover.GlitchCount() << " glitches, peak="
            << setprecision(5) << remover.ProgramPeak() 
            << ", headroom=" << setprecision(1) << remover.HeadroomDecibels() << " dB." 
            << endl;

        cout << remover.FormatGlitchGraph() << endl;
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

                PrintProgramSummary(writer.OutFileName(), remover);                    
                remover.ResetProgram();

                cout << "Splitting program at " << TimeStamp(position + boundary) << endl;
                writer.StartNewFile(OutProgramFileName(outFilePrefix, hour));
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
            programPosition += blockLength;
        }

        remover.Flush();
        PrintProgramSummary(writer.OutFileName(), remover);                    
        remover.ResetProgram();
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

    bool GlitchRemover::IsGlitch(float prevPeak, float newPeak, float otherPeak)
    {
        const float Threshold = 0.27f;
        const float BadJump = 0.08f;

        return 
            (newPeak > sampleLimit) ||
            ((newPeak > Threshold) && (newPeak - prevPeak > BadJump) && (newPeak - otherPeak > BadJump));
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

        float newPeakLeft  = PeakValue(chunk.left);
        float newPeakRight = PeakValue(chunk.right);

        if (IsGlitch(prevPeakLeft, newPeakLeft, newPeakRight) || IsGlitch(prevPeakRight, newPeakRight, newPeakLeft))
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
                prevPeakLeft = newPeakLeft;
                prevPeakRight = newPeakRight;
            }
        }
        else
        {
            // This chunk is not glitchy.

            prevPeakLeft = newPeakLeft;
            prevPeakRight = newPeakRight;

            if (chunklist.empty())
            {
                // We are continuing a good section. This is the most common case. (We hope.)
                WriteChunk(lastGoodChunk);
                lastGoodChunk = chunk;
            }
            else
            {
                // We just found the end of a glitch that is short enough to remove.
                //cout << "Removing glitch at " << TimeStamp(glitchStartSample) << endl;
                ++glitchCount;
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
        Chunk remainder(chunk.left, chunk.right, CrossFadeSamples, chunk.Length() - CrossFadeSamples, chunk.position);
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

    int GlitchRemover::ChunkListSampleCount() const
    {
        int sum = 0;
        for (const Chunk &chunk : chunklist)
            sum += chunk.Length();
        return sum;
    }

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

    std::string GlitchGraph::Format() const
    {
        using namespace std;

        size_t duration = tally.size();

        int maxHeight = 0;
        for (size_t minute=0; minute < duration; ++minute)
        {
            int height = Height(minute);
            if (height > maxHeight)
                maxHeight = height;
        }

        int horizontal = duration;
        int vertical = maxHeight + 1;       // always include horizontal time axis "0----+----1---"

        vector<char> image(vertical * horizontal, ' ');

        for (size_t minute=0; minute < duration; ++minute)
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
            for (int x = 0; x < horizontal; ++x)
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
            return std::min(HEIGHT_LIMIT, tally.at(minute));

        return 0;
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
