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
    }

    void Project::Load(const char *filename)
    {
        using namespace std;
        using namespace tinyxml2;

        XMLDocument doc;
        XMLError status = doc.LoadFile(filename);
        if (status != XML_SUCCESS)
            throw Error(string("Cannot open Audacity project file ") + filename);

        cout << "Loaded xml: " << filename << endl;

        XMLElement *root = doc.RootElement();
        for (XMLElement *trackElem = root->FirstChildElement("wavetrack"); trackElem; trackElem = trackElem->NextSiblingElement("wavetrack"))  
        {
            cout << "Found wavetrack" << endl;
            channelList.push_back(WaveTrack());
            WaveTrack &track = channelList[channelList.size() - 1];
            track.Parse(trackElem);
        }
    }
}

int main(int argc, const char *argv[])
{
    using namespace std;
    using namespace unglitch;

    if (argc != 2)
    {
        cerr << "USAGE: unglitch infile.aup" << endl;
        return 1;
    }

    const char *inAudacityProjectFileName = argv[1];
    try
    {
        Project project;
        project.Load(inAudacityProjectFileName);
        cout << "Loaded project" << endl;
    }
    catch (const Error &error)
    {
        cerr << "ERROR: " << error.Message() << endl;
        return 1;
    }

    return 0;
}
