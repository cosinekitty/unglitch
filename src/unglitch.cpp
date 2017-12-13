/*
    unglitch.cpp  -  Don Cross  -  13 December 2017.

    Removes certain glitches that appear when I record radio programs.
*/

#include <iostream>
#include "tinyxml2.h"

int main(int argc, const char *argv[])
{
    using namespace std;
    using namespace tinyxml2;

    int rc = 1;

    if (argc != 2)
    {
        cerr << "USAGE: unglitch infile.aup" << endl;
        return 1;
    }

    const char *inAudacityProjectFileName = argv[1];
    XMLDocument doc;
    XMLError status = doc.LoadFile(inAudacityProjectFileName);
    if (status != XML_SUCCESS)
    {
        cerr << "ERROR " << status << " trying to open Audacity project file '" << inAudacityProjectFileName << "'" << endl;
        return 1;
    }

    cout << "Loaded xml: " << inAudacityProjectFileName << endl;

    XMLElement *root = doc.RootElement();
    XMLElement *trackElem = root->FirstChildElement("wavetrack");
    while (trackElem)
    {
        cout << "Found wavetrack" << endl;
        trackElem = trackElem->NextSiblingElement("wavetrack");
    }

    return rc;
}
