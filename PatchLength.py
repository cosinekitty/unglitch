#!/usr/bin/env python
#
#   PatchLength.py
#
import sys

def ReadField(file):
    data = file.read(4)
    return ord(data[0]) | (ord(data[1]) << 8) | (ord(data[2]) << 16) | (ord(data[3]) << 24)

def EncodeInt(x):
    return chr(x&0xff) + chr((x>>8)&0xff) + chr((x>>16)&0xff) + chr((x>>24)&0xff)

def Patch(inFileName, secondsToSkip):
    with open(inFileName, 'r+b') as file:
        prefix = file.read(4)
        if prefix != 'dns.':
            raise Exception('File is not .au format!')
        dataOffset = ReadField(file)
        file.seek(12)
        encoding = ReadField(file)
        if encoding != 6:
            raise Exception('Unsupported .au format {0}'.format(encoding))
        rate = ReadField(file)
        if rate != 44100:
            raise Exception('Unsupported sampling rate {0}'.format(rate))
        channels = ReadField(file)
        if channels != 2:
            raise Exception('Unsupported number of channels: {0}'.format(channels))
        
        # Adjust the data offset by the specified number of samples.
        newDataOffset = dataOffset + int(round(4 * channels * secondsToSkip * rate))
        print 'old data offset = {0}'.format(dataOffset)
        print 'new data offset = {0}'.format(newDataOffset)
        file.seek(4)
        file.write(EncodeInt(newDataOffset))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'USAGE: PatchLength.py infile.au seconds_to_skip'
        sys.exit(1)

    inFileName = sys.argv[1]
    secondsToSkip = float(sys.argv[2])
    Patch(inFileName, secondsToSkip)