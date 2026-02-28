// Copyright 2025 Giorgio Gamba

#include <portaudio.h>
#include <fftw3.h>
#include <sndfile.h>

struct StreamCallbackData
{
	int startIndex;
	int sprectrogramSize;
};

// Represents the loaded
struct FileCallbackData
{
    int startIndex;
    int sprectrogramSize;
    SNDFILE* file;
    SF_INFO info;
};
