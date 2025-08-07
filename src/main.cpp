// Copyright 2025 Giorgio Gamba

#include <stdlib.h>
#include <iostream>
#include <cmath>

#include <portaudio.h>

void checkError(const PaError& error)
{
    if (error != paNoError)
    {
        std::cout << "[PortAudio Error] " << Pa_GetErrorText(error) << std::endl;
        exit(EXIT_FAILURE);
    }
}

/* Performs the action following the streaming capture */
int paCallback ( const void* inputBuffer
                , void* outputBuffer
                , unsigned long framesPerBuffer
                , const PaStreamCallbackTimeInfo* timeInfo
                , PaStreamCallbackFlags statusFlags
                , void* userData )
{
    const float* in = (float*) inputBuffer;
    (void) outputBuffer;

    constexpr int displaySize = 50;
    std::cout << "\r";

    // Gets greatest volume in the buffer for L and R channels
    float volume_L = 0;
    float volume_R = 0;
    // Parse over elements considering them in couples L, R
    for (unsigned long i = 0; i < framesPerBuffer * 2; i += 2)
    {
        volume_L = std::max(volume_L, std::abs(in[i]));
        volume_R = std::max(volume_R, std::abs(in[i+1]));
    }

    // Prints resulting volume bars
    for (int i = 0; i < displaySize; ++i)
    {
        const float bar = i / (float) displaySize;
        const std::string volumeBar = (bar <= volume_L && bar <= volume_R) ? "|" : " ";
        std::cout << volumeBar;
    }

    fflush(stdout);

    return 0;
}

int main()
{
    // Every time a Pa operation is performed, its value is checked in order
    // to spot problems
    PaError error;

    error = Pa_Initialize();
    checkError(error);

    const int numDevices = Pa_GetDeviceCount();
    std::cout << "Number of devices: " << numDevices << std::endl;

    for (int i = 0; i < numDevices; ++i)
    {
        const PaDeviceInfo* deviceInfo = Pa_GetDeviceInfo(i);
        std::cout << "Name: " << deviceInfo->name << std::endl;
        std::cout << "Sample rate: " << deviceInfo->defaultSampleRate << std::endl;
        std::cout << "Low input latency: " << deviceInfo->defaultLowInputLatency << std::endl;
        std::cout << "---------"  << std::endl;
    }

    // Configures stream
    constexpr int deviceIdx = 0;

    std::cout << "Input stream configuration..." << std::endl;
    PaStreamParameters inStreamParams;
    memset(&inStreamParams, 0, sizeof(inStreamParams));
    inStreamParams.channelCount = 2;
    inStreamParams.device = deviceIdx;
    inStreamParams.hostApiSpecificStreamInfo = nullptr;
    inStreamParams.sampleFormat = paFloat32;
    inStreamParams.suggestedLatency = Pa_GetDeviceInfo(deviceIdx)->defaultLowInputLatency;

    // Open stream
    PaStream* stream;
    error = Pa_OpenStream(&stream, &inStreamParams, nullptr, 44100.0, 512, paNoFlag, paCallback, nullptr);
    checkError(error);

    std::cout << "Start streaming..." << std::endl;
    error = Pa_StartStream(stream);
    checkError(error);

    // Make the application continue execution until stopping with Ctrl-C
    // #TODO add ctrlC capture for proper streaming closing
    while (true)
    {
        Pa_Sleep(1);
    }

    // Close stream
    error = Pa_StopStream(stream);
    checkError(error);

    error = Pa_CloseStream(stream);
    checkError(error);

    error = Pa_Terminate();
    checkError(error);

    return EXIT_SUCCESS;
}
