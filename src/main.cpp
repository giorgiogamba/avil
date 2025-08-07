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

int paCallback ( const void* inputBuffer
                , void* outputBuffer
                , unsigned long framesPerBuffer
                , const PaStreamCallbackTimeInfo* timeInfo
                , PaStreamCallbackFlags statusFlags
                , void* userData )
{
    float* in = (float*) inputBuffer;
    (void) outputBuffer;

    int displaySize = 100;
    std::cout << "\r";

    float volume_L = 0;
    float volume_R = 0;

    // Parse over elements considering them in couples L, R
    for (unsigned long i = 0; i < framesPerBuffer * 2; i += 2)
    {
        volume_L = std::max(volume_L, std::abs(in[i]));
        volume_R = std::max(volume_R, std::abs(in[i+1]));
    }

    for (int i = 0; i < displaySize; ++i)
    {
        const float bar = i / (float) displaySize;
        if (bar <= volume_L && bar <= volume_R)
        {
            std::cout << "||";
        }
        else if (bar <= volume_L)
        {
            std::cout << "'";
        }
        else if (bar <= volume_R)
        {
            std::cout << ",";
        }
        else
        {
            std::cout << " ";
        }
    }

    fflush(stdout);

    return 0;
}

int main()
{
    PaError error;
    error = Pa_Initialize();
    checkError(error);

    const int numDevices = Pa_GetDeviceCount();
    std::cout << "Number of devices: " << numDevices << std::endl;

    const PaDeviceInfo* deviceInfo = nullptr;
    for (int i = 0; i < numDevices; ++i)
    {
        deviceInfo = Pa_GetDeviceInfo(i);
        /* #TODO print info about device */
    }

    // Creates stream

    const int deviceIdx = 0;

    PaStreamParameters inStreamParams;
    memset(&inStreamParams, 0, sizeof(inStreamParams));
    inStreamParams.channelCount = 2;
    inStreamParams.device = deviceIdx;
    inStreamParams.hostApiSpecificStreamInfo = nullptr;
    inStreamParams.sampleFormat = paFloat32;
    inStreamParams.suggestedLatency = Pa_GetDeviceInfo(deviceIdx)->defaultLowInputLatency;

    PaStreamParameters outStreamParams;
    memset(&outStreamParams, 0, sizeof(outStreamParams));
    outStreamParams.channelCount = 2;
    outStreamParams.device = deviceIdx;
    outStreamParams.hostApiSpecificStreamInfo = nullptr;
    outStreamParams.sampleFormat = paFloat32;
    outStreamParams.suggestedLatency = Pa_GetDeviceInfo(deviceIdx)->defaultLowInputLatency;

    PaStream* stream;
    error = Pa_OpenStream(&stream, &inStreamParams, &outStreamParams, 44100, 512, paNoFlag, paCallback, nullptr);
    checkError(error);

    error = Pa_StartStream(stream);
    checkError(error);

    Pa_Sleep(10 * 1000);

    error = Pa_StopStream(stream);
    checkError(error);

    error = Pa_CloseStream(stream);
    checkError(error);

    error = Pa_Terminate();
    checkError(error);

    return EXIT_SUCCESS;
}
