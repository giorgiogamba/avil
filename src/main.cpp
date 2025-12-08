// Copyright 2025 Giorgio Gamba

#include <stdlib.h>
#include <iostream>
#include <cmath>

#include <portaudio.h>
#include <fftw3.h>
#include <sndfile.h>

#pragma region CONSTANTS

constexpr double SAMPLE_RATE{44100.0};
constexpr int FRAMES_PER_BUFFER{512};
constexpr int NUM_CHANNELS{2};

constexpr int DISPLAY_SIZE{50};

// Defines spectrogram's boundaries
constexpr int SPECTROGRAM_FREQ_START{20};
constexpr int SPECTROGRAM_FREQ_END{20 * 1000};

#pragma endregion

#pragma region Typed

struct StreamCallbackData
{
	double* in;
	double* out;
	fftw_plan p;
	int startIndex;
	int sprectrogramSize;
};
static StreamCallbackData* spectrogramData;

// Represents the loaded
struct FileCallbackData
{
    SNDFILE* file;
    SF_INFO info;
    	double* in;
	double* out;
	fftw_plan p;
	int startIndex;
	int sprectrogramSize;
};
static FileCallbackData* fileSpectrogramData;

#pragma endregion

void checkError(const PaError& error)
{
    if (error != paNoError)
    {
        std::cout << "[PortAudio Error] " << Pa_GetErrorText(error) << std::endl;
        exit(EXIT_FAILURE);
    }
}

void printVolumeGraph(const float* in, unsigned long framesPerBuffer)
{
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
    for (int i = 0; i < DISPLAY_SIZE; ++i)
    {
        const float bar = i / static_cast<float>(DISPLAY_SIZE);
        const std::string volumeBar = (bar <= volume_L && bar <= volume_R) ? "|" : " ";
        std::cout << volumeBar;
    }
}

void printFrequencyGraph(const float* in, unsigned long framesPerBuffer, void* userData)
{
    StreamCallbackData* streamData = reinterpret_cast<StreamCallbackData*>(userData);

    // Collects data for Fourier transform display
    for (auto i{0}; i < framesPerBuffer; ++i)
    {
        // Audio channels samples are set in a matrix where row i contains i-th samples form all the audio channels
        // We own just 2 channels, so we get the even samples
        streamData->in[i] = in[i * NUM_CHANNELS];
    }

    // Executes fourier transform on the stremed data
    fftw_execute(streamData->p);

    for (int i{0}; i < DISPLAY_SIZE; ++i)
    {
        const double step = i / static_cast<double>(DISPLAY_SIZE);
        const auto outIndex = static_cast<int>(streamData->startIndex + step * streamData->sprectrogramSize);
        const auto freq = streamData->out[outIndex];

        if (freq < 0.125)
        {
            printf("▁");
        }
        else if (freq < 0.25)
        {
            printf("▂");
        }
        else if (freq < 0.375)
        {
            printf("▃");
        }
        else if (freq < 0.5)
        {
            printf("▄");
        }
        else if (freq < 0.625)
        {
            printf("▅");
        }
        else if (freq < 0.75)
        {
            printf("▆");
        }
        else if (freq < 0.925)
        {
            printf("▇");
        }
        else
        {
            printf("█");
        }
    }
}

/* Performs the action following the streaming capture */
int streamCallback ( const void* inputBuffer
                , void* outputBuffer
                , unsigned long framesPerBuffer
                , const PaStreamCallbackTimeInfo* timeInfo
                , PaStreamCallbackFlags statusFlags
                , void* userData )
{
    (void) outputBuffer;
    float* out = nullptr;
    out = (float*) outputBuffer;
    FileCallbackData* data = (FileCallbackData*)userData;
    memset(out, 0, sizeof(float) * framesPerBuffer * data->info.channels);

    // Transposes read data into the output buffer
    sf_count_t numRead;
    numRead = sf_read_float(data->file, out, framesPerBuffer * data->info.channels);

    printVolumeGraph(out, framesPerBuffer);
    // const float* in = (float*) inputBuffer; // #TODO move to static cast
    // printFrequencyGraph(in, framesPerBuffer, userData);
    fflush(stdout);

    const bool fileEnded = numRead < framesPerBuffer;
    return fileEnded ? paComplete : paContinue;
}

int main(int argc, const char* argv[])
{
    // Reads cli arguments
    if (argc < 2)
    {
        std::cout << "Error while initializing application. Not enough arguments provided\n";
        return EXIT_FAILURE;
    }

    const char* filePath = argv[1];
    
    FileCallbackData data{};
    data.file = sf_open(filePath, SFM_READ, &data.info);
    if (sf_error(data.file) != SF_ERR_NO_ERROR)
    {
        std::cout << "An error occured while opening file " << filePath << "\n";
        std::cout << sf_strerror(data.file);
        return EXIT_FAILURE;
    }

    // Every time a Pa operation is performed, its value is checked in order
    // to spot problems
    PaError error;

    error = Pa_Initialize();
    checkError(error);

    PaStream* stream;
    constexpr int noInput{0};
    error = Pa_OpenDefaultStream(&stream, noInput, data.info.channels, paFloat32, data.info.samplerate, FRAMES_PER_BUFFER, streamCallback, &data);
    checkError(error);

    std::cout << "Starting stream...\n";
    error = Pa_StartStream(stream);
    checkError(error);

    // Plays until the file continues
    while (Pa_IsStreamActive(stream))
    {
        Pa_Sleep(100);
    }

	// spectrogramData = (StreamCallbackData*)malloc(sizeof(StreamCallbackData));
	// spectrogramData->in = (double*)malloc(FRAMES_PER_BUFFER * sizeof(double));
	// spectrogramData->out = (double*)malloc(FRAMES_PER_BUFFER * sizeof(double));

	// if (spectrogramData->in == nullptr || spectrogramData->out == nullptr)
    // {
    //     std::cout << "Could not allocate spectrogram data\n";
    //     exit(EXIT_FAILURE);
    // }

    // const double sampleRatio = FRAMES_PER_BUFFER / static_cast<double>(SAMPLE_RATE);
    // spectrogramData->startIndex = std::ceil(sampleRatio * SPECTROGRAM_FREQ_START); // rounds to biggernumber for very small numbers
    // constexpr double DEF_SIZE{FRAMES_PER_BUFFER / 2.0};
    // spectrogramData->sprectrogramSize = std::min(std::ceil(sampleRatio * SPECTROGRAM_FREQ_END), DEF_SIZE) - spectrogramData->startIndex;

    // // Defines the Fourier transform. Data need to remain the same as long as this profile is chosen
    // spectrogramData->p = fftw_plan_r2r_1d(FRAMES_PER_BUFFER, spectrogramData->in, spectrogramData->out, FFTW_R2HC, FFTW_ESTIMATE);

    // const int numDevices = Pa_GetDeviceCount();
    // std::cout << "Number of devices: " << numDevices << std::endl;

    // for (int i = 0; i < numDevices; ++i)
    // {
    //     const PaDeviceInfo* deviceInfo = Pa_GetDeviceInfo(i);
    //     std::cout << "Name: " << deviceInfo->name << std::endl;
    //     std::cout << "Sample rate: " << deviceInfo->defaultSampleRate << std::endl;
    //     std::cout << "Low input latency: " << deviceInfo->defaultLowInputLatency << std::endl;
    //     std::cout << "---------"  << std::endl;
    // }

    // // Configures stream
    // constexpr int deviceIdx = 0;

    // std::cout << "Input stream configuration..." << std::endl;
    // PaStreamParameters inStreamParams;
    // memset(&inStreamParams, 0, sizeof(inStreamParams));
    // inStreamParams.channelCount = NUM_CHANNELS;
    // inStreamParams.device = deviceIdx;
    // inStreamParams.hostApiSpecificStreamInfo = nullptr;
    // inStreamParams.sampleFormat = paFloat32;
    // inStreamParams.suggestedLatency = Pa_GetDeviceInfo(deviceIdx)->defaultLowInputLatency;

    // // Open stream
    // PaStream* stream;
    // error = Pa_OpenStream(&stream, &inStreamParams, nullptr, 44100.0, 512, paNoFlag, streamCallback, spectrogramData);
    // checkError(error);

    // std::cout << "Start streaming..." << std::endl;
    // error = Pa_StartStream(stream);
    // checkError(error);

    // // Make the application continue execution until stopping with Ctrl-C
    // // #TODO add ctrlC capture for proper streaming closing
    // while (true)
    // {
    //     Pa_Sleep(1);
    // }

    // // Close stream
    // error = Pa_StopStream(stream);
    // checkError(error);

    error = Pa_CloseStream(stream);
    checkError(error);

    error = Pa_Terminate();
    checkError(error);

    // // Releases memory for fourier transform
    // fftw_destroy_plan(spectrogramData->p);
    // fftw_free(spectrogramData->in);
    // fftw_free(spectrogramData->out);
    // fftw_free(spectrogramData);

    return EXIT_SUCCESS;
}
