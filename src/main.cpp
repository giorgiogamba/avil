// Copyright 2025 Giorgio Gamba

#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "types.h"
#include "constants.h"

//static StreamCallbackData* spectrogramData;
static FileCallbackData* fileSpectrogramData;

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

void printFileFrequencyGraph(const float* in, unsigned long framesPerBuffer, void* userData)
{
    FileCallbackData* streamData = reinterpret_cast<FileCallbackData*>(userData);

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
    const sf_count_t numRead = sf_read_float(data->file, out, framesPerBuffer * data->info.channels);

    printVolumeGraph(out, framesPerBuffer);
    printFileFrequencyGraph(out, framesPerBuffer, data);
    fflush(stdout);

    const bool fileEnded = numRead < framesPerBuffer;
    return fileEnded ? paComplete : paContinue;
}

void listAvailableDevices()
{
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
}

void releaseResources(FileCallbackData* fileSpectrogramData)
{
    if (!fileSpectrogramData)
        return;

    fftw_destroy_plan(fileSpectrogramData->p);
    fftw_free(fileSpectrogramData->in);
    fftw_free(fileSpectrogramData->out);
    fftw_free(fileSpectrogramData);
}

// Applies Hann Window tecnique in order to avoid spectral leakage
void applyHannWindow(double* input, int size)
{
    for (int i = 0; i < size; ++i)
    {
        const double window = 0.5 * (1.0 - cos(2.0 * M_PI * i / (size - 1)));
        input[i] *= window;
    }
}

bool detectFrequencyCutoff(double* fftOutput, int fftSize, int sampleRate)
{
    // Detects if cuts over 16kHz happened. If nyquist is under then we cannot reconstruct the signal 
    constexpr double cutoffFreq = 16000.0;
    const double nyquist = sampleRate / 2.0;
    if (nyquist < cutoffFreq)
    {
        return false;
    }

    // Comparing over fftSize / 2 because using the FFT simmetry property
    const int halfSize = fftSize / 2;

    // Supposed the frequency spectrum is divided in bins, gets the bin were the cutoff frequency appears
    int cutoffBin = static_cast<int>((cutoffFreq * fftSize) / sampleRate);

    const auto getMagnitude = [&](const int bin)
    {
        // Works on the R2HC format, where imaginary part for N is at length-N
        if (bin == 0)
            return fftOutput[0] * fftOutput[0];

        if (bin == halfSize && fftSize % 2 == 0)
            return fftOutput[halfSize] * fftOutput[halfSize];

        const double real = fftOutput[bin];
        const double imag = fftOutput[fftSize - bin];
        return (real * real + imag * imag);
    };

    double highFreqEnergy = 0.0;
    for (int i = cutoffBin; i < halfSize; ++i)
    {
        highFreqEnergy += getMagnitude(i);
    }

    double midFreqEnergy = 0.0;
    int midStart = static_cast<int>((10000.0 * fftSize) / sampleRate);
    for (int i = midStart; i <= (cutoffBin - 1); ++i)
    {
        midFreqEnergy += getMagnitude(i);
    }

    // The frame is too silent
    constexpr double minEnergy = 1e-6;
    if (midFreqEnergy <= minEnergy)
        return false;

    constexpr double minRatio = 0.005;
    const double ratio = highFreqEnergy / midFreqEnergy;
    return (ratio < minRatio) ? true : false;
}

bool isProbablyMP3(SNDFILE* file, const SF_INFO& info, FileCallbackData* data)
{
    // Skips the first few seconds to avoid silence/fades
    constexpr int secondsToSkip = 5;
    sf_seek(file, info.samplerate * info.channels * secondsToSkip, SEEK_SET);

    int cutoffDetections = 0;
    constexpr const int framesToCheck = 100;
    int validFrames  = 0;

    for (int frame = 0; frame < framesToCheck; ++frame)
    {
        // Reads and puts into buffer
        std::vector<float> buffer(FRAMES_PER_BUFFER * info.channels);
        const sf_count_t numRead = sf_readf_float(file, buffer.data(), FRAMES_PER_BUFFER); // sf_read_float advances automatically
        if (numRead < FRAMES_PER_BUFFER)
            break;

        // Downmixes to mono and computes colume (RMS)
        double sumSq = 0.0;
        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            data->in[i] = 0.f;
            for (int c = 0; c < info.channels; ++c)
            {
                data->in[i] += buffer[i * info.channels + c];
            }

            data->in[i] /= info.channels;
            sumSq += data->in[i] * data->in[i];
        }

        // Avoids frequency analysis on bins that are too silent
        const double rms = sqrt(sumSq / FRAMES_PER_BUFFER);
        if (rms < 0.005)
            continue;

        validFrames ++;

        applyHannWindow(data->in, FRAMES_PER_BUFFER);

        // Executes Fourier transform
        fftw_execute(data->p);

        if (detectFrequencyCutoff(data->out, FRAMES_PER_BUFFER, info.samplerate))
        {
            cutoffDetections++;
        }
    }

    // resets file position
    sf_seek(file, 0, SEEK_SET);

    if (validFrames == 0)
        return false;

    const double ratio = (double) cutoffDetections / validFrames;
    std::cout << "#CUTOFFS/TOTAL = " << cutoffDetections << "/" << validFrames << "\n"; 
    return ratio > 0.7;
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
        std::cout << sf_strerror(data.file) << "\n";
        return EXIT_FAILURE;
    }

    // Every time a Pa operation is performed, its value is checked to spot problems
    PaError error;

    error = Pa_Initialize();
    checkError(error);

    fileSpectrogramData = (FileCallbackData*)malloc(sizeof(FileCallbackData));
	fileSpectrogramData->in = (double*)fftw_malloc(FRAMES_PER_BUFFER * sizeof(double));
	fileSpectrogramData->out = (double*)fftw_malloc(FRAMES_PER_BUFFER * sizeof(double));

	if (!fileSpectrogramData->in || !fileSpectrogramData->out)
    {
        std::cout << "Could not allocate spectrogram data\n";
        exit(EXIT_FAILURE);
    }

    constexpr double sampleRatio = FRAMES_PER_BUFFER / static_cast<double>(SAMPLE_RATE);
    fileSpectrogramData->startIndex = std::ceil(sampleRatio * SPECTROGRAM_FREQ_START);
    constexpr double DEF_SIZE{FRAMES_PER_BUFFER / 2.0};
    fileSpectrogramData->sprectrogramSize = std::min(std::ceil(sampleRatio * SPECTROGRAM_FREQ_END), DEF_SIZE) - fileSpectrogramData->startIndex;

    // Defines the Fourier transform. Data need to remain the same as long as this profile is chosen
    fileSpectrogramData->p = fftw_plan_r2r_1d(FRAMES_PER_BUFFER, fileSpectrogramData->in, fileSpectrogramData->out, FFTW_R2HC, FFTW_ESTIMATE);

    PaStream* stream;
    fileSpectrogramData->file = data.file;
    fileSpectrogramData->info = data.info;
    error = Pa_OpenDefaultStream(&stream, 0, data.info.channels, paFloat32, data.info.samplerate, FRAMES_PER_BUFFER, streamCallback, fileSpectrogramData);
    checkError(error);

    if (isProbablyMP3(fileSpectrogramData->file, fileSpectrogramData->info, fileSpectrogramData))
    {
        std::cout << "WARNING\n";
        Pa_CloseStream(stream);
        Pa_Terminate();
        sf_close(data.file);
        releaseResources(fileSpectrogramData);
        return EXIT_FAILURE;
    }

    std::cout << "Starting stream, stop it with Ctrl+C...\n";
    error = Pa_StartStream(stream);
    checkError(error);

    // Plays until the file continues
    while (Pa_IsStreamActive(stream))
    {
        Pa_Sleep(100);
    }

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

    error = Pa_StopStream(stream);
    checkError(error);

    error = Pa_CloseStream(stream);
    checkError(error);

    error = Pa_Terminate();
    checkError(error);

    releaseResources(fileSpectrogramData);

    return EXIT_SUCCESS;
}
