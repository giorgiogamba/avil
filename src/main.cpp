// Copyright 2025 Giorgio Gamba

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "types.h"
#include "constants.h"
#include "fft.hpp"

static StreamCallbackData* spectrogramData;
static FileCallbackData* fileSpectrogramData;

static std::vector<float> magnitudes(DISPLAY_SIZE, 0.f);
static std::mutex mtxMagnitudes;
static float magnitudeL = 0.f;
static float magnitudeR = 0.f;

/**
 * @brief Analyzes the passes error and eventually prints its content
 * @param error The error to be analyzeds
 */
void checkError(const PaError& error)
{
    if (error != paNoError)
    {
        std::cout << "[PortAudio Error] " << Pa_GetErrorText(error) << std::endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Prints a cmd-line volume graph representing the volume of the passed input
 * 
 * @param in Data to be displayed
 * @param framesPerBuffer Number of frames to be analyzed from in buffer
 */
void updateVolumeGraph(const float* in, unsigned long framesPerBuffer)
{
    // Parses over elements considering them in couples L, R
    // and gets the greates value of the buffer
    float volumeL = 0, volumeR = 0;
    for (unsigned long i = 0; i < framesPerBuffer * 2; i += 2)
    {
        volumeL = std::max(volumeL, std::abs(in[i]));
        volumeR = std::max(volumeR, std::abs(in[i+1]));
    }

    // Prints resulting volume bars
    std::unique_lock<std::mutex> lock(mtxMagnitudes);
    magnitudeL = volumeL;
    magnitudeR = volumeR;
}

/**
 * @brief Utility function to retrieve magnitude from the passed fft object at the passed bin
 * 
 * @param fftOutput object representing FFT signal
 * @param bin index to be analuzed
 * @param size Buffers length
 * @return The magnitude of fftOutput of size size at the passed bin
 */
double getMagnitude(double* fftOutput, const int bin, const int size)
{
    // Works on the R2HC format, where imaginary part for N is at length-N
    if (bin == 0)
        return fftOutput[0] * fftOutput[0];

    if (bin == size && size % 2 == 0)
        return fftOutput[bin] * fftOutput[bin];

    const double real = fftOutput[bin];
    const double imag = fftOutput[size - bin];
    return (real * real + imag * imag);
}

/**
 * @brief Prints a cmd-line frequency spectrum of the passed input
 * 
 * @param in input data
 * @param framesPerBuffer length on input data in
 * @param userData data to handle fft
 */
void updateFrequencyGraph(const float* in, unsigned long framesPerBuffer, void* userData)
{
    StreamCallbackData* streamData = reinterpret_cast<StreamCallbackData*>(userData);

    // Collects data for Fourier transform display
    // Audio channels samples are set in a matrix where row i contains i-th samples form all the audio channels
    // We own just 2 channels, so we get the even samples
    for (auto i{0}; i < framesPerBuffer; ++i)
        streamData->in[i] = in[i * NUM_CHANNELS];

    // Executes fourier transform on the stremed data
    fftw_execute(streamData->p);


    std::vector<float> rawMagnitudes(DISPLAY_SIZE);
    double maxRawMagnitude = -__DBL_MAX__;
    for (int i{0}; i < DISPLAY_SIZE; ++i)
    {
        const double step = i / static_cast<double>(DISPLAY_SIZE);
        const auto binIdx = static_cast<int>(streamData->startIndex + step * streamData->sprectrogramSize);
        const auto mag = getMagnitude(streamData->out, binIdx, FRAMES_PER_BUFFER);
        rawMagnitudes[i] = mag;
        maxRawMagnitude = std::max(maxRawMagnitude, mag);
    }

    /**
     * @brief How this value is composed
     * Since magnitudes are complex numbers, they are scaling wiht FRAMES_PER_BUFFER * FRAMES_PER_BUFFER
     * 1e-6 corresponds to an amplitude of 0.001, whihc is -60db, which is a standard noise flooe
     */
    constexpr double MIN_DISPLAY_ENERGY = 1e-6 * FRAMES_PER_BUFFER * FRAMES_PER_BUFFER;

    // Once all the magnitudes are collected, normalize them based on the max value
    std::vector<float> frequenciesMagnitudes(DISPLAY_SIZE);
    const double logMaxRawMagnitude = std::log1p(maxRawMagnitude);
    for (int i{0}; i < DISPLAY_SIZE; ++i)
    {
        // Avoids division with 0
        if (logMaxRawMagnitude <= 0.0 || rawMagnitudes[i] < MIN_DISPLAY_ENERGY)
        {
            frequenciesMagnitudes[i] = 0.f;
            continue;
        }

        const float normalizedMagnitude = static_cast<float>(std::log1p(rawMagnitudes[i]) / logMaxRawMagnitude);
        frequenciesMagnitudes[i] = std::clamp(normalizedMagnitude, 0.f, 1.f);
    }

    std::unique_lock<std::mutex> lock(mtxMagnitudes);
    magnitudes.swap(frequenciesMagnitudes);
}

/**
 * @brief Callback function invoked once that an audio capture of a buffer is completed.
 *        Used to execute data processing on the captured data before sending to the output buffer
 * 
 * @param inputBuffer audio data captured
 * @param outputBuffer audio data transmitted to the output
 * @param framesPerBuffer frames of input and output buffers
 * @param timeInfo 
 * @param statusFlags 
 * @param userData 
 * @return a code representing of the callback ended correctly
 */
int fileStreamCallback  ( const void* inputBuffer
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

    updateVolumeGraph(out, framesPerBuffer);
    updateFrequencyGraph(out, framesPerBuffer, data);

    const bool fileEnded = numRead < framesPerBuffer;
    return fileEnded ? paComplete : paContinue;
}

/**
 * @brief Callback function invoked once that a microphone audio capture of a buffer is completed.
 *        Used to execute data processing on the captured data before sending to the output buffer
 * 
 * @param inputBuffer audio data captured
 * @param outputBuffer audio data transmitted to the output
 * @param framesPerBuffer frames of input and output buffers
 * @param timeInfo 
 * @param statusFlags 
 * @param userData 
 * @return a code representing of the callback ended correctly
 */
int microphoneStreamCallback( const void* inputBuffer
                            , void* outputBuffer
                            , unsigned long framesPerBuffer
                            , const PaStreamCallbackTimeInfo* timeInfo
                            , PaStreamCallbackFlags statusFlags
                            , void* userData )
{
    if (!inputBuffer)
        return paContinue;

    const float* input = (float*) inputBuffer;
    float* output = (float*) outputBuffer;
    StreamCallbackData* data = (StreamCallbackData*)userData;

    std::vector<std::complex<float>> signal(framesPerBuffer);
    for (unsigned long i {0}; i < framesPerBuffer; ++i)
        signal[i] = std::complex<float>(input[i * NUM_CHANNELS], 0.f);

    FFT_AVX2::fft(signal.data(), framesPerBuffer);
    FFT_AVX2::suppressFeedback(signal.data(), framesPerBuffer, 0.3f);
    FFT_AVX2::ifft(signal.data(), framesPerBuffer);

    for (unsigned long i{0}; i < framesPerBuffer; ++i)
    {
        const float processedSample = signal[i].real();
        for (int ch{0}; ch < NUM_CHANNELS; ++ch)
            output[i * NUM_CHANNELS + ch] = processedSample * GAIN_REDUCTION_FACTOR;
    }

    updateVolumeGraph(input, framesPerBuffer);
    updateFrequencyGraph((float*)output, framesPerBuffer, data);

    return paContinue;
}

/**
 * @brief Displays available devices for the current device
 */
void listAvailableDevices()
{
    std::cout << "===========================" << std::endl;
    std::cout << "AVAILABLE DEVICES" << std::endl;
    const int numDevices = Pa_GetDeviceCount();
    std::cout << "Number of devices: " << numDevices << std::endl;

    std::cout << "---------" << std::endl;
    for (int i = 0; i < numDevices; ++i)
    {
        const PaDeviceInfo* deviceInfo = Pa_GetDeviceInfo(i);
        std::cout << "Name: " << deviceInfo->name << std::endl;
        std::cout << "INDEX: " << i << std::endl;
        std::cout << "Sample rate: " << deviceInfo->defaultSampleRate << std::endl;
        std::cout << "Low input latency: " << deviceInfo->defaultLowInputLatency << std::endl;
        std::cout << "---------" << std::endl;
    }

    std::cout << "===========================" << std::endl;
}

/**
 * @brief Frees resources for the passed spectrogram
 * 
 * @param fileSpectrogramData utility struct for file audio playback
 */
void releaseResources(FileCallbackData* fileSpectrogramData)
{
    if (!fileSpectrogramData)
        return;

    fftw_destroy_plan(fileSpectrogramData->p);
    fftw_free(fileSpectrogramData->in);
    fftw_free(fileSpectrogramData->out);
    fftw_free(fileSpectrogramData);
}

/**
 * @brief Frees resources for the passed spectrogram
 * 
 * @param fileSpectrogramData utility struct for microphone audio playback
 */
void releaseStreamResources(StreamCallbackData* streamSpectrogramData)
{
    if (!streamSpectrogramData)
        return;

    fftw_destroy_plan(streamSpectrogramData->p);
    fftw_free(streamSpectrogramData->in);
    fftw_free(streamSpectrogramData->out);
    free(streamSpectrogramData);
}

/**
 * @brief Applies Hann Window tecnique in order to avoid spectral leakage
 * 
 * @param input input data to be handled
 * @param size size of the input data
 */
void applyHannWindow(double* input, int size)
{
    for (int i = 0; i < size; ++i)
    {
        const double window = 0.5 * (1.0 - cos(2.0 * M_PI * i / (size - 1)));
        input[i] *= window;
    }
}

/**
 * @brief Checks if the currently playing file has been upscaled
 * 
 * @param fftOutput the fft to be analyzed
 * @param fftSize size of data to be analyzed
 * @param sampleRate rate of data sampling
 * @return true if the file has been upscaled
 * @return false if the file is in real format
 */
bool detectFrequencyCutoff(double* fftOutput, int fftSize, int sampleRate)
{
    // Detects if cuts over 16kHz happened.
    // If nyquist is under then we cannot reconstruct the signal 
    constexpr double cutoffFreq = 16000.0;
    const double nyquist = sampleRate / 2.0;
    if (nyquist < cutoffFreq)
        return false;

    // Comparing over fftSize / 2 because using the FFT simmetry property
    const int halfSize = fftSize / 2;

    // Supposed the frequency spectrum is divided in bins, 
    // gets the bin were the cutoff frequency appears
    const int cutoffBin = static_cast<int>((cutoffFreq * fftSize) / sampleRate);

    double highFreqEnergy = 0.0;
    for (int i = cutoffBin; i < halfSize; ++i)
        highFreqEnergy += getMagnitude(fftOutput, i, fftSize);

    double midFreqEnergy = 0.0;
    int midStart = static_cast<int>((10000.0 * fftSize) / sampleRate);
    for (int i = midStart; i <= (cutoffBin - 1); ++i)
        midFreqEnergy += getMagnitude(fftOutput, i, fftSize);

    // The frame is too silent
    if (midFreqEnergy <= MIN_ENERGY)
        return false;

    const double ratio = highFreqEnergy / midFreqEnergy;
    return (ratio < MIN_RATIO) ? true : false;
}

/**
 * @brief Checks if the current file has been upscaled
 * 
 * @param file file to check
 * @param info meta informations related to the file
 * @param data callback data used for analysis
 * @return true if the file could have been upsaled
 * @return false if the file could have not been upscaled
 */
bool isProbablyMP3(SNDFILE* file, const SF_INFO& info, FileCallbackData* data)
{
    // Skips the first few seconds to avoid silence/fades
    sf_seek(file, info.samplerate * info.channels * SECONDS_TO_SKIP, SEEK_SET);

    int cutoffDetections = 0;
    int validFrames  = 0;
    for (int frame = 0; frame < FRAMES_TO_CHECK; ++frame)
    {
        // Reads and puts into buffer. sf_read_float advances automatically
        std::vector<float> buffer(FRAMES_PER_BUFFER * info.channels);
        const sf_count_t numRead = sf_readf_float(file, buffer.data(), FRAMES_PER_BUFFER);
        if (numRead < FRAMES_PER_BUFFER)
            break;

        // Downmixes to mono and computes colume (RMS)
        double sumSq = 0.0;
        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            data->in[i] = 0.f;
            for (int c = 0; c < info.channels; ++c)
                data->in[i] += buffer[i * info.channels + c];

            data->in[i] /= info.channels;
            sumSq += data->in[i] * data->in[i];
        }

        // Avoids frequency analysis on bins that are too silent
        if (sqrt(sumSq / FRAMES_PER_BUFFER) < MIN_RATIO)
            continue;

        validFrames ++;

        applyHannWindow(data->in, FRAMES_PER_BUFFER);

        // Executes Fourier transform
        fftw_execute(data->p);

        if (detectFrequencyCutoff(data->out, FRAMES_PER_BUFFER, info.samplerate))
            cutoffDetections++;
    }

    // resets file position
    sf_seek(file, 0, SEEK_SET);

    if (validFrames == 0)
        return false;

    const double ratio = (double) cutoffDetections / validFrames;
    std::cout << "#CUTOFFS/TOTAL = " << cutoffDetections << "/" << validFrames << std::endl; 
    return ratio > 0.7;
}

int main(int argc, const char* argv[])
{
    // The application provides two different features:
    // Microphone audio capturing and upscaled audio files recognition
    // 1. If no arguments are provided, audio capturing is performed,
    // 2. If a audio file path is provided, audio recognition is performed.
    //      2.1. If the file is original, the file is reproduced
    // Both features, when reproducing something, show a cmdline spectrum analyzer and amplitude meter

    PaError error;
    error = Pa_Initialize();
    checkError(error);

    // Reads cli arguments
    if (argc == 1)
    {
        //listAvailableDevices();

        int deviceIdx = Pa_GetDefaultInputDevice();
        if (deviceIdx == paNoDevice)
        {
            std::cout << "[PortAudio Error] No default input device found" << std::endl;
            return EXIT_FAILURE;
        }

        const PaDeviceInfo* deviceInfo = Pa_GetDeviceInfo(deviceIdx);

        // Use minimum of requested channels and device's max channels
        const int inputChannels = std::min(NUM_CHANNELS, deviceInfo->maxInputChannels);
        if (inputChannels <= 0)
        {
            std::cout << "[PortAudio Error] Device has no input channels" << std::endl;
            return EXIT_FAILURE;
        }

        PaStreamParameters inStreamParams;
        memset(&inStreamParams, 0, sizeof(inStreamParams));
        inStreamParams.channelCount = inputChannels;
        inStreamParams.device = deviceIdx;
        inStreamParams.hostApiSpecificStreamInfo = nullptr;
        inStreamParams.sampleFormat = paFloat32;
        inStreamParams.suggestedLatency = deviceInfo->defaultLowInputLatency;

        PaStreamParameters outParams;
        const int outputDevice = Pa_GetDefaultOutputDevice();
        outParams.device = outputDevice;
        const PaDeviceInfo* outDeviceInfo = Pa_GetDeviceInfo(outputDevice);
        const int outputChannels = std::min(NUM_CHANNELS, outDeviceInfo->maxOutputChannels);
        outParams.channelCount = outputChannels;
        outParams.sampleFormat = paFloat32;
        outParams.suggestedLatency = outDeviceInfo->defaultLowOutputLatency;
        outParams.hostApiSpecificStreamInfo = nullptr;

        spectrogramData = (StreamCallbackData*)malloc(sizeof(StreamCallbackData));
        spectrogramData->in = (double*)fftw_malloc(FRAMES_PER_BUFFER * sizeof(double));
        spectrogramData->out = (double*)fftw_malloc(FRAMES_PER_BUFFER * sizeof(double));
        if (!spectrogramData->in || !spectrogramData->out)
        {
            std::cout << "[AVIL ERROR] Could not allocate spectrogram data" << std::endl;
            exit(EXIT_FAILURE);
        }

        constexpr double sampleRatio = FRAMES_PER_BUFFER / static_cast<double>(SAMPLE_RATE);
        spectrogramData->startIndex = std::ceil(sampleRatio * SPECTROGRAM_FREQ_START);
        constexpr double DEF_SIZE{FRAMES_PER_BUFFER / 2.0};
        spectrogramData->sprectrogramSize = std::min(std::ceil(sampleRatio * SPECTROGRAM_FREQ_END), DEF_SIZE) - spectrogramData->startIndex;

        // Defines the Fourier transform. Data need to remain the same as long as this profile is chosen
        spectrogramData->p = fftw_plan_r2r_1d(FRAMES_PER_BUFFER, spectrogramData->in, spectrogramData->out, FFTW_R2HC, FFTW_ESTIMATE);

        PaStream* stream;
        error = Pa_OpenStream(&stream, &inStreamParams, &outParams, 44100.0, 512, paNoFlag, microphoneStreamCallback, spectrogramData);
        checkError(error);
        error = Pa_StartStream(stream);
        checkError(error);

        while (true)
            Pa_Sleep(1);
    }
    else if (argc == 2)
    {
        const char* filePath = argv[1];
        FileCallbackData data{};
        data.file = sf_open(filePath, SFM_READ, &data.info);
        if (sf_error(data.file) != SF_ERR_NO_ERROR)
        {
            std::cout << "[AVIL ERROR] An error occured while opening file " << filePath << ". ("<< sf_strerror(data.file) <<")" << std::endl;
            return EXIT_FAILURE;
        }

        // Validate channel count
        if (data.info.channels <= 0 || data.info.channels > 32)
        {
            std::cout << "[AVIL ERROR] Invalid channel count in file: " << data.info.channels << std::endl;
            sf_close(data.file);
            return EXIT_FAILURE;
        }

        fileSpectrogramData = (FileCallbackData*)malloc(sizeof(FileCallbackData));
        fileSpectrogramData->in = (double*)fftw_malloc(FRAMES_PER_BUFFER * sizeof(double));
        fileSpectrogramData->out = (double*)fftw_malloc(FRAMES_PER_BUFFER * sizeof(double));
        if (!fileSpectrogramData->in || !fileSpectrogramData->out)
        {
            std::cout << "[AVIL ERROR] Could not allocate spectrogram data" << std::endl;
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
        error = Pa_OpenDefaultStream(&stream, 0, data.info.channels, paFloat32, data.info.samplerate, FRAMES_PER_BUFFER, fileStreamCallback, fileSpectrogramData);
        checkError(error);

        if (isProbablyMP3(fileSpectrogramData->file, fileSpectrogramData->info, fileSpectrogramData))
        {
            std::cout << "WARNING. The file is probably a mp3 file converted to wav file" << std::endl;
            Pa_CloseStream(stream);
            Pa_Terminate();
            sf_close(data.file);
            releaseResources(fileSpectrogramData);
            return EXIT_FAILURE;
        }

        error = Pa_StartStream(stream);
        checkError(error);

        // Plays until the file continues
        while (Pa_IsStreamActive(stream))
            Pa_Sleep(100);

        error = Pa_StopStream(stream);
        checkError(error);
        error = Pa_CloseStream(stream);
        checkError(error);
    }
    else
    {
        std::cout << "[AVIL ERROR] Invalid arguments provided. Exiting..." << std::endl;
        return EXIT_FAILURE;
    }

    error = Pa_Terminate();
    checkError(error);

    releaseResources(fileSpectrogramData);
    releaseStreamResources(spectrogramData);

    return EXIT_SUCCESS;
}
