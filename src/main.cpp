// Copyright 2025 Giorgio Gamba

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ftxui/dom/elements.hpp>
#include <ftxui/screen/screen.hpp>
#include <ftxui/component/loop.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/component/component.hpp>

#include "types.h"
#include "constants.h"
#include "fft.hpp"

static StreamCallbackData* spectrogramData;
static FileCallbackData* fileSpectrogramData;

static std::vector<float> magnitudes(DISPLAY_SIZE, 0.f);
static std::mutex mtxMagnitudes;
static float magnitudeL = 0.f;
static float magnitudeR = 0.f;
static std::string upscalingDetectionResultString = "";

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
 * @brief Returns the magnitude of the provided complex number
 * @param c the complex number
 */
inline float getComplexMagnitude(const std::complex<float>& c) noexcept
{
    return c.real() * c.real() + c.imag() * c.imag();
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
 * @brief Prints a cmd-line frequency spectrum of the passed input
 * 
 * @param in input data
 * @param framesPerBuffer length on input data in
 * @param userData data to handle fft
 */
void updateFrequencyGraph(const float* in, unsigned long framesPerBuffer, const double startIndex, const double spectrogramSize)
{
    // Converts data in complex format
    std::vector<std::complex<float>> signal(framesPerBuffer);
    for (unsigned long i {0}; i < framesPerBuffer; ++i)
        signal[i] = std::complex<float>(in[i * NUM_CHANNELS], 0.f);

    FFT_AVX2::fft(signal.data(), framesPerBuffer);

    std::vector<float> rawMagnitudes(DISPLAY_SIZE);
    double maxRawMagnitude = -__DBL_MAX__;
    for (int i{0}; i < DISPLAY_SIZE; ++i)
    {
        const double step = i / static_cast<double>(DISPLAY_SIZE);
        const auto binIdx = static_cast<int>(startIndex + step * spectrogramSize);
        rawMagnitudes[i] = getComplexMagnitude(signal[binIdx]);
        maxRawMagnitude = std::max(maxRawMagnitude, rawMagnitudes[i]);
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
    updateFrequencyGraph(out, framesPerBuffer, data->startIndex, data->sprectrogramSize);

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
    updateFrequencyGraph((float*)output, framesPerBuffer, data->startIndex, data->sprectrogramSize);

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
 * @brief Checks if the currently playing file has been upscaled
 * 
 * @param fftOutput the fft to be analyzed
 * @param fftSize size of data to be analyzed
 * @param sampleRate rate of data sampling
 * @return true if the file has been upscaled
 * @return false if the file is in real format
 */
bool detectFrequencyCutoff(const std::vector<std::complex<float>>& fftOutput, int fftSize, int sampleRate)
{
    // Detects if cuts over 16kHz happened.
    // If nyquist is under then we cannot reconstruct the signal 
    constexpr double cutoffFreq = 16000.0;
    const double nyquist = sampleRate / 2.0;
    if (nyquist < cutoffFreq)
        return false;

    // Never read beyond half the buffer (FFT symmetry) or the buffer size itself
    const int halfSize = std::min(fftSize / 2, static_cast<int>(fftOutput.size()));
    const int cutoffBin = static_cast<int>((cutoffFreq * fftSize) / sampleRate);
    const int midStart  = static_cast<int>((10000.0 * fftSize) / sampleRate);

    // Guard against bad bin calculations exceeding buffer
    if (cutoffBin >= halfSize || midStart >= cutoffBin)
        return false;

    double highFreqEnergy = 0.0;
    for (int i = cutoffBin; i < halfSize; ++i)
        highFreqEnergy += getComplexMagnitude(fftOutput[i]);

    double midFreqEnergy = 0.0;
    const int midStart = static_cast<int>((10000.0 * fftSize) / sampleRate);
    for (int i = midStart; i <= (cutoffBin - 1); ++i)
        midFreqEnergy += getComplexMagnitude(fftOutput[i]);

    // The frame is too silent
    if (midFreqEnergy <= MIN_ENERGY)
        return false;

    const double ratio = highFreqEnergy / midFreqEnergy;
    const bool isOverThreshold = ratio < MIN_RATIO;
    return isOverThreshold ? true : false;
}

void stopStream(PaStream* stream)
{
    PaError error{};
    error = Pa_StopStream(stream);
    checkError(error);
    error = Pa_CloseStream(stream);
    checkError(error);
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
bool isProbablyMP3(SNDFILE* file, const SF_INFO& info)
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

        std::vector<std::complex<float>> signal;
        signal.resize(FRAMES_PER_BUFFER);

        // Downmixes to mono and computes volume (RMS)
        double sumSq = 0.0;
        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            float mono = 0.f;
            for (int c = 0; c < info.channels; ++c)
                mono += buffer[i * info.channels + c];

            mono /= info.channels;
            signal[i] = std::complex<float>(mono, 0.f);
            sumSq += getComplexMagnitude(signal[i]);
        }

        // Avoids frequency analysis on bins that are too silent
        if (sqrt(sumSq / FRAMES_PER_BUFFER) < MIN_RATIO)
            continue;

        validFrames++;

        // Applies Hann Window
        for (int i = 0; i < FRAMES_PER_BUFFER; ++i)
        {
            const double window = 0.5 * (1.0 - cos(2.0 * M_PI * i / (FRAMES_PER_BUFFER - 1)));
            signal[i] = std::complex<float>(signal[i].real() * static_cast<float>(window));
        }

        FFT_AVX2::fft(signal.data(), FRAMES_PER_BUFFER);

        if (detectFrequencyCutoff(signal, FRAMES_PER_BUFFER, info.samplerate))
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

void runTUI(ftxui::ScreenInteractive& screen, std::atomic<bool>* running)
{
    using namespace ftxui;

    // Refreshes TUI at 60 fps
    std::thread refreshThread([&screen, running]()
    {
        while (*running)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(16));
            screen.PostEvent(Event::Custom);
        }
    });

    // Creates the TUI
    auto renderer = Renderer([&] {
        std::lock_guard<std::mutex> lock(mtxMagnitudes);

        // Spectrum bars
        Elements bars;
        for (const float mag : magnitudes)
            // Need to take 1.f - value becuase ftxui doesn't have a
            // bar representation going from bottom to top
            bars.push_back(gaugeDown(1.f - mag) | color(Color::Green) | flex);

        // Y axis labels, evenly spaced top to bottom (0dB at top, -60dB at bottom)
        const auto yAxis = vbox({
            text("0dB")  | color(Color::GrayLight),
            filler(),
            text("-20")  | color(Color::GrayLight),
            filler(),
            text("-40")  | color(Color::GrayLight),
            filler(),
            text("-60")  | color(Color::GrayLight),
        }) | size(WIDTH, EQUAL, 4);

        // X axis frequency labels evenly spaced
        const auto xAxis = hbox({
            text(std::to_string((int)SPECTROGRAM_FREQ_START) + "Hz") | color(Color::GrayLight),
            filler(),
            text(std::to_string((int)(SPECTROGRAM_FREQ_START + (SPECTROGRAM_FREQ_END - SPECTROGRAM_FREQ_START) * 0.25f) / 1000) + "kHz") | color(Color::GrayLight),
            filler(),
            text(std::to_string((int)(SPECTROGRAM_FREQ_START + (SPECTROGRAM_FREQ_END - SPECTROGRAM_FREQ_START) * 0.5f) / 1000) + "kHz") | color(Color::GrayLight),
            filler(),
            text(std::to_string((int)(SPECTROGRAM_FREQ_START + (SPECTROGRAM_FREQ_END - SPECTROGRAM_FREQ_START) * 0.75f) / 1000) + "kHz") | color(Color::GrayLight),
            filler(),
            text(std::to_string((int)SPECTROGRAM_FREQ_END / 1000) + "kHz") | color(Color::GrayLight),
        });

        const auto volumeMeters = hbox({
            text("L ")          | color(Color::Yellow),
            gauge(magnitudeL)   | color(Color::Yellow) | flex,
            text(" R ")         | color(Color::Yellow),
            gauge(magnitudeR)   | color(Color::Yellow) | flex,
        });

        const auto upscalingDetectionResult = upscalingDetectionResultString.empty()
            ? text("")
            : text(upscalingDetectionResultString) | color(upscalingDetectionResultString == "WARNING" ? Color::Red : Color::Green);

        return vbox({
            text("  avil — spectrum visualizer  ") | bold | center | color(Color::Cyan),
            separator(),
            hbox({
                yAxis,
                separator(),
                hbox(bars) | flex,
            }) | flex,
            separator(),
            hbox({text("     "), xAxis | flex}),
            separator(),
            volumeMeters,
            separator(),
            upscalingDetectionResult | center,
            separator(),
            text("Ctrl+C to quit") | dim | center,
        }) | border;
    });

    auto ioComponent = CatchEvent(renderer, [&](Event event)
    {
        if (event == Event::Character('q'))
        {
            *running = false;
            screen.ExitLoopClosure()();
            return true;
        }

        return false;
    });

    screen.Loop(ioComponent);

    *running = false;
    refreshThread.join();

    constexpr const char* CLEAR_TERMINAL_CODE{"\033[2J\033[H"};
    std::cout << CLEAR_TERMINAL_CODE << std::flush;
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

    constexpr double DEF_SIZE{FRAMES_PER_BUFFER / 2.0};

    // Reads cli arguments
    if (argc == 1)
    {
        const int deviceIdx = Pa_GetDefaultInputDevice();
        if (deviceIdx == paNoDevice)
        {
            std::cout << "[PortAudio Error] No default input device found" << std::endl;
            return EXIT_FAILURE;
        }

        const PaDeviceInfo* deviceInfo = Pa_GetDeviceInfo(deviceIdx);
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
        outParams.device = Pa_GetDefaultOutputDevice();
        const PaDeviceInfo* outDeviceInfo = Pa_GetDeviceInfo(outParams.device);
        outParams.channelCount = std::min(NUM_CHANNELS, outDeviceInfo->maxOutputChannels);
        outParams.sampleFormat = paFloat32;
        outParams.suggestedLatency = outDeviceInfo->defaultLowOutputLatency;
        outParams.hostApiSpecificStreamInfo = nullptr;

        spectrogramData = (StreamCallbackData*)malloc(sizeof(StreamCallbackData));

        constexpr double sampleRatio = FRAMES_PER_BUFFER / static_cast<double>(SAMPLE_RATE);
        spectrogramData->startIndex = std::ceil(sampleRatio * SPECTROGRAM_FREQ_START);
        spectrogramData->sprectrogramSize = std::min(std::ceil(sampleRatio * SPECTROGRAM_FREQ_END), DEF_SIZE) - spectrogramData->startIndex;

        PaStream* stream;
        error = Pa_OpenStream(&stream, &inStreamParams, &outParams, 44100.0, 512, paNoFlag, microphoneStreamCallback, spectrogramData);
        checkError(error);
        error = Pa_StartStream(stream);
        checkError(error);

        // Activates microphone player
        std::atomic<bool> running = true;
        auto screen = ftxui::ScreenInteractive::Fullscreen();
        runTUI(screen, &running);
        stopStream(stream);
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
        constexpr double sampleRatio = FRAMES_PER_BUFFER / static_cast<double>(SAMPLE_RATE);
        fileSpectrogramData->startIndex = std::ceil(sampleRatio * SPECTROGRAM_FREQ_START);
        fileSpectrogramData->sprectrogramSize = std::min(std::ceil(sampleRatio * SPECTROGRAM_FREQ_END), DEF_SIZE) - fileSpectrogramData->startIndex;
        fileSpectrogramData->file = data.file;
        fileSpectrogramData->info = data.info;

        if (isProbablyMP3(fileSpectrogramData->file, fileSpectrogramData->info))
        {
            upscalingDetectionResultString += "Probably the file has been upscaled";
            sf_close(data.file);
            return EXIT_FAILURE;
        }
        else
        {
            upscalingDetectionResultString = "File doesn't look upscaled";
        }

        PaStream* stream;
        error = Pa_OpenDefaultStream(&stream, 0, data.info.channels, paFloat32, data.info.samplerate, FRAMES_PER_BUFFER, fileStreamCallback, fileSpectrogramData);
        checkError(error);
        error = Pa_StartStream(stream);
        checkError(error);

        std::atomic<bool> running{true};
        auto screen = ftxui::ScreenInteractive::Fullscreen();

        std::thread filePlayer([&]()
        {
            while (Pa_IsStreamActive(stream) && running)
                Pa_Sleep(100);

            // File ended, closes TUI
            running = false;
            screen.ExitLoopClosure()();
        });

        runTUI(screen, &running);
        filePlayer.join();
        stopStream(stream);
    }
    else
    {
        std::cout << "[AVIL ERROR] Invalid arguments provided. Exiting..." << std::endl;
        return EXIT_FAILURE;
    }

    error = Pa_Terminate();
    checkError(error);

    free(fileSpectrogramData);
    free(spectrogramData);

    return EXIT_SUCCESS;
}
