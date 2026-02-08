// Copyright 2025 Giorgio Gamba

constexpr double SAMPLE_RATE{44100.0};
constexpr int FRAMES_PER_BUFFER{512};
constexpr int NUM_CHANNELS{2};

constexpr int DISPLAY_SIZE{25};

// Defines spectrogram's boundaries
constexpr int SPECTROGRAM_FREQ_START{20};
constexpr int SPECTROGRAM_FREQ_END{20 * 1000};