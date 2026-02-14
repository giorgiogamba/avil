// Copyright 2025 Giorgio Gamba

constexpr double SAMPLE_RATE{44100.0};
constexpr int FRAMES_PER_BUFFER{512};
constexpr int NUM_CHANNELS{2};

constexpr int DISPLAY_SIZE{25};

// Defines spectrogram's boundaries
constexpr int SPECTROGRAM_FREQ_START{20};
constexpr int SPECTROGRAM_FREQ_END{20 * 1000};

constexpr float GAIN_REDUCTION_FACTOR{0.1f};
constexpr double MIN_ENERGY{1e-6};
constexpr double MIN_RATIO{0.005};
constexpr int SECONDS_TO_SKIP{5};
constexpr const int FRAMES_TO_CHECK{100};