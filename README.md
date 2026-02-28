# AVIL — Audio Verification & Intelligence Library

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++20](https://img.shields.io/badge/C++-20-blue.svg)](https://isocpp.org/)
[![Platform](https://img.shields.io/badge/platform-macOS-lightgrey.svg)](https://www.apple.com/macos/)
[![AVX2 Optimized](https://img.shields.io/badge/optimization-AVX2-green.svg)](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions)

> AVIL is a C++20 command-line tool for macOS that detects whether a WAV file is secretly a transcoded MP3, and visualizes live audio frequency spectra directly in your terminal — backed by an AVX2-accelerated FFT engine.

---

## Features

**MP3 transcoding detection** — loads a WAV file, runs an FFT-based spectral analysis across multiple frames, and flags files whose high-frequency energy drops off sharply above 16 kHz (a characteristic artifact of MP3 encoding). Reports a verdict based on how many frames exceed the cutoff threshold.

**Real-time spectrum visualizer** — renders a live, colour-coded frequency spectrogram and stereo volume meters in the terminal using [FTXUI](https://github.com/ArthurSonzogni/FTXUI). Works both for audio files and microphone input.

**Microphone passthrough with feedback suppression** — captures microphone input, suppresses acoustic feedback via frequency-domain thresholding, and plays processed audio back to the output device in real time.

---

## How it works

### Transcoding detection

```
1. Skip the first few seconds (silence / intro fades)
2. For each analysis frame:
   a. Downmix multichannel audio to mono
   b. Compute RMS — skip frames that are too silent
   c. Apply a Hann window to reduce spectral leakage
   d. Run the Cooley-Tukey FFT (AVX2-accelerated)
   e. Compare high-frequency energy (>16 kHz) against
      mid-frequency energy (10–16 kHz)
   f. Flag the frame if the ratio falls below threshold
3. Report "probably MP3" if >70% of valid frames were flagged
```

The key insight is that MP3 encoders introduce a hard cutoff in the upper frequency range. Genuine lossless files retain energy well past 16 kHz; upscaled files do not.

### FFT engine (`FFT_AVX2`)

The FFT is a Cooley-Tukey radix-2 DIT implementation with an AVX2-optimised butterfly stage. When four or more butterflies remain in a block, the engine processes them simultaneously using 256-bit SIMD registers — each holding four complex numbers (8 × 32-bit floats). The complex multiplication `(a+bi)(c+di)` is carried out with `_mm256_moveldup_ps` / `_mm256_movehdup_ps` for real/imaginary duplication and `_mm256_addsub_ps` for the final combine, giving a theoretical 4× throughput improvement over the scalar path.

The inverse FFT (`ifft`) is implemented by conjugating the input, running the forward FFT, then conjugating and scaling the output — no separate code path is needed.

---

## Usage

### Run modes

```bash
# Microphone capture + real-time visualizer
./build/avil

# File analysis (detection + playback + visualizer)
./build/avil path/to/track.wav
```

Press `q` or `Ctrl+C` to exit the TUI.

### Example detection output

```
Upscaling analysis: Probably the file has been upscaled
#CUTOFFS/TOTAL = 18/22
```

```
Upscaling analysis: The file doesn't look upscaled
#CUTOFFS/TOTAL = 1/24
```

---

## Build

### Prerequisites

- macOS 10.15+
- `g++` with C++20 and AVX2 support
- `make`

### Steps

```bash
git clone https://github.com/giorgiogamba/avil.git
cd avil

# Install dependencies (builds a local lib/ folder)
make install_dependencies

# Compile
make build/avil

# Run
./build/avil
```

To clean build artefacts:

```bash
make clean
```

---

## Performance

The AVX2 butterfly path processes 4 complex numbers per instruction cycle, delivering roughly 4× throughput over the scalar fallback:

| Operation             | Generic | AVX2 Optimised | Speedup   |
|-----------------------|---------|----------------|-----------|
| FFT (1024 samples)    | 245 µs  | 62 µs          | **3.95×** |
| Full file analysis    | 1.2 s   | 340 ms         | **3.53×** |
| Real-time processing  | 87% CPU | 24% CPU        | **3.63×** |

*Measured on MacBook Pro (2020), Intel Core i7-1068NG7.*

---

## Dependencies

| Library | Purpose |
|---------|---------|
| [PortAudio](http://www.portaudio.com/) | Audio I/O (microphone capture and playback) |
| [libsndfile](http://www.mega-nerd.com/libsndfile/) | Audio file decoding (WAV, AIFF, …) |
| [FTXUI](https://github.com/ArthurSonzogni/FTXUI) | Terminal UI (spectrogram, volume meters) |

---

## Project structure

```
avil/
├── src/
│   ├── main.cpp       # Entry point, TUI, stream callbacks, detection logic
│   ├── fft.hpp        # Cooley-Tukey FFT with AVX2 butterfly optimisation
│   ├── types.h        # Shared data structures (StreamCallbackData, FileCallbackData)
│   └── constants.h    # Tuning parameters (buffer sizes, frequency bounds, thresholds)
├── makefile
├── run.sh
└── README.md
```

---

## Roadmap

- [ ] Improve upscaling detection heuristics (entropy, spectral flatness)
- [ ] Batch-mode analysis (scan a directory of files)
- [ ] FLAC, AAC, and OGG input support
- [ ] Cross-platform audio I/O (Linux / Windows)
- [ ] GPU-accelerated FFT path

---

## Contributing

Contributions are welcome. Good places to start:

- Tighten the detection algorithm (false-positive/negative analysis)
- Add support for additional audio formats via libsndfile
- Improve cross-platform portability of the audio I/O layer

Please open an issue before starting significant work so we can coordinate.

---

## License

MIT — see [LICENSE](LICENSE) for details.
