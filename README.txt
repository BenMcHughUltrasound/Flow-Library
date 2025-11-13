# Wave_Lib – MATLAB Waveform Processing Library

`Wave_Lib` is a MATLAB class that provides a comprehensive toolkit for **waveform analysis**, including FFT, CWT, signal processing, and flow-rate computations.

It integrates seamlessly with `WaveformObj`, `fftObj`, and `cwtObj` classes and supports both quick-look plotting and advanced batch processing.

---

## 📘 Overview

**Wave_Lib** allows you to:
- Load, save, and convert waveform data.
- Perform FFTs and inverse FFTs.
- Compute Continuous Wavelet Transforms (CWTs).
- Plot time-domain, frequency-domain, and CWT results.
- Apply windowing, smoothing, and interpolation.
- Calculate flow velocities and rates.
- Detect signal arrivals via thresholding.

---

## ⚙️ Requirements

Ensure the following are on your MATLAB path:
- `WaveformObj.m`
- `fftObj.m`
- `cwtObj.m`
- `STANDARDIZE_FIGURE.m`

Tested with MATLAB **R2021b+**.

---

## 🚀 Quick Start

```matlab
% Initialize
lib = Wave_Lib();

% Load waveform data from a CSV file
wave = Wave_Lib.LoadWaveforms("example_data.csv");

% Plot waveform
Wave_Lib.plt(wave(1), 1);

% Perform FFT and plot
fft_result = Wave_Lib.fftObj(wave(1), 1, 1, 2048);
Wave_Lib.pltFFT(fft_result, 2);

% Reconstruct waveform from FFT
wave_recon = Wave_Lib.IFFTObj(fft_result);
Wave_Lib.plt(wave_recon, 3);

% Compute CWT
cwt_result = Wave_Lib.cwtObj(wave(1), 10e3, 1e6);
Wave_Lib.pltcwt(cwt_result(1), 4);
