# FFT on Teensy 4.0 — ARM Cortex-M7 DSP Hardware Acceleration & CMSIS

> **MEng Final Year Project** — Dublin City University, School of Electronic Engineering  
> **Author:** Hamdan Bin Abdul Rahman | **Supervisor:** Conor McArdle | **August 2025**

---

## 📌 Project Overview

This project implements a **hybrid FFT framework** on the **Teensy 4.0 (ARM Cortex-M7 @ 600 MHz)** that dynamically selects the most efficient FFT algorithm based on transform size. It combines ARM's CMSIS-DSP optimized routines for small sizes with custom Radix-2 and Radix-4 implementations for large transforms up to **16,384 points**, enabling real-time high-resolution spectral analysis on a resource-constrained embedded platform.

---

## 🧠 Abstract

High-resolution spectral analysis is essential for audio, communications, and monitoring applications. Standard ARM CMSIS-DSP FFT routines are limited to ~4096 points, restricting frequency resolution. This project extends that capability with a hybrid framework that dynamically selects between CMSIS-DSP, custom Radix-2, and custom Radix-4 algorithms — achieving real-time performance across all tested transform sizes, verified with both synthesised tones and live audio input.

---

## ⚙️ Algorithm Selection Logic

| FFT Size | Algorithm Used |
|----------|---------------|
| ≤ 4096   | ARM CMSIS-DSP `arm_rfft_fast_f32` |
| 8192     | Custom Radix-2 FFT |
| 16384    | Custom Radix-4 FFT |

The selection is made dynamically at runtime based on input size `N`.

---

## 🏗️ Hardware

| Component | Details |
|-----------|---------|
| Microcontroller | Teensy 4.0 |
| Processor | ARM Cortex-M7 @ 600 MHz |
| FPU | Single-precision hardware FPU |
| ADC | 12-bit, A0 pin |
| Audio Input | External microphone sensor module |
| Communication | USB Serial (115200 baud) |

---

## 📊 Performance Results

| FFT Size | Algorithm | Execution Time |
|----------|-----------|---------------|
| 4096 | CMSIS-DSP RFFT | 0.33 ms |
| 8192 | Custom Radix-4 FFT | 2.68 ms |
| 16384 | Custom Radix-4 FFT | 4.14 ms |

All transforms completed well under the **10 ms real-time threshold**.

### Precomputed vs On-the-Fly Twiddle Factors

| FFT Size | On-the-Fly `arm_sin/cos` (ms) | Precomputed (ms) | Speedup |
|----------|-------------------------------|------------------|---------|
| 8192 | 16.95 | 2.20 | ~670% faster |
| 16384 | 26.34 | 4.44 | ~493% faster |


## 🐍 Python Visualization

FFT results are streamed from the Teensy 4.0 to a host PC over USB Serial and plotted using **Matplotlib**.

### Static Sine Wave Test
```bash
python sinewavetest.py
```
Plots a single captured frame of the time-domain signal and FFT magnitude.

### Real-Time Audio Plot
```bash
python realtimeaudioplot.py
```
Continuously updates both the time-domain waveform and FFT magnitude spectrum in real time.

### Requirements
```bash
pip install pyserial numpy matplotlib
```

> ⚠️ Update the serial port in the Python scripts to match your system:
> - **Mac:** `/dev/tty.usbmodem######`
> - **Windows:** `COM3` (or whichever port Teensy appears on)

---

## 🔑 Key Features

- **Hybrid algorithm selection** — automatically picks the best FFT routine for the given transform size
- **Precomputed twiddle factors** — eliminates runtime `sin()`/`cos()` calls for up to 670% speedup
- **Hann windowing** — reduces spectral leakage for real-world signals
- **Parabolic interpolation** — sub-bin frequency estimation for improved peak accuracy
- **Auto-gain control** — handles quiet or low-amplitude input signals
- **Real-time serial streaming** — time and frequency domain data sent to PC for visualization
- **RAM-aware design** — memory usage profiled across all FFT sizes on Teensy 4.0

---

## 🧪 Test Results — 250 Hz Sine Wave

**16384-point Custom Radix-4 FFT (Serial Output):**

Used Custom Radix-4 FFT
Custom FFT - Bin: 410
Re: -4714.67, Im: -1533.45, Mag: 4957.77
Normalized Magnitude: 0.61
Generation Time (ms): 5.63
FFT Time (ms): 4.14
PEAK_FREQ (Hz): 250.16, PEAK_MAG: 0.6052

**Verification:**
- Expected bin: `(250 × 16384) / 10000 = 410` ✅
- Expected magnitude: `√(4714.67² + 1533.45²) = 4957.77` ✅
- Peak frequency: `410 × 10000 / 16384 = 250.16 Hz` ✅

---

## 💾 RAM Usage

| FFT Size | Algorithm | RAM Used | RAM Free |
|----------|-----------|----------|----------|
| 4096 | CMSIS-DSP | ~140 KB | ~351 KB |
| 8192 | Custom Radix-4 | ~265 KB | ~227 KB |
| 16384 | Custom Radix-4 | ~445 KB | ~46 KB |

---

## 📚 References

- Cooley & Tukey, "An algorithm for the machine calculation of complex Fourier series," *Math. Comput.*, 1965.
- ARM CMSIS-DSP Library — [github.com/ARM-software/CMSIS-DSP](https://github.com/ARM-software/CMSIS-DSP)
- Frigo & Johnson, "The design and implementation of FFTW3," *Proc. IEEE*, 2005.

---

## 📄 License

This project was developed as part of an MEng dissertation at Dublin City University. The custom FFT implementations are available for academic and educational use.

---

*Dublin City University — School of Electronic Engineering — MEng Electronic and Computer Engineering*


