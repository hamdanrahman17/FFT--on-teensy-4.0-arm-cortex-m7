/*
  Teensy Hybrid FFT (No Window, No Precomputed Twiddles)
  - CMSIS RFFT for sizes <= 4096
  - Custom radix-2 / radix-4 for larger sizes
  - Inline twiddles computed using arm_cos_f32()/arm_sin_f32()
  - Test-sine injection mode included for verification
 */
#include <Arduino.h>
#include <arm_math.h>  // CMSIS-DSP

// ----------------------- Config -----------------------
#define FFT_SIZE       16384  // choose 4096 or 8192 or 16384 (power of two)
#define CMSIS_MAX      4096

#define SAMPLING_FREQ  10000.0f // Hz
#define TEST_MODE      1        // 0 = real ADC sampling, 1 = inject test sine
#define TEST_FREQ_HZ   250.0f  // test sine freq used when TEST_MODE == 1
#define TEST_AMP       0.8f

// ----------------------- Buffers -----------------------
float32_t inputBuffer[FFT_SIZE];
float32_t magnitude[FFT_SIZE / 2];
float32_t fftBuffer[2 * FFT_SIZE]; // complex buffer [Re0,Im0,Re1,Im1,...]

arm_rfft_fast_instance_f32 cmsisFFT;

// ----------------------- Utilities -----------------------
bool isPowerOfTwo(int n) {
  return (n > 0) && ((n & (n - 1)) == 0);
}
bool isPowerOfFour(int n) {
  return isPowerOfTwo(n) && ((__builtin_ctz(n) % 2) == 0);
}

// ----------------------- Sampling -----------------------
void injectTestSine(float32_t* buffer, int size, float freq, float amp) {
  for (int i = 0; i < size; i++) {
    buffer[i] = amp * sinf(2.0f * PI * freq * ((float)i) / SAMPLING_FREQ);
  }
}





// ----------------------- Bit reversal helpers -----------------------
void bitReverse(float32_t* data, int n) {
  int j = 0;
  for (int i = 0; i < n; i++) {
    if (i < j) {
      float32_t tRe = data[2*i];
      float32_t tIm = data[2*i + 1];
      data[2*i]     = data[2*j];
      data[2*i + 1] = data[2*j + 1];
      data[2*j]     = tRe;
      data[2*j + 1] = tIm;
    }
    int m = n >> 1;
    while (j >= m && m >= 2) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
}

int reverseBase4(int x, int log4n) {
  int result = 0;
  for (int i = 0; i < log4n; i++) {
    result <<= 2;
    result |= (x & 0x3);
    x >>= 2;
  }
  return result;
}

void bitReverseRadix4(float32_t* data, int n) {
  int log4n = 0;
  int temp = n;
  while (temp > 1) {
    temp >>= 2;
    log4n++;
  }
  for (int i = 0; i < n; i++) {
    int j = reverseBase4(i, log4n);
    if (i < j) {
      float32_t tRe = data[2*i];
      float32_t tIm = data[2*i + 1];
      data[2*i]     = data[2*j];
      data[2*i + 1] = data[2*j + 1];
      data[2*j]     = tRe;
      data[2*j + 1] = tIm;
    }
  }
}

// ----------------------- Radix-4 FFT (inline twiddles via CMSIS) -----------------------
void fftRadix4(float32_t* data, int n) {
  bitReverseRadix4(data, n);

  for (int m = 4; m <= n; m <<= 2) {
    int groupSize = m;                // current DFT size for this stage
    int butterflies = m / 4;          // number of butterflies per group
    int step = n / m;                 // twiddle step in terms of bins

    for (int k = 0; k < n; k += groupSize) {
      for (int j = 0; j < butterflies; j++) {
        // compute twiddle angles
        float32_t angle1 = -2.0f * PI * ((float)j * (float)step) / (float)n;
        float32_t angle2 = -2.0f * PI * ((float)(2*j) * (float)step) / (float)n;
        float32_t angle3 = -2.0f * PI * ((float)(3*j) * (float)step) / (float)n;

        // compute twiddles using CMSIS functions (they return float32_t)
        float32_t w1Re = arm_cos_f32(angle1);
        float32_t w1Im = arm_sin_f32(angle1);
        float32_t w2Re = arm_cos_f32(angle2);
        float32_t w2Im = arm_sin_f32(angle2);
        float32_t w3Re = arm_cos_f32(angle3);
        float32_t w3Im = arm_sin_f32(angle3);

        int a = 2 * (k + j);
        int b = 2 * (k + j + butterflies);
        int c = 2 * (k + j + 2 * butterflies);
        int d = 2 * (k + j + 3 * butterflies);

        float32_t aRe = data[a];
        float32_t aIm = data[a + 1];
        float32_t bRe = data[b] * w1Re - data[b + 1] * w1Im;
        float32_t bIm = data[b] * w1Im + data[b + 1] * w1Re;
        float32_t cRe = data[c] * w2Re - data[c + 1] * w2Im;
        float32_t cIm = data[c] * w2Im + data[c + 1] * w2Re;
        float32_t dRe = data[d] * w3Re - data[d + 1] * w3Im;
        float32_t dIm = data[d] * w3Im + data[d + 1] * w3Re;

        float32_t t0Re = aRe + cRe;
        float32_t t0Im = aIm + cIm;
        float32_t t1Re = aRe - cRe;
        float32_t t1Im = aIm - cIm;
        float32_t t2Re = bRe + dRe;
        float32_t t2Im = bIm + dIm;
        float32_t t3Re = bIm - dIm;
        float32_t t3Im = dRe - bRe;

        data[a]     = t0Re + t2Re;
        data[a + 1] = t0Im + t2Im;
        data[b]     = t1Re + t3Re;
        data[b + 1] = t1Im + t3Im;
        data[c]     = t0Re - t2Re;
        data[c + 1] = t0Im - t2Im;
        data[d]     = t1Re - t3Re;
        data[d + 1] = t1Im - t3Im;
      }
    }
  }
}

// ----------------------- Radix-2 FFT (inline twiddles via CMSIS) -----------------------
void fftRadix2(float32_t* data, int n) {
  bitReverse(data, n);
  int stages = 0;
  {
    int tmp = n;
    while (tmp > 1) { tmp >>= 1; stages++; }
  }

  for (int s = 1; s <= stages; s++) {
    int m = 1 << s;            // block size
    int halfm = m >> 1;
    int step = n / m;

    for (int k = 0; k < n; k += m) {
      for (int j = 0; j < halfm; j++) {
        float32_t angle = -2.0f * PI * ((float)j * (float)step) / (float)n;
        float32_t wRe = arm_cos_f32(angle);
        float32_t wIm = arm_sin_f32(angle);

        int i0 = 2 * (k + j);
        int i1 = 2 * (k + j + halfm);

        float32_t tRe = wRe * data[i1] - wIm * data[i1 + 1];
        float32_t tIm = wRe * data[i1 + 1] + wIm * data[i1];

        float32_t uRe = data[i0];
        float32_t uIm = data[i0 + 1];

        data[i0]     = uRe + tRe;
        data[i0 + 1] = uIm + tIm;
        data[i1]     = uRe - tRe;
        data[i1 + 1] = uIm - tIm;
      }
    }
  }
}

// ----------------------- Hybrid runner -----------------------
void runFFT(float32_t* input, int size) {
  if (size <= CMSIS_MAX) {
    // Use CMSIS RFFT (real input -> complex packed output)
    // Prepare output buffer sized for complex pairs (arm_rfft_fast_f32 expects output len == size)
    static float32_t cmsisOut[CMSIS_MAX]; // enough for <= CMSIS_MAX
    arm_rfft_fast_f32(&cmsisFFT, input, cmsisOut, 0);

    for (int i = 0; i < size / 2; i++) {
      float32_t re = cmsisOut[2*i];
      float32_t im = cmsisOut[2*i + 1];
      magnitude[i] = re * re + im * im;
    }

    // find peak bin
    int peakBin = 1;
    float32_t peakMag = magnitude[1];
    for (int i = 2; i < size / 2; i++) {
      if (magnitude[i] > peakMag) { peakMag = magnitude[i]; peakBin = i; }
    }
    float32_t Re = cmsisOut[2*peakBin];
    float32_t Im = cmsisOut[2*peakBin + 1];
    float32_t mag = sqrtf(Re*Re + Im*Im);
    Serial.print("CMSIS RFFT - peak bin: "); Serial.print(peakBin);
    Serial.print("  freq(Hz): "); Serial.print((peakBin * SAMPLING_FREQ) / (float)size, 2);
    Serial.print("  normMag: "); Serial.println(mag * (2.0f / size), 6);
  } else {
    // Convert real input to complex (interleaved)
    for (int i = 0; i < size; i++) {
      fftBuffer[2*i]     = input[i];
      fftBuffer[2*i + 1] = 0.0f;
    }

    if (!isPowerOfTwo(size)) {
      Serial.println("ERROR: FFT_SIZE must be power of two!");
      return;
    }

    if (isPowerOfFour(size)) {
      fftRadix4(fftBuffer, size);
      Serial.println("Used custom Radix-4 FFT");
    } else {
      fftRadix2(fftBuffer, size);
      Serial.println("Used custom Radix-2 FFT");
    }

    for (int i = 0; i < size / 2; i++) {
      float32_t re = fftBuffer[2*i];
      float32_t im = fftBuffer[2*i + 1];
      magnitude[i] = re*re + im*im;
    }

    int peakBin = 1;
    float32_t peakMag = magnitude[1];
    for (int i = 2; i < size / 2; i++) {
      if (magnitude[i] > peakMag) { peakMag = magnitude[i]; peakBin = i; }
    }
    float32_t Re = fftBuffer[2*peakBin];
    float32_t Im = fftBuffer[2*peakBin + 1];
    float32_t mag = sqrtf(Re*Re + Im*Im);
    Serial.print("Custom FFT - peak bin: "); Serial.print(peakBin);
    Serial.print("  freq(Hz): "); Serial.print((peakBin * SAMPLING_FREQ) / (float)size, 2);
    Serial.print("  normMag: "); Serial.println(mag * (2.0f / size), 6);
  }
}

// ----------------------- Peak detection with parabolic interpolation -----------------------
void detectPeakFreq(int size) {
  // Finds peak in magnitude[] (ignores DC and Nyquist) and does parabolic interpolation
  uint32_t startIdx = 1;
  uint32_t endIdx = (size / 2) - 1;
  float32_t peakVal = 0.0f;
  uint32_t peakIdx = 0;
  for (uint32_t i = startIdx; i < endIdx; i++) {
    if (magnitude[i] > peakVal) { peakVal = magnitude[i]; peakIdx = i; }
  }
  float32_t refined = (float32_t)peakIdx;
  if (peakIdx > 0 && peakIdx < (uint32_t)((size/2) - 1)) {
    float32_t magL = magnitude[peakIdx - 1];
    float32_t magR = magnitude[peakIdx + 1];
    float32_t magC = magnitude[peakIdx];
    float32_t denom = magL - 2.0f * magC + magR;
    if (fabsf(denom) > 1e-12f) {
      float32_t delta = 0.5f * (magL - magR) / denom;
      refined += delta;
    }
  }
  float32_t freq = refined * SAMPLING_FREQ / (float)size;
  float32_t normMag = sqrtf(peakVal) * (2.0f / (float)size);
  Serial.print("PEAK_FREQ(Hz): "); Serial.print(freq, 2);
  Serial.print("  PEAK_MAG: "); Serial.println(normMag, 6);
}

// ----------------------- Setup / Loop -----------------------
void setup() {
  Serial.begin(115200);
  while (!Serial) { } // wait for serial console
  analogReadResolution(12);
  analogReference(3.3); // or EXTERNAL / INTERNAL depending on your board

  if (arm_rfft_fast_init_f32(&cmsisFFT, CMSIS_MAX) != ARM_MATH_SUCCESS) {
    Serial.println("CMSIS RFFT init failed!");
    while (1) { delay(1000); }
  }

  Serial.println("Hybrid FFT ready (No Window, inline CMSIS twiddles).");
}

void loop() {
  Serial.println("\n--- New frame ---");

  uint32_t t0 = micros();
  injectTestSine(inputBuffer, FFT_SIZE, TEST_FREQ_HZ, TEST_AMP);
  uint32_t t1 = micros();
  Serial.print("Sampling time (ms): ");
  Serial.println((t1 - t0) / 1000.0f, 3);

  uint32_t t2 = micros();
  runFFT(inputBuffer, FFT_SIZE);
  uint32_t t3 = micros();
  Serial.println("DATA_START");

  // Send time domain
  for (int i = 0; i < FFT_SIZE; i++) {
    Serial.print(inputBuffer[i], 6);
    if (i < FFT_SIZE - 1) Serial.print(",");
  }
  Serial.println();

  // Send frequency domain (linear magnitude)
  for (int i = 0; i < FFT_SIZE / 2; i++) {
    Serial.print(sqrtf(magnitude[i]), 6); // linear magnitude
    if (i < FFT_SIZE / 2 - 1) Serial.print(",");
  }
  Serial.println();

  Serial.println("DATA_END");

  Serial.print("FFT time (ms): ");
  Serial.println((t3 - t2) / 1000.0f, 3);

  detectPeakFreq(FFT_SIZE);

  delay(500); // small gap between frames
}


