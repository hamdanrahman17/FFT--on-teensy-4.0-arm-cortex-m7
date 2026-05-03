#include <Arduino.h>
#include <arm_math.h>  // CMSIS DSP

#define FFT_SIZE 8192// Set 4096 or 16384
#define CMSIS_MAX 4096

#define SAMPLING_FREQ 10000.0f  // Hz
#define TEST_FREQ 250.0f
#define AMP 0.8f

float32_t inputBuffer[FFT_SIZE];
float32_t window[FFT_SIZE];
float32_t magnitude[FFT_SIZE / 2];
float32_t fftBuffer[2 * FFT_SIZE];  // [Re0, Im0, Re1, Im1, ...]

bool isPowerOfTwo(int n) {
  return (n > 0) && ((n & (n - 1)) == 0);
}

bool isPowerOfFour(int n) {
  return isPowerOfTwo(n) && (__builtin_ctz(n) % 2 == 0);
}


arm_rfft_fast_instance_f32 cmsisFFT;

// ===============================
// Hann Window
// ===============================
void generateHannWindow() {
  for (int i = 0; i < FFT_SIZE; i++) {
    window[i] = 0.5f * (1.0f - cosf(2.0f * PI * i / (FFT_SIZE - 1)));
  }
}

void applyHann(float32_t* buffer, int size) {
  for (int i = 0; i < size; i++) {
    buffer[i] *= window[i];
  }
}

// ===============================
// Inject Test Sine
// ===============================
void injectTestSine(float32_t* buffer, int size, float freq, float amp) {
  for (int i = 0; i < size; i++) {
    buffer[i] = amp * sinf(2.0f * PI * freq * i / SAMPLING_FREQ);
  }
}

// ===============================
// Bit-Reversal for FFT
// ===============================
void bitReverse(float32_t* data, int n) {
  int j = 0;
  for (int i = 0; i < n; i++) {
    if (i < j) {
      float32_t tempRe = data[2 * i];
      float32_t tempIm = data[2 * i + 1];
      data[2 * i] = data[2 * j];
      data[2 * i + 1] = data[2 * j + 1];
      data[2 * j] = tempRe;
      data[2 * j + 1] = tempIm;
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
        result |= (x & 0x3);  // extract 2 bits (one base-4 digit)
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
            float32_t tempRe = data[2 * i];
            float32_t tempIm = data[2 * i + 1];
            data[2 * i] = data[2 * j];
            data[2 * i + 1] = data[2 * j + 1];
            data[2 * j] = tempRe;
            data[2 * j + 1] = tempIm;
        }
    }
}
// ===============================
// Radix-4 FFT
// ===============================
void fftRadix4(float32_t* data, int n) {
   bitReverseRadix4(data, n);

  for (int m = 4; m <= n; m <<= 2) {
    float32_t theta = -2.0f * PI / m;
    float32_t w1Re = cosf(theta);
    float32_t w1Im = sinf(theta);
    float32_t w2Re = cosf(2 * theta);
    float32_t w2Im = sinf(2 * theta);
    float32_t w3Re = cosf(3 * theta);
    float32_t w3Im = sinf(3 * theta);

    for (int k = 0; k < n; k += m) {
      float32_t wRe = 1.0f, wIm = 0.0f;
      float32_t w2r = 1.0f, w2i = 0.0f;
      float32_t w3r = 1.0f, w3i = 0.0f;

      for (int j = 0; j < m / 4; j++) {
        int a = 2 * (k + j);
        int b = 2 * (k + j + m / 4);
        int c = 2 * (k + j + m / 2);
        int d = 2 * (k + j + 3 * m / 4);

        float32_t aRe = data[a];
        float32_t aIm = data[a + 1];
        float32_t bRe = data[b] * wRe - data[b + 1] * wIm;
        float32_t bIm = data[b] * wIm + data[b + 1] * wRe;
        float32_t cRe = data[c] * w2r - data[c + 1] * w2i;
        float32_t cIm = data[c] * w2i + data[c + 1] * w2r;
        float32_t dRe = data[d] * w3r - data[d + 1] * w3i;
        float32_t dIm = data[d] * w3i + data[d + 1] * w3r;

        float32_t t0Re = aRe + cRe;
        float32_t t0Im = aIm + cIm;
        float32_t t1Re = aRe - cRe;
        float32_t t1Im = aIm - cIm;
        float32_t t2Re = bRe + dRe;
        float32_t t2Im = bIm + dIm;
        float32_t t3Re = bIm - dIm;
        float32_t t3Im = dRe - bRe;

        data[a] = t0Re + t2Re;
        data[a + 1] = t0Im + t2Im;
        data[b] = t1Re + t3Re;
        data[b + 1] = t1Im + t3Im;
        data[c] = t0Re - t2Re;
        data[c + 1] = t0Im - t2Im;
        data[d] = t1Re - t3Re;
        data[d + 1] = t1Im - t3Im;

        float wr = wRe * w1Re - wIm * w1Im;
        float wi = wRe * w1Im + wIm * w1Re;
        wRe = wr; wIm = wi;

        float w2t = w2r * w2Re - w2i * w2Im;
        w2i = w2r * w2Im + w2i * w2Re;
        w2r = w2t;

        float w3t = w3r * w3Re - w3i * w3Im;
        w3i = w3r * w3Im + w3i * w3Re;
        w3r = w3t;
      }
    }
  }
}
void fftRadix2(float32_t* data, int n) {
  bitReverse(data, n);  // Assuming you already have this implemented

  for (int s = 1; s <= log2(n); s++) {
    int m = 1 << s;
    float32_t wmRe = cosf(-2.0f * PI / m);
    float32_t wmIm = sinf(-2.0f * PI / m);

    for (int k = 0; k < n; k += m) {
      float32_t wRe = 1.0f;
      float32_t wIm = 0.0f;

      for (int j = 0; j < m / 2; j++) {
        int i0 = 2 * (k + j);
        int i1 = 2 * (k + j + m / 2);

        float32_t tRe = wRe * data[i1] - wIm * data[i1 + 1];
        float32_t tIm = wRe * data[i1 + 1] + wIm * data[i1];

        float32_t uRe = data[i0];
        float32_t uIm = data[i0 + 1];

        data[i0] = uRe + tRe;
        data[i0 + 1] = uIm + tIm;
        data[i1] = uRe - tRe;
        data[i1 + 1] = uIm - tIm;

        float32_t tempRe = wRe * wmRe - wIm * wmIm;
        wIm = wRe * wmIm + wIm * wmRe;
        wRe = tempRe;
      }
    }
  }
}

// ===============================
// Hybrid FFT Execution
// ===============================
void runFFT(float32_t* input, int size) {
  if (size <= CMSIS_MAX) {
    float32_t fftOut[CMSIS_MAX];
    arm_rfft_fast_f32(&cmsisFFT, input, fftOut, 0);

    for (int i = 0; i < size / 2; i++) {
      float32_t re = fftOut[2 * i];
      float32_t im = fftOut[2 * i + 1];
      magnitude[i] = re * re + im * im;
    }

    // Output CMSIS result
    int expectedBin = round((TEST_FREQ * size) / SAMPLING_FREQ);
    float Re = fftOut[2 * expectedBin];
    float Im = fftOut[2 * expectedBin + 1];
    float mag = sqrtf(Re * Re + Im * Im);

    Serial.print("CMSIS RFFT - Bin: "); Serial.println(expectedBin);
    Serial.print("Re: "); Serial.print(Re);
    Serial.print(", Im: "); Serial.print(Im);
    Serial.print(", Mag: "); Serial.println(mag);
    Serial.print("Normalized Magnitude: ");
    Serial.println(mag * (2.0f / size));

  } else {
    // Convert real input to complex format
    for (int i = 0; i < size; i++) {
      fftBuffer[2 * i] = input[i];
      fftBuffer[2 * i + 1] = 0.0f;
    }

    if (!isPowerOfTwo(size)) {
      Serial.println("ERROR: FFT size must be a power of 2!");
      return;
    }

    if (isPowerOfFour(size)) {
      fftRadix4(fftBuffer, size);
      Serial.println("Used Custom Radix-4 FFT");
    } else {
      fftRadix2(fftBuffer, size);
      Serial.println("Used Custom Radix-2 FFT");
    }

    // Compute magnitude
    for (int i = 0; i < size / 2; i++) {
      float32_t re = fftBuffer[2 * i];
      float32_t im = fftBuffer[2 * i + 1];
      magnitude[i] = re * re + im * im;
    }

    // Output result
    int expectedBin = round((TEST_FREQ * size) / SAMPLING_FREQ);
    float Re = fftBuffer[2 * expectedBin];
    float Im = fftBuffer[2 * expectedBin + 1];
    float mag = sqrtf(Re * Re + Im * Im);
    float normMag = mag * (2.0f / size);

    Serial.print("Custom FFT - Bin: "); Serial.println(expectedBin);
    Serial.print("Re: "); Serial.print(Re);
    Serial.print(", Im: "); Serial.print(Im);
    Serial.print(", Mag: "); Serial.println(mag);
    Serial.print("Normalized Magnitude: ");Serial.println(normMag);
  }
}




// ===============================
// Peak Detection
// ===============================
void detectPeakFreq(int size) {
  float32_t peakVal = 0.0f;
  uint32_t peakIdx = 0;

  // Ignore DC (bin 0) and Nyquist (bin size/2 - 1) to avoid noise/artifacts
  uint32_t startIdx = 1;
  uint32_t endIdx = (size / 2) - 1;

  for (uint32_t i = startIdx; i < endIdx; i++) {
    if (magnitude[i] > peakVal) {
      peakVal = magnitude[i];
      peakIdx = i;
    }
  }

  // Parabolic interpolation for more precise frequency estimate
  float32_t refinedIdx = (float32_t)peakIdx;

    if (peakIdx > 0 && peakIdx < (uint32_t)((size / 2) - 1))  { {
    float32_t magL = magnitude[peakIdx - 1];
    float32_t magR = magnitude[peakIdx + 1];
    float32_t magC = magnitude[peakIdx];

    float32_t denominator = magL - 2.0f * magC + magR;
    if (fabsf(denominator) > 1e-6f) {
      float32_t delta = 0.5f * (magL - magR) / denominator;
      refinedIdx += delta;  // Peak is refined to sub-bin accuracy
    }
  }

  float32_t freq;
 if (size <= CMSIS_MAX) {
    // CMSIS RFFT
    freq = refinedIdx * SAMPLING_FREQ / size;
 } else {
    // Radix-4 FFT
    freq = refinedIdx * SAMPLING_FREQ /  size;
 }

  Serial.print("PEAK_FREQ (Hz): ");
  Serial.print(freq, 2);
  Serial.print(", PEAK_MAG: ");
  float32_t normMag = sqrtf(peakVal) * (2.0f / size);
  Serial.println(normMag, 4);
 }

}
// ===============================
// Setup
// ===============================
void setup() {
  Serial.begin(115200);
  while (!Serial);
  analogReadResolution(12);
  analogReference(3.3);

  if (arm_rfft_fast_init_f32(&cmsisFFT, CMSIS_MAX) != ARM_MATH_SUCCESS) {
    Serial.println("CMSIS FFT Init failed!");
    while (1);
  }

  generateHannWindow();
  Serial.println("Teensy Hybrid FFT Ready");
}

// ===============================
// Loop
// ===============================
void loop() {
  Serial.println("\nSampling...");
  uint32_t t0 = micros();
  injectTestSine(inputBuffer, FFT_SIZE, TEST_FREQ, AMP);
  uint32_t t1 = micros();

  applyHann(inputBuffer, FFT_SIZE);

  Serial.println("Running FFT...");
  uint32_t t2 = micros();
  runFFT(inputBuffer, FFT_SIZE);
  uint32_t t3 = micros();
    // Send data to Python
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


  Serial.print("Sampling Time (ms): ");
  Serial.println((t1 - t0) / 1000.0f, 2);
  Serial.print("FFT Time (ms): ");
  Serial.println((t3 - t2) / 1000.0f, 2);

  detectPeakFreq(FFT_SIZE);

  delay(5000);
} 