#include <Arduino.h>
#include <arm_math.h>  // CMSIS DSP

#define FFT_SIZE 16384// Set 4096 or 16384
#define CMSIS_MAX 4096
#define ANALOG_PIN A0
#define SAMPLING_FREQ 10000 // Hz
#define AMP 1.0f

float32_t inputBuffer[FFT_SIZE];
float32_t window[CMSIS_MAX];
float32_t magnitude[FFT_SIZE / 2];
float32_t fftBuffer[2 * FFT_SIZE]; 
const float PEAK_THRESHOLD = 0.0001f;// [Re0, Im0, Re1, Im1, ...]
float32_t twiddleRe[FFT_SIZE / 2];
float32_t twiddleIm[FFT_SIZE / 2];

arm_rfft_fast_instance_f32 cmsisFFT;

bool isPowerOfTwo(int n) {
  return (n > 0) && ((n & (n - 1)) == 0);
}

bool isPowerOfFour(int n) {
  return isPowerOfTwo(n) && (__builtin_ctz(n) % 2 == 0);
}

// ===============================
// Hann Window
// ===============================
void generateHannWindow() {
  for (int i = 0; i < CMSIS_MAX; i++) {
    window[i] = 0.5f * (1.0f - cosf(2.0f * PI * i / (CMSIS_MAX - 1)));
  }
}

void applyHann(float32_t* buffer, int size) {
  int winSize = min(size, CMSIS_MAX);
  for (int i = 0; i < winSize; i++) {
    buffer[i] *= window[i];
  }
}


void sampleAudio(float32_t* buffer, int size) {
  uint32_t period = 1000000 / SAMPLING_FREQ;
  float32_t sum = 0.0f;

  // Step 1: Capture raw samples
  for (int i = 0; i < size; i++) {
    float32_t val = (float32_t)analogRead(ANALOG_PIN);
    buffer[i] = val;
    sum += val;
    delayMicroseconds(period);
  }

  // Step 2: Remove DC offset
  float32_t mean = sum / size;
  for (int i = 0; i < size; i++) {
    buffer[i] = (buffer[i] - mean) / 2048.0f; // normalize to -1..+1
  }

  // Step 3: Add a small noise floor in silence
  const float NOISE_LEVEL = 0.002f; // ~ -54 dB
  for (int i = 0; i < size; i++) {
    if (fabsf(buffer[i]) < NOISE_LEVEL) {
      buffer[i] += ((float)rand() / RAND_MAX - 0.5f) * NOISE_LEVEL;
    }
  }

  // Step 4: Optional auto-gain for low signals
  const float AUTO_GAIN = 2.0f; // 2x boost
  for (int i = 0; i < size; i++) {
    buffer[i] *= AUTO_GAIN;
    if (buffer[i] > 1.0f) buffer[i] = 1.0f;     // clamp
    if (buffer[i] < -1.0f) buffer[i] = -1.0f;   // clamp
  }
}

  // PRECOMPUTED TWIDDILE FACTORS//
void computeTwiddles(int fftSize) {
    for (int k = 0; k < fftSize / 2; k++) {
        float32_t angle = -2.0f * PI * k / fftSize;
        twiddleRe[k] = cosf(angle);
        twiddleIm[k] = sinf(angle);
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
        int step = FFT_SIZE / m;

        for (int k = 0; k < n; k += m) {
            for (int j = 0; j < m / 4; j++) {
                float32_t w1Re = twiddleRe[j * step];
                float32_t w1Im = twiddleIm[j * step];
                float32_t w2Re = twiddleRe[(2 * j) * step];
                float32_t w2Im = twiddleIm[(2 * j) * step];
                float32_t w3Re = twiddleRe[(3 * j) * step];
                float32_t w3Im = twiddleIm[(3 * j) * step];

                int a = 2 * (k + j);
                int b = 2 * (k + j + m / 4);
                int c = 2 * (k + j + m / 2);
                int d = 2 * (k + j + 3 * m / 4);

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

void fftRadix2(float32_t* data, int n) {
    bitReverse(data, n);
    int stages = log2(n);

    for (int s = 1; s <= stages; s++) {
        int m = 1 << s;             // FFT size at this stage
        int step = FFT_SIZE / m;    // Twiddle factor step size

        for (int k = 0; k < n; k += m) {
            for (int j = 0; j < m / 2; j++) {
                float32_t wRe = twiddleRe[j * step];
                float32_t wIm = twiddleIm[j * step];

                int i0 = 2 * (k + j);
                int i1 = 2 * (k + j + m / 2);

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
    for (int i = 0; i < size / 2; i++) {
      float32_t re = fftOut[2 * i];
      float32_t im = fftOut[2 * i + 1];
      magnitude[i] = re * re + im * im;
    }

    // Find peak bin dynamically
    int peakBin = 1;
    float32_t peakMag = magnitude[1];
    for (int i = 2; i < size / 2; i++) {
      if (magnitude[i] > peakMag) {
        peakMag = magnitude[i];
        peakBin = i;
      }
    }
    float Re = fftOut[2 * peakBin];
    float Im = fftOut[2 * peakBin + 1];
    float mag = sqrtf(Re * Re + Im * Im);

    Serial.print("CMSIS RFFT - Bin: "); Serial.println(peakBin);
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
    for (int i = 0; i < size / 2; i++) {
      float32_t re = fftBuffer[2 * i];
      float32_t im = fftBuffer[2 * i + 1];
      magnitude[i] = re * re + im * im;
    }

    // Find peak bin dynamically
    int peakBin = 1;
    float32_t peakMag = magnitude[1];
    for (int i = 2; i < size / 2; i++) {
      if (magnitude[i] > peakMag) {
        peakMag = magnitude[i];
        peakBin = i;
      }
    }
    float Re = fftBuffer[2 * peakBin];
    float Im = fftBuffer[2 * peakBin + 1];
    float mag = sqrtf(Re * Re + Im * Im);
    float normMag = mag * (2.0f / size);

    Serial.print("Custom FFT - Bin: "); Serial.println(peakBin);
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
  computeTwiddles(FFT_SIZE);
  while (!Serial);
  analogReadResolution(12);
  analogReference(3.3);
  

  if (arm_rfft_fast_init_f32(&cmsisFFT, CMSIS_MAX) != ARM_MATH_SUCCESS) {
    Serial.println("CMSIS FFT Init failed!");
    while (1);
  }

  //generateHannWindow();
  Serial.println("Teensy Hybrid FFT Ready");
}

// ===============================
// Loop
// ===============================
void loop() {
  Serial.println("\nSampling...");
  uint32_t t0 = micros();
  //injectTestSine(inputBuffer, FFT_SIZE, TEST_FREQ, AMP);
  sampleAudio(inputBuffer, FFT_SIZE);
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

  delay(10);
}