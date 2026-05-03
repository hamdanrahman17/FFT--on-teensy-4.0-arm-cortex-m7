import serial
import numpy as np
import matplotlib.pyplot as plt

# Match this to Arduino
FFT_SIZE = 4096
FS = 10000.0  # Sampling frequency

ser = serial.Serial("/dev/tty.usbmodem160990301", 115200)

time_data = []
freq_data = []

while True:
    line = ser.readline().decode().strip()

    if line == "DATA_START":
        # Read time-domain data
        td_line = ser.readline().decode().strip()
        time_data = np.array([float(x) for x in td_line.split(",")])

        # Read frequency-domain data
        fd_line = ser.readline().decode().strip()
        freq_data = np.array([float(x) for x in fd_line.split(",")])

        # Read footer
        end_line = ser.readline().decode().strip()
        if end_line == "DATA_END":
            # Plot
            t = np.arange(FFT_SIZE) / FS
            freqs = np.linspace(0, FS/2, FFT_SIZE//2, endpoint=False)

            plt.figure(figsize=(12,5))
            plt.subplot(1, 2, 1)
            plt.plot(t, time_data)
            plt.title("Input Signal")
            plt.xlabel("Time [s]")
            plt.ylabel("Amplitude")
            plt.grid(True)

            plt.subplot(1, 2, 2)
            plt.plot(freqs, freq_data)
            plt.title("FFT Magnitude")
            plt.xlabel("Frequency [Hz]")
            plt.ylabel("Magnitude")
            plt.grid(True)

            plt.tight_layout()
            plt.show()
