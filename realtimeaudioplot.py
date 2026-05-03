import serial
import numpy as np
import matplotlib.pyplot as plt

# Match to Arduino
FFT_SIZE = 16384
FS = 10000.0  # Sampling frequency

ser = serial.Serial("/dev/tty.usbmodem160990301", 115200)  # change port

# Prepare figure
plt.ion()
fig, (ax_time, ax_freq) = plt.subplots(1, 2, figsize=(12, 5))

t = np.arange(FFT_SIZE) / FS
freqs = np.linspace(0, FS/2, FFT_SIZE//2, endpoint=False)

line_time, = ax_time.plot(t, np.zeros(FFT_SIZE))
ax_time.set_title("Input Signal")
ax_time.set_xlabel("Time [s]")
ax_time.set_ylabel("Amplitude")
ax_time.grid(True)

line_freq, = ax_freq.plot(freqs, np.zeros(FFT_SIZE//2))
ax_freq.set_title("FFT Magnitude")
ax_freq.set_xlabel("Frequency [Hz]")
ax_freq.set_ylabel("Magnitude")
ax_freq.grid(True)

# Loop
while True:
    line = ser.readline().decode().strip()

    if line == "DATA_START":
        # Time-domain data
        td_line = ser.readline().decode().strip()
        time_data = np.array([float(x) for x in td_line.split(",")])

        # Frequency-domain data
        fd_line = ser.readline().decode().strip()
        freq_data = np.array([float(x) for x in fd_line.split(",")])

        # End marker
        end_line = ser.readline().decode().strip()
        if end_line == "DATA_END":
            # Update plots
            line_time.set_ydata(time_data)
            line_freq.set_ydata(freq_data)

            ax_time.set_ylim(np.min(time_data)*1.1, np.max(time_data)*1.1)
            ax_freq.set_ylim(0, np.max(freq_data)*1.1)

            fig.canvas.draw()
            fig.canvas.flush_events()
