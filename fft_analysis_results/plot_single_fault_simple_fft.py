from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

base = Path(__file__).resolve().parent.parent
xlsx = base / 'data.xlsx'
xl = pd.ExcelFile(xlsx)
sheet = xl.sheet_names[0]
df = pd.read_excel(xlsx, sheet_name=sheet)

t = df.iloc[:, 0].to_numpy(float)
x = df.iloc[:, 1].to_numpy(float)

dt = t[1] - t[0]
n = len(x)
fs = 1 / dt

X = np.fft.fft(x)
freq = np.fft.fftfreq(n, dt)

pos = freq > 0
freq_pos = freq[pos]
amp_pos = 2 * np.abs(X[pos]) / n

peak = np.argmax(amp_pos)
print('sheet:', sheet)
print('fs:', fs)
print('peak frequency:', freq_pos[peak])
print('peak amplitude:', amp_pos[peak])

plt.figure(figsize=(11, 7))

plt.subplot(2, 1, 1)
plt.plot(t, x)
plt.xlabel('Time / s')
plt.ylabel('x(t)')
plt.title('Time domain')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(freq_pos, amp_pos)
plt.axvline(freq_pos[peak], color='r', linestyle='--', label=f'peak = {freq_pos[peak]:.6g} Hz')
plt.xlabel('Frequency / Hz')
plt.ylabel('Amplitude')
plt.title('Simple FFT amplitude spectrum')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
