from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

base = Path(__file__).resolve().parent.parent
xlsx = base / 'data.xlsx'
xl = pd.ExcelFile(xlsx)
sheet = next((s for s in xl.sheet_names if '单' in s or 'Դ' in s), xl.sheet_names[0])
df = pd.read_excel(xlsx, sheet_name=sheet)

t = pd.to_numeric(df.iloc[:, 0], errors='coerce').to_numpy(float)
x = pd.to_numeric(df.iloc[:, 1], errors='coerce').to_numpy(float)
ok = np.isfinite(t) & np.isfinite(x)
t, x = t[ok], x[ok]
order = np.argsort(t)
t, x = t[order], x[order]

dt = np.median(np.diff(t))
n = len(x)
y = x - np.mean(x)
win = np.hanning(n)
X = np.fft.rfft(y * win)
freq = np.fft.rfftfreq(n, dt)
amp = 2 * np.abs(X) / win.sum()
peak = np.argmax(amp[1:]) + 1

fig, axes = plt.subplots(2, 1, figsize=(11, 7))

axes[0].plot(t, x, linewidth=0.8)
axes[0].set_title(f'{sheet} - time domain')
axes[0].set_xlabel('Time / s')
axes[0].set_ylabel('x(t)')
axes[0].grid(True, alpha=0.3)

axes[1].plot(freq[1:], amp[1:], linewidth=0.8)
axes[1].axvline(freq[peak], color='r', linestyle='--', label=f'peak = {freq[peak]:.6g} Hz')
axes[1].set_title(f'{sheet} - frequency domain')
axes[1].set_xlabel('Frequency / Hz')
axes[1].set_ylabel('Amplitude')
axes[1].grid(True, alpha=0.3)
axes[1].legend()

plt.tight_layout()
plt.show()
