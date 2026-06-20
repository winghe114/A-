import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def fit_at_freq(t, x, freq):
    tc = t - t.mean()
    w = 2 * np.pi * freq
    a = np.column_stack([np.ones_like(t), tc, np.sin(w * t), np.cos(w * t)])
    coef, *_ = np.linalg.lstsq(a, x, rcond=None)
    fit = a @ coef
    sse = float(np.sum((x - fit) ** 2))
    offset, slope, b, c = coef
    amp = float(np.hypot(b, c))
    phase = float(np.arctan2(c, b))
    return sse, offset, slope, amp, phase


def refine_freq(t, x, f0, df, min_freq):
    center = f0
    half_width = max(3 * df, f0 * 0.02)
    best = None
    for _ in range(4):
        lo = max(min_freq, center - half_width)
        hi = center + half_width
        for f in np.linspace(lo, hi, 81):
            row = fit_at_freq(t, x, float(f))
            if best is None or row[0] < best[0]:
                best = (row[0], float(f), *row[1:])
        center = best[1]
        half_width /= 8
    return best


def estimate_sheet(df, sheet_name, min_freq=0.01, top_n=8):
    t = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(float)
    x = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(float)
    ok = np.isfinite(t) & np.isfinite(x)
    t, x = t[ok], x[ok]
    order = np.argsort(t)
    t, x = t[order], x[order]

    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    n = len(x)
    tc = t - t.mean()
    trend = np.polyval(np.polyfit(tc, x, 1), tc)
    y = x - trend
    win = np.hanning(n)
    spec = np.fft.rfft(y * win)
    freqs = np.fft.rfftfreq(n, dt)
    fft_amp = 2 * np.abs(spec) / win.sum()
    dfreq = freqs[1] - freqs[0]

    mask = freqs >= min_freq
    search_freqs = freqs[mask]
    search_amp = fft_amp[mask]
    local_max = np.r_[False, (search_amp[1:-1] > search_amp[:-2]) & (search_amp[1:-1] >= search_amp[2:]), False]
    peak_idx = np.flatnonzero(local_max)
    if len(peak_idx) == 0:
        peak_idx = np.arange(len(search_amp))
    ranked = peak_idx[np.argsort(search_amp[peak_idx])[-top_n:]][::-1]

    rows = []
    best = None
    for idx in ranked:
        f0 = float(search_freqs[idx])
        sse, freq, offset, slope, amp, phase = refine_freq(t, x, f0, dfreq, min_freq)
        residual_rms = float(np.sqrt(sse / n))
        signal_rms = amp / np.sqrt(2)
        snr_db = float(20 * np.log10(signal_rms / residual_rms)) if residual_rms > 0 else np.inf
        phase = float(np.arctan2(np.sin(phase), np.cos(phase)))
        row = {
            "sheet": sheet_name,
            "samples": n,
            "dt": dt,
            "fs_hz": fs,
            "freq_hz": freq,
            "amplitude": amp,
            "phase_rad": phase,
            "phase_deg": float(np.degrees(phase)),
            "offset": float(offset),
            "slope": float(slope),
            "residual_rms": residual_rms,
            "snr_db": snr_db,
            "nearest_fft_freq_hz": f0,
            "nearest_fft_amp": float(search_amp[idx]),
        }
        rows.append(row)
        if best is None or row["residual_rms"] < best["residual_rms"]:
            best = row

    return best, pd.DataFrame(rows), freqs, fft_amp


def save_plot(sheet_name, freqs, amp, best):
    import matplotlib.pyplot as plt

    safe_name = "".join(ch if ch.isalnum() else "_" for ch in sheet_name)
    plt.figure(figsize=(10, 5))
    plt.semilogy(freqs[1:], amp[1:])
    plt.axvline(best["freq_hz"], color="r", linestyle="--", label=f"{best['freq_hz']:.6g} Hz")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude")
    plt.title(f"FFT spectrum - {sheet_name}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"fft_spectrum_{safe_name}.png", dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="FFT analysis for weak sinusoidal signal in data.xlsx")
    parser.add_argument("xlsx", nargs="?", default="data.xlsx")
    parser.add_argument("--min-freq", type=float, default=0.01)
    parser.add_argument("--top-n", type=int, default=8)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    path = Path(args.xlsx)
    all_sheets = pd.read_excel(path, sheet_name=None)
    all_candidates = []
    best_rows = []

    for sheet_name, df in all_sheets.items():
        best, candidates, freqs, amp = estimate_sheet(df, sheet_name, args.min_freq, args.top_n)
        best_rows.append(best)
        all_candidates.append(candidates)
        print(f"\n[{sheet_name}]")
        print(f"sampling interval dt = {best['dt']:.10g} s, fs = {best['fs_hz']:.10g} Hz, samples = {best['samples']}")
        print(f"frequency = {best['freq_hz']:.10g} Hz")
        print(f"amplitude = {best['amplitude']:.10g}")
        print(f"initial phase = {best['phase_rad']:.10g} rad = {best['phase_deg']:.10g} deg")
        print(f"residual RMS = {best['residual_rms']:.10g}, SNR = {best['snr_db']:.3f} dB")
        if not args.no_plot:
            save_plot(sheet_name, freqs, amp, best)

    pd.DataFrame(best_rows).to_csv("fft_best_results.csv", index=False, encoding="utf-8-sig")
    pd.concat(all_candidates, ignore_index=True).to_csv("fft_candidates.csv", index=False, encoding="utf-8-sig")
    print("\nSaved: fft_best_results.csv, fft_candidates.csv")
    if not args.no_plot:
        print("Saved: fft_spectrum_*.png")


if __name__ == "__main__":
    main()
