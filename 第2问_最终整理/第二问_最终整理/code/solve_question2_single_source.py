import json
import math
from pathlib import Path

import numpy as np
import openpyxl
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "work" / "xlsx" / "data.xlsx"
OUTPUT_DIR = ROOT / "outputs"
OUT = OUTPUT_DIR / "question2_single_source_fit.json"
TIME_FIG = OUTPUT_DIR / "q2_time_domain_comparison_paper_cn.png"
FREQ_FIG = OUTPUT_DIR / "q2_frequency_domain_validation_paper_cn.png"
RESIDUAL_FIG = OUTPUT_DIR / "q2_residual_signal_paper_cn.png"
FONT_REGULAR = Path(r"C:\Windows\Fonts\times.ttf")
FONT_BOLD = Path(r"C:\Windows\Fonts\timesbd.ttf")
FONT_CN_REGULAR = Path(r"C:\Windows\Fonts\simsun.ttc")
FONT_CN_BOLD = Path(r"C:\Windows\Fonts\simhei.ttf")


def load_single_source():
    wb = openpyxl.load_workbook(DATA, data_only=True, read_only=True)
    ws = wb.worksheets[0]
    t, x = [], []
    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[0] is None or row[1] is None:
            continue
        t.append(float(row[0]))
        x.append(float(row[1]))
    return np.asarray(t, dtype=float), np.asarray(x, dtype=float)


def fit_sine(t, x, f0):
    omega = 2 * np.pi * f0
    design = np.column_stack([
        np.sin(omega * t),
        np.cos(omega * t),
        np.ones_like(t),
    ])
    coef, _, _, _ = np.linalg.lstsq(design, x, rcond=None)
    a, b, c = [float(v) for v in coef]
    yhat = design @ coef
    residual = x - yhat
    amp = math.sqrt(a * a + b * b)
    phase = math.atan2(b, a)
    phase_mod = (phase + np.pi) % (2 * np.pi) - np.pi
    return {
        "a_sin": a,
        "b_cos": b,
        "offset": c,
        "amplitude": amp,
        "phase_rad": float(phase_mod),
        "phase_deg": float(phase_mod * 180 / np.pi),
        "fit": yhat,
        "residual": residual,
    }


def metrics(x, yhat, residual, k_params):
    n = len(x)
    sst = float(np.sum((x - np.mean(x)) ** 2))
    sse = float(np.sum(residual ** 2))
    mse = sse / n
    rmse = math.sqrt(mse)
    mae = float(np.mean(np.abs(residual)))
    r2 = 1 - sse / sst
    noise_std = float(np.std(residual, ddof=k_params))
    signal_rms = float(np.sqrt(np.mean((yhat - np.mean(yhat)) ** 2)))
    residual_rms = float(np.sqrt(np.mean(residual ** 2)))
    return {
        "n": int(n),
        "sse": sse,
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "r2": float(r2),
        "signal_rms": signal_rms,
        "residual_rms": residual_rms,
        "estimated_snr_db": float(20 * np.log10(signal_rms / (residual_rms + 1e-300))),
        "residual_mean": float(np.mean(residual)),
        "residual_std": noise_std,
    }


def spectrum(signal, fs):
    y = signal - np.mean(signal)
    win = np.hanning(len(y))
    freq = np.fft.rfftfreq(len(y), 1 / fs)
    power = np.abs(np.fft.rfft(y * win)) ** 2
    return freq, power


def nearest_index(freq, target):
    return int(np.argmin(np.abs(freq - target)))


def frequency_metrics(raw_x, yhat, residual, fs, f0):
    freq_raw, power_raw = spectrum(raw_x, fs)
    freq_fit, power_fit = spectrum(yhat, fs)
    freq_res, power_res = spectrum(residual, fs)
    k_raw = nearest_index(freq_raw, f0)
    k_fit = nearest_index(freq_fit, f0)
    k_res = nearest_index(freq_res, f0)

    band = 0.05
    fit_band = (freq_fit >= f0 - band) & (freq_fit <= f0 + band)
    fit_energy_ratio = float(np.sum(power_fit[fit_band]) / (np.sum(power_fit) + 1e-300))
    peak_suppression_ratio = float(power_raw[k_raw] / (power_res[k_res] + 1e-300))

    return {
        "target_frequency_hz": float(f0),
        "raw_peak_power_at_f0": float(power_raw[k_raw]),
        "recovered_peak_power_at_f0": float(power_fit[k_fit]),
        "residual_peak_power_at_f0": float(power_res[k_res]),
        "peak_suppression_ratio": peak_suppression_ratio,
        "recovered_energy_concentration_ratio": fit_energy_ratio,
    }


def residual_spectrum_peak(freq, power):
    mask = (freq > 0.05) & (freq < freq[-1])
    k = np.where(mask)[0][np.argmax(power[mask])]
    return {
        "peak_hz": float(freq[k]),
        "peak_power": float(power[k]),
    }


def load_font(size, bold=False):
    font_path = FONT_BOLD if bold else FONT_REGULAR
    return ImageFont.truetype(str(font_path), size=size)


def load_cn_font(size, bold=False):
    font_path = FONT_CN_BOLD if bold else FONT_CN_REGULAR
    return ImageFont.truetype(str(font_path), size=size)


def format_tick(value):
    if abs(value) < 1e-10:
        value = 0.0
    if abs(value - round(value)) < 1e-10:
        return str(int(round(value)))
    if abs(value) >= 1:
        return f"{value:.1f}"
    return f"{value:.2f}"


def make_ticks(vmin, vmax, count=6):
    if abs(vmax - vmin) < 1e-12:
        vmax = vmin + 1.0
    return np.linspace(vmin, vmax, count)


def data_rect(panel_rect):
    left, top, right, bottom = panel_rect
    return (left + 110, top + 55, right - 35, bottom - 85)


def map_x(x, rect, x_min, x_max):
    left, _, right, _ = rect
    width = max(1.0, right - left)
    return left + (float(x) - x_min) / (x_max - x_min + 1e-300) * width


def map_y(y, rect, y_min, y_max):
    _, top, _, bottom = rect
    height = max(1.0, bottom - top)
    return bottom - (float(y) - y_min) / (y_max - y_min + 1e-300) * height


def scale_points(xs, ys, rect, x_min=None, x_max=None, y_min=None, y_max=None):
    x_min = float(xs[0])
    x_max = float(xs[-1]) if x_max is None else float(x_max)
    x_min = float(xs[0]) if x_min is None else float(x_min)
    if y_min is None:
        y_min = float(np.min(ys))
    if y_max is None:
        y_max = float(np.max(ys))
    if abs(y_max - y_min) < 1e-12:
        y_max = y_min + 1.0
    pts = []
    for x, y in zip(xs, ys):
        pts.append((map_x(x, rect, x_min, x_max), map_y(y, rect, y_min, y_max)))
    return pts


def draw_panel(draw, panel_rect, title, x_label, y_label, x_ticks, y_ticks, x_min, x_max, y_min, y_max):
    left, top, right, bottom = panel_rect
    plot_rect = data_rect(panel_rect)
    p_left, p_top, p_right, p_bottom = plot_rect

    title_font = load_cn_font(28, bold=True)
    label_font = load_cn_font(24)
    tick_font = load_font(22)

    draw.rectangle([left, top, right, bottom], fill="white")
    draw.rectangle([p_left, p_top, p_right, p_bottom], outline="black", width=2)

    for xt in x_ticks:
        px = map_x(xt, plot_rect, x_min, x_max)
        draw.line([(px, p_bottom), (px, p_bottom + 10)], fill="black", width=2)
        draw.line([(px, p_top), (px, p_bottom)], fill=(220, 220, 220), width=1)
        text = format_tick(float(xt))
        bbox = draw.textbbox((0, 0), text, font=tick_font)
        draw.text((px - (bbox[2] - bbox[0]) / 2, p_bottom + 16), text, fill="black", font=tick_font)

    for yt in y_ticks:
        py = map_y(yt, plot_rect, y_min, y_max)
        draw.line([(p_left - 10, py), (p_left, py)], fill="black", width=2)
        draw.line([(p_left, py), (p_right, py)], fill=(220, 220, 220), width=1)
        text = format_tick(float(yt))
        bbox = draw.textbbox((0, 0), text, font=tick_font)
        draw.text((p_left - 20 - (bbox[2] - bbox[0]), py - (bbox[3] - bbox[1]) / 2), text, fill="black", font=tick_font)

    title_bbox = draw.textbbox((0, 0), title, font=title_font)
    draw.text((left + (right - left - (title_bbox[2] - title_bbox[0])) / 2, top + 8), title, fill="black", font=title_font)

    xlabel_bbox = draw.textbbox((0, 0), x_label, font=label_font)
    draw.text((p_left + (p_right - p_left - (xlabel_bbox[2] - xlabel_bbox[0])) / 2, bottom - 48), x_label, fill="black", font=label_font)

    ylabel_bbox = draw.textbbox((0, 0), y_label, font=label_font)
    draw.text((left + 10, p_top - 36), y_label, fill="black", font=label_font)
    return plot_rect


def draw_series(draw, rect, xs, ys, color, width=2, x_min=None, x_max=None, y_min=None, y_max=None):
    pts = scale_points(xs, ys, rect, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
    draw.line(pts, fill=color, width=width)


def draw_legend(draw, items, anchor):
    x0, y0 = anchor
    font = load_cn_font(22)
    line_len = 36
    row_h = 34
    text_width = max(draw.textbbox((0, 0), label, font=font)[2] for _, label in items)
    box_w = 30 + line_len + 16 + text_width + 20
    box_h = 18 + row_h * len(items)
    draw.rectangle([x0, y0, x0 + box_w, y0 + box_h], fill="white", outline="black", width=1)
    for i, (color, label) in enumerate(items):
        y = y0 + 16 + i * row_h + 10
        draw.line([(x0 + 16, y), (x0 + 16 + line_len, y)], fill=color, width=4)
        draw.text((x0 + 16 + line_len + 16, y - 13), label, fill="black", font=font)


def draw_vertical_marker(draw, plot_rect, x_value, x_min, x_max, label, color):
    px = map_x(x_value, plot_rect, x_min, x_max)
    _, top, _, bottom = plot_rect
    step = 12
    y = top
    while y < bottom:
        draw.line([(px, y), (px, min(y + step // 2, bottom))], fill=color, width=2)
        y += step
    font = load_font(22, bold=True)
    draw.text((px + 8, top + 8), label, fill=color, font=font)


def nice_limits(y_min, y_max, pad_ratio=0.08):
    span = y_max - y_min
    if span < 1e-12:
        span = 1.0
    pad = span * pad_ratio
    return y_min - pad, y_max + pad


def create_time_domain_figure(t, raw_x, recovered_s, fs):
    n_zoom = min(len(t), int(round(5 * fs)))
    xs = t[:n_zoom]
    raw = raw_x[:n_zoom]
    rec = recovered_s[:n_zoom]

    img = Image.new("RGB", (1600, 950), "white")
    draw = ImageDraw.Draw(img)
    panel = (70, 50, 1540, 890)
    x_min, x_max = float(xs[0]), float(xs[-1])
    y_min, y_max = nice_limits(float(min(np.min(raw), np.min(rec))), float(max(np.max(raw), np.max(rec))))
    plot_rect = draw_panel(
        draw,
        panel,
        "观测信号与恢复信号时域对比",
        "时间 / s",
        "幅值",
        x_ticks=np.arange(0, 5.1, 0.5),
        y_ticks=make_ticks(y_min, y_max, 7),
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
    )
    draw_series(draw, plot_rect, xs, raw, (110, 110, 110), width=2, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
    draw_series(draw, plot_rect, xs, rec, (200, 30, 30), width=3, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
    draw_legend(draw, [((110, 110, 110), "观测信号 x(t)"), ((200, 30, 30), "恢复信号 s(t)")], (1135, 105))
    img.save(TIME_FIG, dpi=(300, 300))


def create_frequency_domain_figure(freq_raw, power_raw, freq_fit, power_fit, freq_res, power_res, f0):
    mask = freq_raw <= 10.0
    xs = freq_raw[mask]
    raw = power_raw[mask]
    fit = power_fit[mask]
    res = power_res[mask]

    raw_n = raw / (np.max(raw) + 1e-300)
    fit_n = fit / (np.max(fit) + 1e-300)
    res_n = res / (np.max(res) + 1e-300)

    img = Image.new("RGB", (1600, 1180), "white")
    draw = ImageDraw.Draw(img)
    rect_top = (70, 45, 1540, 545)
    rect_bottom = (70, 610, 1540, 1110)
    x_ticks = np.arange(0, 10.1, 1.0)
    y_ticks = np.linspace(0, 1.0, 6)
    top_plot = draw_panel(
        draw,
        rect_top,
        "原始频谱与恢复频谱对比",
        "频率 / Hz",
        "归一化功率",
        x_ticks=x_ticks,
        y_ticks=y_ticks,
        x_min=0.0,
        x_max=10.0,
        y_min=0.0,
        y_max=1.0,
    )
    bottom_plot = draw_panel(
        draw,
        rect_bottom,
        "残差信号频谱",
        "频率 / Hz",
        "归一化功率",
        x_ticks=x_ticks,
        y_ticks=y_ticks,
        x_min=0.0,
        x_max=10.0,
        y_min=0.0,
        y_max=1.0,
    )
    draw_series(draw, top_plot, xs, raw_n, (110, 110, 110), width=2, x_min=0.0, x_max=10.0, y_min=0.0, y_max=1.0)
    draw_series(draw, top_plot, xs, fit_n, (200, 30, 30), width=3, x_min=0.0, x_max=10.0, y_min=0.0, y_max=1.0)
    draw_series(draw, bottom_plot, xs, res_n, (30, 100, 200), width=2, x_min=0.0, x_max=10.0, y_min=0.0, y_max=1.0)
    draw_vertical_marker(draw, top_plot, f0, 0.0, 10.0, f"f0 = {f0:.2f} Hz", (0, 0, 0))
    draw_vertical_marker(draw, bottom_plot, f0, 0.0, 10.0, f"f0 = {f0:.2f} Hz", (0, 0, 0))
    draw_legend(draw, [((110, 110, 110), "原始信号频谱"), ((200, 30, 30), "恢复信号频谱")], (1140, 110))
    draw_legend(draw, [((30, 100, 200), "残差信号频谱")], (1180, 680))
    img.save(FREQ_FIG, dpi=(300, 300))


def create_residual_figure(t, residual, fs):
    n_zoom = min(len(t), int(round(5 * fs)))
    xs = t[:n_zoom]
    ys = residual[:n_zoom]

    img = Image.new("RGB", (1600, 950), "white")
    draw = ImageDraw.Draw(img)
    panel = (70, 50, 1540, 890)
    x_min, x_max = float(xs[0]), float(xs[-1])
    y_min, y_max = nice_limits(float(np.min(ys)), float(np.max(ys)))
    plot_rect = draw_panel(
        draw,
        panel,
        "恢复周期分量后的残差信号",
        "时间 / s",
        "残差幅值",
        x_ticks=np.arange(0, 5.1, 0.5),
        y_ticks=make_ticks(y_min, y_max, 7),
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
    )
    draw_series(draw, plot_rect, xs, ys, (30, 100, 200), width=2, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
    draw_legend(draw, [((30, 100, 200), "残差信号 r(t)")], (1180, 105))
    img.save(RESIDUAL_FIG, dpi=(300, 300))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    t, x = load_single_source()
    fs = 1.0 / np.median(np.diff(t))
    f0_from_q1 = 2.0
    initial = fit_sine(t, x, f0_from_q1)
    yhat = initial["fit"]
    residual = initial["residual"]
    recovered_s = yhat - initial["offset"]

    m = metrics(x, yhat, residual, 3)
    freq_raw, power_raw = spectrum(x, fs)
    freq_fit, power_fit = spectrum(recovered_s, fs)
    freq_res, power_res = spectrum(residual, fs)
    fmetrics = frequency_metrics(x, recovered_s, residual, fs, f0_from_q1)

    create_time_domain_figure(t, x, recovered_s, fs)
    create_frequency_domain_figure(freq_raw, power_raw, freq_fit, power_fit, freq_res, power_res, f0_from_q1)
    create_residual_figure(t, residual, fs)

    result = {
        "question": "Q2 single-source waveform recovery",
        "data": {
            "samples": int(len(x)),
            "fs_hz": round(float(fs), 6),
            "duration_s": round(float(t[-1] - t[0]), 6),
        },
        "model": {
            "fitted_linear_form": "x(t)=a sin(2*pi*f0*t)+b cos(2*pi*f0*t)+c+e(t)",
            "recovered_signal_form": "s_hat(t)=A sin(2*pi*f0*t+phi)",
            "phase_relation": "A=sqrt(a^2+b^2), phi=atan2(b,a)",
            "f0_hz": f0_from_q1,
        },
        "estimated_parameters": {
            "a_sin": round(initial["a_sin"], 10),
            "b_cos": round(initial["b_cos"], 10),
            "offset_c": round(initial["offset"], 10),
            "amplitude_A": round(initial["amplitude"], 10),
            "phase_phi_rad": round(initial["phase_rad"], 10),
            "phase_phi_deg": round(initial["phase_deg"], 6),
        },
        "recovered_signal": {
            "expression": "s_hat(t)=0.0351903 sin(4*pi*t + 1.57555)",
            "with_offset_expression": "x_hat(t)=0.0351903 sin(4*pi*t + 1.57555)+4.379e-5",
        },
        "observable_error_metrics": {
            "note": "The noiseless true s(t) samples are not provided; these metrics are residual errors against observed x(t), not direct true-signal errors.",
            "rmse": round(m["rmse"], 9),
            "mae": round(m["mae"], 9),
            "r2": round(m["r2"], 9),
            "estimated_snr_db": round(m["estimated_snr_db"], 6),
            "residual_mean": round(m["residual_mean"], 12),
            "residual_std": round(m["residual_std"], 9),
        },
        "frequency_domain_evaluation": {
            "target_frequency_hz": round(fmetrics["target_frequency_hz"], 6),
            "raw_peak_power_at_f0": round(fmetrics["raw_peak_power_at_f0"], 6),
            "recovered_peak_power_at_f0": round(fmetrics["recovered_peak_power_at_f0"], 6),
            "residual_peak_power_at_f0": round(fmetrics["residual_peak_power_at_f0"], 6),
            "peak_suppression_ratio": round(fmetrics["peak_suppression_ratio"], 6),
            "recovered_energy_concentration_ratio": round(fmetrics["recovered_energy_concentration_ratio"], 9),
            "residual_main_peak_hz": round(residual_spectrum_peak(freq_res, power_res)["peak_hz"], 6),
        },
        "figures": {
            "time_domain": str(TIME_FIG),
            "frequency_domain": str(FREQ_FIG),
            "residual_time_domain": str(RESIDUAL_FIG),
        },
        "true_error_handling": {
            "status": "not directly observable from the provided attachment",
            "reason": "the problem states the form of s(t), but does not provide the true A, phi, or noiseless s(t) sequence",
            "recommended_paper_statement": "Use least-squares residuals, peak suppression at f0, and recovered-signal energy concentration as observable recovery evidence; use same-SNR simulation if a direct true-signal error metric is required.",
        },
    }
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
