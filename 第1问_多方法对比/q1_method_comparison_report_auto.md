# 第1问多方法频率检测结果分析

采样频率为 100.000000 Hz，观测时长为 400.000000 s，全长理论频率分辨率为 0.002500000000 Hz。

普通全长 FFT 的最大频点为 1.999950001252 Hz；本文以全长投影/GLRT 的 2.000001615205 Hz 作为主参考结果。

## 方法对比表

| 方法 | 频率/Hz | 能量占比 | 与投影法差值/Hz | 评价 |
|---|---:|---:|---:|---|
| Plain full-length FFT bin | 1.99995000125 | 0.0582149274411 | -5.16139536115e-05 | FFT bin only; has grid bias. |
| Quadratic interpolated FFT | 2.00012431776 | 0.0578483223584 | 0.000122702549699 | Low-cost correction of FFT grid bias. |
| Zero-padded FFT x16 | 2.00000210225 | 0.0582940223597 | 4.87047161268e-07 | Densifies the spectrum display but does not add information. |
| Full-length projection / GLRT | 2.00000161521 | 0.0582940294129 | 0 | Best main method; continuous-frequency coherent detection. |
| Dense cos-sin correlation scan | 2.000001615 | 0.0582940294129 | -2.03383088149e-10 | Equivalent to matched filtering; peak score 24.7675. |
| AR(24) prewhitened projection | 2.00000385968 | 0.058293879738 | 2.24447710062e-06 | Tests robustness after reducing colored noise correlation. |
| SSA denoise rank-2 projection | 2.00000536258 | 0.0582936122231 | 3.74737082742e-06 | Data-adaptive low-rank denoising, window 400. |
| Adaptive stochastic resonance + projection | 1.99995336981 | 0.0582249115061 | -4.82453980915e-05 | Nonlinear enhancement; best gain 0.2733, output SNR 28.922 dB. |
| Autocorrelation order-1 + projection | 2.00000532873 | 0.0582936197251 | 3.71352232342e-06 | Time-domain periodicity enhancement; autocorr-domain fraction 0.951217, local SNR 58.755 dB. |
| Autocorrelation order-2 + projection | 2.00000763812 | 0.0582929517864 | 6.02291025764e-06 | Time-domain periodicity enhancement; autocorr-domain fraction 0.732558, local SNR 82.471 dB. |
| Autocorrelation order-3 + projection | 2.00001135741 | 0.0582912100695 | 9.74220897065e-06 | Time-domain periodicity enhancement; autocorr-domain fraction 0.586827, local SNR 79.121 dB. |
| MUSIC pseudospectrum | 1.806640625 | 4.91058895854e-05 | -0.193360990204 | Subspace high-resolution spectrum; toolbox dependent. |

## 结论

1. 全长投影/GLRT 与密集相关扫描本质一致，均利用 cos/sin 双参考信号进行相干积累，结果最稳定，适合作为第一问主方法。
2. 插值 FFT 和零填充 FFT 能修正普通 FFT 的频点栅栏误差，但统计判决和参数恢复能力不如投影/GLRT 完整。
3. AR 预白化可检验有色噪声影响。若其频率结果与主方法一致，可作为鲁棒性证据。
4. SSA 去噪可作为数据驱动预处理，但可能改变幅值结构，适合作辅助展示，不建议替代主频率估计。
5. 随机共振输出也能检测到接近 2 Hz 的频率，但增强优势有限。当前最佳输入增益 0.273333，输出投影占比 0.0649052，局部 SNR 28.9224 dB，适合作非线性辅助验证。

综合建议：第一问主线采用“全长投影搜索 + GLRT/相关检测”；插值 FFT、AR 预白化、SSA 和随机共振作为对比或鲁棒性分析。
