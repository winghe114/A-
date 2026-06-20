# 全长投影模型脚本说明

- 脚本：`run_full_length_projection_model.m`
- 用途：生成第 1 问主模型的全长正弦投影频率估计、相关检测曲线和对比图。
- 输入：`../../data.xlsx` 第 1 个 sheet。
- 输出图片：`fig_projection_energy_curve.png`、`fig_fft_vs_projection.png`、`fig_correlation_glrt_curve.png`。
- 输出数据：`full_length_projection_results.csv`、`projection_search_curve.csv`、`correlation_glrt_curve.csv`。
- 正文引用：5.1 正文引用投影能量曲线和 FFT/投影对比图。
- 简略代码逻辑：先用全长 FFT 粗定位，再在局部频带上最大化正弦投影能量，并输出频率、幅值、相位和相关检测结果。
