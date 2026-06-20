# 近频超分辨实验脚本说明

- 脚本：`run_q3_superresolution_experiment.m`
- 用途：比较传统全长 FFT 双峰拾取与局部双频变量投影精修在近频分辨中的表现。
- 输入：`../../data.xlsx` 第 2 个 sheet，仅用于取得采样时间轴。
- 输出图片：`fig_q3_superresolution_fft_vs_varpro.png`。
- 输出数据：`q3_superresolution_comparison.csv`、`q3_superresolution_fft_success_matrix.csv`、`q3_superresolution_varpro_success_matrix.csv`、`q3_superresolution_summary.txt`。
- 正文引用：5.3.5 正文引用本组 FFT 与变量投影近频分辨率对比图。
- 简略代码逻辑：构造两个真实频率为 \(10\) Hz 与 \(10+\Delta f\) 的仿真信号，扫描 \(\Delta f\) 和 SNR，统计两频率是否被不同估计频率正确匹配。
