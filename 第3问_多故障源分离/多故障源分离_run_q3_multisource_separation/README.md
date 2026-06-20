# 多故障源分离脚本说明

- 脚本：`run_q3_multisource_separation.m`
- 用途：实现第 3 问主模型，即全长正弦投影-OMP、多频联合最小二乘、稳健阈值源数判定和近频热力图辅助仿真。
- 输入：`../../data.xlsx` 第 2 个 sheet。
- 输出图片：`fig_q3_original_spectrum_identified.png`、`fig_q3_omp_residual_projection_spectra.png`、`fig_q3_time_reconstruction_short.png`、`fig_q3_separated_components.png`、`fig_q3_frequency_before_after.png`、`fig_q3_near_frequency_resolution_heatmap.png`。
- 输出数据：`q3_multisource_results.csv`、`q3_omp_history.csv`、`q3_residual_projection_spectra.csv`、`q3_summary_metrics.csv`、`q3_summary.txt`、`q3_threshold_sensitivity.csv`、`q3_simulation_resolution_results.csv`、`q3_resolution_success_matrix.csv`。
- 正文引用：5.3 正文引用本组除近频热力图外的主要分离图。
- 简略代码逻辑：在残差投影谱中逐次寻找最显著频率，每加入一个频率后联合重拟合所有分量，并用 median+MAD 阈值自动停止。
