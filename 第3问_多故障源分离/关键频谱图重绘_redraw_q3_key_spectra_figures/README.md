# 关键频谱图重绘脚本说明

- 脚本：`redraw_q3_key_spectra_figures.m`
- 用途：只根据已有 CSV 快速重绘第 3 问两张关键频谱图，避免重新运行完整主模型和近频仿真。
- 输入：`../多故障源分离_run_q3_multisource_separation/` 中的 `q3_residual_projection_spectra.csv`、`q3_omp_history.csv`、`q3_multisource_results.csv`。
- 输出图片：覆盖主分离组中的 `fig_q3_original_spectrum_identified.png` 和 `fig_q3_omp_residual_projection_spectra.png`。
- 输出数据：无新增数据文件。
- 正文引用：正文引用的图片实体位于主分离组，本组只保存重绘工具。
- 简略代码逻辑：读取主分离结果表，绘制 \(J_1(f)\) 投影能量谱和每步只标当前新峰的 OMP 残差投影谱。
