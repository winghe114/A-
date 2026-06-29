# MWS 随机共振脚本说明

- 脚本：`run_mws_stochastic_resonance.m`
- 用途：作为第 1 问辅助增强实验，尝试 MWS 势函数随机共振模型。
- 输入：`../../data.xlsx` 第 1 个 sheet。
- 输出图片：`fig_mws_potential.png`、`fig_mws_time_spectrum_comparison.png`、`fig_mws_parameter_scan.png`、`fig_mws_vs_classic_sr.png`。
- 输出数据：`mws_sr_parameter_scan_results.csv`、`mws_sr_best_output.csv`、`mws_sr_summary.txt`。
- 正文引用：5.1.1 正文引用 `fig_mws_potential.png` 和 `fig_mws_time_spectrum_comparison.png`，其余图表作为辅助对比材料保留。
- 简略代码逻辑：对 MWS 随机共振参数进行有限网格搜索，选取投影响应较好的输出，并与经典随机共振分支进行对比。
