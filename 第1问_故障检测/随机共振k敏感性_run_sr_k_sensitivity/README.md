# 随机共振 k 敏感性脚本说明

- 脚本：`run_sr_k_sensitivity.m`
- 用途：分析随机共振输入增益 \(k\) 对目标频率响应的影响。
- 输入：`../../data.xlsx` 第 1 个 sheet。
- 输出图片：`fig_sr_k_sensitivity.png`、`fig_sr_k_best_cases.png`。
- 输出数据：`sr_k_sensitivity_results.csv`、`sr_k_sensitivity_summary.csv`。
- 正文引用：当前 5.1 正文未直接引用本组图片，作为随机共振参数敏感性材料保留。
- 简略代码逻辑：固定随机共振基础参数，扫描输入增益 \(k\)，统计目标频率投影、输出 RMS 和增益指标。
