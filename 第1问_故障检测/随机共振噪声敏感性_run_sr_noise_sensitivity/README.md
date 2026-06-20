# 随机共振噪声敏感性脚本说明

- 脚本：`run_sr_noise_sensitivity.m`
- 用途：分析额外加入高斯噪声时随机共振增强效果是否改善。
- 输入：`../../data.xlsx` 第 1 个 sheet。
- 输出图片：`fig_sr_noise_sensitivity.png`、`fig_sr_best_added_noise_preview.png`。
- 输出数据：`sr_noise_sensitivity_results.csv`。
- 正文引用：当前 5.1 正文未直接引用本组图片，作为随机共振噪声敏感性材料保留。
- 简略代码逻辑：扫描不同人工噪声强度，比较随机共振输出在目标频率附近的投影响应和信噪指标。
