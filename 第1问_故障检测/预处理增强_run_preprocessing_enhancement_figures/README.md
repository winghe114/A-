# 预处理增强脚本说明

- 脚本：`run_preprocessing_enhancement_figures.m`
- 用途：生成 5.1.1 中随机共振增强和一阶自相关增强的正式图表与数据。
- 输入：`../../data.xlsx` 第 1 个 sheet。
- 输出图片：`fig_sr_bistable_potential.png`、`fig_sr_time_spectrum_comparison.png`、`fig_acf_time_spectrum_comparison.png`。
- 输出数据：`stochastic_resonance_output.csv`、`autocorrelation_order1.csv`。
- 正文引用：5.1 正文引用本组 3 张图片。
- 简略代码逻辑：读取单故障源信号，去均值并标准化，分别生成随机共振输出和一阶自相关信号，再绘制处理前后时域与频域对比。
