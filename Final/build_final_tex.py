from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_generated(name):
    text = (ROOT / "_generated" / name).read_text(encoding="utf-8")
    replacements = {
        "figures/MWS随机共振_run_mws_stochastic_resonance/fig_mws_potential.png": "figures/q1_mws_potential.png",
        "figures/MWS随机共振_run_mws_stochastic_resonance/fig_mws_time_spectrum_comparison.png": "figures/q1_mws_time_spectrum_comparison.png",
        "figures/预处理增强_run_preprocessing_enhancement_figures/fig_acf_time_spectrum_comparison.png": "figures/q1_acf_time_spectrum_comparison.png",
        "figures/全长投影模型_run_full_length_projection_model/fig_projection_energy_curve.png": "figures/q1_projection_energy_curve.png",
        "figures/全长投影模型_run_full_length_projection_model/fig_fft_vs_projection.png": "figures/q1_fft_vs_projection.png",
        "figures/多故障源分离_run_q3_multisource_separation/fig_q3_original_spectrum_identified.png": "figures/q3_original_spectrum_identified.png",
        "figures/多故障源分离_run_q3_multisource_separation/fig_q3_omp_residual_projection_spectra.png": "figures/q3_omp_residual_projection_spectra.png",
        "figures/多故障源分离_run_q3_multisource_separation/fig_q3_time_reconstruction_short.png": "figures/q3_time_reconstruction_short.png",
        "figures/多故障源分离_run_q3_multisource_separation/fig_q3_frequency_before_after.png": "figures/q3_frequency_before_after.png",
        "figures/近频超分辨实验_run_q3_superresolution_experiment/fig_q3_superresolution_fft_vs_varpro.png": "figures/q3_superresolution_fft_vs_varpro.png",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


PREAMBLE = r"""% !TeX program = xelatex
\documentclass[UTF8,a4paper,12pt]{ctexart}

\usepackage{amsmath,amssymb}
\usepackage{array}
\usepackage{booktabs}
\usepackage{caption}
\usepackage{float}
\usepackage{geometry}
\usepackage{graphicx}
\usepackage{listings}
\usepackage{longtable}
\usepackage{setspace}
\usepackage{tabularx}
\usepackage{tocloft}
\usepackage{xcolor}

\setmainfont{Times New Roman}
\setCJKmainfont{SimSun}[AutoFakeBold=2.5]
\setCJKsansfont{SimHei}
\setCJKmonofont{NSimSun}

\geometry{
  a4paper,
  top=2.54cm,
  bottom=2.54cm,
  left=2.70cm,
  right=2.70cm
}

\setlength{\parindent}{2em}
\setlength{\parskip}{0pt}
\linespread{1.25}
\raggedbottom

\captionsetup{font=small,labelsep=quad}
\renewcommand{\figurename}{图}
\renewcommand{\tablename}{表}
\numberwithin{figure}{section}
\numberwithin{table}{section}
\renewcommand{\cftsecleader}{\cftdotfill{\cftdotsep}}
\lstset{
  basicstyle=\ttfamily\small,
  breaklines=true,
  columns=fullflexible,
  frame=single,
  keepspaces=true,
  numbers=left,
  numberstyle=\tiny,
  showstringspaces=false,
  keywordstyle=\bfseries\color{blue!50!black},
  commentstyle=\color{gray!80!black}
}

\newcommand{\fillline}[1]{\underline{\makebox[#1][l]{}}}
\newcommand{\fronttitle}[1]{\begin{center}{\songti\zihao{-2}#1}\end{center}}
"""

FRONT_AND_MID = r"""
\begin{document}
\pagestyle{plain}

% ==================== 论文模板第 1 页：承诺书 ====================
\thispagestyle{empty}

\begin{center}
  {\songti\zihao{-2}2026年NUDT-CSU数学建模联赛}

  \vspace{0.8cm}
  {\songti\zihao{2}承\quad 诺\quad 书}
\end{center}

\vspace{0.8cm}

我们完全明白，在竞赛开始后参赛队员不能以任何方式（包括电话、电子邮件、网上咨询等）与本队以外的任何人（包括指导教师）研究、讨论与赛题有关的问题。

我们知道，抄袭别人的成果和买卖论文、代码都是违反竞赛规则的，如果引用别人的成果或其它公开的资料（包括网上查到的资料），必须按照规定的参考文献的表述方式在正文引用处和参考文献中明确列出；按要求使用AI开展竞赛。

我们郑重承诺，严格遵守竞赛规则，以保证竞赛的公正、公平性。如有违反竞赛规则的行为，我们愿意承担由此引起的一切后果。

\vspace{1.2cm}

\noindent 所属学校（学校全称）：\fillline{9.5cm}

\vspace{0.6cm}

\noindent 参赛队员：

\vspace{0.4cm}

\noindent 队员1姓名：\fillline{4.0cm}\quad 学号：\fillline{3.2cm}

\vspace{0.4cm}

\noindent 队员2姓名：\fillline{4.0cm}\quad 学号：\fillline{3.2cm}

\vspace{0.4cm}

\noindent 队员3姓名：\fillline{4.0cm}\quad 学号：\fillline{3.2cm}

\vfill

\begin{flushright}
日期：\quad 2026 年 6 月 19 日
\end{flushright}

\clearpage

% ==================== 论文模板第 2 页：题目与摘要页 ====================
\thispagestyle{empty}

\noindent 附件2：

\vspace{0.8cm}

\fronttitle{2026年NUDT-CSU数学建模联赛}

\vspace{1.2cm}

\noindent\hspace*{2em}题\quad 目：\fillline{10.0cm}

\vspace{0.8cm}

\noindent 摘\quad 要：

\vspace{11.2cm}

\noindent 关键词：

\clearpage

\tableofcontents
\clearpage

% ==================== 正文写作区 ====================
\title{论文标题}
\author{作者}
\date{}
\maketitle

\section*{摘要}
% 摘要正文

\section*{关键词}
% 关键词

\section{问题重述}
\subsection{问题背景}
机械设备在长期高速、高负荷运行过程中，齿轮、轴承等关键零部件可能因磨损、疲劳或局部破损而产生故障。当齿轮发生破损时，故障部位会在转动过程中周期性地激发冲击振动，从而在振动加速度信号中形成具有固定特征频率的微弱周期分量。由于机械电机运转、传动系统振动以及外部环境干扰共同作用，实际采集信号往往包含强背景噪声，故障信号能量较弱，容易被噪声淹没，导致直接采用常规时域波形观察或普通频谱分析难以准确判断故障是否存在。

\subsection{问题提出}
题目给出了机械设备振动加速度数据 data.xlsx，其中存在一个微弱正弦周期信号，要求结合信号处理、统计分析及优化理论，建立数学模型来解决以下问题：

\textbf{问题一}：要求在单源故障情形下，利用“单源故障”表中的数据，建立能够从强噪声背景中检测微弱周期信号的模型，准确估计故障特征频率，并说明所用方法相较于传统傅里叶变换的优势。

\textbf{问题二}：要求在已获得故障特征频率的基础上，对原始微弱周期信号进行波形恢复，估计其振幅和初相，给出恢复信号表达式，并从时域、频域及定量误差角度评价恢复效果。

\textbf{问题三}：进一步考虑多处故障同时存在的情形，要求基于前两问模型扩展出多频率分量的自动分离与识别方法，估计各故障源对应的频率、振幅和初相，并通过数值实验讨论故障频率接近时的辨识能力。

\textbf{问题四}：则从检测系统设计角度出发，在最多布置 3 个加速度传感器的条件下，综合考虑信号传播衰减、不同测点灵敏度和噪声空间差异，建立传感器布局优化模型，以提高系统对微弱故障信号的检测概率并降低误报率。

\section{问题分析}
本题本质上是强噪声背景下的微弱周期信号参数估计问题。由于故障信号具有较明确的周期结构，而噪声成分随机性较强，可以将建模重点放在“周期性提取”和“信号参数估计”上。整体思路为：首先对原始振动数据时域进行去均值、标准化等预处理；其次利用频域粗定位和时域相关检测相结合的方法确定故障特征频率；然后基于已估计频率建立参数化正弦模型，通过最小二乘等方法恢复波形；最后将单频模型推广到多频叠加模型，并进一步结合多传感器检测性能指标完成布局优化。

\begin{figure}[H]
\centering
\includegraphics[width=0.95\textwidth]{figures/mermaid-diagram-2026-06-21-161711.png}
\caption{整体建模思路}
\end{figure}

\subsection{问题一的分析}
问题一要求在单源故障条件下识别被强噪声淹没的故障特征频率。普通傅里叶变换虽然能够提供频谱信息，但其频率估计受采样间隔、频率栅格和噪声峰值影响较大，当故障分量能量较弱时，直接寻找频谱最大峰可能产生偏差。因此，可先利用全长 FFT 对可能的故障频带进行粗略定位，再在候选频率附近构造正弦、余弦基函数，对原始信号进行相关投影搜索。

从统计检测角度看，若某一频率处存在周期分量，则观测信号在该频率对应的正交正弦基上的投影能量应显著高于邻近频率。由此可将故障频率估计转化为连续频率范围内投影能量最大化问题。同时，可结合随机共振、自相关增强等方法作为辅助验证，以增强周期结构的可见性和检测可靠性。

\subsection{问题二的分析}
问题二在问题一检测出故障特征频率的基础上，进一步恢复单源故障信号的波形。由于题目给出的故障信号可近似表示为单一正弦周期信号，因此在频率已知或已高精度估计的条件下，未知量主要为振幅和初相。直接估计振幅和初相具有一定非线性，可将正弦信号展开为同频率的正弦项与余弦项线性组合，将问题转化为线性最小二乘参数估计问题，再由系数换算得到振幅和初相，从而得到恢复后的故障信号表达式。

\subsection{问题三的分析}
问题三将单源故障扩展为多源故障，即观测信号中可能同时包含多个不同频率、不同振幅和不同初相的微弱周期分量。此时，如果直接套用单频检测方法，强一些的故障分量可能掩盖弱分量，频率接近时还可能出现谱峰混叠。因此，需要在前两问模型基础上建立多频分量的自动识别与分离方法。

可将多源故障信号表示为若干个正弦分量的叠加，并采用“检测-估计-剔除-再检测”的迭代策略：先通过投影能量谱或稀疏频谱搜索找出最显著的候选频率，再对已选频率进行联合最小二乘拟合，估计各分量的振幅和初相；随后从观测信号中扣除已解释的周期成分，对残差信号继续搜索潜在故障频率，直到新增分量的检测统计量不再显著。为避免误检，还需设定模型阶数选择或显著性判别准则。

对于频率十分接近的多个故障源，模型的辨识能力主要受采样时长、噪声强度、频率间隔和分量幅值差异影响。可通过构造含有已知频率、振幅和相位的仿真信号，在不同频率间隔和信噪比条件下重复实验，统计频率估计误差和分离成功率，从而分析模型的分辨极限。

\subsection{问题四的分析}
问题四关注检测系统层面的优化。前三问主要从算法角度提高弱故障信号的识别能力，而实际检测效果还与传感器安装位置密切相关。不同测点到故障源的传播路径不同，信号衰减程度、结构传递特性和噪声水平也不同，因此同一故障在不同传感器上表现出的信噪比可能存在明显差异。合理的传感器布局应尽可能选择对故障信号敏感、噪声较小且信息互补性较强的位置。

建模时，可将每个候选测点对应的检测性能量化为信号增益、噪声方差、局部信噪比或单点检测概率，并考虑不同传感器之间的相关性。由于最多只能布置 3 个传感器，问题可转化为有限候选测点集合上的组合优化问题：目标是在误报率受控的条件下最大化系统综合检测概率，或等价地最大化融合后的检测统计量信噪比。多传感器信息融合可采用加权投影能量、似然比融合或投票判决等方式，其中权重可由各测点的噪声水平和故障灵敏度确定。

为验证布局方案的有效性，可基于前三问得到的故障检测模型设计仿真实验，在不同信噪比、不同故障位置和不同噪声空间分布条件下，对比最优布局、随机布局和单传感器布局的检测概率与误报率。若最优布局在低信噪比下仍能保持较高检测概率且误报率较低，则说明该布局能够提升系统鲁棒性和可靠性。

\section{模型假设}
\begin{enumerate}
  \item 假设附件中给出的振动加速度数据真实、可靠，采样时间和采样值准确，不存在明显缺失值、重复采样点或人为录入错误。
  \item 假设单源故障情形下，机械设备中仅存在一个主要故障齿轮，其故障产生的特征信号可近似表示为单一正弦周期信号：
  \[
  s(t)=A\sin(2\pi f_0t+\varphi_0).
  \]
  \item 假设采集信号由故障周期信号和背景噪声线性叠加而成，即
  \[
  x(t)=s(t)+n(t),
  \]
  不考虑传感器非线性失真、信号饱和和采集系统量化误差对结果的影响。
  \item 假设在观测时间内，故障特征频率 $f_0$、振幅 $A$ 和初相 $\varphi_0$ 保持不变，即故障信号具有平稳周期性。
  \item 假设背景噪声 $n(t)$ 在观测时间内统计特性近似平稳，且与故障周期信号不相关。因此，通过长时间相干积累、相关检测和最小二乘拟合可以抑制随机噪声影响。
  \item 假设单源故障信号在全观测区间内持续存在，而非短时突发或间歇出现。因此可以使用全长数据进行频率搜索和波形恢复。
\end{enumerate}

\section{符号说明}
\begin{longtable}{>{\centering\arraybackslash}p{3.0cm}p{7.9cm}>{\centering\arraybackslash}p{2.8cm}}
\caption{主要符号说明}\\
\toprule
符号 & 含义 & 单位/说明\\
\midrule
\endfirsthead
\toprule
符号 & 含义 & 单位/说明\\
\midrule
\endhead
$t$ & 时间变量 & s\\
$t_i$ & 第 $i$ 个采样时刻 & s\\
$\Delta t$ & 采样时间间隔 & s\\
$f_s$ & 采样频率 & Hz\\
$N$ & 采样点数 & 1\\
$x_i, x(t_i)$ & 原始观测振动加速度信号 & 附件数据量纲\\
$y_i$ & 去均值后的观测信号 & 附件数据量纲\\
$u_i$ & 去均值并标准化后的信号 & 无量纲\\
$s(t)$ & 故障周期信号 & 附件数据量纲\\
$\hat{s}(t)$ & 恢复得到的故障周期信号 & 附件数据量纲\\
$n(t)$ & 背景噪声信号 & 附件数据量纲\\
$r(t_i)$ & 恢复后的残差信号 & 附件数据量纲\\
$A,\hat A$ & 故障周期信号振幅及其估计值 & 附件数据量纲\\
$f_0,\hat f_0$ & 故障特征频率及其估计值 & Hz\\
$\varphi_0,\hat\varphi_0$ & 故障信号初相及其估计值 & rad\\
$a,b,c$ & 正弦、余弦展开系数和直流偏置 & 待估参数\\
$\mathbf H(f)$ & 候选频率 $f$ 下的正弦/余弦设计矩阵 & 无\\
$\boldsymbol\theta$ & 线性模型参数向量 & 无\\
$J(f)$ & 候选频率 $f$ 下的投影能量 & 信号平方量纲\\
$\varepsilon_i$ & 随机误差项或未建模扰动 & 附件数据量纲\\
$\mathrm{RMSE}$ & 残差均方根误差 & 附件数据量纲\\
$\mathrm{MAE}$ & 平均绝对误差 & 附件数据量纲\\
$R^2$ & 拟合解释率或恢复能量占比 & 无\\
$\mathrm{SNR}$ & 信噪比 & dB\\
$K,\hat K$ & 多故障源个数及其估计值 & 1\\
$f_k,A_k,\varphi_k$ & 第 $k$ 个故障源的频率、振幅和初相 & Hz/量纲/rad\\
$\mathbf r^{(m)}$ & OMP 第 $m$ 次迭代后的残差 & 向量\\
$\boldsymbol\Psi(\mathbf f)$ & 多频正弦/余弦联合设计矩阵 & 无\\
\bottomrule
\end{longtable}

\section{模型的建立与求解}
"""

Q2 = r"""
\subsection{问题二模型的建立与求解}

第一问已得到单源故障特征频率 \(f_0\approx2.00\ \mathrm{Hz}\)。第二问的任务是在频率已知的条件下估计幅值 \(A\) 和初相 \(\varphi_0\)，恢复微弱周期分量
\[
s(t)=A\sin(2\pi f_0t+\varphi_0).
\]
由于附件未给出无噪声真实信号，本文采用时域波形、频域能量集中性以及可观测残差共同评价恢复效果。

\subsubsection{数据预处理}

沿用第一问对“单源故障”工作表的读取结果，数据包含 \(N=40001\) 个采样点，采样频率为 \(f_s=100.00\ \mathrm{Hz}\)。为保留幅值尺度，参数估计直接使用原始观测序列 \(x(t_i)\)，并在模型中显式加入常数偏置项 \(c\)，由最小二乘同时估计。

\subsubsection{正弦最小二乘恢复模型的建立}

在 \(f_0\) 已知时，将正弦信号按同频正弦、余弦基展开：
\[
a=A\cos\varphi_0,\qquad b=A\sin\varphi_0,
\]
\[
s(t)=a\sin(2\pi f_0t)+b\cos(2\pi f_0t).
\]
考虑直流偏置和随机扰动，观测模型为
\[
x(t_i)=a\sin(2\pi f_0t_i)+b\cos(2\pi f_0t_i)+c+\varepsilon_i.
\]

令
\[
\mathbf{x}=\begin{bmatrix}x(t_1)\\x(t_2)\\\vdots\\x(t_N)\end{bmatrix},\quad
\boldsymbol{\theta}=\begin{bmatrix}a\\b\\c\end{bmatrix},\quad
\boldsymbol{\varepsilon}=\begin{bmatrix}\varepsilon_1\\\varepsilon_2\\\vdots\\\varepsilon_N\end{bmatrix},
\]
并构造设计矩阵
\[
\mathbf H=\begin{bmatrix}
\sin(2\pi f_0t_1)&\cos(2\pi f_0t_1)&1\\
\sin(2\pi f_0t_2)&\cos(2\pi f_0t_2)&1\\
\vdots&\vdots&\vdots\\
\sin(2\pi f_0t_N)&\cos(2\pi f_0t_N)&1
\end{bmatrix}.
\]
则
\[
\mathbf{x}=\mathbf H\boldsymbol{\theta}+\boldsymbol{\varepsilon}.
\]
在零均值高斯误差假设下，极大似然估计等价于最小化残差平方和，因此取目标函数
\[
J(\boldsymbol{\theta})=\|\mathbf{x}-\mathbf H\boldsymbol{\theta}\|_2^2.
\]

\subsubsection{模型的求解}

对 \(J(\boldsymbol{\theta})\) 求导并令其为零，得到正规方程
\[
\mathbf H^T\mathbf H\boldsymbol{\theta}=\mathbf H^T\mathbf{x}.
\]
当 \(\mathbf H^T\mathbf H\) 可逆时，线性参数的解析解为
\[
\hat{\boldsymbol{\theta}}=(\mathbf H^T\mathbf H)^{-1}\mathbf H^T\mathbf{x}.
\]
再由 \(\hat a,\hat b\) 反推幅值和初相：
\[
\hat A=\sqrt{\hat a^2+\hat b^2},
\qquad
\hat\varphi_0=\operatorname{atan2}(\hat b,\hat a),\qquad
\hat c=\hat\theta_3.
\]

将 \(f_0=2.000\ \mathrm{Hz}\) 代入设计矩阵，对 \(N=40001\) 个采样点执行最小二乘求解，得到：
\begin{table}[H]
\centering
\caption{问题二线性参数估计结果}
\begin{tabular}{ccc}
\toprule
参数 & 符号 & 估计值\\
\midrule
正弦系数 & $\hat a$ & $-0.0001672284$\\
余弦系数 & $\hat b$ & $0.0351899101$\\
直流偏置 & $\hat c$ & $4.37905\times10^{-5}$\\
\bottomrule
\end{tabular}
\end{table}

由此得到
\[
\hat A=0.0351903,
\qquad
\hat\varphi_0=1.57555\ \mathrm{rad}.
\]
因此，恢复的故障周期信号为
\[
\hat s(t)=0.0351903\sin(4\pi t+1.57555).
\]
含直流偏置的完整拟合信号为
\[
\hat x(t)=0.0351903\sin(4\pi t+1.57555)+4.379\times10^{-5}.
\]

\subsubsection{结果分析}

\paragraph{（1）时域恢复效果}

将原始含噪观测信号与恢复周期信号在前 5 s 时间范围内进行对比。

\begin{figure}[H]
\centering
\includegraphics[width=0.82\textwidth]{figures/q2_time_recovery.png}
\caption{观测信号与恢复信号的时域对比}
\end{figure}

图中灰色曲线为原始含噪观测信号，红色曲线为恢复周期信号。恢复结果呈稳定的 2 Hz 正弦结构，说明模型提取出了目标周期分量。

\paragraph{（2）频域一致性验证}

为验证恢复分量是否对应第一问得到的故障频率，本文对原始信号、恢复信号和残差信号分别进行频谱分析，并比较 2 Hz 处的峰值功率。

\begin{figure}[H]
\centering
\includegraphics[width=0.82\textwidth]{figures/q2_frequency_validation.png}
\caption{原始频谱、恢复频谱与残差频谱的频域一致性验证}
\end{figure}

恢复信号在 \(f_0=2.00\ \mathrm{Hz}\) 处保持显著主峰，而残差信号在目标频率处的峰值明显减弱，说明 2 Hz 成分主要被提取到恢复信号中。

\begin{table}[H]
\centering
\caption{不同信号在 2 Hz 处的峰值功率对比}
\begin{tabular}{cc}
\toprule
频谱 & 2 Hz 处峰值功率\\
\midrule
原始信号直接 FFT & 130729.20\\
恢复信号 & 123771.90\\
残差信号 & 117.68\\
\bottomrule
\end{tabular}
\end{table}

由表中结果可得目标频率处的峰值抑制比为
\[
\frac{P_{\mathrm{raw}}(f_0)}{P_{\mathrm{res}}(f_0)}=\frac{130729.20}{117.68}\approx1110.86.
\]
该比值说明，原始信号中 2 Hz 处的功率在残差信号中被削弱超过 1100 倍；残差频谱主峰位于 \(35.83\ \mathrm{Hz}\)，目标频率已不再占主导。

\paragraph{（3）定量评价指标}

附件只给出了含噪观测信号，未给出无噪声真实信号，因此无法直接计算逐点真实误差。本文采用可观测残差作为辅助误差评价：
\[
r(t_i)=x(t_i)-\hat{x}(t_i).
\]
残差主要反映恢复周期分量后剩余的背景扰动，并不等同于真实恢复误差。

\begin{figure}[H]
\centering
\includegraphics[width=0.82\textwidth]{figures/q2_residual_signal.png}
\caption{恢复周期分量后的残差信号}
\end{figure}

残差信号仍呈随机波动，说明原始观测中除目标 2 Hz 周期分量外还含有大量背景扰动；结合频域结果，其主导频率已不再是 2 Hz，可作为目标分量被剥离的证据。

计算得到残差评价指标如下：
\begin{table}[H]
\centering
\caption{恢复结果的可观测残差评价指标}
\begin{tabular}{ccc}
\toprule
指标 & 公式 & 数值\\
\midrule
RMSE & $\sqrt{\frac{1}{N}\sum r_i^2}$ & 0.100013683\\
MAE & $\frac{1}{N}\sum |r_i|$ & 0.079703345\\
决定系数 $R^2$ & $1-\frac{\sum r_i^2}{\sum (x_i-\bar x)^2}$ & 0.058293952\\
估计 SNR & $10\log_{10}(P_s/P_r)$ & $-12.082919\ \mathrm{dB}$\\
\bottomrule
\end{tabular}
\end{table}

其中，\(R^2\) 较低并不表示目标分量恢复失败，因为本题恢复的是强噪声中的微弱周期成分，而不是完整拟合全部含噪观测。本文以频域一致性和目标频率抑制比作为主要依据，残差指标作为辅助评价。

\subsubsection{本节小结}

本问在第一问已知故障频率 \(f_0=2.00\ \mathrm{Hz}\) 的基础上，建立含直流偏置的正弦最小二乘恢复模型，求得信号幅值和初相。

最终恢复得到的单源故障周期信号为
\[
\hat s(t)=0.0351903\sin(4\pi t+1.57555).
\]
时域结果表明恢复信号具有稳定的 2 Hz 正弦结构；频域结果表明恢复信号能量集中于目标频率，且残差中 2 Hz 分量被削弱超过 1100 倍。因此，该方法能够有效恢复强噪声中的单源故障微弱周期信号。
"""

BACK = r"""
\subsection{问题四模型的建立与求解}
% 问题四模型建立与求解待补充。

\section{模型的分析与检验}
\subsection{灵敏度分析}
\subsection{误差分析}
\subsection{稳定性检验}

\section{模型的评价、改进与推广}
\subsection{模型的优点}
\subsection{模型的缺点}
\subsection{模型的改进}
\subsection{模型的推广}

\section{参考文献}

\appendix
\section{核心程序代码}

\subsection{问题一：MWS 随机共振核心代码}
\begin{lstlisting}[language=Matlab]
[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
y = xRaw(:) - mean(xRaw, 'omitnan');
u = y / std(y, 0, 'omitnan');
fs = 1 / median(diff(t));
f0 = 2.000001615205;

baseParams = struct();
baseParams.a = 1.0;
baseParams.downsampleFactor = 2;
baseParams.keepFraction = 0.10;
baseParams.initialState = 0.05;
baseParams.stateLimit = 15;
baseParams.smoothingWindow = max(5, round(fs/baseParams.downsampleFactor/f0/8));

bList = [0.8, 1.2, 1.8, 2.5, 3.5, 5.5];
v0List = [0.8, 1.5, 2.5, 4.0, 8.0, 16.0, 24.0];
rList = [0.25, 0.50, 0.80, 1.20, 1.60];
cList = [0.12, 0.20, 0.35, 0.60, 0.90];
kList = [0.04, 0.08, 0.14, 0.22, 0.32, 0.48, 0.70];

rowIdx = 0;
for ib = 1:numel(bList)
    for iv = 1:numel(v0List)
        for ir = 1:numel(rList)
            for ic = 1:numel(cList)
                for ik = 1:numel(kList)
                    p = baseParams;
                    p.b = bList(ib); p.V0 = v0List(iv);
                    p.R = rList(ir); p.C = cList(ic); p.k = kList(ik);

                    [tm, xm, stableFlag] = simulateMwsSr(t, u, p);
                    if ~stableFlag || numel(xm) < 200
                        continue;
                    end
                    projFrac = coherentProjectionFraction(tm, xm, f0);
                    [localSnrDb, peakAmp] = localFrequencySnr(tm, xm, f0);
                    roughFreq = localPeakFrequency(tm, xm, f0, 0.35);

                    rowIdx = rowIdx + 1;
                    score(rowIdx, :) = [projFrac, localSnrDb, peakAmp, roughFreq];
                    paramSet(rowIdx) = p; %#ok<SAGROW>
                    outputSet{rowIdx} = xm; timeSet{rowIdx} = tm; %#ok<SAGROW>
                end
            end
        end
    end
end

[~, bestIdx] = max(0.42*zscore(score(:,1)) + 0.43*zscore(score(:,2)) ...
    - 0.15*zscore(abs(score(:,4)-f0)));
bestParams = paramSet(bestIdx);
bestTime = timeSet{bestIdx};
bestOutput = outputSet{bestIdx};

function [tout, xout, stableFlag] = simulateMwsSr(t, xNorm, p)
    ds = max(1, round(p.downsampleFactor));
    td = t(1:ds:end);
    u = xNorm(1:ds:end);
    if p.smoothingWindow > 1
        u = movmean(u, p.smoothingWindow);
    end
    dt = median(diff(td));
    x = p.initialState;
    keepStart = max(1, floor(numel(td) * p.keepFraction));
    xout = zeros(numel(td)-keepStart, 1);
    tout = zeros(numel(td)-keepStart, 1);
    outIdx = 0;
    stableFlag = true;
    for n = 1:numel(td)-1
        inputValue = p.k * u(n);
        k1 = mwsDeriv(x, inputValue, p);
        k2 = mwsDeriv(x + 0.5*dt*k1, inputValue, p);
        k3 = mwsDeriv(x + 0.5*dt*k2, inputValue, p);
        k4 = mwsDeriv(x + dt*k3, inputValue, p);
        x = x + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
        if ~isfinite(x)
            stableFlag = false;
            break;
        end
        x = min(max(x, -p.stateLimit), p.stateLimit);
        if n >= keepStart
            outIdx = outIdx + 1;
            tout(outIdx) = td(n);
            xout(outIdx) = x;
        end
    end
    tout = tout(1:outIdx);
    xout = xout(1:outIdx) - mean(xout(1:outIdx), 'omitnan');
end

function dx = mwsDeriv(x, inputValue, p)
    expArg = min(max((abs(x) - p.R) / p.C, -60), 60);
    q = exp(expArg);
    wsForce = (p.V0 / p.C) * sign(x) * q / (1 + q)^2;
    dx = -p.b * x - wsForce + inputValue;
end
\end{lstlisting}

\subsection{问题一：一阶自相关增强核心代码}
\begin{lstlisting}[language=Matlab]
[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
y = xRaw(:) - mean(xRaw, 'omitnan');
u = y / std(y, 0, 'omitnan');
fs = 1 / median(diff(t));

maxLagSeconds = 200;
[acfLag, acfSignal, acfRaw] = firstOrderAutocorrelation(u, fs, maxLagSeconds);

function [lag, acfSignal, acfRaw] = firstOrderAutocorrelation(u, fs, maxLagSeconds)
    u = u(:) - mean(u(:), 'omitnan');
    n = numel(u);
    nfft = 2 ^ nextpow2(2 * n - 1);
    U = fft(u, nfft);
    r = real(ifft(U .* conj(U)));
    r = r(1:n) / n;
    r = r(2:end);

    keepCount = min(numel(r), max(10, floor(maxLagSeconds * fs)));
    acfRaw = r(1:keepCount);
    lag = (1:keepCount)' / fs;
    acfSignal = detrend(acfRaw, 'linear');

    acfSignal = acfSignal - mean(acfSignal, 'omitnan');
    scale = std(acfSignal, 0, 'omitnan');
    if ~isfinite(scale) || scale <= eps
        scale = 1;
    end
    acfSignal = acfSignal / scale;
end
\end{lstlisting}

\subsection{问题一：全长正弦投影测频核心代码}
\begin{lstlisting}[language=Matlab]
% read data and remove DC component
[t, xRaw] = readSingleSourceData(dataFile);
t = t(:);
y = xRaw(:) - mean(xRaw, 'omitnan');
fs = 1 / median(diff(t));
totalEnergy = sum(y.^2);

% FFT coarse localization
[freq, ~, powerRaw] = fullLengthSpectrum(y, fs);
valid = freq > 0.05 & freq < min(49.5, fs/2 - 0.1);
validIdx = find(valid);
[~, localMax] = max(powerRaw(validIdx));
fFft = freq(validIdx(localMax));

% full-length sinusoidal projection search
halfWidth = 0.04;
fGrid = linspace(fFft-halfWidth, fFft+halfWidth, 6001)';
projEnergy = zeros(size(fGrid));
for k = 1:numel(fGrid)
    [~, fit, ~] = fitSinusoidAtFrequency(t, y, fGrid(k));
    projEnergy(k) = sum(fit.^2);
end
[~, bestGridIdx] = max(projEnergy);
fCoarse = fGrid(bestGridIdx);
gridStep = fGrid(2) - fGrid(1);
objective = @(f) -projectionEnergyAtFrequency(t, y, f);
fHat = fminbnd(objective, fCoarse-5*gridStep, fCoarse+5*gridStep);

% estimate component and residual at the final frequency
[theta, fit, residual] = fitSinusoidAtFrequency(t, y, fHat);
J = sum(fit.^2);
energyFraction = J / totalEnergy;
\end{lstlisting}

\subsection{问题二：单源波形恢复核心代码}
\begin{lstlisting}[language=Matlab]
% f0 is obtained from Question 1
f0 = 2.000001615205;
[t, x] = readSingleSourceData(dataFile);
t = t(:);
x = x(:);

% linearized sinusoidal model with offset
H = [sin(2*pi*f0*t), cos(2*pi*f0*t), ones(size(t))];
theta = H \ x;
a = theta(1);
b = theta(2);
c = theta(3);

% recover amplitude, phase and signal
Ahat = sqrt(a^2 + b^2);
phihat = atan2(b, a);
xFit = H * theta;
sHat = xFit - c;
residual = x - xFit;

% observable evaluation
rmse = sqrt(mean(residual.^2));
mae = mean(abs(residual));
r2 = 1 - sum(residual.^2) / sum((x - mean(x)).^2);
fs = 1 / median(diff(t));
[freqRaw, pRaw] = singleSidedPower(x, fs);
[freqRes, pRes] = singleSidedPower(residual, fs);
[~, k0] = min(abs(freqRaw - f0));
[~, kr] = min(abs(freqRes - f0));
suppressionRatio = pRaw(k0) / pRes(kr);
\end{lstlisting}

\subsection{问题三：投影 OMP 多源分离核心代码}
\begin{lstlisting}[language=Matlab]
[t, xRaw] = readMultiSourceData(dataFile);
t = t(:);
y = xRaw(:) - mean(xRaw, 'omitnan');
fs = 1 / median(diff(t));
settings = defaultOmpSettings(fs, t(end)-t(1));

selectedFreqs = [];
residual = y;
for iter = 1:settings.maxSources
    [gridFreq, gridEnergy] = gridProjectionSpectrum(residual, fs);
    valid = gridFreq >= settings.fMinHz & gridFreq <= settings.fMaxHz;
    for j = 1:numel(selectedFreqs)
        valid = valid & abs(gridFreq-selectedFreqs(j)) > settings.minSeparationHz;
    end

    [noiseMedian, noiseMadSigma, threshold] = ...
        robustNoiseThreshold(gridEnergy(valid), settings.thresholdMultiplier);
    [peakEnergy, peakIdx] = max(gridEnergy(valid));
    validFreq = gridFreq(valid);
    coarseFreq = validFreq(peakIdx);
    peakToMedian = peakEnergy / max(noiseMedian, eps);

    if peakEnergy < threshold || peakToMedian < settings.minPeakToMedianRatio
        break;
    end

    refinedFreq = refineFrequency(t, residual, coarseFreq, ...
        settings.refineHalfWidthHz, settings);
    selectedFreqs = sort([selectedFreqs; refinedFreq]);
    selectedFreqs = coordinateRefineFrequencies(t, y, selectedFreqs, settings);
    [theta, fitAll, residual, componentMatrix] = ...
        fitMultiSinusoid(t, y, selectedFreqs);
end

% convert linear coefficients to amplitude and phase
for k = 1:numel(selectedFreqs)
    a = theta(2*k-1);
    b = theta(2*k);
    Ahat(k) = sqrt(a^2 + b^2); %#ok<AGROW>
    phihat(k) = atan2(a, b);   %#ok<AGROW>
end
\end{lstlisting}

\subsection{问题四：三传感器布局优化核心代码}
\begin{lstlisting}[language=Matlab]
sensors = make_candidate_sensors('M12_baseline');
sources = make_fixed_sources(3, f0);
layouts = enumerate_layouts(numel(sensors), 3);

bestObj = -inf;
bestLayout = layouts{1};
for k = 1:numel(layouts)
    S = layouts{k};
    lambdas = zeros(numel(sources), 1);
    for j = 1:numel(sources)
        lambdas(j) = layout_lambda(sensors, sources(j), S, ...
            baseAlpha, baseBeta, baseAmp, sigma0, L);
    end
    objective = baseOmega * min(lambdas) + ...
        (1-baseOmega) * mean(lambdas);
    if objective > bestObj
        bestObj = objective;
        bestLayout = S;
    end
end

eta = chi2_threshold_even(2*numel(bestLayout), pfaMax);
simBest = simulate_layout_pd(sensors, sources, bestLayout, ...
    baseAlpha, baseBeta, baseAmp, sigma0, fs, t, eta, nMCDetect);

function lam = layout_lambda(sensors, source, layout, alpha, beta, amp, sigma0, L)
    lam = 0;
    for idx = layout
        d = norm(sensors(idx).point - source.point);
        g = exp(-alpha*d) / (1 + beta*d);
        sigmai = sigma0 * sqrt(sensors(idx).nu);
        ai = sensors(idx).kappa * g * amp * source.ampScale;
        lam = lam + L * ai^2 / (2 * sigmai^2);
    end
end
\end{lstlisting}

\end{document}
"""


def main():
    q1 = read_generated("5_1.tex")
    q3 = read_generated("5_3.tex")
    full = "\n".join([PREAMBLE, FRONT_AND_MID, q1, Q2, q3, BACK])
    (ROOT / "Final.tex").write_text(full, encoding="utf-8")
    print(f"wrote {ROOT / 'Final.tex'}: {len(full.splitlines())} lines, {len(full)} chars")


if __name__ == "__main__":
    main()
