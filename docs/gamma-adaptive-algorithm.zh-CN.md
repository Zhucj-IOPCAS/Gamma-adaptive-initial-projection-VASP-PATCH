# Gamma-Adaptive Initial-Projection Algorithm

## Overview

Gamma-adaptive 方法使用 Gamma 点 Bloch 波函数自动构造 Wannier 初始投影

## Contents

- [1. Gamma 点波函数的 Wannier 初猜](#1-gamma-点波函数的-wannier-初猜)
  - [1.1 Disentanglement](#11-disentanglement选择最优-gamma-子空间)
  - [1.2 Jacobi rotation](#12-jacobi-rotationgamma-初猜的-foster-boys-局域化)
  - [1.3 Locate center](#13-locate-centergamma-初猜的轨道中心)
  - [1.4 Phase fix](#14-phase-fixgamma-初猜的整体相位规范)
- [2. Gamma-adaptive .amn 矩阵](#2-gamma-adaptive-amn-矩阵)

## 1. Gamma 点波函数的 Wannier 初猜

### 1.1 Disentanglement：选择最优 Gamma 子空间

在 Gamma 点的 active Bloch 本征态中选择一个 $J$ 维子空间，其中 $J$ 是 Wannier 轨道数。

$$\ket{\phi_{n\Gamma}^{(S)}}=\sum_{m=1}^{N}\ket{\psi_{m\Gamma}} S_{mn}, n=1,\ldots,J.$$

其中，

$$S\in\mathbb C^{N\times J}, S^\dagger S=I_J.$$

矩阵 $S$ 的列给出所选子空间的一组正交基，相应的投影算符为

$$P=SS^\dagger.$$

Disentanglement 的目标是确定投影算符 $P$。

#### 1.1.1 Gamma 点波函数的 overlap

Gamma 点虽然只有一个独立的晶体动量，但周期位置算符的信息仍可由一组连接 Gamma 点及其倒格矢像的 overlap 表示。对每个邻接倒空间向量 $\mathbf b$，定义

$$M^{(b)}_{mn}=\bra{\psi_{m\Gamma}}e^{-i\mathbf b\cdot\mathbf r}\ket{\psi_{n\Gamma}}.$$

每个邻接向量对应一个几何权重 $w_b$。在完整 active Bloch 空间中，

$$M^{(b)}\in\mathbb C^{N\times N}.$$

若选择矩阵为 $S$，则该 $J$ 维子空间中的 overlap 为

$$M_S^{(b)}=S^\dagger M^{(b)}S\in\mathbb C^{J\times J}.$$

#### 1.1.2 Outer window 与 frozen window

Outer window 定义允许参与子空间选择的候选 Bloch 态集合

$$\mathcal O=\lbrace m\mid E_{\mathrm{win}}^{\min}\leq\epsilon_m\leq E_{\mathrm{win}}^{\max}\rbrace.$$

目标子空间只能从 $\mathcal O$ 中选取，因此对所有 $m\notin\mathcal O$，必须有

$$S_{mn}=0.$$

Frozen window 定义必须完整保留在目标子空间中的 Bloch 本征态集合

$$\mathcal F=\lbrace m\in\mathcal O\mid E_{\mathrm{froz}}^{\min}\leq\epsilon_m\leq E_{\mathrm{froz}}^{\max}\rbrace.$$

记：

- $N$：active Bloch 态数；
- $B=|\mathcal O|$：outer-window Bloch 态数；
- $J$：目标 Wannier 轨道数；
- $M=|\mathcal F|$：frozen Bloch 态数；
- $\mathcal G=\mathcal O\setminus\mathcal F$：non-frozen outer states 组成的子空间；
- $L=J-M$：需要从 $\mathcal G$ 中选出的自由方向数。

如果没有 frozen window，则 $\mathcal F=\varnothing$ 且 $M=0$。各维数满足

$$M\leq J\leq B\leq N.$$

Frozen 的含义是子空间包含关系。对每个 $f\in\mathcal F$，令 $\ket{e_f}$ 表示该 Bloch 本征态在 active-band 基底中的单位坐标向量，则要求

$$P\ket{e_f}=\ket{e_f}.$$

设 frozen band indices 为 $f_1,\ldots,f_M$。把相应的单位坐标列组成

$$F=(\ket{e_{f_1}},\ldots,\ket{e_{f_M}})\in\mathbb C^{N\times M}, F^\dagger F=I_M.$$

Frozen 子空间的投影算符为

$$P_F=FF^\dagger.$$

目标空间还需要从 $\mathcal G$ 中选择

$$L=J-M$$

个正交方向。将这些方向组成

$$C\in\mathbb C^{N\times L}, C^\dagger C=I_L, F^\dagger C=0.$$

并要求 $C$ 在 $\mathcal G$ 以外的行全部为零。完整的子空间矩阵写为

$$S=(F,C).$$

相应投影算符为

$$P=SS^\dagger=FF^\dagger+CC^\dagger.$$


初始子空间包含全部 frozen states，并从 $\mathcal G$ 中选取能量最低的 $L$ 个 Bloch states 作为初始自由方向。

#### 1.1.3 Gauge-invariant spread $\Omega_I$

对任意候选子空间 $S$，定义

$$\Omega_I[S]=\sum_b w_b\left[J-\sum_{m,n=1}^{J}|(M_S^{(b)})_{mn}|^2\right].$$

利用 $P=SS^\dagger$，也可以写成

$$\Omega_I[P]=\sum_b w_b\left\lbrace J-\mathrm{Tr}\left[PM^{(b)}PM^{(b)\dagger}\right]\right\rbrace.$$

$\Omega_I$ 度量目标子空间在周期 overlap 映射下向外部空间的泄漏。Disentanglement 通过最小化 $\Omega_I$ 选择目标子空间。

#### 1.1.4 SMV 自洽迭代

定义

$$Z[P]=\sum_b w_bM^{(b)}PM^{(b)\dagger}.$$

则 gauge-invariant spread 可以写成

$$\Omega_I[P]=J\sum_b w_b-\mathrm{Tr}\left[PZ[P]\right].$$

因此，最小化 $\Omega_I[P]$ 等价于最大化

$$\Phi[P]=\mathrm{Tr}\left[PZ[P]\right].$$

由于 $Z[P]$ 依赖于待求投影算符 $P$，该问题通过自洽迭代求解。给定上一轮的投影算符 $P^{(i-1)}$，构造

$$Z^{(i)}=\sum_b w_bM^{(b)}P^{(i-1)}M^{(b)\dagger}.$$

在本轮迭代中，$Z^{(i)}$ 被视为固定矩阵。以下将其简记为 $Z$，并通过最大化

$$\mathrm{Tr}\left[PZ\right]$$

确定新的目标子空间。Frozen 子空间的投影算符可以写成

$$P=P_F+CC^\dagger.$$

代入目标函数得到

$$\mathrm{Tr}\left[PZ\right]=\mathrm{Tr}\left[P_FZ\right]+\mathrm{Tr}\left[C^\dagger ZC\right].$$

第一项由 frozen 子空间决定，因此只需在 non-frozen outer space $\mathcal G=\mathcal O\setminus\mathcal F$ 中最大化第二项。令 $R_{\mathcal G}$ 的列为 $\mathcal G$ 中的单位坐标向量，并写成

$$C=R_{\mathcal G}V, V^\dagger V=I_L.$$

将 $Z$ 限制到 non-frozen outer space，得到

$$Z_{\mathcal G}=R_{\mathcal G}^\dagger ZR_{\mathcal G}.$$

目标函数相应化为

$$\mathrm{Tr}\left[C^\dagger ZC\right]=\mathrm{Tr}\left[V^\dagger Z_{\mathcal G}V\right].$$

根据 Ky Fan 最大值原理，该目标函数的最大值由 $Z_{\mathcal G}$ 最大的 $L$ 个本征值决定。设

$$Z_{\mathcal G}v_a=\lambda_av_a, \lambda_1\geq\lambda_2\geq\cdots,$$

则

$$\max_{V^\dagger V=I_L}\mathrm{Tr}\left[V^\dagger Z_{\mathcal G}V\right]=\sum_{a=1}^{L}\lambda_a.$$

因此

$$V=(v_1,\ldots,v_L).$$

由此迭代自由子空间、完整子空间及其投影算符：

$$C=R_{\mathcal G}V, S=(F,C), P=SS^\dagger.$$

再用迭代后的 $P$ 构造下一轮 $Z$，重复上述过程。每轮迭代后监控 gauge-invariant spread 的变化

$$\Delta\Omega_I^{(i)}=|\Omega_I^{(i)}-\Omega_I^{(i-1)}|,$$

以及投影算符的变化

$$\Delta P^{(i)}=\|P^{(i)}-P^{(i-1)}\|_F.$$

当 $\Delta\Omega_I^{(i)}$ 和 $\Delta P^{(i)}$ 同时稳定时，得到最终子空间，并将 Gamma 周期 overlap 投影到选定的 $J$ 维子空间：

$$M_{\mathrm{S}}^{(b)}=S^\dagger M^{(b)}S.$$

### 1.2 Jacobi rotation：Gamma 初猜的 Foster-Boys 局域化

将 disentanglement 得到的 $N_w=J$ 个正交态记为 $\ket{\psi_{m\Gamma}}$。Jacobi rotation 在这个固定子空间内构造 Gamma 点局域轨道

$$\ket{\phi_{\alpha}^{\Gamma}}=\sum_{m=1}^{N_w}U_{m\alpha}^{(\Gamma)}\ket{\psi_{m\Gamma}}.$$

#### 1.2.1 Gamma 点 Foster-Boys 泛函

Jacobi rotation 前的 overlap 为

$$M_{mn}^{(0,\mathbf b)}=\bra{\psi_{m\Gamma}}e^{-i\mathbf b\cdot\mathbf r}\ket{\psi_{n\Gamma}}.$$

经过 $U^{(\Gamma)}$ 旋转后，

$$M^{(\mathbf b)}=U^{(\Gamma)\dagger}M^{(0,\mathbf b)}U^{(\Gamma)}.$$

Gamma 点周期 Foster-Boys spread 为

$$\Omega_{\Gamma}=\sum_{\mathbf b}w_{\mathbf b}\sum_{n=1}^{N_w}\left(1-|M_{nn}^{(\mathbf b)}|^2\right).$$

因此，最小化 $\Omega_{\Gamma}$ 等价于最大化

$$F=\sum_{\mathbf b}w_{\mathbf b}\sum_{n=1}^{N_w}|M_{nn}^{(\mathbf b)}|^2.$$

#### 1.2.2 二轨道 $SU(2)$ 旋转

Jacobi 方法依次选择一对轨道，并在对应的二维子空间内施加

$$M'^{(\mathbf b)}=V^\dagger M^{(\mathbf b)}V, V\in SU(2).$$

将这一对轨道对应的 $2\times2$ overlap block 写成

$$M^{(\mathbf b)}=a_0^{(\mathbf b)}I+\mathbf a^{(\mathbf b)}\cdot\boldsymbol{\sigma},$$

其中

$$a_0^{(\mathbf b)}=\frac{1}{2}\left(M_{11}^{(\mathbf b)}+M_{22}^{(\mathbf b)}\right),$$

$$a_1^{(\mathbf b)}=\frac{1}{2}\left(M_{12}^{(\mathbf b)}+M_{21}^{(\mathbf b)}\right),$$

$$a_2^{(\mathbf b)}=\frac{i}{2}\left(M_{12}^{(\mathbf b)}-M_{21}^{(\mathbf b)}\right),$$

$$a_3^{(\mathbf b)}=\frac{1}{2}\left(M_{11}^{(\mathbf b)}-M_{22}^{(\mathbf b)}\right).$$

$SU(2)$ 旋转可以写成

$$V=e^{i\frac{\theta}{2}\mathbf n\cdot\boldsymbol{\sigma}}=\cos\frac{\theta}{2}I+i\sin\frac{\theta}{2}\mathbf n\cdot\boldsymbol{\sigma}.$$

它在 Pauli 向量空间中对应一个 $SO(3)$ 旋转 $R$：

$$V^\dagger\left(\mathbf a\cdot\boldsymbol{\sigma}\right)V=(R\mathbf a)\cdot\boldsymbol{\sigma}.$$

相应的三维旋转为

$$R_{\mathbf n}(\theta)\mathbf a=\mathbf a\cos\theta+(\mathbf n\times\mathbf a)\sin\theta+\mathbf n(\mathbf n\cdot\mathbf a)(1-\cos\theta).$$

#### 1.2.3 确定最优旋转方向

旋转后的两个对角元满足

$$|M_{11}'^{(\mathbf b)}|^2+|M_{22}'^{(\mathbf b)}|^2=2|a_0^{(\mathbf b)}|^2+2|a_3'^{(\mathbf b)}|^2.$$

由于 $a_0^{(\mathbf b)}$ 在相似变换下不变，轨道对的最优旋转等价于最大化

$$\sum_{\mathbf b}w_{\mathbf b}|a_3'^{(\mathbf b)}|^2.$$

定义

$$\mathbf u=R^T\mathbf e_3,$$

则

$$a_3'^{(\mathbf b)}=\mathbf u^T\mathbf a^{(\mathbf b)}.$$

因此

$$\sum_{\mathbf b}w_{\mathbf b}|a_3'^{(\mathbf b)}|^2=\mathbf u^TW\mathbf u,$$

其中 $W$ 是 $3\times3$ 实对称矩阵，

$$W_{\alpha\beta}=\sum_{\mathbf b}w_{\mathbf b}\mathrm{Re}\left[a_{\alpha}^{(\mathbf b)}a_{\beta}^{(\mathbf b)\ast}\right].$$

根据 Rayleigh-Ritz 变分原理，

$$\max_{|\mathbf u|=1}\mathbf u^TW\mathbf u=\lambda_{\max},$$

最优方向 $\mathbf u_{\max}$ 是 $W$ 最大本征值对应的归一化本征向量：

$$W\mathbf u_{\max}=\lambda_{\max}\mathbf u_{\max}.$$

#### 1.2.4 累计幺正变换

由 $\mathbf u_{\max}=(u_x,u_y,u_z)$ 定义

$$c=\sqrt{\frac{1+u_z}{2}},z=\frac{u_x+iu_y}{2c},$$

则 $V$ 为：

$$V=\begin{pmatrix}c & -z^\ast \\ 
z & c\end{pmatrix}.$$

对所有 $\mathbf b$ 迭代该轨道对的 overlap：

$$M^{(\mathbf b)}\rightarrow V^\dagger M^{(\mathbf b)}V.$$

同时累计 Gamma 点 Jacobi 变换：

$$U^{(\Gamma)}\rightarrow U^{(\Gamma)}V.$$

依次旋转所有轨道对并重复上述过程，直至 $\Omega_\Gamma$ 稳定。最终与 disentanglement 变换 $S$ 合并：

$$U_{\mathrm{total}}^{(\Gamma)}=SU^{(\Gamma)}.$$

### 1.3 Locate center：Gamma 初猜的轨道中心

对每个 Gamma 点局域轨道，其中心定义为能够用最小 Fermi 包络包含给定比例电荷的位置。

#### 1.3.1 实空间电荷密度

轨道 $\alpha$ 的电荷密度为

$$\rho_\alpha(\mathbf r)=\sum_s|\phi_{\alpha s}^{\Gamma}(\mathbf r)|^2,$$

其中标量波函数只有一个分量，spinor 波函数则对两个自旋分量求和。

#### 1.3.2 Fermi 包络与最小包络半径

对于每个候选包络中心 $\mathbf r_c$，构造以 $\mathbf r_c$ 为中心的 Wigner-Seitz 原胞。以下 $\mathbf r$ 均表示该 WS 原胞内相对于包络中心的坐标，因此 $|\mathbf r|$ 是到包络中心的距离。

Fermi 包络定义为

$$f_R(\mathbf r)=\left[1+\exp\left(\frac{|\mathbf r|-R}{\sigma}\right)\right]^{-1}.$$

包络内的电荷比例为

$$q_\alpha(\mathbf r_c,R)=\frac{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_R(\mathbf r)d\mathbf r}{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)d\mathbf r}.$$

给定目标电荷比例 $\eta$，对每个候选中心确定满足该条件的最小包络半径：

$$R_\alpha(\mathbf r_c)=\min\lbrace R\mid q_\alpha(\mathbf r_c,R)\geq\eta\rbrace.$$

中心搜索从多个候选位置分别开始，包括轨道电荷密度最大的位置和均匀分布在周期原胞中的位置。首先将候选位置按照包含目标电荷比例所需的最小包络半径排序。随后从包络半径较小的候选位置开始进行完整的自洽中心迭代。

#### 1.3.3 包络中心的自洽迭代

包络内的电荷中心相对于包络中心的位移为

$$\mathbf s_\alpha(\mathbf r_c)=\frac{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_{R_\alpha}(\mathbf r)\mathbf r d\mathbf r}{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_{R_\alpha}(\mathbf r)d\mathbf r}.$$

将包络中心向包络内的电荷中心移动：

$$\mathbf r_{\mathrm{trial}}=\mathbf r_c+\gamma\mathbf s_\alpha(\mathbf r_c).$$

在试探中心重新构造 WS 原胞，并重新计算最小包络半径和包络中心位移。包络中心位移减小或包络半径减小时接受试探中心，重复上述过程直至

$$\left|\mathbf s_\alpha(\mathbf r_c)\right|<\epsilon_c.$$

对每个收敛中心，计算 Fermi 包络内的电荷二阶矩：

$$\mu_{2,\alpha}(\mathbf r_c)=\frac{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_{R_\alpha}(\mathbf r)|\mathbf r|^2d\mathbf r}{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_{R_\alpha}(\mathbf r)d\mathbf r}.$$

最终选择包络半径最小的自洽中心；当包络半径接近时，选择二阶矩更小的中心。

### 1.4 Phase fix：Gamma 初猜的整体相位规范

#### 1.4.1 标量 wannier 轨道的相位规范

对于标量轨道，取波函数振幅最大的位置

$$\mathbf r_\alpha^{\max}=\mathrm\ast{argmax}_{\mathbf r}|\phi_\alpha^\Gamma(\mathbf r)|,$$

并使用相位修正

$$p_\alpha=\frac{\phi_\alpha^\Gamma(\mathbf r_\alpha^{\max})^\ast}{|\phi_\alpha^\Gamma(\mathbf r_\alpha^{\max})|}.$$

相应地迭代 Gamma 点变换矩阵：

$$\widetilde U_{\mathrm{total},m\alpha}^{(\Gamma)}=p_\alpha U_{\mathrm{total},m\alpha}^{(\Gamma)}.$$

#### 1.4.2 Spinor wannier 轨道的相位规范

对于 spinor 轨道，写成

$$\Phi_{\alpha,z}^\Gamma(\mathbf r)=\begin{pmatrix}a_{\alpha,z}(\mathbf r)\\
b_{\alpha,z}(\mathbf r)\end{pmatrix}.$$

所有轨道使用同一个公共自旋基底 $V$：

$$\Phi_{\alpha,V}^\Gamma(\mathbf r)=V^\dagger\Phi_{\alpha,z}^\Gamma(\mathbf r)=\begin{pmatrix}a_{\alpha,V}(\mathbf r)\\
b_{\alpha,V}(\mathbf r)\end{pmatrix}.$$

在该基底中，选取对称参考 spinor

$$\chi_V(\theta,\phi)=e^{-i\pi/4}\begin{pmatrix}e^{-i\phi/2}\cos(\theta/2)\\
e^{i\phi/2}\sin(\theta/2)\end{pmatrix}.$$

spinor function表示为

$$\Phi_{\alpha,V}^\Gamma(\mathbf r)=\sqrt{\rho_\alpha(\mathbf r)}e^{i[\gamma_{\alpha,V}(\mathbf r)+\pi/4]}\chi_V(\theta_{\alpha,V}(\mathbf r),\phi_{\alpha,V}(\mathbf r)).$$

其中

$$\rho_\alpha(\mathbf r)=|a_{\alpha,V}(\mathbf r)|^2+|b_{\alpha,V}(\mathbf r)|^2.$$

当两个分量均非零时，局域相位为

$$\gamma_{\alpha,V}(\mathbf r)=\frac{\arg a_{\alpha,V}(\mathbf r)+\arg b_{\alpha,V}(\mathbf r)}{2}\pmod\pi.$$

当 $a_{\alpha,V}$ 或 $b_{\alpha,V}$ 接近零时，相应的 $\arg a_{\alpha,V}$ 或 $\arg b_{\alpha,V}$ 失去意义。因此需要选择公共自旋主轴，使主要电荷分布尽量远离新自旋基底的两个极点。

定义局域自旋密度

$$\mathbf s_\alpha(\mathbf r)=\Phi_{\alpha,z}^{\Gamma\dagger}(\mathbf r)\boldsymbol\sigma\Phi_{\alpha,z}^\Gamma(\mathbf r).$$

若新基底的 $z$ 轴方向为单位向量 $\mathbf n$，则

$$\cos\theta_{\alpha,V}(\mathbf r)=\frac{\mathbf n\cdot\mathbf s_\alpha(\mathbf r)}{\rho_\alpha(\mathbf r)}.$$

两个极点对应 $|\cos\theta_{\alpha,V}|=1$。因此通过最小化全部轨道的平均极向权重选择 $\mathbf n$：

$$\mathcal P(\mathbf n)=\frac{1}{N_w}\sum_{\alpha=1}^{N_w}\frac{\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)\cos^2\theta_{\alpha,V}(\mathbf r)d\mathbf r}{\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)d\mathbf r}.$$

该泛函可以写成

$$\mathcal P(\mathbf n)=\mathbf n^TK\mathbf n,$$

其中

$$K_{ij}=\frac{1}{N_w}\sum_{\alpha=1}^{N_w}\frac{\int_{\Omega_\alpha}s_{\alpha i}(\mathbf r)s_{\alpha j}(\mathbf r)/\rho_\alpha(\mathbf r)d\mathbf r}{\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)d\mathbf r}.$$

根据 Rayleigh-Ritz 变分原理，最优主轴是 $K$ 最小本征值对应的归一化本征向量：

$$K\mathbf n=\kappa_{\min}\mathbf n.$$

由 $\mathbf n=(n_x,n_y,n_z)$ 构造公共自旋基底

$$c=\sqrt{\frac{1+n_z}{2}},$$

$$z=\frac{n_x+in_y}{2c},$$

$$V=\begin{pmatrix}c&-z^\ast\\
z&c\end{pmatrix}.$$

在确定公共主轴后，使用 doubled phase

$$D_\alpha(\mathbf r)=e^{2i[\gamma_{\alpha,V}(\mathbf r)+\pi/4]}.$$

在北极区域，使用

$$D_\alpha^N(\mathbf r)=i\left[\frac{a_{\alpha,V}(\mathbf r)}{|a_{\alpha,V}(\mathbf r)|}\right]^2,$$

在南极区域，使用其时间反演对应形式

$$D_\alpha^S(\mathbf r)=-i\left[\frac{b_{\alpha,V}(\mathbf r)}{|b_{\alpha,V}(\mathbf r)|}\right]^2.$$

接下来为整条轨道选择一个共同相位 $e^{-i\vartheta_\alpha}$，使修正后轨道相对于参考 spinor 的实部平方在 WS 原胞内最大：

$$F_\alpha(\vartheta_\alpha)=\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)\left[\mathrm{Re}\left(e^{-i\vartheta_\alpha}e^{i[\gamma_{\alpha,V}(\mathbf r)+\pi/4]}\right)\right]^2d\mathbf r.$$

定义 doubled-phase 积分

$$C_\alpha=\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)D_\alpha(\mathbf r)d\mathbf r.$$

则

$$F_\alpha(\vartheta_\alpha)=\frac{1}{2}\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)d\mathbf r+\frac{1}{2}\mathrm{Re}\left[e^{-2i\vartheta_\alpha}C_\alpha\right].$$

因此 $F_\alpha$ 在

$$\vartheta_\alpha=\frac{1}{2}\arg C_\alpha\pmod\pi$$

时达到最大。最终相位修正为

$$p_\alpha=e^{-i\vartheta_\alpha}=e^{-i\arg(C_\alpha)/2},$$

并迭代

$$\widetilde U_{\mathrm{total},m\alpha}^{(\Gamma)}=p_\alpha U_{\mathrm{total},m\alpha}^{(\Gamma)}.$$

对称参考 spinor 中的 $e^{-i\pi/4}$ 保证

$$\Theta\chi_V(\theta,\phi)=\chi_V(\pi-\theta,\phi+\pi).$$

因此对于 Kramers pair，即使不进行配对，相位修正后自动满足

$$\widetilde\Phi_{\bar\alpha}^\Gamma=\pm\Theta\widetilde\Phi_\alpha^\Gamma.$$

## 2. Gamma-adaptive .amn 矩阵

### 2.1 Gamma 点与全 $\mathbf k$ 点波函数的 overlap

对于 Wannier 轨道 $\alpha$，以其中心 $\boldsymbol\tau_\alpha$ 构造 Wigner--Seitz 原胞 $\Omega_{\mathrm{WS}}(\boldsymbol\tau_\alpha)$。定义全 $\mathbf k$ Bloch 态与 Gamma 点 Bloch 态在该 WS 原胞内的 overlap：

$$I_{nm}^{(\mathbf k;\alpha)}=\int_{\Omega_{\mathrm{WS}}(\boldsymbol\tau_\alpha)}d^3r\,\psi_{n\mathbf k}^\ast(\mathbf r)\psi_{m\Gamma}(\mathbf r).$$

对于 spinor 波函数，内积包含两个 spinor 分量之和：

$$I_{nm}^{(\mathbf k;\alpha)}=\sum_s\int_{\Omega_{\mathrm{WS}}(\boldsymbol\tau_\alpha)}d^3r\,\psi_{n\mathbf k s}^\ast(\mathbf r)\psi_{m\Gamma s}(\mathbf r).$$

由于$\mathbf k$ 点波函数与 Gamma 点波函数的 Bloch phase 在原胞内不能相互消去，因此该 integrand
$$\psi_{n\mathbf k}^\ast(\mathbf r)\psi_{m\Gamma}(\mathbf r)=e^{-i\mathbf k\cdot\mathbf r}u_{n\mathbf k}^\ast(\mathbf r)u_{m\Gamma}(\mathbf r).$$

 不具有晶格周期性。该 overlap 不能直接由周期平面波系数的内积得到，而需要在选定的实空间 WS 原胞内显式积分。

### 2.2 构造 .amn 矩阵

经过 disentanglement、Jacobi rotation 和 phase fix 后，Gamma-adaptive Wannier 初猜为

$$\ket{g_{\alpha\mathbf 0}}=\sum_{m=1}^{N}\widetilde U_{\mathrm{total},m\alpha}^{(\Gamma)}\ket{\psi_{m\Gamma}}.$$

.amn 矩阵定义为全 $\mathbf k$ Bloch 态与该 Gamma-adaptive 初猜的 overlap：

$$A_{n\alpha}^{(\mathbf k)}=\braket{\psi_{n\mathbf k}|g_{\alpha\mathbf 0}}.$$

代入 Gamma 点初猜，得到

$$A_{n\alpha}^{(\mathbf k)}=\sum_{m=1}^{N}I_{nm}^{(\mathbf k;\alpha)}\widetilde U_{\mathrm{total},m\alpha}^{(\Gamma)}.$$
