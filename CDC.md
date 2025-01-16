# 股票和劳动力市场协同整合时集体确定缴款计划的最优策略

[toc]



## Introduction

在集体确定缴费（Collective Defined Contribution）养老金计划中，基金的资产集中由金融机构管理，福利取决于基金的财务状况。集体养老金计划通过强制参与能够实现代际之间的风险转移，相应的风险由当前和未来的几代人共同承担从而改善福利。

本文研究了一个CDC养老金计划中的随机最优问题，在此模型中假定缴费率是固定的，福利支付取决于最终的薪资水平。养老金基金可以投资于一种风险资产和一种无风险资产组成的金融市场。我们将劳动收入过程描述为几何随机游走，其中的漂移想依赖于当前股息与当前劳动收入的比率，其中股息过程遵循几何布朗运动。劳动收入被假定为总劳动收入与成员特质性冲击（shocks）的乘积，且其恒定增长率是未知的。我们通过一个连续时间的二状态隐马尔科夫链进行建模。此外，我们假设养老金成员分为在职和退休。在职成员依然在职并向养老金基金缴费而退休成员从养老金基金中获得福利。每个参与者在年龄$a_0$时加入养老金计划并在年龄$a_1$​时退休。在处理具有常数相对风险厌恶(CRRA)效用函数的最优投资问题时，传统的猜解不易用于求解相关的HJB方程，因此我们在模型不确定性的框架下进行研究，假设基金经历的目标是寻找最优的比例投资和福利替代政策，以最大化社会福利和终端盈余财富。通过求解HJB方程得到了最优资产配置和福利支付政策的显式解。

与现有的CDC基金方案最优设计相比，我们的工作考虑了劳动和股市之间的协整关系。在现实中，劳动收入与金融市场密切相关。劳动收入的随机性和金融市场的回报是基金经理面临的不确定性来源。因此，我们通过假设劳动收入和股息的对数差值遵循均值回归过程，建立了两者之间的真实关系。劳动收入和风险资产的价格趋势具有协整性，意味着它们在长期内具有相同的趋势。当协整关系较强时，退休收入的波动性大于基金组合回报的波动性。我们的模型通过考虑劳动市场和股市之间的关系，更加贴近现实的金融市场。此外，我们还考虑了个体性劳动收入冲击的不确定性。

因此，除了让劳动收入冲击的对数遵循算术布朗运动，我们的模型进一步提出劳动冲击的增长率可能随时间随机变化。我们通过一个连续时间二态隐马尔可夫链来建模这种不确定性，这增加了解HJB方程的数学难度。此外，本文首次考虑了个体性劳动收入冲击的不确定性，并研究了劳动和股市之间协整关系对CDC养老金系统中福利替代比例的影响。这些假设使得我们的模型在寿险领域中更加实用。

本文其余部分的组织结构如下：第2节介绍了资产和劳动收入的公式化以及相关假设；第3节推导了最优投资策略和福利替代率的显式表达式，并展示了一个应用实例；第4节通过数值实例来说明我们的结果；第5节提供了结论性意见。

##  模型的公式化

设$T>0$为有限时间区间，定义三个标准布朗运动$Z_D=\{Z_D(t),t\geq0\},Z_1=\{ Z_1(t),t\geq0 \},Z_L=\{Z_L(t),t\geq0  \}$.它们定义在概率空间$(\Omega,\mathcal{F},\mathbb{P})$上，令$\mathbb{F}=\{\mathcal{F}_t,t\geq0\}$是由布朗运动生成的增广滤波。假设三个标准布朗运动都是独立的。本文假设在一个过滤完整的概率空间中，对于$p\geq1,$定义
$$
\begin{align}
&L^p_{\mathcal{F}_t}(\Omega;\mathbb{R})=\left\{ X:\Omega\rightarrow\mathbb{R}\big\vert X(t)是\mathcal{F}_t可测的，\mathbb{E}[|X(t)|^p]<\infin \right\},\\
&L^2_\mathbb{F}(s,t;\mathbb{R})=\left\{ X:[s,t]\times \Omega\rightarrow\mathbb{R}\big|X(t)是\mathbb{F}适应的，\mathbb{E}\left[ \int_s^t|X(\nu)|^2d\nu \right]<\infin \right\},\\
&L^p_\mathbb{F}(\Omega;L^2(s,t;\mathbb{R}))=\left\{ X:[s,t]\times\Omega\rightarrow \mathbb{R}\big|X(t)是\mathbb{F}适应的，\mathbb{E}\left[\left( \int_s^t|X(\nu)|^2d\nu\right)^p \right]<\infin  \right\},\\
&L^p_{\mathbb{F}}(\Omega;C([s,t];\mathbb{R}))=\left\{ X:[s,t]\times\Omega\rightarrow\mathbb{R}\big|X(t)是有界的\mathbb{F}适应且有连续路径，\mathbb{E}\left[ \sup_{\nu\in[s,t]}|X(\nu)|^p \right]<\infin \right\}.
\end{align}
$$

### 金融市场

在我们的模型中我们假设基金经理可以在时间区间$[0,T]$内对有一种无风险资产和一种风险资产组成的金融市场进行投资。风险资产向投资者支付持续的股息流，令$D(t)$表示风险资产在时间$t$的股息过程。股息的动态过程表示为：
$$
\begin{cases}
\frac{dD(t)}{D(t)}=g_Ddt+\sigma dZ_D(t),\quad t\in[0,T]\\
D(0)=d_0>0
\end{cases}
$$
其中$g_D,\sigma$分别是股息的增长率和波动率。令定价核$M(t)$的动态表示为：
$$
\frac{dM(t)}{M(t)}=-rdt-\lambda_m dZ_D(t)
$$
其中$r>0$表示常数的无风险利率，$\lambda_m$表示常数的风险价格。

> [!Note]
>
> 定价核(pricing kernel, AKA Stochastic Discount Factor, SDF)是一个用来将未来现金流折现到当前时刻的工具，反映了时间和不确定性对未来现金流的影响，通常用于描述金融资产价格的动态和资产的风险溢价。

令$X(t)$为风险资产的价格过程。资产价格可以通过以下方式描述，即股息的贴现总和：
$$
X(t)=\int^\infin_t\mathbb{E}_t\left[ \frac{M(s)}{M(t)}D(s) \right]ds.
$$
因此，我们可以推导出$D(t)$和$X(t)$成正比，即
$$
X(t)=\frac{D(t)}{r+\lambda_m\sigma-g_D},
$$
推导如下图

<img src="https://s2.loli.net/2025/01/14/EYsFTLcvHl7aPSx.jpg" alt="xtdt" style="zoom:55%;" />

由于股息和股票价格是常数，于是风险资产的动态过程服从下面的几何布朗运动：
$$
\begin{cases}
\frac{dX(t)}{X(t)}=g_Ddt+\sigma dZ_D(t),\quad t\in[0,T],\\
X(0)=x_0>0.
\end{cases}
$$
令$S(t)$为超额盈余过程，随着时间$t$变化可以描述为
$$
\begin{cases}
\frac{dS(t)}{S(t)}=\frac{dX(t)+D(t)dt}{X(t)}=\mu dt+\sigma dZ_D(t),\quad t\in[0,T],\\
S(0)=s_0>0.\tag{2.1}\label{2.1}
\end{cases}
$$
其中，期望回报由定价核的定义得到$\mu=r+\lambda_m\sigma.$无风险资产$S_0(t)$的动态过程如下
$$
\begin{cases}
\frac{dS_0(t)}{S_0(t)}=rdt,\\
S_0(0)=s_{00}>0.\label{2.2}\tag{2.2}
\end{cases}
$$

### 劳动收入

假设劳动收入$L(t)$是总劳动收入$L_1(t)$和成员的个体性冲击$L_2(t)$的乘积，在对数化模型中有
$$
l(t)=l_1(t)+l_2(t),\label{2.3}\tag{2.3}
$$
文章通过让总劳动收入与股市之间的对数差值服从均值回归过程，来模拟总劳动收入和股市的协整性。令差值$y(t)$满足
$$
y(t)=\log L_1(t)-\log D(t)-\lambda,\label{2.4}\tag{2.4}
$$
其中正的常数$\lambda$总劳动收入与股息的长期对数比率。差值$y(t)$的动态过程由下面的方程描述
$$
\begin{cases}
dy(t)=-ky(t)dt+v_LdZ_L(t)-v_DdZ_D(t),\quad t\in[0,T],\\
y(0)=y_0.
\end{cases}
$$
其中$k$决定了变量$y(t)$向长期均值回归的速度并捕捉了总劳动收入与股息之间的协整性（即当$k=0$时二者不存在协整性）。$v_L,v_D$分别是条件波动率，$Z_L(t)$是与总劳动收入不确定性相关的SBM.

个体性冲击的对数过程由下面的式子描述
$$
\begin{cases}
dl_2(t)=\left( \alpha(t)-\frac{v_1^2}{2} \right)dt+v_1dZ_1(t),\quad t\in[0,T],\\
l_2(0)=l_{20}.
\end{cases}
$$
其中$\alpha(t)$是增长率,$v_1$是相应的波动率，$Z_1(t)$是与$Z_L(t),Z_D(t)$相互独立的SBM.

### 养老金系统

本文考虑的养老金保险计划中，成员分为两组，在职成员是指向养老基金缴纳费用的工作成员，而退休成员从养老基金中获取福利。所有成员假设从$a_0$开始加入计划，直到退休年龄$a_1$为止，且死亡年龄为$a_2.$我们还假设生存函数$s(x),$且$s(a_0)=1$​.

用$n(t)$表示在时间$t$时年龄为$a_0$的新成员加入养老金计划的密度。$n(t)$是一个非负函数表示新成员加入养老金计划的密度不可能为负。然后在时间$t$时年龄为$x$的人数为
$$
n(t-(x-a_0))s(x),\quad t>0.
$$
其中$t-(x-a_0)$可能为负，这意味着年龄为$x$的个体是在$x-a_0$年前加入计划的。当年龄为$x$的个体尚未加入养老金计划时，$n(t-(x-a_0))=0.$

劳动力市场中在职成员和退休成员的总人数分别用$M_1(t),M_2(t)$表示：
$$
M_1(t)=\int_{a_0}^{a_1}n(t-(x-a_0))s(x)dx,\\
M_2(t)=\int_{a_1}^{a_2}n(t-(x-a_0))s(x)dx.
$$
在职成员的数量决定了总的缴费率，假设$C_0$时时间$0$时的即时缴费率，$\eta_1$是缴费的指数增长率，因此养老基金在时间$t$时的总缴费率为:
$$
C(t)=\int_{a_0}^{a_1}n(t-x+a_0)s(x)C_0e^{\eta_1t}dx.
$$
退休成员的数量决定了养老基金的总福利、薪资结构以及初始年度养老金的支付率。初始养老金支付率被假定为退休时最终薪资的一定比例。对于在时间$t$退休的成员，初始养老金的支付额为$b(t)L(t),$其中$b(t)$是替代率可以视为控制变量。对于在时间$t$,年龄为$x$的成员，最终薪资是$x-a_1$年前的薪水。为了确定$x$岁的退休成员在时间$t$的年度养老金支付率，引入一个新的量$FL(x,t)$表示该成员在退休$x-a_1$年后的假定薪资，这个量定义为
$$
FL(x,t)=L(t)e^{-\eta_0(x-a_1)},\; t\geq0,\; x\geq a_1.
$$
其中$L(t)$表示在时间$t$退休成员的薪资（但实际上这时候已经退休了），假定薪资通过指数增长率$\eta_0$进行确定性地向后推算。这个方法与成员退休时的实际薪资不同，尤其是在$\eta_0>0,x>a_1$时，实际薪资与假定薪资之间的差异随着年龄增大而增大。我们假设对假定的薪资率应用一个调整因子，因此时间$t$时年龄为$x$的成员的养老金支付率为：
$$
B(x,t)=b(t)L(t)e^{-\eta_0(x-a_1)}.
$$
所有退休成员在时间$t$时的实际总退休福利（涵盖从$a_1\sim a_2$的退休成员）可以通过下式得到
$$
B(t)=\int_{a_1}^{a_2}n(t-x+a_0)s(x)B(x,t)dx=F(t)b(t)L(t),
$$
其中$F(t)$是一个正的函数，定义为
$$
F(t)=\int_{a_1}^{a_2}n(t-x+a_0)s(x)e^{-\eta_0(x-a_1)}dx.
$$
为了增加养老金的福利，基金经理将在缴费和福利之间的盈余部分进行动态投资。在本文的模型中，金融市场由一个风险资产和一个无风险资产组成，设$\pi(t)$是时间$t$投资于风险资产的比例。则养老金基金的财富过程为
$$
\begin{cases}
dW(t)=\pi(t)W(t)\frac{dX(t)+D(t)dt}{X(t)}+(1-\pi(t))W(t)\frac{dS_0(t)}{S_0(t)}+C(t)dt-B(t)dt,\\
w(0)=w_0>0.\tag{2.5}\label{2.5}
\end{cases}
$$

>[!Note]
>
>- $W(t)$是养老金基金的财富过程
>- $\pi(t)$是投资于风险资产的比例
>- $X(t)$是风险资产的价格过程
>- $D(t)$是风险资产的股息支付
>- $S_0(t)$是无风险资产的价格过程
>- $C(t)$是缴费项
>- $B(t)$是福利支付项

接下来我们定义可接受的策略和本文研究的主要问题。

>**Definition 2.1**
>
>对于任何固定$t\in[0,T],$策略对$(\pi(t),b(t))$被称为可接受的，如果它满足以下条件：
>
>- 投资策略和初始福利支付政策$(\pi(t),b(t))$是$\mathcal{F}_t$适应的，以使得SED$(\ref{3.4})$存在唯一解$W_{\pi,b}(t).$
>- $\pi(t)\in L_{\mathbb{F}}^2(0,T;\mathbb{R}^+)$且$b(t)\in L_{\mathbb{F}}^2(0,T;\mathbb{R}^+)$对所有$t>0$成立。

>**Problem 2.1**
>
>对于初始状态$(t,W_t),$养老金基金经理的目标函数是最大化：
>$$
>J(t,w,l,y)=\mathbb{E}_{\pi,b}\left[ \int_t^Te^{-rs}U(b(s)F(s)L(S))ds+\lambda_1e^{-rT}U(W(T)) \right],\tag{2.6}\label{2.6}
>$$
>其中$\lambda_1$是一个非负常数表示对终期财富带来的效用的权重。$\mathbb{E}_{\pi,b}$是在概率测度$\mathbb{P}$下给定$W(t)=w,L(t)=l,y(t)=y$时的条件期望。于是这个问题的价值函数就是
>$$
>V(t,w,l,y)=\sup_{(\pi,b)\in\mathcal{A}}J(t,w,l,y),
>$$
>其中$\mathcal{A}$是一组控制对的集合。

## 养老金基金计划的最优策略

本节通过标准的动态规划方法研究这个随机最优控制问题。当劳动市场和股票市场是协整关系时得出了最优控制的显式表达式。本文考虑两种劳动收入过程的情况，具有未知增长率的特质性冲击以及不考虑特质性冲击。

### 具有特质性冲击的最优策略

这一小节中假设员工面临关于劳动收入的不确定性，通过考虑具有未知增长率的特质性冲击来捕捉这种不确定性。本文使用一个连续时间的二状态的隐马尔科夫链来描述劳动收入冲击的内在不确定性和动态特性。具体来说，$\alpha(t)$是未知的增长率，使用连续时间的二状态隐马尔可夫链在$(\Omega,\mathcal{F},\mathbb{P})$上建模，并在$\alpha_1,\alpha_2$之间变化，其中$\alpha_1>\alpha_2.$增长率$\alpha(t)$可能取高值$\alpha_1$或低值$\alpha_2.$在一个小时间间隔$\Delta t$内，增长率$\alpha(t)$在时刻$t$取$\alpha_1$的概率是$1-p_1\Delta t$,保持$\alpha_2$的概率为$1-p_2\Delta t$，其中$p_1,p_2$是二状态隐马尔科夫链的转移强度，从$\alpha_1$转移到$\alpha_2$的概率为$p_1\Delta t,$从$\alpha_2$转移到$\alpha_1$的概率为$p_2\Delta t.$对于任何的$t$设$P(t)$为条件概率，表示给定观测信息$\mathcal{F}_t,$增长率$\alpha(t)$取$\alpha_1$的概率，
$$
P(t)=\mathbb{P}(\alpha(t)=\alpha_1|\mathcal{F}_t).
$$
对数化特质性冲击的预期增长率$\mu_0(t)$是两种可能增长率的加权平均值，给定为
$$
\mu_0(t)=P(t)\alpha_1+(1-P(t))\alpha_2=\alpha_2+\beta P(t),\label{31}\tag{3.1}
$$
其中$\beta=\alpha_1-\alpha_2.$对于一个小时间间隔$\Delta t,$对数特质性劳动冲击的总变化是$l_2(t+\Delta t)-l_2(t),$其预期变化为$\mu_0(t)\Delta t-\frac{v_1^2}{2}$.对相应波动率进行归一化得到
$$
dZ_1(t)=\frac{1}{v_1\Delta t}\left(l_2(t+\Delta t)-l_2(t)-\mu_0(t)\Delta t-\frac{v_1^2}{2}\right).\tag{3.2}\label{3.2}
$$
将式子$(\ref{3.2})$带入式子$(\ref{31})$就得到
$$
dl_2(t)=(\alpha_2+\beta P(t)-\frac{v_1^2}{2})dt+v_1dZ_1(t).
$$
我们还可以得到$P(t)$的动态过程
$$
dP(t)=(p_2-(p_1+p_2)P(t))dt+v_1^{-1}\beta P(t)(1-P(t))dZ_1(t).
$$
由式$(\ref{2.3}),(\ref{2.4})$可知对数劳动收入可以写作
$$
l(t)=y(t)+d(t)+\lambda+l_2(t),
$$
其中$d(t)=\log D(t)$且其动态过程可以写作
$$
\begin{cases}
d\,d(t)=(g_D-\frac{\sigma^2}{2})dt+\sigma dZ_D(t),\quad t\in[0,T],\\
d(0)=d_0.
\end{cases}
$$
应用Ito公式，对数化劳动收入就导出为
$$
d\,l(t)=[-ky(t)+g_D-\frac{\sigma^2}{2}+\lambda+\beta P(t)+\alpha_2-\frac{v_1^2}{2}]dt\\
+(\sigma-v_D)dZ_D(t)+v_LdZ_L(t)+v_1dZ_1(t).
$$
使用Ito公式就可得到
$$
\begin{cases}
&\frac{dL(t)}{L(t)}&=[-ky(t)+g_D-\frac{\sigma^2}{2}+\lambda+\beta P(t)+\alpha_2+\frac{v_L^2}{2}+\frac{1}{2}(\sigma-v_D)^2]dt\\
& &+(\sigma-v_D)dZ_D(t)+v_LdZ_L(t)+v_1dZ_1(t),\quad t\in[0,T].\\
&L(0)&=l_0.\label{3.3}\tag{3.3}
\end{cases}
$$
结合式子$(\ref{2.1}),(\ref{2.2})$​就可以重写
$$
\begin{cases}
\frac{dW(t)}{W(t)}=[(\mu-r)\pi+r+\frac{C(t)}{W(t)}-\frac{F(t)L(t)b(t)}{W(t)}]dt+\pi\sigma dZ_D(t),\\
W(0)=w_0>0,\tag{3.4}\label{3.4}
\end{cases}
$$
其中$L(t),P(t)$的动态过程也已经给出。

随后我们用$P$表示$P(t)$,$L(t)$记作$l,W(t)$记作$w,$稍作符号上的简化，对于初始状态$(t,W_t)$,基金经理的目标是最大化以下的目标函数
$$
J(t,w,l,y,P)=\mathbb{E}_{\pi,b}\left[ \int_t^Te^{-rs}U(b(s)F(s)L(s))ds+\lambda_1e^{-rT}U(W(T))\right],\label{3.6}\tag{3.6}
$$
该问题的价值函数由下面的公式给出
$$
V(t,w,l,y,P)=\sup_{(\pi,b)\in\mathcal{A}}J(t,w,l,y,P),\tag{3.7}\label{3.7}
$$
其中$\mathcal{A}$是控制对的集合。假设效用函数表达式如下
$$
U(w)=-\frac{1}{\gamma}e^{-\gamma w},
$$
其中$\gamma>0$是常数的绝对风险厌恶系数。基金经理的目标是最大化福利和最终的剩余财富。为了简便起见，我们定义
$$
\begin{align}
C^{1,2,2,2}&([0,T]\times\mathbb{R}^+\times\mathbb{R}^+\times\mathbb{R}^+\times\mathbb{R}^+)\\
=&\{V(t,w,l,y,P)|V(t,\cdot,\cdot,\cdot,\cdot) \text{ is once continuously differentiable on }[0,T],\\
&V(\cdot,w,l,y,P)\text{ is twice continuously differentiable for }\\
&w\in \mathbb{R}^+,l\in\mathbb{R}^+,y\in\mathbb{R}^+,P\in\mathbb{R}^+.\}
\end{align}
$$
使用标准方法即可得到下面的HJB方程满足$V(t,w,l,y,P)\in C^{1,2,2,2} $


$$
0 = \sup_{\pi, b} \left\{ 
\begin{align}

    V_t &+ w V_w \left[ r + (\mu - r) \pi + \frac{C(t)}{w}-\frac{F(t)lb}{w} \right] + \frac{1}{2} w^2 V_{ww} \sigma^2 \pi^2 \\
    &+ lv_L \left[ -ky + g_D + \lambda - \frac{1}{2} \sigma^2 + \beta P + \alpha_2 + \frac{1}{2} v^2_L +\frac{1}{2}(\sigma-v_D)^2\right] \\
    &+ \frac{1}{2} l^2 V_{ll} \left[ (\sigma - v_D)^2 + v_1^2+v_L^2 \right]-kyV_y+\frac{1}{2}V_{yy}(v_D^2+v_L^2)+wlV_{wl}\pi\sigma(\sigma-v_D) \\
    &- ky V_y + \frac{1}{2} V_{yy} (v^2_D + v^2_L) + w l V_{wl} \sigma (\sigma - v_D) \\
    &- w V_{wy} \pi \sigma v_D + l V_{ly} \left[ v^2_L - v_D (\sigma - v_D) + lP V_{lP} \beta(1 - P) + V_P \left[p_2 - (p_1 + p_2) P \right] \right] \\
   &+ \frac{1}{2v_1^2}  V_{P P} P^2 (1 - P)^2 \beta^2 - e^{-rt} \frac{1}{\gamma} e^{-\gamma F(t) lb
   
\end{align}
\right\}
$$


下面是推导过程及结果

<img src="https://s2.loli.net/2025/01/15/9fDl5RJTwaBO3Yb.jpg" alt="HJB方程的推导（变分法）" style="zoom:45%;" /><img src="https://s2.loli.net/2025/01/15/zROLChYXqIPeAmK.jpg" alt="HJB推导" style="zoom:50%;" />

上面最后一个式子记为$(3.8)$,再加上边界条件
$$
V(T,w,l,y,P)=-e^{-rT}\frac{\lambda_1}{\gamma}e^{-\gamma w}.\tag{3.9}\label{3.9}
$$


下面的定理可以阐述这个随机控制问题的结果

>**Theorem 3.1**
>
>对于任意的$t\in[0,T],$最优投资策略和福利调整策略分别由下式给出
>$$
>\begin{align}
>&\pi^*(t,w,l,y,P)=\frac{\mu-r}{\gamma f_1(t)\sigma^2 w},\tag{3.10}\label{3.10}\\
>&b^*(t,w,l,y,P)=\frac{\ln \lambda_1+\ln f_1(t)-\gamma f_1(t)w-\gamma f_2(t)-\gamma f_5(P)}{-\gamma F(t)l},\label{3.11}\tag{3.11}
>\end{align}
>$$
>且最终相应的价值函数就是
>$$
>V(t,w,l,y,P)=-\frac{\lambda_1}{\gamma}e^{-\gamma[f_1(t)w+f_2(t)+f_5(P)]-rt},
>$$
>其中
>$$
>\begin{align}
>f_1(t)&=\left[ e^{-\int_t^Trds}+\int_t^Te^{-\int_t^srdu}ds \right]^{-1},\\
>f_2(t)&=\int_t^Te^{-\int_t^sf_1(u)du}\times\left[ f_1(s)\left( C(s)-\frac{1-\ln f_1(s)-\ln\lambda_1}{\gamma} \right)+\frac{1}{2}\frac{(\mu-r)^2}{\gamma\sigma^2}+\frac{r}{\gamma} \right]ds-f_5(P(T)).
>\end{align}
>$$
>$f_5(P)$满足下面的常微分方程，其中$0<P<1,$
>$$
>f_1(t)f_5(P)=f_5'(P)[p_2-(p_1+p_2)P]-\frac{1}{2v^2}[\gamma f'_5(P)-f''_5(P)]P^2(1-P)^2\beta^2,
>$$
>其边界条件是
>$$
>\begin{cases}
>f_1(t)f_5(0)=p_2f'_5(0),\\
>f_(t)f_5(1)=-p_1f'_5(1).
>\end{cases}
>$$

>**Proof.**
>
>由HJB方程$(3.8)$可以导出两个控制的一阶条件为
>$$
>\begin{align}
>&0=wV_w(\mu-r)+w^2V_{ww}\sigma^2\pi+wlV_{wl}\sigma(\sigma-v_D)-wV_{wy}\sigma v_D,\\
>&0=-V_{w}F(t)l+F(t)le^{-rt}e^{-\gamma F(t)lb}.
>\end{align}
>$$
>因此最优投资策略和福利调整策略由下式给出
>$$
>\begin{align}
>&\pi^*=\frac{V_w\sigma v_D-V_w(\mu-r)-lV_{wl}\sigma(\sigma-v_D)}{wV_{ww}\sigma^2},\\
>&b^*=\frac{\ln V_w+rt}{-\gamma F(t)l}.
>\end{align}
>\tag{3.13}\label{3.13}
>$$
>猜测解的形式为
>$$
>V(t,w,l,y,P)=-\frac{\lambda_1}{\gamma}\exp\{ -\gamma[f_1(t)w+f_2(t)+f_3(t,l)+f_4(t)y+f_5(P)]-rt \},
>$$
>且有边界条件$f_1(T)=1,f_2(T)=-f_5(P(T)),f_3(T)=f_4(T)=0.$





## 数值分析

本节探究金融市场中的参数对于最优投资分配和福利策略的影响。

首先假设$\mu_1(x)=A+Bc^x$是年龄$x$岁时的死亡力(Force of Mortality)，于是存活函数就可以写为
$$
\begin{align}
s(x)&=e^{-\int_0^{x-a_0}\mu(a_0+s)ds}\\
&=e^{-A(x-a_0)-\frac{B}{\ln c}(c^x-c^{a_0})},\;a_0\leq x\leq a_1.
\end{align}
$$
基础参数设置为$A=0.00022,B=2.7\times10^{-6},c=1.124,a_0=30,a_1=65,a_2=100,n=10$.

