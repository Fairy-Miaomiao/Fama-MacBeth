# Fama-Macbeth回归

## 一、理论概述

Fama-Macbeth回归是1973年Fama和Macbeth为验证CAPM模型而提出的一种因子统计方法，该模型现如今被广泛用被广泛用于计量经济学的面板数据分析，而在金融领域在用于多因子模型的回归检验，用于估计各类模型中的因子暴露和因子收益（风险溢价）。

Fama-Macbeth与传统的截面回归类似，本质上也与是一个两阶段回归，不同的是它用了巧妙的方法解决了截面相关性的问题，从而得出更加无偏，相合的估计。

**时间序列回归**

Fama-Macbeth模型与传统截面回归相同，第一步都是做时间序列回归。在因子分析框架中，时间序列回归是为了获得个股在因子上的暴露。如果模型中的因子是 portfolio returns（即使用投资组合收益率作为因子，例如Fama-French三因子模型中的SMB，HML和市场因子），那么可以通过时间序列回归（time-series regression）来分析$E[R_i]$和$\beta_i$在截面上的关系。

令$f_t$为因子组合在t期的收益率，$R_{it}$为个股i在t期的收益率，用$f_t$对每只股票的$R_{it}$回归，即可得到每只股票的全样本因子暴露$\beta_i$。

$R_{it}=\alpha_i+\beta_if_t+\epsilon_{it},t=1,2,...,T,\forall i$ （1）

也可滚动计算某个时间段的因子暴露$\beta _{it}$，体现个股随市场的变化设置时间段长度为period

$R_{ik}=\alpha_i+\beta _{it}f_k+\epsilon_{ik},k=t-period,2,...,t,\forall i$ （2）

**截面回归**

传统截面回归的第一步是通过时间序列回归得到个股暴露，这一步与Fama-Macbeth回归相同，而第二步回归体现了传统截面回归和Fama-Macbeth的最大不同。

传统截面回归：

在时序回归中回归式在时间序列上取均值，在$E[\epsilon]=0$的假设下可以得出：

$E[R_i]=\alpha _i+\beta _iE[f]$ （3）

上式正是个股的期望收益与因子暴露在截面上的关系，截距$\alpha_i$为个股的错误定价。

那么便可通过截面回归找到因子的期望收益率$E[f]$，方法是最小化个股定价错误$\alpha _i$的平方和。对个股的的收益在时序上取均值得到个股期望收益$E[R_i]$，用全样本的个股因子暴露对个股期望收益做无截距回归。

$E[R_i]=\beta _i \lambda+\alpha _i$ （4）

回归残差$\alpha _i$为个股的错误定价，$\lambda$为因子的期望收益率。

截面回归最大的缺陷在于忽略了截面上的残差相关性，使得OLS给出的标准误存在巨大的低估。

**Fama-Macbeth回归**

与截面回归相同，Fama-Macbeth回归第一步是通过时间序列回归得到因子暴露值，不同的是，第二步中，Fama-Macbeth在每个t上都做了一次无截距截面回归：

$R_{it}=\beta _i \lambda _t+\alpha _{it},i=1,2,...,N, \forall t$ (5)

上式中的$\beta_i$为全样本$\beta$，当然若使用滚动回归数据，也可以在不同截面的回归上使用对应时期的$\beta _{i,t}$。

Fama-Macbeth回归相当于在每个t上做一次独立的截面回归，这T次回归的参数取均值作为回归的估计值：

$\hat{\lambda}=\frac{1}{T}\sum_{t=1}^T\hat{\lambda}_t, \hat{\alpha}_i=\frac{1}{T} \sum_{t=1}^T\hat{\alpha}_{it}$ (6)

上述方法的巧妙之处在于它把 T 期的回归结果当作 T 个独立的样本。参数的 standard errors 刻画的是样本统计量在不同样本间是如何变化的。在传统的截面回归中，我们只进行一次回归，得到$\lambda$和$\alpha _i$的一个样本估计。而在Fama-Macbeth截面回归中，把T期样本点独立处理，得到T个$\lambda$和$\alpha _i$的样本估计。

若使用全样本因子暴露$\beta _i$进行估计,截面回归和Fama-Macbeth的估计结果相同，当使用滚动窗口进行估计时（Fama and MacBeth (1973)中作者使用了滚动窗口），截面回归和Fama-Macbeth回归会得到完全不同的估计结果。

Fama-Macbeth回归很好的解决了截面相关性的问题，但对于时间序列上的相关性仍然无力。



## 二、实现

为简单说明Fama-Macbeth两阶段回归的主要步骤，以下用投资组合数据估计一个简单的 CAPM 模型。数据主要使用了「25 Portfolios Formed on Size and Book-to-Market」 中的 25 个投资组合 1926.7-2020.10 期间的月度收益率，和「Fama/French 3 Factors」 中的无风险收益、市场超额收益数据。主要变量如下：