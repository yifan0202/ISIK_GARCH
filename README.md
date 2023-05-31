# ISIK_GARCH
A project of VaR prediction with models including GARCH, SK-GARCH and ISIK-GARCH.
## Abstract implied higher moments with Edgeworth expansion
$$C_{approx}=C+IS_TA_3+(IK_T-3)A_4+IS_T3A_6$$

$$C=S_0e^{\delta\sigma}\Phi\left(d\right)-Ke^{-rT}\Phi\left(d-\sigma\right)$$

$$A_3=\frac{1}{6}S_0e^{\delta\sigma}\sigma\left[\left(2d-d\right)\phi\left(d\right)+\sigma^{2\Phi\left(d\right)}\right]$$

$$A_4=\frac{1}{24}S_0e^{\delta\sigma}\sigma\left[\left(d^2-1-3\sigma\left(d-\sigma\right)\right)\phi\left(d\right)+\sigma^3\Phi\left(d\right)\right]$$

$$A_6=\frac{1}{72}S_0e^{\delta\sigma}\sigma\left[\sigma^5\Phi\left(d\right)+\left(3-6d^2+d^4+5\sigma\left(d-\left(d-\sigma\right)\left(\sigma d-2\right)-\left(d-\sigma\right)^3\right)\right)\phi\left(d\right)\right]$$

$$d=\frac{\log{\left(\frac{S_0}{K}\right)}+\mu}{\sigma}+\sigma$$

$$\delta=\frac{\mu-rT+\frac{\sigma^2}{2}}{\sigma}$$
## GARCH models
### GARCH
$$r_t=\mu+\sqrt h_tz_t$$

$$h_t=\omega+\alpha_1h_{t-1}+\alpha_2h_tz_t^2$$

The maximum likelihood is $log(L_{rt})=-1/2[log2π+logh_t+z_t^2]$
### SK-GARCH
$$s_t=\beta_0+\beta_1s_{t-1}+\beta_2z_{t-1}$$

$$\ k_t=\gamma_0+\gamma_1k_{t-1}+\theta_2\left|z_{t-1}\right|$$
### ISIK-GARCH
The implied skewness and kurtosis satisfy the following equations:
$$s_t=\beta_0+\beta_1s_{t-1}+\beta_2IS_{t-1},\ \ IS_t=\delta_0+\delta_1s_t+\delta_2z_t+\eta_t$$
$$\ k_t=\gamma_0+\gamma_1k_{t-1}+\gamma_2IS_{t-1},\ \ IK_t=\theta_0+\theta_1k_t+\theta_2\left|z_t\right|+\vartheta_t$$

## VaR forecast
$$\int_{-\infty}^{-VaR_{q,t+1}}{f\left(r_{t+1}\middle| I_t\right)dr_{t+1}}=\alpha$$
