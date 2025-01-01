# Molecular Dynamics Simulation
## Velocity Verlet
$$
{\displaystyle 
{\begin{aligned}
\mathbf {x} (t+\Delta t)&=\mathbf {x} (t)+\mathbf {v} (t)\,\Delta t+{\tfrac {1}{2}}\,\mathbf {a} (t)\Delta t^{2},\\
\mathbf {v} (t+\Delta t)&=\mathbf {v} (t)+{\frac {\mathbf {a} (t)+\mathbf {a} (t+\Delta t)}{2}}\Delta t.
\end{aligned}}}
$$
## Lennard-Jones Potential
### Potential
$$
 U(r)=4\epsilon \left[\left({\frac {\sigma }{r}}\right)^{p}-\left({\frac {\sigma }{r}}\right)^{q}\right]　
$$
### Force
$$
F(r)=-{\frac {d}{dr}}U(r)=4\epsilon \left(p{\frac {{\sigma }^{p}}{r^{p+1}}}-q{\frac {{\sigma }^{q}}{r^{q+1}}}\right)
$$
### Dimensionless 
|Property|Symbol|Reduced Form|
|-|-|-|
|長さ|$$r^{*}$$|$$\frac {r}{\sigma}$$|
|時間|$$t^{*}$$|$$t{\sqrt {\frac {\varepsilon }{m\sigma ^{2}}}}$$|
|温度|$$T^{*}$$|$$\frac {k_{B}T}{\varepsilon }$$|
|力|$$F^{*}$$|$${\frac {F\sigma }{\varepsilon }}$$|
|エネルギー|$$U^{*}$$|$$\frac {U}{\varepsilon }$$|
|圧力|$$p^{*}$$|$$\frac {p\sigma ^{3}}{\varepsilon }$$|
|密度|$$\rho ^{*}$$|$$\rho \sigma ^{3}$$|
|表面張力|$$\gamma ^{*}$$|$$\frac {\gamma \sigma ^{2}}{\varepsilon }$$|
### Patameter
|  | $\sigma$[nm] | $\epsilon$[J] | $\epsilon/k_{B}$[K] | m [kg] |
|:---:|:---:|:---:|:---:| :---:|
| Ne | 0.274 | $0.50 \times 10^{-21}$ | 36.2 |  
| Ar | 0.340 | $1.67 \times 10^{-21}$ | 121 |$6.634 \times 10^{-26}$|   
| Kr | 0.365 | $2.25 \times 10^{-21}$ | 163 |
| Xe | 0.398 | $3.20 \times 10^{-21}$ | 232 |   