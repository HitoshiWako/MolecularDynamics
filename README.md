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
## 蒸発源
任意の確率分布を作るには,作りたい確率密度関数を$f(x)$とする。$f(x)$の累積分布関数を$F(x)$としたとき一様分布乱数$R$を用いれば求めたい確率密度をもつ乱数$X$は
$$
X = F^{-1}(R) 
$$
にて計算される 。
[任意の確率密度を持つ乱数を作る #Python](https://qiita.com/kaityo256/items/95b4a3f61c963f08c899)
今回、蒸発源から余弦則に従い粒子が放出されているとする。この際の確率密度関数は
$$
f(x) = \cos x
$$
この$f(X)$について累積分布関数$F(x)$は区間が$-\pi/2\sim\pi/2$であるので
$$
F(x)=\frac{1}{2}\int_{-\pi/2}^{x}\cos\theta d\theta = \frac{\sin x + 1}{2}
$$
ここで$1/2$は出力範囲を0〜1の範囲にする係数である。逆関数は
$$
F^{-1}(R) = \sin^{-1}(2R-1)
$$
となる。
## 温度と速度
熱運動速度は以下のような関係式で表せる
$$
v = \sqrt{\frac{\alpha k_{B}T}{m}}
$$
ここで$\alpha$は自由度とする
無次元化した速度は
$$
v^{*} = v\sqrt{\frac{m}{\epsilon}}
$$
なので温度$T$と$v^{*}$の関係式は
$$
v^{*} = \sqrt{\frac{m}{\epsilon}}\sqrt{\frac{\alpha k_{B}T}{m}}=\sqrt{\frac{\alpha k_{B}T}{\epsilon}}=\sqrt{\alpha T^{*}}
$$
## 系の温度
[コンピュータ利用による薄膜作製シミュレーション](https://www.newglass.jp/mag/TITL/maghtml/27-pdf/+27-p315.pdf)
系の温度は
$$
T = \frac{2}{\alpha k_{B}N}\frac{1}{2}\sum_{i}^{N}m_{i}\cdot v_{i}^2
$$