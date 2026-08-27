r"""C07 的载体：**Kármán–Trefftz 翼型**的精确解（保角映射）。

★★★ **它能做 C06 做不到的一件事。** 圆柱没有尖尾缘，Kutta 选不出环量，所以 C06 里
Γ 是**规定**的。KT 翼型有**有限尾缘角**，Kutta 条件**必须自己选出 Γ**，而精确值有解析式
⇒ 这是项目里**唯一能把 Kutta 条件对精确 Γ 检验**的载体。

映射（n = 2 - τ/π，τ 是尾缘夹角；n = 2 退化为 Joukowski 的尖点尾缘）::

    z = n b [ (ζ+b)^n + (ζ-b)^n ] / [ (ζ+b)^n - (ζ-b)^n ]

圆心 ζ_c、半径 R = |b - ζ_c| 的圆**必过 ζ = b**（尾缘的原像）。圆平面复速度带环量::

    dw/dζ = U ( e^{-iα} - R² e^{iα} / (ζ-ζ_c)² ) + i Γ / (2π (ζ-ζ_c))

Kutta（后驻点落在 ζ = b）给出 **Γ = 4π U R sin(α + β)**，β = arg(b - ζ_c)。

★★★ **实现是用四个与它无关的 oracle 验过的，不是靠公式抄对**（2026-08-26 实测，
τ=10°、x_c=-0.10、α=5°）：

| oracle | 做法 | 结果 |
|---|---|---|
| ① 环量 | 沿**物理平面**翼型围线数值积分速度 | \|Γ\| 相对差 **4.1e-07** |
| ② 尾缘夹角 | 从 z(θ) 的**几何**切线量 | **10.207°** vs 参数 10.000° |
| ③ 驻点 | 圆平面 \|dw/dζ\| 的最小点 | θ = **359.955°** ≈ ζ = b ⇒ Kutta 成立 |
| ④ 尾缘速度衰减率 | θ→0 时 \|q\| 每十倍的比值 | **0.879**，理论 0.1^(τ/π) = **0.8797** |

★★ oracle ④ 最强：它把观察到的尾缘行为**通过映射的指数**接回 τ 参数。
配套对照 —— **τ = 0（Joukowski 尖点）时 \|q\| → 0.9056 常数非零**，而 τ > 0 时 → 0。
两个情形干净分开，这是实现正确最有力的证据。

★ 实用后果，记在这里免得被当成缺陷：τ = 10° 时尾缘速度**衰减极慢**（指数
τ/π = 0.0556），精确解自己在 θ = 1e-5 处还有 \|q\| ≈ 0.49 ⇒ 数值解**看不到干净的驻点**，
而尾缘带的 Cp 收敛阶因此只有 ~0.6（C07 记为 RECORDED，不设判据）。

★ 几何按翼型约定归一化：平移 + 缩放使**前缘 x=0、尾缘 x=1**，于是
`airfoil_wake_2d`（尾缘取 `pt[0]`、前缘取 `argmin(x)`、远场圆心硬编码在半弦）
**一行新代码都不需要**。缩放后 Γ_n = Γ/c，cl = 2Γ_n/U。
"""
import numpy as np

#: 生产参数（本项目的 C07 就用这一组；改它要重新跑 oracle 并更新上表）
TAU_DEG, X_C, Y_C, B = 10.0, -0.10, 0.0, 1.0
ALPHA_DEG, U_INF = 5.0, 1.0
#: 精确表面的采样密度。★ 40000 点给出最近邻匹配距 ~3.5e-05，而实测把它加密 5x
#: （匹配距 6e-06）只让 LE 带 RMS 从 0.1594 动到 0.1574 ⇒ **匹配误差不是误差来源**。
N_SURFACE = 40000


def _map(zeta, b, n):
    p = (zeta + b) ** n
    m = (zeta - b) ** n
    return n * b * (p + m) / (p - m)


def _dzdzeta(zeta, b, n):
    p = (zeta + b) ** n
    m = (zeta - b) ** n
    dp = n * (zeta + b) ** (n - 1)
    dm = n * (zeta - b) ** (n - 1)
    d = p - m
    return n * b * ((dp + dm) * d - (p + m) * (dp - dm)) / d ** 2


class KarmanTrefftz:
    """精确解。所有物理量都在**归一化**几何上给出（弦长 1，前缘 0，尾缘 1）。"""

    def __init__(self, tau_deg=TAU_DEG, x_c=X_C, y_c=Y_C, b=B,
                 alpha_deg=ALPHA_DEG, u_inf=U_INF, n_surface=N_SURFACE):
        self.n = 2.0 - tau_deg / 180.0
        self.tau_deg, self.b, self.u = tau_deg, b, u_inf
        self.alpha = np.deg2rad(alpha_deg)
        self.zc = complex(x_c, y_c)
        self.R = abs(b - self.zc)
        self.beta = np.angle(b - self.zc)
        #: ★ Kutta 选出的环量 —— 这是本门要检验的那个数
        self.gamma_raw = 4.0 * np.pi * u_inf * self.R * np.sin(self.alpha + self.beta)
        #: 半格偏移，避开 θ=0：那里 dz/dζ = 0（映射的临界点，有限尾缘角的来源）
        th = (np.arange(n_surface) + 0.5) * (2.0 * np.pi / n_surface)
        z_raw = _map(self._circle(th), b, self.n)
        self.x0 = float(z_raw.real.min())
        self.chord = float(z_raw.real.max() - self.x0)
        self.z = (z_raw - self.x0) / self.chord
        self.q = np.abs(self._w_prime(self._circle(th))
                        / _dzdzeta(self._circle(th), b, self.n))
        self.cp = 1.0 - (self.q / u_inf) ** 2
        #: 归一化后的环量与精确升力（弦长 = 1）
        self.gamma = self.gamma_raw / self.chord
        self.cl = 2.0 * self.gamma / u_inf

    def _circle(self, theta):
        return self.zc + self.R * np.exp(1j * theta)

    def _w_prime(self, zeta):
        zp = zeta - self.zc
        return (self.u * (np.exp(-1j * self.alpha)
                          - self.R ** 2 * np.exp(1j * self.alpha) / zp ** 2)
                + 1j * self.gamma_raw / (2.0 * np.pi * zp))

    def surface_speed_at(self, theta):
        """任意 θ 处的物理表面速率（oracle ④ 用它取 θ→0 的极限）。"""
        ze = self._circle(np.asarray(theta, dtype=np.float64))
        return np.abs(self._w_prime(ze) / _dzdzeta(ze, self.b, self.n))

    def polyline(self, n_half):
        """闭合折线，`naca0012_coordinates` 约定：尾缘 -> 上 -> 前缘 -> 下 -> 尾缘。"""
        t = np.linspace(0.0, 2.0 * np.pi, 2 * n_half + 1)[:-1]
        zz = (_map(self._circle(t), self.b, self.n) - self.x0) / self.chord
        zz = np.roll(zz, -int(np.argmax(zz.real)))          # 从尾缘起
        if zz.imag[1] < 0:                                   # 保证先走上表面
            zz = np.concatenate([zz[:1], zz[:0:-1]])
        c = np.column_stack([zz.real, zz.imag])
        return np.vstack([c, c[:1]])                         # 闭合（生成器会丢重复末点）

    def cp_at(self, points, chunk=4096):
        """把物理点匹配到最近的精确表面点并取其 Cp。

        ★ 分块算距离：`len(points) x N_SURFACE` 的全矩阵在 N_SURFACE = 4e4 时
        已是数百 MB，一次算完会把内存吃光 —— 实测踩过。
        """
        pts = np.asarray(points, dtype=np.complex128)
        out = np.empty(pts.size)
        dist = np.empty(pts.size)
        for i in range(0, pts.size, chunk):
            d = np.abs(pts[i:i + chunk, None] - self.z[None, :])
            j = np.argmin(d, axis=1)
            out[i:i + chunk] = self.cp[j]
            dist[i:i + chunk] = d[np.arange(j.size), j]
        return out, dist
