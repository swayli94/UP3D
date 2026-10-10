r"""C08 的载体：**Ringleb 流** —— 唯一 2-D 跨声速精确解，且是**全速势方程本身**的解。

★★★ **它能给而 C05 给不了的东西。** 喷管里有激波 ⇒ 数值解**本该**与等熵精确解在激波处
不同，所以那里只能判**激波位置**。Ringleb 的超声速区是**无激波、光滑**的 ⇒ 人工密度在那里
加的每一分耗散**都是纯误差**。**这是 `upwind_c` / `m_crit` 第一次能对真值定价。**

hodograph 解（滞止量归一，a0 = rho0 = 1，q <= k）::

    a = sqrt(1 - (g-1)/2 q^2),  rho = a^(2/(g-1))
    J = 1/a + 1/(3a^3) + 1/(5a^5) - 1/2 ln((1+a)/(1-a))
    x = 1/(2 rho) (2/k^2 - 1/q^2) - J/2
    y = +- sqrt(1 - (q/k)^2) / (k rho q)

k 标记流线，q 沿流线变化；域取两条流线为壁、两条等速线为进出口的 **U 形通道**：
沿一条流线从上支加速到鼻子（q = k，y = 0）再沿下支减速回出口 ⇒
**进出口亚声速，超声速口袋嵌在中间的转弯处**。生产参数 K_IN = 1.1（内壁，鼻子 M = 1.264）、
K_OUT = 0.9（外壁，鼻子 M = 0.983）、Q0 = 0.45（进出口 M = 0.459）—— 峰值 Mach 刻意压在
**真实翼型激波的量级**，与喷管同一条理由。

★★★ **两个真实的符号错误是 oracle 抓出来的，写在这里因为第二个的教训更值钱：**

| # | 错误 | 抓它的 oracle | 症状 |
|---|---|---|---|
| ① | 速度角写成 `+asin(q/k)` | **PDE 本身**（∇×u、∇·(ρu)） | `\|∇×u\| ~ O(1)`；实测切向角是 **-asin(q/k)** |
| ② | 下支流向反了 | ★ **PDE 看不见**；靠**鼻子处速度连续性** | 现约定跳 1.9~2.1，取负后跳 2.8e-05 |

★★★ **② 的教训**：`∇×u = 0` 与 `∇·(ρu) = 0` 都是**齐次条件，对速度整体变号不变** ⇒
**PDE oracle 对这个符号是瞎的**，两个符号都"通过"。判别条件必须是**界面连续性**。
**一个 oracle 通过，不等于它有能力分辨你正在选的那件事。**

★★ 另外五个 oracle（全部实测）：两套密度约定一致 3.7e-16 · Mach 一致 ·
φ 路径无关 8e-09 · 逆映射往返 9.1e-15 · **∇φ = u** 6.9e-07（把求积点向鼻子端聚束之前是
4.07e-05 —— 那是 Dirichlet 数据的精度地板）。

★ 归一化换算（实测精确到 3.7e-16）：取参考 Ringleb 速度 `QS_INF`，则
`M_inf = sqrt(qs^2 / (1 - (g-1)/2 qs^2))`，而求解器的 `q_hat = q / QS_INF`。
"""
import os
import numpy as np
G = 1.4

def a_of(q):   return np.sqrt(1.0 - 0.5*(G-1.0)*q**2)
def rho_of(q): return a_of(q)**(2.0/(G-1.0))
def J_of(q):
    a = a_of(q)
    return 1.0/a + 1.0/(3.0*a**3) + 1.0/(5.0*a**5) - 0.5*np.log((1.0+a)/(1.0-a))

def xy(q, k, upper=True):
    """(q, k) -> 物理坐标。k 是流线参数，q <= k。"""
    r = rho_of(q)
    x = (0.5/r)*(2.0/k**2 - 1.0/q**2) - 0.5*J_of(q)
    s = np.sqrt(np.clip(1.0 - (q/k)**2, 0.0, None))
    y = s/(k*r*q)
    return x, (y if upper else -y)

def velocity(q, k, upper=True):
    """速度矢量：|u| = q，方向沿流线。sinθ = q/k。"""
    #: ★ 实测：k=const 流线的切向角精确等于 **-asin(q/k)**（上支）。
    #: 第一版写成 +asin(q/k)，PDE oracle 立刻以 |∇×u| ~ O(1) 判它错。
    sin_t = -q/k
    cos_t = np.sqrt(np.clip(1.0 - (q/k)**2, 0.0, None))
    #: ★★★ 下支必须**整体取负**：沿一条流线从上支加速到鼻子(q=k)，再沿下支**减速**回出口
    #: ⇒ 下支的流向是 q 减小。判别条件是**鼻子处的速度连续性**（现约定跳 1.9~2.1，
    #: 取负后跳 2.7e-05），而**不是** PDE —— ∇×u=0 与 ∇·(ρu)=0 都是齐次条件，
    #: 对整体变号不变，所以 PDE oracle 对这个符号是瞎的。
    return (q*cos_t, q*sin_t) if upper else (-q*cos_t, q*sin_t)


K_IN, K_OUT, Q0 = 1.1, 0.9, 0.45      # 内壁 k=1.1（鼻子 M=1.264）/ 外壁 k=0.9（M=0.983）

def wall_points(k, n):
    """一条 k=const 流线的全支：上支 q:Q0->k、鼻子、下支 q:k->Q0。"""
    q = Q0 + (k - Q0)*(1.0 - np.cos(np.linspace(0, np.pi/2, n)))   # 向鼻子聚束
    up = np.array([xy(qq, k, True) for qq in q[:-1]])
    nose = np.array([[xy(k*(1-1e-14), k, True)[0], 0.0]])
    lo = np.array([xy(qq, k, False) for qq in q[:-1][::-1]])
    return np.vstack([up, nose, lo])

def cap_points(k_a, k_b, upper, n):
    """等速线 q=Q0 上从 k_a 到 k_b（进/出口）。"""
    ks = np.linspace(k_a, k_b, n)
    return np.array([xy(Q0, kk, upper) for kk in ks])

def _line_int(f, vel, t0, t1, n, cluster_at_t1=False):
    """★ `cluster_at_t1`：把求积点向 t1 端聚束。沿 k=const 流线积到鼻子(q=k)时被积函数
    变刚性（cosθ→0，曲线掉头），均匀点在下支给出 ~4e-05 的 ∇φ=u 误差 —— 那是 Dirichlet
    数据的精度地板，必须压掉。"""
    if cluster_at_t1:
        u = np.linspace(0.0, 1.0, n)
        ts = t1 - (t1 - t0)*(1.0 - u)**2
    else:
        ts = np.linspace(t0, t1, n)
    P = np.array([f(t) for t in ts]); U = np.array([vel(t) for t in ts])
    d = np.diff(P, axis=0); m = 0.5*(U[:-1] + U[1:])
    return float(np.sum(m[:, 0]*d[:, 0] + m[:, 1]*d[:, 1]))


def phi_exact(q, k, upper, n=4001):
    """势，**沿流向连续构造**：参考点取 (q=Q0, k=K_OUT) 的上支进口。

    上支：沿 q=Q0 从 K_OUT 走到 k，再沿 k 从 Q0 走到 q。
    下支：先走到该 k 的鼻子（上支，q->k），再沿下支从鼻子降到 q。
    ★ 下支不能独立以 (Q0, K_OUT, lower) 起算 —— 那是另一个物理点，φ 会在 y=0 断开。
    """
    eps = 1e-12
    a = _line_int(lambda t: np.array(xy(Q0, t, True)),
                  lambda t: np.array(velocity(Q0, t, True)), K_OUT, k, n)
    if upper:
        b = _line_int(lambda t: np.array(xy(t, k, True)),
                      lambda t: np.array(velocity(t, k, True)), Q0, q, n)
        return a + b
    b = _line_int(lambda t: np.array(xy(t, k, True)),
                  lambda t: np.array(velocity(t, k, True)), Q0, k*(1-eps), n,
                  cluster_at_t1=True)
    c = -_line_int(lambda t: np.array(xy(t, k, False)),
                   lambda t: np.array(velocity(t, k, False)), q, k*(1-eps), n,
                   cluster_at_t1=True)
    return a + b + c


class PhiTable:
    """φ 的预积分表：φ(q,k,upper) = A(k) + B(q,k)，两段都用累积梯形。

    ★ 直接对每个节点做 4001 点求积要 3000 万次求值（实测一次 h=0.06 的准备就 ~2 分钟）。
    本表把它降到一次性 O(nq·nk)，插值误差**实测**（见 C08 门里的断言）。
    """
    def __init__(self, nq=601, nk=201, upper=True):
        self.upper = upper
        self.ks = np.linspace(K_OUT, K_IN, nk)
        # A(k)：沿 q=Q0 从 K_OUT 到 k
        P = np.array([xy(Q0, kk, upper) for kk in self.ks])
        U = np.array([velocity(Q0, kk, upper) for kk in self.ks])
        d = np.diff(P, axis=0); m = 0.5*(U[:-1] + U[1:])
        self.A = np.concatenate([[0.0], np.cumsum(m[:,0]*d[:,0] + m[:,1]*d[:,1])])
        # B(q,k)：每条 k 上沿 q 从 Q0 累积。q 上界随 k 变 ⇒ 用归一化参数 t = (q-Q0)/(k-Q0)
        self.ts = np.linspace(0.0, 1.0, nq)
        self.B = np.empty((nk, nq))
        for i, kk in enumerate(self.ks):
            q = Q0 + self.ts*(kk - Q0)
            q[-1] = kk*(1 - 1e-12)
            P = np.array([xy(qq, kk, upper) for qq in q])
            U = np.array([velocity(qq, kk, upper) for qq in q])
            d = np.diff(P, axis=0); m = 0.5*(U[:-1] + U[1:])
            self.B[i] = np.concatenate([[0.0], np.cumsum(m[:,0]*d[:,0] + m[:,1]*d[:,1])])

    def __call__(self, q, k):
        t = np.clip((q - Q0)/(k - Q0), 0.0, 1.0)
        i = np.clip(np.searchsorted(self.ks, k) - 1, 0, len(self.ks)-2)
        wk = (k - self.ks[i])/(self.ks[i+1] - self.ks[i])
        j = np.clip(np.searchsorted(self.ts, t) - 1, 0, len(self.ts)-2)
        wt = (t - self.ts[j])/(self.ts[j+1] - self.ts[j])
        B = ((1-wk)*((1-wt)*self.B[i, j] + wt*self.B[i, j+1])
             + wk*((1-wt)*self.B[i+1, j] + wt*self.B[i+1, j+1]))
        A = (1-wk)*self.A[i] + wk*self.A[i+1]
        return A + B


def invert(x, y, q_lo=0.2, q_hi=1.45):
    """(x,y) -> (q,k)。由 x 方程解出 s² = ρq²(x+J/2)+1/2（s=q/k），代入 y 方程二分。"""
    ay = abs(y)
    def resid(q):
        r = rho_of(q); s2 = min(max(r*q*q*(x + 0.5*J_of(q)) + 0.5, 0.0), 1.0)
        return s2*(1.0 - s2)/(q**4 * r*r) - ay*ay
    lo, hi = q_lo, q_hi; flo = resid(lo)
    for _ in range(200):
        mid = 0.5*(lo + hi); fm = resid(mid)
        if flo*fm <= 0.0: hi = mid
        else: lo, flo = mid, fm
        if hi - lo < 1e-14: break
    q = 0.5*(lo + hi); r = rho_of(q)
    s = np.sqrt(min(max(r*q*q*(x + 0.5*J_of(q)) + 0.5, 0.0), 1.0))
    return q, (q/s if s > 0 else np.inf)


def ringleb_mesh_2d(h, n_wall=200, n_cap=40):
    """gmsh 平面网格：两条流线为壁，两条等速线为进出口。"""
    import gmsh, sys as _s
    _s.path.insert(0, "/home/lrz/codes/UP3D")
    from pyfp3d.meshgen.planar import _collect_2d
    wi = wall_points(K_IN, n_wall); wo = wall_points(K_OUT, n_wall)
    inl = cap_points(K_OUT, K_IN, True, n_cap); out = cap_points(K_IN, K_OUT, False, n_cap)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("ringleb"); geo = gmsh.model.geo
        def spline(P, first=None, last=None):
            pts = []
            for i, (x, y) in enumerate(P):
                if i == 0 and first is not None: pts.append(first); continue
                if i == len(P)-1 and last is not None: pts.append(last); continue
                pts.append(geo.addPoint(float(x), float(y), 0.0))
            return geo.addSpline(pts), pts[0], pts[-1]
        c_in, p0, p1 = spline(inl)
        c_wi, _, p2 = spline(wi, first=p1)
        c_out, _, p3 = spline(out, first=p2)
        c_wo, _, _ = spline(wo[::-1], first=p3, last=p0)
        geo.addPlaneSurface([geo.addCurveLoop([c_in, c_wi, c_out, c_wo])])
        geo.synchronize()
        for o, v in (("Mesh.MeshSizeMin", h), ("Mesh.MeshSizeMax", h),
                     ("Mesh.MeshSizeExtendFromBoundary", 0),
                     ("Mesh.MeshSizeFromPoints", 0), ("Mesh.MeshSizeFromCurvature", 0)):
            gmsh.option.setNumber(o, v)
        gmsh.model.mesh.generate(2)
        return _collect_2d({"wall": [c_wi, c_wo], "inlet": [c_in], "outlet": [c_out]})
    finally:
        gmsh.finalize()


#: 参考态：Ringleb 速度 0.7 ⇒ M∞ = 0.737046（两套密度公式实测一致到 3.7e-16）
QS_INF = 0.7
M_INF = float(np.sqrt(QS_INF**2 / (1.0 - 0.5*(G - 1.0)*QS_INF**2)))


def build_case(h):
    """网格 + 每节点的精确 (q, k) + 精确 φ（项目归一化）。"""
    import sys as _s
    _s.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pyfp3d.meshgen.extrude import extrude_single_layer
    p2, tri, eg = ringleb_mesh_2d(h)
    mesh = extrude_single_layer(p2, tri, eg, None, dz=h, name="ringleb")
    N = mesh.nodes
    qk = np.array([invert(x, y) for x, y in N[:, :2]])
    up = N[:, 1] >= 0
    phi = np.array([phi_exact(q, k, bool(u), 801)
                    for q, k, u in zip(qk[:, 0], qk[:, 1], up)]) / QS_INF
    ctr = N[mesh.elements].mean(axis=1)
    qkc = np.array([invert(x, y) for x, y in ctr[:, :2]])
    return dict(mesh=mesh, phi_exact=phi, q_node=qk[:, 0], k_node=qk[:, 1],
                q_cell=qkc[:, 0], k_cell=qkc[:, 1], ctr=ctr,
                qhat_cell=qkc[:, 0]/QS_INF,
                mach_cell=qkc[:, 0]/a_of(qkc[:, 0]))
