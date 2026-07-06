"""
ONERA M6 半翼三维非结构网格生成（无边界层）
几何数据来源: NASA GRC NPARC Alliance Validation Archive
  https://www.grc.nasa.gov/www/wind/valid/m6wing/m6wing.html
翼型坐标: foilmod.txt (零后缘厚度修正版)
"""
import gmsh
import math

# =====================================================================
# 1. ONERA-D 翼型坐标 (foilmod.txt, 零后缘厚度)
# =====================================================================
FOILMOD_UPPER = [
    (0.0000000,0.0000000),(0.0000165,0.0006914),(0.0000696,0.0014416),
    (0.0001675,0.0022554),(0.0003232,0.0031382),(0.0005508,0.0040959),
    (0.0008657,0.0051343),(0.0012868,0.0062598),(0.0018364,0.0074784),
    (0.0025441,0.0087958),(0.0034428,0.0102163),(0.0045704,0.0117419),
    (0.0059751,0.0133708),(0.0077112,0.0150951),(0.0098413,0.0168984),
    (0.0124479,0.0187537),(0.0156171,0.0206220),(0.0194609,0.0224545),
    (0.0241067,0.0242004),(0.0297008,0.0258245),(0.0364261,0.0273317),
    (0.0444852,0.0287912),(0.0541248,0.0303278),(0.0656303,0.0320138),
    (0.0793366,0.0338372),(0.0956354,0.0357742),(0.1149796,0.0377923),
    (0.1378963,0.0398522),(0.1649976,0.0419089),(0.1919327,0.0436214),
    (0.2187096,0.0450507),(0.2453310,0.0462358),(0.2717978,0.0471987),
    (0.2981113,0.0479494),(0.3242726,0.0484902),(0.3502830,0.0488183),
    (0.3761446,0.0489296),(0.4018567,0.0488202),(0.4274223,0.0484833),
    (0.4528441,0.0479351),(0.4781197,0.0471661),(0.5032514,0.0461903),
    (0.5282426,0.0450209),(0.5530937,0.0436741),(0.5778043,0.0421684),
    (0.6023757,0.0405241),(0.6268104,0.0387613),(0.6511093,0.0368990),
    (0.6752726,0.0349542),(0.6993027,0.0329402),(0.7231995,0.0308662),
    (0.7469658,0.0287365),(0.7705998,0.0265505),(0.7941055,0.0243027),
    (0.8174828,0.0219842),(0.8407324,0.0195838),(0.8638564,0.0170915),
    (0.8868235,0.0145051),(0.9061905,0.0121952),(0.9225336,0.0101138),
    (0.9363346,0.0083265),(0.9479946,0.0068038),(0.9578511,0.0055144),
    (0.9661860,0.0044240),(0.9732361,0.0035015),(0.9792020,0.0027211),
    (0.9842508,0.0020606),(0.9885252,0.0015014),(0.9921438,0.0010280),
    (0.9952080,0.0006271),(0.9978030,0.0002876),(1.0000000,0.0000000),
]

# 闭合翼型: 上表面(LE->TE) + 下表面(TE->LE, z取负, 去掉首尾避免重复)
AIRFOIL_CLOSED = list(FOILMOD_UPPER) + [
    (x, -z) for x, z in reversed(FOILMOD_UPPER[1:-1])
]

# =====================================================================
# 2. M6 翼平面参数 (官网表2)
# =====================================================================
B_SEMI   = 1.1963
TAPER    = 0.562
MAC      = 0.64607
LE_SWEEP = math.radians(30.0)
C_ROOT   = MAC / (2.0/3.0 * (1 + TAPER + TAPER**2) / (1 + TAPER))
C_TIP    = TAPER * C_ROOT

ETAS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]

def chord_at(eta):
    return C_ROOT * (1.0 - (1.0 - TAPER) * eta)

def xle_at(eta):
    return eta * B_SEMI * math.tan(LE_SWEEP)

# =====================================================================
# 3. gmsh 初始化
# =====================================================================
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("onera_m6")

# =====================================================================
# 4. 构造各展向站位的闭合截面
# =====================================================================
face_tags = []
n_upper = len(FOILMOD_UPPER)  # 72

for eta in ETAS:
    y      = eta * B_SEMI
    chord  = chord_at(eta)
    xle    = xle_at(eta)

    # 创建所有点
    pts = []
    for x_n, z_n in AIRFOIL_CLOSED:
        x = xle + x_n * chord
        z = z_n * chord
        pts.append(gmsh.model.occ.addPoint(x, y, z))

    # --- 修正: 确保两条样条共享端点 ---
    # 上样条: LE(pts[0]) → TE(pts[n_upper-1])
    upper_pts = pts[:n_upper]

    # 下样条: TE(pts[n_upper-1]) → 下表面各点 → LE(pts[0])
    # 关键: 首部补 TE 点, 尾部补 LE 点, 与上样条端点完全一致
    lower_pts = [pts[n_upper - 1]] + pts[n_upper:] + [pts[0]]

    spl_up = gmsh.model.occ.addSpline(upper_pts)
    spl_lo = gmsh.model.occ.addSpline(lower_pts)
    cl = gmsh.model.occ.addCurveLoop([spl_up, spl_lo])
    face_tags.append(gmsh.model.occ.addPlaneSurface([cl]))

gmsh.model.occ.synchronize()

# =====================================================================
# 5. 放样生成翼面实体
# =====================================================================
thru = gmsh.model.occ.addThruSections(face_tags, makeSolid=True, makeRuled=False)
gmsh.model.occ.synchronize()
wing_vols = [v[1] for v in gmsh.model.getEntities(3)]
print(f"翼面实体: {wing_vols}")

# =====================================================================
# 6. 远场流域 + 布尔减
# =====================================================================
Lx = 10.0 * C_ROOT
Ly = B_SEMI + 5.0 * C_ROOT
Lz = 10.0 * C_ROOT
x0 = -5.0 * C_ROOT
y0 = 0.0
z0 = -5.0 * C_ROOT

box = gmsh.model.occ.addBox(x0, y0, z0, Lx, Ly, Lz)
cut_dim_tags = [(3, v) for v in wing_vols]
gmsh.model.occ.cut([(3, box)], cut_dim_tags, removeObject=True, removeTool=False)
gmsh.model.occ.synchronize()

fluid_vols = [v[1] for v in gmsh.model.getEntities(3)]
print(f"流体域实体: {fluid_vols}")

# =====================================================================
# 7. 网格尺寸控制
# =====================================================================
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.005)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.8)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)

# 翼面附近加密
all_surfs = gmsh.model.getEntities(2)
wing_surf_tags = []
for dim, tag in all_surfs:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    if xmin >= -0.01 and ymin >= -0.01 and ymax <= B_SEMI + 0.01 \
       and abs(zmin) < 0.1 and abs(zmax) < 0.1:
        wing_surf_tags.append(tag)

if wing_surf_tags:
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", wing_surf_tags)
    gmsh.model.mesh.field.setNumber(dist, "Sampling", 100)

    thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(thr, "InField", dist)
    gmsh.model.mesh.field.setNumber(thr, "SizeMin", 0.008)
    gmsh.model.mesh.field.setNumber(thr, "SizeMax", 0.4)
    gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.05)
    gmsh.model.mesh.field.setNumber(thr, "DistMax", 1.5)
    gmsh.model.mesh.field.setAsBackgroundMesh(thr)

# =====================================================================
# 8. 非结构网格算法
# =====================================================================
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.Smoothing", 10)
gmsh.option.setNumber("Mesh.Format", 1)

# =====================================================================
# 9. 生成 + 输出
# =====================================================================
gmsh.model.mesh.generate(3)
gmsh.write("m6_wing_unstructured.msh")

nodes = gmsh.model.mesh.getNodes()[0]
print(f"节点数: {len(nodes)}")
tet = gmsh.model.mesh.getElementsByType(4)
print(f"四面体数: {len(tet[0])}")

# gmsh.fltk.run()
gmsh.finalize()
