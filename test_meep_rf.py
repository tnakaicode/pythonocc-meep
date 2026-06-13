import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. 電圧信号（時間波形 V(t)）の定義
# ==========================================
def custom_voltage_signal(t):
    # 【修正】パルスではなく、ずっとパタパタと入れ替わる高周波電圧（交流）にし続けます
    freq = 0.15          # 高周波の周波数
    return np.sin(2 * np.pi * freq * t)

# ==========================================
# 2. 計算領域と解像度
# ==========================================
resolution = 25
cell_x = 16.0
cell_y = 6.0
cell = mp.Vector3(cell_x, cell_y, 0)

# ==========================================
# 3. 任意の電極形状の定義（壁から少し離す）
# ==========================================
geometry = []

# ① 上部電極（反時計回りに修正：左上 → 右上 → 右下 → 左下 の順）
top_vertices = [
    mp.Vector3(-1.5,  2.5, 0),  # ① 根元 左 (左上)
    mp.Vector3( 1.5,  2.5, 0),  # ② 根元 右 (右上)
    mp.Vector3( 0.3,  0.8, 0),  # ③ 先端 右 (右下)
    mp.Vector3(-0.3,  0.8, 0)   # ④ 先端 左 (左下)
]
geometry.append(mp.Prism(vertices=top_vertices, height=mp.inf, material=mp.metal))

# ② 下部電極（こちらは元から反時計回りなのでそのままでOK：左下 → 右下 → 右上 → 左上 の順）
bottom_vertices = [
    mp.Vector3(-1.5, -2.5, 0),  # ① 根元 左 (左下)
    mp.Vector3( 1.5, -2.5, 0),  # ② 根元 右 (右下)
    mp.Vector3( 0.3, -0.8, 0),  # ③ 先端 右 (右上)
    mp.Vector3(-0.3, -0.8, 0)   # ④ 先端 左 (左上)
]
geometry.append(mp.Prism(vertices=bottom_vertices, height=mp.inf, material=mp.metal))


# ==========================================
# 4. 給電点（Feed Port）の配置（隙間に配置）
# ==========================================
# 上壁（GND）と電極の隙間（Y=2.5 〜 3.0 の間）に電圧信号を印加します。
src = mp.Source(
    mp.CustomSource(src_func=custom_voltage_signal),
    component=mp.Ey,
    center=mp.Vector3(0, 2.75, 0),  # 隙間の中心
    size=mp.Vector3(1.0, 0.5, 0)    # 隙間を埋めるサイズ
)

# 左右のみPML（上下の端は自動的に完全導体の導波管壁になります）
pml_layers = [mp.PML(1.0, direction=mp.X)]


# ==========================================
# 5. シミュレーション実行と可視化
# ==========================================
sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=[src],
    resolution=resolution
)

# 信号が電極を伝わり、導波管内に広がっていく様子を計算
sim.run(until=100)

# ==========================================
# 5. 電界（Ex, Ey）データの抽出と絶対値の計算
# ==========================================
# 【修正①】EyだけでなくExも取得して絶対値を計算します。
# 電界ベクトルの絶対値 Magnitude |E| = sqrt(Ex^2 + Ey^2) を計算
ey_data = sim.get_array(component=mp.Ey)
ex_data = sim.get_array(component=mp.Ex)
e_mag = np.sqrt(ey_data**2 + ex_data**2)

# ==============================================================================
# 変更後：等高線プロットと電極の同時可視化
# ==============================================================================
plt.figure(figsize=(10, 4.5))

# 1. 先に等高線（contourf）を描画して空間の電界強度分布をプロット
# （※データの配列サイズに合わせた座標格子 X, Y を使用してください）
X, Y = np.meshgrid(np.linspace(-8, 8, ex_data.shape[0]), np.linspace(-3, 3, ex_data.shape[1]), indexing='ij')
contour = plt.contourf(X, Y, e_mag, levels=10, cmap='jet', vmin=0, vmax=0.64)

# 2. その上から Meep の幾何構造（電極）を重ね書き
# plot2Dを後から呼び出すことで、電極形状（白い領域）が最前面に描画されます
sim.plot2D(fill_quantities=False, intensities=False)

# 3. 必要に応じて電極の輪郭線（PECマトリクスなど）を黒線で強調する場合
# （pec_flag が導体部=1, 空間=0 のバイナリデータの場合）
plt.contour(X, Y, e_mag, levels=[0.5], colors='black', linewidths=1.5)

plt.colorbar(contour, label="Electric Field Magnitude (|E|)")
plt.title("Realistic RF Feed with Custom Electrode Shape & Contourf Plot")
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-8, 8)
plt.ylim(-3, 3)
plt.savefig("test_meep_rf.png")
