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

# ① 上部電極（根元を Y=2.5 にして、上壁 Y=3.0 との間に隙間を作ります）
top_vertices = [
    mp.Vector3(-1.5,  2.5, 0),  # 根元 左
    mp.Vector3(-0.3,  0.8, 0),  # 先端 左
    mp.Vector3( 0.3,  0.8, 0),  # 先端 右
    mp.Vector3( 1.5,  2.5, 0)   # 根元 右
]
geometry.append(mp.Prism(vertices=top_vertices, height=mp.inf, material=mp.metal))

# ② 下部電極（こちらも下壁 Y=-3.0 との間に隙間を作ります）
bottom_vertices = [
    mp.Vector3(-1.5, -2.5, 0),  # 根元 左
    mp.Vector3( 1.5, -2.5, 0),  # 根元 右
    mp.Vector3( 0.3, -0.8, 0),  # 先端 右
    mp.Vector3(-0.3, -0.8, 0)   # 先端 左
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

# プロット
plt.figure(figsize=(10, 5))
sim.plot2D()  # 幾何構造の描画（PEC内部は黒くなります）

# 電界強度の分布を重ね書き
# 【修正②】カラーマップを電界強度集中が見やすい 'magma' に変更します。
# 【修正③】絶対値なので vmin=0 に設定します。
plt.contourf(e_mag.T, cmap='magma')

plt.colorbar(label="Electric Field Magnitude (|E|)")
plt.title("Realistic RF Feed with Custom Electrode Shape & Magnitude Plot")
plt.xlabel("X")
plt.ylabel("Y")
plt.savefig("test_meep_rf.png")
