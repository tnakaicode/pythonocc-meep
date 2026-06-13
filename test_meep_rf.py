import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. 電圧信号（時間波形 V(t)）の定義
# ==========================================
# Pythonの関数として、任意の入力電圧のタイムチャートを定義できます。
def custom_voltage_signal(t):
    # 例：ゆっくり立ち上がって消える、高周波のガウシアンパルス信号
    # （実験データのCSVなどから補間して波形を作ることも可能です）
    freq = 0.15          # 高周波の周波数
    pulse_center = 40.0  # ピーク時間
    pulse_width = 12.0   # パルス幅
    
    # 交流成分 × 時間エンベロープ
    return np.sin(2 * np.pi * freq * t) * np.exp(-((t - pulse_center) / pulse_width)**2)


# ==========================================
# 2. 計算領域と解像度
# ==========================================
resolution = 25
cell_x = 16.0
cell_y = 6.0
cell = mp.Vector3(cell_x, cell_y, 0)


# ==========================================
# 3. 任意の電極形状（多角形）の定義
# ==========================================
# mp.Prism を使うと、頂点（vertices）を結んで好きな形状の金属を作れます。
# ここでは「先端に電界が集中するテーパー形状」を作ります。

geometry = []

# ① 上部電極（導波管の上壁 Y=3.0 から中央へ突き出る）
top_vertices = [
    mp.Vector3(-1.5,  3.0, 0),  # 根元 左
    mp.Vector3(-0.3,  0.8, 0),  # 先端 左（細くなっている）
    mp.Vector3( 0.3,  0.8, 0),  # 先端 右
    mp.Vector3( 1.5,  3.0, 0)   # 根元 右
]
geometry.append(mp.Prism(vertices=top_vertices, height=mp.inf, material=mp.metal))

# ② 下部電極（下壁 Y=-3.0 から突き出る対向電極）
bottom_vertices = [
    mp.Vector3(-1.5, -3.0, 0),  # 根元 左
    mp.Vector3( 1.5, -3.0, 0),  # 根元 右
    mp.Vector3( 0.3, -0.8, 0),  # 先端 右
    mp.Vector3(-0.3, -0.8, 0)   # 先端 左
]
geometry.append(mp.Prism(vertices=bottom_vertices, height=mp.inf, material=mp.metal))


# ==========================================
# 4. 給電点（Feed Port）の配置
# ==========================================
# 【最重要】光源を隙間ではなく、上部電極の「根元の壁際」に小さく配置します。
# これにより、高周波電流がまず金属（mp.metal）の表面を伝わり、
# 電極形状に応じて自然に歪みながら、導波管の空間へ放射されます（プローブ給電の模倣）。

src = mp.Source(
    mp.CustomSource(src_func=custom_voltage_signal), # 上で定義したV(t)関数を指定
    component=mp.Ey,
    center=mp.Vector3(0, 2.9, 0),   # 電極の根元（上壁のすぐ下）
    size=mp.Vector3(0.4, 0.1, 0)    # 給電ピンのような小さなサイズ
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

# 電界（Ey）の取得
ey_data = sim.get_array(component=mp.Ey)

# プロット
plt.figure(figsize=(10, 5))
sim.plot2D() # 幾何構造の描画（テーパー形状が確認できます）

# 電界分布を重ね書き（エッジ付近の強い電界が見えるよう vmin/vmax を少し絞ります）
plt.imshow(ey_data.T, cmap='RdBu', vmin=-0.2, vmax=0.2,
           extent=[-cell_x/2, cell_x/2, -cell_y/2, cell_y/2], origin='lower', alpha=0.6)
plt.colorbar(label="Electric Field (Ey)")
plt.title("Realistic RF Feed with Custom Electrode Shape & Signal")
plt.xlabel("X")
plt.ylabel("Y")
plt.savefig("test_meep_rf.png")
