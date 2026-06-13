import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. パラメータの設定
# ==========================================
resolution = 20  # メッシュ細かさ
cell_x = 16.0    # 導波管の長さ（X方向）
cell_y = 6.0     # 導波管の幅（Y方向）
cell = mp.Vector3(cell_x, cell_y, 0)

w_electrode = 1.5   # 電極の横幅
h_electrode = 2.0   # 電極の縦幅
gap = 0.8           # 電極間のギャップ（ここに電圧がかかる）

# ==========================================
# 2. 幾何構造（上下の電極）の定義
# ==========================================
geometry = []

# 上部電極
electrode_top = mp.Block(
    center=mp.Vector3(0, (gap / 2) + (h_electrode / 2), 0),
    size=mp.Vector3(w_electrode, h_electrode, mp.inf),
    material=mp.metal
)
geometry.append(electrode_top)

# 下部電極
electrode_bottom = mp.Block(
    center=mp.Vector3(0, -((gap / 2) + (h_electrode / 2)), 0),
    size=mp.Vector3(w_electrode, h_electrode, mp.inf),
    material=mp.metal
)
geometry.append(electrode_bottom)

# ==========================================
# 3. 光源（高周波電圧の印加）と境界条件
# ==========================================
# 【ポイント①】電極の「隙間」に光源を配置。
# 上下の電極間に電圧をかけるため、電場成分は Y 方向（mp.Ey）にします。
fcen = 0.2  # 印加する高周波の周波数
src = mp.Source(
    mp.ContinuousSource(frequency=fcen),
    component=mp.Ey,                     # Y方向の電場 ＝ 上下電極間の電圧
    center=mp.Vector3(0, 0, 0),          # 電極の真ん中の隙間
    size=mp.Vector3(w_electrode, gap, 0) # 隙間のサイズ全体に給電
)

# 【ポイント②】PML（吸収壁）を「左右のみ」に設定。
# これにより、シミュレーション領域の「上下の端」は自動的に完全導体（金属壁）になり、
# 左右に無限に続く金属導波管（Waveguide）が作られます。
pml_layers = [mp.PML(1.0, direction=mp.X)]

# ==========================================
# 4. シミュレーションの実行
# ==========================================
sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=[src],
    resolution=resolution
)

# 高周波が十分に伝播するまで走らせる
sim.run(until=150)

# ==========================================
# 5. 電界（Ey）データの抽出と可視化
# ==========================================
# 今回は Ey（上下方向の電場）を抽出します
ey_data = sim.get_array(component=mp.Ey)

# プロット
plt.figure(figsize=(10, 5))
sim.plot2D()  # 構造の表示

# 電界分布を重ね合わせ
plt.imshow(ey_data.T, cmap='RdBu', vmin=-0.5, vmax=0.5,
           extent=[-cell_x/2, cell_x/2, -cell_y/2, cell_y/2], origin='lower', alpha=0.6)
plt.colorbar(label="Electric Field (Ey)")
plt.title("RF Voltage Application & Waveguide Excitation")
plt.xlabel("X")
plt.ylabel("Y")
plt.savefig("test_meep_rf.png")
