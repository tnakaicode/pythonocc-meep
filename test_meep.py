import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. パラメータの設定
# ==========================================
resolution = 20  # メッシュ細かさ (ピクセル/単位長さ)
cell_x = 16.0    # シミュレーション領域のX幅
cell_y = 8.0     # シミュレーション領域のY幅
cell = mp.Vector3(cell_x, cell_y, 0)

w_waveguide = 1.0   # 導波路の幅
w_electrode = 2.0   # 電極の横幅
h_electrode = 1.5   # 電極の縦幅
gap = 0.3           # 導波路と電極の隙間

# ==========================================
# 2. 幾何構造（導波路と電極）の定義
# ==========================================
geometry = []

# ① 誘電体導波路 (X方向に無限に伸びるコア)
waveguide = mp.Block(
    center=mp.Vector3(0, 0, 0),
    size=mp.Vector3(mp.inf, w_waveguide, mp.inf),
    material=mp.Medium(index=3.45)  # 例：シリコンなど
)
geometry.append(waveguide)

# ② 上部電極 (金属: mp.metal を指定)
electrode_top = mp.Block(
    center=mp.Vector3(0, (w_waveguide / 2) + gap + (h_electrode / 2), 0),
    size=mp.Vector3(w_electrode, h_electrode, mp.inf),
    material=mp.metal  # 完全導体(PEC)として振る舞います
)
geometry.append(electrode_top)

# ③ 下部電極
electrode_bottom = mp.Block(
    center=mp.Vector3(0, -((w_waveguide / 2) + gap + (h_electrode / 2)), 0),
    size=mp.Vector3(w_electrode, h_electrode, mp.inf),
    material=mp.metal
)
geometry.append(electrode_bottom)

# ==========================================
# 3. 光源（Source）と境界条件（PML）
# ==========================================
# 左端から連続波（Continuous Source）を入力
fcen = 0.15  # 中心周波数
src = mp.Source(
    mp.ContinuousSource(frequency=fcen),
    component=mp.Ez,
    center=mp.Vector3(-cell_x / 2 + 1.0, 0, 0),
    size=mp.Vector3(0, w_waveguide * 2, 0)
)

# 周囲を吸収境界（PML）で囲む
pml_layers = [mp.PML(1.0)]

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

# 指定した時間（ここでは200ステップ分）走らせる
sim.run(until=200)

# ==========================================
# 5. 電界データの抽出と可視化
# ==========================================
# Ez電界の2次元配列を取得
ez_data = sim.get_ez_array()

# プロット
plt.figure(figsize=(10, 5))
# 構造のレイアウトを背景に重ねる
sim.plot2D()
# 電界分布をカラーマップで表示 (シミュレーション終了時点の瞬間値)
plt.imshow(ez_data.T, cmap='RdBu', vspace=np.min(ez_data), vmax=np.max(ez_data),
           extent=[-cell_x/2, cell_x/2, -cell_y/2, cell_y/2], origin='lower', alpha=0.6)
plt.colorbar(label="Electric Field (Ez)")
plt.title("E-Field Distribution with Electrodes")
plt.show()