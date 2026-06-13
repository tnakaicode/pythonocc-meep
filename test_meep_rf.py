import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ==========================================
# 1. 電圧信号（時間波形 V(t)）の定義
# ==========================================
frequency = 0.15          # ソースの基本周波数
voltage_amplitude = 1.0   # 電圧振幅
waveform_type = "sine"   # 任意波形: "sine", "square", "triangle", "gaussian_pulse", またはユーザー定義

# 任意波形をここで定義します。必要に応じてこの関数を編集してください。
def voltage_waveform(t):
    if waveform_type == "sine":
        return voltage_amplitude * np.sin(2 * np.pi * frequency * t)
    elif waveform_type == "square":
        return voltage_amplitude * np.sign(np.sin(2 * np.pi * frequency * t))
    elif waveform_type == "triangle":
        return voltage_amplitude * (2 / np.pi) * np.arcsin(np.sin(2 * np.pi * frequency * t))
    elif waveform_type == "gaussian_pulse":
        return voltage_amplitude * np.exp(-((t - 50) / 10) ** 2) * np.sin(2 * np.pi * frequency * t)
    else:
        # カスタム波形をここに直接定義してください。
        return voltage_amplitude * np.sin(2 * np.pi * frequency * t)

# ==========================================
# 2. 計算領域と解像度
# ==========================================
resolution = 25
cell_x = 16.0
cell_y = 8.0
cell = mp.Vector3(cell_x, cell_y, 0)

# ==========================================
# 3. 任意の電極形状の定義（壁から少し離す）
# ==========================================
geometry = []

# ① 上部電極
# 反時計回りに頂点を指定することで幾何形状が正しく認識されます。
top_vertices = [
    mp.Vector3(-1.5,  2.5, 0),
    mp.Vector3( 1.5,  2.5, 0),
    mp.Vector3( 0.3,  0.8, 0),
    mp.Vector3(-0.3,  0.8, 0)
]
geometry.append(mp.Prism(vertices=top_vertices, height=mp.inf, material=mp.metal))

# ② 下部電極
bottom_vertices = [
    mp.Vector3(-1.5, -2.5, 0),
    mp.Vector3( 1.5, -2.5, 0),
    mp.Vector3( 0.3, -0.8, 0),
    mp.Vector3(-0.3, -0.8, 0)
]
geometry.append(mp.Prism(vertices=bottom_vertices, height=mp.inf, material=mp.metal))


# ==========================================
# 4. 給電点（Feed Port）の配置（隙間に配置）
# ==========================================
# 上部電極と下部電極の隙間にカスタム波形ソースを置きます。
# `voltage_waveform` の内容を変更すれば任意波形を入力できます。
src = mp.Source(
    mp.CustomSource(src_func=voltage_waveform, center_frequency=frequency, fwidth=0.01),
    component=mp.Ey,
    center=mp.Vector3(0, 2.75, 0),
    size=mp.Vector3(1.0, 0.4, 0)
)

# 全周をPMLで囲み、不要な反射を低減します。
pml_layers = [mp.PML(1.0)]


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
sim.run(until=200)

# ==========================================
# 6. 電界（Ex, Ey）データの抽出と絶対値の計算
# ==========================================
ex_data = sim.get_array(component=mp.Ex)
ey_data = sim.get_array(component=mp.Ey)
e_mag = np.sqrt(ex_data**2 + ey_data**2)

# 座標グリッドを再構築
nx, ny = ex_data.shape
x_coords = np.linspace(-cell_x / 2, cell_x / 2, nx)
y_coords = np.linspace(-cell_y / 2, cell_y / 2, ny)
X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')

# ===============================================================================
# 7. 結果の表示
# ===============================================================================
plt.figure(figsize=(10, 5))
contour = plt.contourf(X, Y, e_mag, levels=40, cmap='jet')
plt.colorbar(contour, label="Electric Field Magnitude |E|")
plt.title("Electric Field Magnitude (|E|) with Metal Electrodes")
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-cell_x / 2, cell_x / 2)
plt.ylim(-cell_y / 2, cell_y / 2)

plt.gca().add_patch(patches.Polygon(
    [(v.x, v.y) for v in top_vertices], closed=True,
    facecolor='none', edgecolor='white', linewidth=1.5
))
plt.gca().add_patch(patches.Polygon(
    [(v.x, v.y) for v in bottom_vertices], closed=True,
    facecolor='none', edgecolor='white', linewidth=1.5
))

plt.savefig("test_meep_rf.png", dpi=200)
plt.close()
print("Saved: test_meep_rf.png")
