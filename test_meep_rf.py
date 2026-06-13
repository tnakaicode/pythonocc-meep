import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

tempdir = "/mnt/win_temp/"

# ==========================================
# 1. 電圧信号（時間波形 V(t)）の定義
# ==========================================
frequency = 0.15          # ソースの基本周波数
voltage_amplitude = 1.0   # 電圧振幅
waveform_type = "sine"   # 任意波形: "sine", "square", "triangle", "gaussian_pulse", "custom"
ramp_time = 20.0         # 立ち上がり時間（過渡応答抑制用）

# 任意波形をここで定義します。必要に応じてこの関数を編集してください。
def voltage_waveform(t):
    if t < 0:
        return 0.0
    envelope = 1.0 - np.exp(-t / ramp_time)
    if waveform_type == "sine":
        value = voltage_amplitude * np.sin(2 * np.pi * frequency * t)
    elif waveform_type == "square":
        value = voltage_amplitude * np.sign(np.sin(2 * np.pi * frequency * t))
    elif waveform_type == "triangle":
        value = voltage_amplitude * (2 / np.pi) * np.arcsin(np.sin(2 * np.pi * frequency * t))
    elif waveform_type == "gaussian_pulse":
        value = voltage_amplitude * np.exp(-((t - 50) / 10) ** 2) * np.sin(2 * np.pi * frequency * t)
    elif waveform_type == "custom":
        # ここに任意のユーザー波形を実装してください。
        value = voltage_amplitude * (np.sin(2 * np.pi * frequency * t) + 0.5 * np.sin(4 * np.pi * frequency * t))
    else:
        value = voltage_amplitude * np.sin(2 * np.pi * frequency * t)
    return envelope * value

# ==========================================
# 2. 計算領域と解像度
# ==========================================
resolution = 30
cell_x = 20.0
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
# 4. 導波管壁としての上下金属ブロック追加
# ==========================================
wall_thickness = 0.5
geometry.append(mp.Block(
    center=mp.Vector3(0, cell_y / 2 - wall_thickness / 2, 0),
    size=mp.Vector3(cell_x, wall_thickness, mp.inf),
    material=mp.metal,
))
geometry.append(mp.Block(
    center=mp.Vector3(0, -cell_y / 2 + wall_thickness / 2, 0),
    size=mp.Vector3(cell_x, wall_thickness, mp.inf),
    material=mp.metal,
))

# ==========================================
# 5. 給電点（Feed Port）の配置（隙間に配置）
# ==========================================
# 上部電極と下部電極の隙間にカスタム波形ソースを置きます。
# `voltage_waveform` の内容を変更すれば任意波形を入力できます。
# 本来の給電点は上部と下部電極の間のギャップ領域です。
src = mp.Source(
    mp.CustomSource(src_func=voltage_waveform, center_frequency=frequency, fwidth=0.01),
    component=mp.Ey,
    center=mp.Vector3(0, 0.0, 0),
    size=mp.Vector3(3.0, 0.4, 0)
)

# 左右のみPMLにして、上下は金属導波管壁とします。
pml_layers = [mp.PML(2.0, direction=mp.X)]


# ==========================================
# 6. シミュレーション実行と可視化（PVD・VTK一括出力）
# ==========================================

step_counter = 0
vtk_initialized = False
pec_flag = None
nx, ny = 0, 0
dx, dy = 0, 0

# 💡 時間とファイル名のペアを記録するリスト
time_history = []

def output_vtk_step(sim):
    global step_counter, vtk_initialized, pec_flag, nx, ny, dx, dy, time_history
    
    ex_data = sim.get_array(center=mp.Vector3(), size=sim.cell_size, component=mp.Ex)
    ey_data = sim.get_array(center=mp.Vector3(), size=sim.cell_size, component=mp.Ey)
    e_mag = np.sqrt(ex_data**2 + ey_data**2)
    
    if not vtk_initialized:
        nx, ny = e_mag.shape
        x_coords = np.linspace(-cell_x / 2, cell_x / 2, nx)
        y_coords = np.linspace(-cell_y / 2, cell_y / 2, ny)
        X_local, Y_local = np.meshgrid(x_coords, y_coords, indexing='ij')
        
        pec_flag = np.zeros_like(e_mag)
        wall_mask = (Y_local >= cell_y / 2 - wall_thickness) | (Y_local <= -cell_y / 2 + wall_thickness)
        pq_top = Path([(v.x, v.y) for v in top_vertices])
        pq_bottom = Path([(v.x, v.y) for v in bottom_vertices])
        points = np.vstack((X_local.ravel(), Y_local.ravel())).T
        mask_top = pq_top.contains_points(points).reshape(X_local.shape)
        mask_bottom = pq_bottom.contains_points(points).reshape(X_local.shape)
        pec_flag[mask_top | mask_bottom | wall_mask] = 1.0
        
        dx = cell_x / (nx - 1) if nx > 1 else 1.0
        dy = cell_y / (ny - 1) if ny > 1 else 1.0
        vtk_initialized = True

    current_time = sim.meep_time()
    metadata = f"time:{current_time:.2f} wave:{waveform_type} freq:{frequency} res:{resolution}"
    metadata = metadata[:255]
    
    # 💡 後でPVDファイルに紐付けるため、ファイル名単体と時間を記録
    filename_local = f"test_meep_rf_{step_counter:06d}.vtk"
    vtk_path = f"{tempdir}{filename_local}"
    time_history.append((current_time, filename_local))
    
    with open(vtk_path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{metadata}\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write("DIMENSIONS {} {} 1\n".format(nx, ny))
        f.write("ORIGIN {} {} 0.0\n".format(-cell_x / 2, -cell_y / 2))
        f.write("SPACING {} {} 1.0\n".format(dx, dy))
        f.write("POINT_DATA {}\n".format(nx * ny))
        
        f.write("SCALARS E_magnitude float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(ny):
            for i in range(nx):
                f.write(f"{float(e_mag[i, j]):.6f}\n")
                
        f.write("SCALARS Ex float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(ny):
            for i in range(nx):
                f.write(f"{float(ex_data[i, j]):.6f}\n")
                
        f.write("SCALARS Ey float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(ny):
            for i in range(nx):
                f.write(f"{float(ey_data[i, j]):.6f}\n")
                
        f.write("SCALARS Geometry_PEC float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(ny):
            for i in range(nx):
                f.write(f"{float(pec_flag[i, j]):.1f}\n")
                
    print(f"Saved: {vtk_path} (t={current_time:.2f})")
    step_counter += 1


# --- シミュレーション実行 ---
sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=[src],
    resolution=resolution
)

sim.run(
    mp.at_every(2.0, output_vtk_step),
    until=200
)

# ==========================================
# 💡 【追加】すべてのVTKを1つに束ねる PVD ファイルの作成
# ==========================================
pvd_path = f"{tempdir}test_meep_rf.pvd"
with open(pvd_path, "w") as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
    f.write('  <Collection>\n')
    for t, fname in time_history:
        f.write(f'    <DataSet timestep="{t}" group="" part="0" file="{fname}"/>\n')
    f.write('  </Collection>\n')
    f.write('<</VTKFile>\n')

print(f"\n✨ [完了] 全データを束ねた親玉ファイルを作成しました: {pvd_path}")
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

# 金属領域をマスクして描画すると、導体内部の非物理的な数値表示を避けられます。
polygon_top = Path([(v.x, v.y) for v in top_vertices])
polygon_bottom = Path([(v.x, v.y) for v in bottom_vertices])
points = np.vstack((X.ravel(), Y.ravel())).T
mask = polygon_top.contains_points(points) | polygon_bottom.contains_points(points)
mask = mask.reshape(X.shape)
e_mag_plot = np.ma.masked_where(mask, e_mag)

# ===============================================================================
# 7. 結果の表示
# ===============================================================================
plt.figure(figsize=(10, 5))
contour = plt.contourf(X, Y, e_mag_plot, levels=40, cmap='jet')
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

plt.savefig(tempdir + "test_meep_rf.png", dpi=200)
plt.close()
print("Saved: test_meep_rf.png")

# ==========================================
# 8. 入力波形のプロット出力
# ==========================================
time_points = np.linspace(0.0, 10.0 / frequency, 2000)
voltage_values = np.array([voltage_waveform(t) for t in time_points])
plt.figure(figsize=(8, 3))
plt.plot(time_points, voltage_values, color="tab:blue")
plt.title("Input Voltage Waveform")
plt.xlabel("Time")
plt.ylabel("Voltage")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(tempdir + "test_meep_rf_waveform.png", dpi=200)
plt.close()
print("Saved: test_meep_rf_waveform.png")

# ==========================================
# 9. 電極形状・導波管壁の図として出力
# ==========================================
plt.figure(figsize=(8, 4))
ax = plt.gca()
# 電極と壁の描画
ax.add_patch(patches.Polygon(
    [(v.x, v.y) for v in top_vertices], closed=True,
    facecolor="lightgray", edgecolor="black", linewidth=1.5, label="Top Electrode"
))
ax.add_patch(patches.Polygon(
    [(v.x, v.y) for v in bottom_vertices], closed=True,
    facecolor="lightgray", edgecolor="black", linewidth=1.5, label="Bottom Electrode"
))
ax.add_patch(patches.Rectangle(
    (-cell_x / 2, cell_y / 2 - wall_thickness), cell_x, wall_thickness,
    facecolor="silver", edgecolor="black", linewidth=1.0, label="Waveguide Wall"
))
ax.add_patch(patches.Rectangle(
    (-cell_x / 2, -cell_y / 2), cell_x, wall_thickness,
    facecolor="silver", edgecolor="black", linewidth=1.0
))
ax.add_patch(patches.Rectangle(
    (-1.5, -0.2), 3.0, 0.4,
    facecolor="none", edgecolor="red", linewidth=1.5, linestyle="--", label="Source Region"
))
ax.set_title("Electrode Geometry and Waveguide Walls")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_xlim(-cell_x / 2, cell_x / 2)
ax.set_ylim(-cell_y / 2, cell_y / 2)
ax.set_aspect("equal", adjustable="box")
ax.legend(loc="upper right")
plt.grid(True, linestyle="--", alpha=0.3)
plt.tight_layout()
plt.savefig(tempdir + "test_meep_rf_geometry.png", dpi=200)
plt.close()
print("Saved: test_meep_rf_geometry.png")

# ==========================================
# 10. VTK 出力(フィールド・幾何形状・入力条件付き)
# ==========================================
pec_flag = np.zeros_like(e_mag)
wall_mask = (Y >= cell_y / 2 - wall_thickness) | (Y <= -cell_y / 2 + wall_thickness)
pq_top = Path([(v.x, v.y) for v in top_vertices])
pq_bottom = Path([(v.x, v.y) for v in bottom_vertices])
points = np.vstack((X.ravel(), Y.ravel())).T
mask_top = pq_top.contains_points(points).reshape(X.shape)
mask_bottom = pq_bottom.contains_points(points).reshape(X.shape)
pec_flag[mask_top | mask_bottom | wall_mask] = 1.0

# 💡 メタデータは2行目の「タイトル識別子」に1行で埋め込む（最大256文字）
metadata = f"wave:{waveform_type} freq:{frequency} volt:{voltage_amplitude} ramp:{ramp_time} res:{resolution}cx:{cell_x} cy:{cell_y}"
metadata = metadata[:255] # 安全のため256文字未満にカット

vtk_path = tempdir + "test_meep_rf.vtk"
with open(vtk_path, "w") as f:
    f.write("# vtk DataFile Version 3.0\n")         # 1行目: 決まり文句
    f.write(f"{metadata}\n")                        # 2行目: タイトル（ここに情報を集約）
    f.write("ASCII\n")                              # 3行目: データ形式
    f.write("DATASET STRUCTURED_POINTS\n")          # 4行目: ⚠️ ASCIIの直後は必ずこれ！
    f.write("DIMENSIONS {} {} 1\n".format(nx, ny))
    f.write("ORIGIN {} {} 0.0\n".format(-cell_x / 2, -cell_y / 2))
    
    # 💡 POINT（格子点）のデータなので、間隔は (n - 1) で割るのがVTKの正確な仕様です
    dx = cell_x / (nx - 1) if nx > 1 else 1.0
    dy = cell_y / (ny - 1) if ny > 1 else 1.0
    f.write("SPACING {} {} 1.0\n".format(dx, dy))
    
    f.write("POINT_DATA {}\n".format(nx * ny))
    f.write("SCALARS E_magnitude float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for j in range(ny):
        for i in range(nx):
            f.write(f"{float(e_mag[i, j]):.6f}\n")
    f.write("SCALARS Ex float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for j in range(ny):
        for i in range(nx):
            f.write(f"{float(ex_data[i, j]):.6f}\n")
    f.write("SCALARS Ey float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for j in range(ny):
        for i in range(nx):
            f.write(f"{float(ey_data[i, j]):.6f}\n")
    f.write("SCALARS Geometry_PEC float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for j in range(ny):
        for i in range(nx):
            f.write(f"{float(pec_flag[i, j]):.1f}\n")
            
print(f"Saved: {vtk_path}")
