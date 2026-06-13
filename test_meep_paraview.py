import meep as mp
import numpy as np

# ==============================================================================
# 1. 入力条件・計算モデルパラメータの定義
# ==============================================================================
cell_x, cell_y = 16.0, 6.0
resolution = 25
end_time = 100.0
frequency = 0.15

cell = mp.Vector3(cell_x, cell_y, 0)
pml_layers = [mp.PML(thickness=1.0, direction=mp.X)]

# 幾何構造（完全導体: PECテーパー電極）の配置
geometry = [
    mp.Prism(vertices=[
        mp.Vector3(-1.5,  2.5, 0), mp.Vector3(-0.3,  0.8, 0),
        mp.Vector3( 0.3,  0.8, 0), mp.Vector3( 1.5,  2.5, 0)
    ], height=mp.inf, material=mp.metal),
    mp.Prism(vertices=[
        mp.Vector3(-1.5, -2.5, 0), mp.Vector3( 1.5, -2.5, 0),
        mp.Vector3( 0.3, -0.8, 0), mp.Vector3(-0.3, -0.8, 0)
    ], height=mp.inf, material=mp.metal)
]

def custom_voltage_signal(t):
    return np.sin(2 * np.pi * frequency * t)

sources = [mp.Source(
    mp.CustomSource(src_func=custom_voltage_signal),
    component=mp.Ey,
    center=mp.Vector3(0, 2.75, 0),
    size=mp.Vector3(1.0, 0.5, 0)
)]

# ==============================================================================
# 2. シミュレーションの実行
# ==============================================================================
sim = mp.Simulation(cell_size=cell, boundary_layers=pml_layers,
                    geometry=geometry, sources=sources, resolution=resolution)

print("--- FDTDシミュレーションを開始します ---")
sim.run(until=end_time)
print("--- シミュレーションが正常に終了しました ---")

# ==============================================================================
# 3. データの抽出と物理量の計算
# ==============================================================================
# 空間内のEx成分、Ey成分の電界データを取得
ex_data = sim.get_array(component=mp.Ex)
ey_data = sim.get_array(component=mp.Ey)

# 電界強度絶対値 |E| = sqrt(Ex^2 + Ey^2)
e_mag = np.sqrt(ex_data**2 + ey_data**2)

# 幾何形状マトリクス（PEC領域のフラグ化）
# Meepの物性値（誘電率など）をグリッドごとに取得するか、
# もしくは電界が完全に0（金属内部）の領域を幾何形状フラグとして判定します。
# ここでは物理的に電界が完全にシールドされている領域（金属内部）を1、それ以外を0とします。
# より正確には、各点で極めて小さい、または初期状態から変化がない領域を判定します。
# 確実な方法として、電極内部のマスクを座標から再計算します。
nx, ny = e_mag.shape
x_coords = np.linspace(-cell_x/2, cell_x/2, nx)
y_coords = np.linspace(-cell_y/2, cell_y/2, ny)
X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')

# 上部電極のマスク
mask_top = (Y >= 0.8) & (Y <= 2.5) & (np.abs(X) <= (0.3 + (1.2/1.7)*(Y - 0.8)))
mask_top_ext = (Y > 2.5) & (np.abs(X) <= 1.5) 
# 下部電極のマスク
mask_bottom = (Y <= -0.8) & (Y >= -2.5) & (np.abs(X) <= (0.3 + (1.2/1.7)*(-Y - 0.8)))
mask_bottom_ext = (Y < -2.5) & (np.abs(X) <= 1.5)

pec_flag = np.zeros_like(e_mag)
pec_flag[mask_top | mask_top_ext | mask_bottom | mask_bottom_ext] = 1.0

# ==============================================================================
# 4. ParaView用 VTKレガシーフォーマット (.vtk) への書き出し
# ==============================================================================
output_vtk_path = "test_meep.vtk"
dx = cell_x / nx
dy = cell_y / ny

print(f"ParaView用VTKファイル '{output_vtk_path}' を出力しています...")

with open(output_vtk_path, "w") as f:
    # VTK ヘッダー情報
    f.write("# vtk DataFile Version 3.0\n")
    f.write("2D FDTD Analysis Results - Meep Simulation\n")
    f.write("ASCII\n")
    f.write("DATASET STRUCTURED_POINTS\n")
    
    # グリッドサイズ定義 (2Dデータは Z=1 として3次元格子として定義)
    f.write(f"DIMENSIONS {nx} {ny} 1\n")
    f.write(f"ORIGIN {-cell_x/2} {-cell_y/2} 0.0\n")
    f.write(f"SPACING {dx} {dy} 1.0\n")
    
    # データポイントの総数
    num_points = nx * ny
    f.write(f"POINT_DATA {num_points}\n")
    
    # データ1: 電界強度絶対値 |E|
    f.write("SCALARS E_magnitude float 1\n")
    f.write("LOOKUP_TABLE default\n")
    # VTKフォーマットは Y軸 -> X軸 の順で平坦化（C言語順、ただし3Dインデックス順に注意）
    # 通常、STRUCTURED_POINTSは Xが一番早く動き、次にY、次にZ
    # したがって、(X, Y) グリッドにおいて、Yをループの外、Xを内側にするか、フラット化の軸順を合わせる
    # numpyで indexing='ij' の場合、e_mag[x_idx, y_idx] となるので、
    # 期待される順序 (x0 y0, x1 y0, x2 y0...) にするには、転置してフラット化、または順序通りループ
    for j in range(ny):
        for i in range(nx):
            f.write(f"{e_mag[i, j]:.6f}\n")
            
    # データ2: 幾何形状（PECフラグ: 導体は1.0, 空間は0.0）
    f.write("SCALARS Geometry_PEC float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for j in range(ny):
        for i in range(nx):
            f.write(f"{pec_flag[i, j]:.1f}\n")

print(f"[ファイル出力成功] ParaView用ファイル: '{output_vtk_path}'")
print("ParaViewを起動し、このファイルを読み込んで 'E_magnitude' や 'Geometry_PEC' を切り替えて可視化してください。")
