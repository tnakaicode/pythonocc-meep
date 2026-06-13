import numpy as np
import matplotlib.pyplot as plt

# this try except attempts to import freecad (lowercase) which is the conda
# package name for FreeCAD (mixed case) upon import the conda package appends
# the sys path for Conda installed FreeCAD, consequently FreeCAD can then be
# found by subsequent import statements through out the code base
try:
    import freecad
except ImportError:
    pass

from OCC.Core.gp import gp_Pnt

#from OCC.Display.SimpleGui import init_display

from OCC.Display.OCCViewer import Viewer3d
# 形状作成のためのモジュール（例としてBoxを作成）
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

# 1. 画面を表示しないレンダラーを初期化
display = Viewer3d()
display.Create()
display.SetModeShaded()

# 2. 3Dモデル（形状）を作成してディスプレイに登録
# (ここはご自身の作成した形状に置き換えてください)
my_box = BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape()
display.DisplayShape(my_box, update=True)

# 3. カメラ位置をモデル全体に合わせる
display.FitAll()

# 4. 画像を保存
display.View.Dump("occ.png")
print("画像を保存しました: occ.png")
