---
title: MEEP
---

- <https://flex.phys.tohoku.ac.jp/~maru/>
  - <https://flex.phys.tohoku.ac.jp/~maru/drive-open/>

```bash
wget -r -np -l 0 https://flex.phys.tohoku.ac.jp/~maru/drive-open/
```

install meep

```bash
conda create -c conda-forge -n meep python=3.9 pymeep=*=mpi_mpich* 
#conda install r-nloptr
pip install meshio[all]
pip install autograd bempp-cl necpp scikit-fdiff[interactive,numba] docopt
pip install opencv-python opencv-contrib-python
pip install reportlab python-pptx docx2pdf markdown
pip install moviepy
pip install nlopt
pip install PyMieScatt

conda install -c conda-forge pymeep=*=mpi_mpich*

conda install -c conda-forge pythonocc-core=7.7.2
conda install -c conda-forge freecad
pip install "git+https://github.com/tpaviot/pythonocc-utils.git"

conda env export > meep_occt.yaml
```

## meep

MEEPで用いられるMaxwell方程式は

単位長さaをユーザーが設定する。
真空の誘電率&mathjax{ε_0};が1。
真空の透磁率&mathjax{μ_0};が1。
など通常のMaxwell方程式と異なります。

MEEPで用いられるMaxwell方程式は次のような式変形の最後の形です。

（ドキュメント参考）

通常の以下の形のMaxwell方程式を考えます。

$$ \begin{aligned}
    \nabla\cdot\left(\varepsilon_0\varepsilon_r\boldsymbol{E}\right)=\rho \\
    \nabla\cdot\left(\mu_0 \boldsymbol{H}\right)= 0\\
    \nabla\times\boldsymbol{E}= - \frac{\partial}{\partial t}\left(\mu_0 \boldsymbol{H}\right)\\
    \nabla\times\boldsymbol{H}=\boldsymbol{J} + \frac{\partial}{\partial t}\left(\varepsilon_0\varepsilon_r\boldsymbol{E}\right)
\end{aligned} $$

次のように変形します。

$$ \begin{aligned}
    &(a\nabla)\cdot \left(1+i\frac{(\frac{a\sigma}{c\varepsilon_0\varepsilon_r})}{(\frac{a\omega}{c})}\right) (\varepsilon_r\sqrt{\varepsilon_0}\boldsymbol{E})=0\\
    &(a\nabla)\cdot(\sqrt{\mu_0}\boldsymbol{H})=0\\
    &(a\nabla)\times(\sqrt\varepsilon_0\boldsymbol{E})=-\frac{\partial}{\partial(\frac{ct}{a})}(\sqrt{\mu_0}\boldsymbol{H})\\
    &(a\nabla)\times(\sqrt{\mu_0}\boldsymbol{H})=\frac{\partial}{\partial(\frac{ct}{a})}\left(1+i\frac{(\frac{a\sigma}{c\varepsilon_0\varepsilon_r})}{(\frac{a\omega}{c})}\right)(\varepsilon_r\sqrt{\varepsilon_0}\boldsymbol{E})
\end{aligned} $$

$a$はシミュレーションで使用する単位長さです。$c$は光速です。

以下のように&mathjax{x',y',z',t',ω',E',H',σ'};を定めます。

$$ \begin{aligned}
    & x'=\frac{x}{a},y'=\frac{y}{a},z'=\frac{z}{a}\\
    &\boldsymbol{E}'=\sqrt{\varepsilon_0}\boldsymbol{E},\boldsymbol{H}'=\sqrt{\varepsilon_0}\boldsymbol{H}\\
    &t'=\frac{ct}{a},\omega'=\frac{a\omega}{c},\sigma'=\frac{a\sigma}{c\varepsilon_0\varepsilon_r}\\
\end{aligned} $$

このとき、Maxwell方程式は次のようになります。

$$ \begin{aligned}
    &\nabla'\cdot\left(1+i\frac{\sigma'}{\omega'}\right)\varepsilon_r\boldsymbol{E}'=0\\
    &\nabla'\cdot\boldsymbol{H}'=0\\
    &\nabla'\times\boldsymbol{E}'=-\frac{\partial}{\partial t'}\boldsymbol{H}'\\
    &\nabla'\times\boldsymbol{H}'=\frac{\partial}{\partial t'}\left(1+i\frac{\sigma'}{\omega'}\right)\varepsilon_r\boldsymbol{E}'\\
\end{aligned} $$

この形のMaxwell方程式をMEEPで使っています。このとき、真空の誘電率、透磁率はどちらも1です。

また、単位長さは$a[m]$です。MEEP中の単位時間$a/c$に真空中の光は単位長さaだけ進みます。

この光の速度をSI単位系に直すと、&mathjax{a / (a/c) = c}; となり、矛盾していないことが確認できますね。

MEEPの単位系では光は単位時間&mathjax{c/a};に単位長さ&mathjax{a};だけ進むことから、周波数 &mathjax{f [c/a]}; の真空中の光の波長は &mathjax{1/f [a]}; です。

&mathjax{x',y',z',t',ω',σ'};はすべて無次元です。

金属が存在しない場合(&mathjax{σ=0};)は、&mathjax{a};の具体的な値を決めずにシミュレーションを行ったあとに、&mathjax{a};の値を決めてスケール変換することで&mathjax{a};の値に沿った正しいシミュレーション結果を得られます。

つまり、一回のシミュレーションで任意の単位長さaに対するシミュレーション結果が実質的に得られます。

一方、金属が存在する場合(σ!=0)では、MEEPのシミュレーションに用いる伝導率が単位長さaの値に依存するので、事前にaの値を定めておく必要があります。

電場E'と磁場H'のMEEP中での単位は等しくなっており、平面波の場合には電場と磁場の値が等しくなっています。

確かめてみましょう。平面波の場合、インピーダンス$Z=\sqrt{\mu/\varepsilon}$を用いると磁場$H'=\sqrt{\mu}H ={\sqrt{\mu}/Z}E=\sqrt{\varepsilon}E=E'$となります。
$H'=E'$であることが確認できましたね。

角振動数ω'（MEEP単位系）で 比誘電率を &mathjax{ε=ε_{re} + i ε_{im}}; と定めたいときは

$$\epsilon_{r} = \epsilon_{re}$$

$$\sigma' = \omega' \frac{\epsilon_{im}}{ \epsilon_{re}}$$

とすればいいことも上の式から分かります。

## meep Simulation Class

All `Simulation` attributes are described in further detail below.
          In brackets after each variable is the type of value that it should hold.
          The classes, complex datatypes like `GeometricObject`, are described in a later subsection.
          The basic datatypes, like `integer`, `boolean`, `complex`, and `string` are defined by Python.
          `Vector3` is a `meep` class.

**`geometry` [ list of `GeometricObject` class ]** —
          Specifies the geometric objects making up the structure being simulated.
          When objects overlap, later objects in the list take precedence.
          Defaults to no objects (empty list).

**`geometry_center` [ `Vector3` class ]** —
          Specifies the coordinates of the center of the cell.
          Defaults to (0, 0, 0), but changing this allows you to shift the coordinate system used in Meep
          (for example, to put the origin at the corner).
          Passing `geometry_center=c` is equivalent to adding the `c` vector to the coordinates of every other object in the simulation,
          i.e. `c` becomes the new origin that other objects are defined with respect to.

**`sources` [ list of `Source` class ]** —
          Specifies the current sources to be present in the simulation. Defaults to none (empty list).

**`symmetries` [ list of `Symmetry` class ]** — Specifies the spatial symmetries
          (mirror or rotation) to exploit in the simulation. Defaults to none (empty
          list). The symmetries must be obeyed by *both* the structure and the sources.
          See also [Exploiting Symmetry](Exploiting_Symmetry.md).

**`boundary_layers` [ list of `PML` class ]** — Specifies the
          [PML](Perfectly_Matched_Layer.md) absorbing boundary layers to use. Defaults to
          none.

**`cell_size` [ `Vector3` ]** — Specifies the size of the cell which is centered
          on the origin of the coordinate system. Any sizes of 0 imply a
          reduced-dimensionality calculation. Strictly speaking, the dielectric function
          is taken to be uniform along that dimension. A 2d calculation is especially
          optimized. See `dimensions` below. **Note:** because Maxwell's equations are
          scale invariant, you can use any units of distance you want to specify the cell
          size: nanometers, microns, centimeters, etc. However, it is usually convenient
          to pick some characteristic lengthscale of your problem and set that length to 1.
          See also [Units](Introduction.md#units-in-meep). Required argument (no default).

**`default_material` [`Medium` class ]** — Holds the default material that is
          used for points not in any object of the geometry list. Defaults to `air` (ε=1).
          This can also be a NumPy array that defines a dielectric function much like
          `epsilon_input_file` below (see below). If you want to use a material function
          as the default material, use the `material_function` keyword argument (below).

**`material_function` [ function ]** — A Python function that takes a `Vector3`
          and returns a `Medium`. See also [Material Function](#medium).
          Defaults to `None`.

**`epsilon_func` [ function ]** — A Python function that takes a `Vector3` and
          returns the dielectric constant at that point. See also [Material
          Function](#medium). Defaults to `None`.

**`epsilon_input_file` [`string`]** — If this string is not empty (the default),
          then it should be the name of an HDF5 file whose first/only dataset defines a
          scalar, real-valued, frequency-independent dielectric function over some
          discrete grid. Alternatively, the dataset name can be specified explicitly if
          the string is in the form "filename:dataset". This dielectric function is then
          used in place of the ε property of `default_material` (i.e. where there are no
          `geometry` objects). The grid of the epsilon file dataset need *not* match the
          computational grid; it is scaled and/or linearly interpolated as needed to map
          the file onto the cell. The structure is warped if the proportions of the grids
          do not match. **Note:** the file contents only override the ε property of the
          `default_material`, whereas other properties (μ, susceptibilities,
          nonlinearities, etc.) of `default_material` are still used.

**`dimensions` [`integer`]** — Explicitly specifies the dimensionality of the
          simulation, if the value is less than 3. If the value is 3 (the default), then
          the dimensions are automatically reduced to 2 if possible when `cell_size` in
          the $z$ direction is `0`. If `dimensions` is the special value of `CYLINDRICAL`,
          then cylindrical coordinates are used and the $x$ and $z$ dimensions are
          interpreted as $r$ and $z$, respectively. If `dimensions` is 1, then the cell
          must be along the $z$ direction and only $E_x$ and $H_y$ field components are
          permitted. If `dimensions` is 2, then the cell must be in the $xy$ plane.

**`m` [`number`]** — For `CYLINDRICAL` simulations, specifies that the angular
          $\\phi$ dependence of the fields is of the form $e^{im\\phi}$ (default is `m=0`).
          If the simulation cell includes the origin $r=0$, then `m` must be an integer.

**`accurate_fields_near_cylorigin` [`boolean`]** — For `CYLINDRICAL` simulations
          with |*m*| &gt; 1, compute more accurate fields near the origin $r=0$ at the
          expense of requiring a smaller Courant factor. Empirically, when this option is
          set to `True`, a Courant factor of roughly $\\min[0.5, 1 / (|m| + 0.5)]$ or
          smaller seems to be needed. Default is `False`, in which case the $D_r$, $D_z$,
          and $B_r$ fields within |*m*| pixels of the origin are forced to zero, which
          usually ensures stability with the default Courant factor of 0.5, at the expense
          of slowing convergence of the fields near $r=0$.

**`resolution` [`number`]** — Specifies the computational grid resolution in
          pixels per distance unit. Required argument. No default.

**`k_point` [`False` or `Vector3`]** — If `False` (the default), then the
          boundaries are perfect metallic (zero electric field). If a `Vector3`, then the
          boundaries are Bloch-periodic: the fields at one side are
          $\\exp(i\\mathbf{k}\\cdot\\mathbf{R})$ times the fields at the other side, separated
          by the lattice vector $\\mathbf{R}$. A non-zero `Vector3` will produce complex
          fields. The `k_point` vector is specified in Cartesian coordinates in units of
          2π/distance. Note: this is *different* from [MPB](https://mpb.readthedocs.io),
          equivalent to taking MPB's `k_points` through its function
          `reciprocal->cartesian`.

**`kz_2d` [`"complex"`, `"real/imag"`, or `"3d"`]** — A 2d cell (i.e.,
          `dimensions=2`) combined with a `k_point` that has a *non-zero* component in $z$
          would normally result in a 3d simulation with complex fields. However, by
          default (`kz_2d="complex"`), Meep will use a 2d computational cell in which
          $k_z$ is incorporated as an additional term in Maxwell's equations, which still
          results in complex fields but greatly improved performance. Setting `kz_2d="3d"`
          will instead use a 3d cell that is one pixel thick (with Bloch-periodic boundary
          conditions), which is considerably more expensive. The third possibility,
          `kz_2d="real/imag"`, saves an additional factor of two by storing some field
          components as purely real and some as purely imaginary in a "real" field, but
          this option requires some care to use. See [2d Cell with Out-of-Plane
          Wavevector](2d_Cell_Special_kz.md).

**`ensure_periodicity` [`boolean`]** — If `True` (the default) *and* if the
          boundary conditions are periodic (`k_point` is not `False`), then the geometric
          objects are automatically repeated periodically according to the lattice vectors
          which define the size of the cell.

**`eps_averaging` [`boolean`]** — If `True` (the default), then [subpixel
          averaging](Subpixel_Smoothing.md) is used when initializing the dielectric
          function. For simulations involving a [material function](#medium),
          `eps_averaging` is `False` (the default) and must be
          [enabled](Subpixel_Smoothing.md#enabling-averaging-for-material-function) in
          which case the input variables `subpixel_maxeval` (default 10<sup>4</sup>) and
          `subpixel_tol` (default 10<sup>-4</sup>) specify the maximum number of function
          evaluations and the integration tolerance for the adaptive numerical
          integration. Increasing/decreasing these, respectively, will cause a more
          accurate but slower computation of the average ε with diminishing returns for
          the actual FDTD error. Disabling subpixel averaging will lead to [staircasing
          effects and irregular
          convergence](Subpixel_Smoothing.md#what-happens-when-subpixel-smoothing-is-disabled).

**`force_complex_fields` [`boolean`]** — By default, Meep runs its simulations
          with purely real fields whenever possible. It uses complex fields which require
          twice the memory and computation if the `k_point` is non-zero or if `m` is
          non-zero. However, by setting `force_complex_fields` to `True`, Meep will always
          use complex fields.

**`force_all_components` [`boolean`]** — By default, in a 2d simulation Meep
          uses only the field components that might excited by your current sources:
          either the in-plane $(E_x,E_y,H_z)$ or out-of-plane $(H_x,H_y,E_z)$ polarization,
          depending on the source.  (Both polarizations are excited if you use multiple source
          polarizations, or if an anisotropic medium is present that couples the two
          polarizations.)   In rare cases (primarily for combining results of multiple
          simulations with differing polarizations), you might want to force it to
          simulate all fields, even those that remain zero throughout the simulation, by
          setting `force_all_components` to `True`.

**`filename_prefix` [`string`]** — A string prepended to all output filenames
          (e.g., for HDF5 files). If `None` (the default), then Meep constructs a default
          prefix based on the current Python filename ".py" replaced by "-" (e.g. `foo.py`
          uses a `"foo-"` prefix). You can get this prefix by calling `get_filename_prefix`.

**`Courant` [`number`]** — Specify the
          [Courant factor](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)
          $S$ which relates the time step size to the spatial discretization: $cΔ t = SΔ x$.
          Default is 0.5. For numerical stability, the Courant factor must be *at
          most* $n_\\textrm{min}/\\sqrt{\\textrm{# dimensions}}$, where $n_\\textrm{min}$ is
          the minimum refractive index (usually 1), and in practice $S$ should be slightly
          smaller.

**`output_volume` [`Volume` class ]** — Specifies the default region of space
          that is output by the HDF5 output functions (below); see also the `Volume` class
          which manages `meep::volume*` objects. Default is `None`, which means that the
          whole cell is output. Normally, you should use the `in_volume(...)` function to
          modify the output volume instead of setting `output_volume` directly.

**`output_single_precision` [`boolean`]** — Meep performs its computations in
          [double precision](https://en.wikipedia.org/wiki/double_precision), and by
          default its output HDF5 files are in the same format. However, by setting this
          variable to `True` (default is `False`) you can instead output in [single
          precision](https://en.wikipedia.org/wiki/single_precision) which saves a factor
          of two in space.

**`progress_interval` [`number`]** — Time interval (seconds) after which Meep
          prints a progress message. Default is 4 seconds.

**`extra_materials` [ list of `Medium` class ]** — By default, Meep turns off
          support for material dispersion ([susceptibilities](#susceptibility) or
          [conductivity](Materials.md#conductivity-and-complex)) or nonlinearities if none
          of the objects in `geometry` have materials with these properties &mdash; since
          they are not needed, it is faster to omit their calculation. This doesn't work,
          however, if you use a `material_function`: materials via a user-specified
          function of position instead of just geometric objects. If your material
          function only returns a nonlinear material, for example, Meep won't notice this
          unless you tell it explicitly via `extra_materials`. `extra_materials` is a list
          of materials that Meep should look for in the cell in addition to any materials
          that are specified by geometric objects. You should list any materials other
          than scalar dielectrics that are returned by `material_function` here.

**`load_structure` [`string`]** — If not empty, Meep will load the structure
          file specified by this string. The file must have been created by
          `mp.dump_structure`. Defaults to an empty string. See [Load and Dump
          Structure](#load-and-dump-structure) for more information.

**`chunk_layout` [`string` or `Simulation` instance or `BinaryPartition` class]** —
          This will cause the `Simulation` to use the chunk layout described by either
          (1) an `.h5` file (created using `Simulation.dump_chunk_layout`), (2) another
          `Simulation` instance, or (3) a [`BinaryPartition`](#binarypartition) class object.
          For more information, see [Load and Dump Structure](#load-and-dump-structure) and
          [Parallel Meep/User-Specified Cell Partition](Parallel_Meep.md#user-specified-cell-partition).
          The following require a bit more understanding of the inner workings of Meep to use. See also [SWIG Wrappers](#swig-wrappers).

**`structure` [`meep::structure*`]** — Pointer to the current structure being
          simulated; initialized by `_init_structure` which is called automatically by
          `init_sim()` which is called automatically by any of the [run
          functions](#run-functions). The structure initialization is handled by the
          `Simulation` class, and most users will not need to call `_init_structure`.

**`fields` [`meep::fields*`]** — Pointer to the current fields being simulated;
          initialized by `init_sim()` which is called automatically by any of the [run
          functions](#run-functions).

**`num_chunks` [`integer`]** — Minimum number of "chunks" (subregions) to divide
          the structure/fields into. Overrides the default value determined by
          the number of processors, PML layers, etcetera. Mainly useful for debugging.

**`split_chunks_evenly` [`boolean`]** — When `True` (the default), the work per
          [chunk](Chunks_and_Symmetry.md) is not taken into account when splitting chunks
          up for multiple processors. The cell is simply split up into equal chunks (with
          the exception of PML regions, which must be on their own chunk). When `False`,
          Meep attempts to allocate an equal amount of work to each processor, which can
          increase the performance of [parallel simulations](Parallel_Meep.md).

## Import/Export

- DGSII file
  - <https://meep.readthedocs.io/en/latest/Python_User_Interface/#gdsii-support>
  - <https://meep.readthedocs.io/en/latest/Python_Tutorials/GDSII_Import/>
