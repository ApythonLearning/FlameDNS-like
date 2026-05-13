# FlameDNS-like

本仓库用于准备、运行和后处理真实的 **PeleLMeX 低马赫数二维轴对称 H2/O2 预混球形火焰算例**。

重要原则：
- Cantera 只用于计算 1D 预混火焰 profile 和参考量，例如 `S_L`、`delta_T`、`Tb`、`rho_u`、`rho_b`。
- 真实二维结果必须来自 PeleLMeX 运行生成的 AMReX plotfile，例如 `plt00020/`、`plt00040/`。
- 后处理必须读取真实 plotfile 中的温度场，而不是读取 synthetic 火焰半径。

官方参考：

- PeleLMeX 文档：https://amrex-combustion.github.io/PeleLMeX/
- PeleLMeX 源码：https://github.com/AMReX-Combustion/PeleLMeX
- PeleLMeX 输入参数说明：https://amrex-combustion.github.io/PeleLMeX/manual/html/LMeXControls.html

## 1. 推荐环境

建议使用 Linux 或 WSL2 Ubuntu。原生 Windows 上编译 PeleLMeX/AMReX 通常比较麻烦。

安装基础依赖：

```bash
sudo apt update
sudo apt install -y git build-essential gfortran cmake mpich libopenmpi-dev python3 python3-pip
```

## 2. 安装 PeleLMeX / AMReX

下载 PeleLMeX 及 submodules：

```bash
git clone --recursive --shallow-submodules --single-branch https://github.com/AMReX-Combustion/PeleLMeX.git
cd PeleLMeX
```

设置环境变量：

```bash
export PELE_HOME=$PWD
export PELE_PHYSICS_HOME=$PELE_HOME/Submodules/PelePhysics
export AMREX_HOME=$PELE_PHYSICS_HOME/Submodules/amrex
export AMREX_HYDRO_HOME=$PELE_HOME/Submodules/AMReX-Hydro
export SUNDIALS_HOME=$PELE_PHYSICS_HOME/Submodules/sundials
```

建议写入 `~/.bashrc`：

```bash
echo "export PELE_HOME=$PELE_HOME" >> ~/.bashrc
echo "export PELE_PHYSICS_HOME=$PELE_PHYSICS_HOME" >> ~/.bashrc
echo "export AMREX_HOME=$AMREX_HOME" >> ~/.bashrc
echo "export AMREX_HYDRO_HOME=$AMREX_HYDRO_HOME" >> ~/.bashrc
echo "export SUNDIALS_HOME=$SUNDIALS_HOME" >> ~/.bashrc
```

先编译官方示例，确认 PeleLMeX 环境可用：

```bash
cd $PELE_HOME/Exec/RegTests/FlameSheet
make TPL
make -j 8
mpiexec -np 2 ./PeleLMeX2d.gnu.MPI.ex input.2d-regt
```

如果官方示例不能运行，先修复 PeleLMeX/AMReX 环境，不要继续跑本仓库算例。

## 3. 安装本仓库 Python 依赖

进入本仓库：

```bash
cd /mnt/e/Projects/CombustionSim/FlameDNS-like/spherical_flame_dns_like
python3 -m pip install -r requirements.txt
```

Windows PowerShell 路径：

```powershell
cd E:\Projects\CombustionSim\FlameDNS-like\spherical_flame_dns_like
python -m pip install -r requirements.txt
```

如果 Windows/conda 环境中 `h5py` 出现 DLL import 错误，优先使用 conda-forge：

```bash
conda install -c conda-forge cantera h5py hdf5 yt pandas matplotlib pyyaml
```

检查依赖：

```bash
python3 scripts/check_dependencies.py
```

需要确认：

- Cantera 可 import
- `h5py` 可 import
- `yt` 可 import
- `PELE_HOME` 已设置
- `AMREX_HOME` 已设置
- PeleLMeX 可在 case 目录中编译运行

## 4. PeleLMeX 工作流总览

PeleLMeX 的工作流可以理解为下面这条链：

```text
Cantera 1D flame
    -> 生成 H2/O2 profile 和参考尺度
    -> 复制/生成 PeleLMeX Exec case
    -> 修改 inputs.axisym_h2o2
    -> pelelmex_prob.H/.cpp 初始化二维 RZ 场和高温火核
    -> make TPL && make
    -> 运行 PeleLMeX
    -> 生成 pltXXXXX/ plotfile
    -> 用 yt 读取 plotfile 后处理
    -> 输出 flame_radii.csv 和图像
```

PeleLMeX case 目录的关键文件：

```text
AxisymmetricH2O2/
  GNUmakefile
  inputs.axisym_h2o2
  pelelmex_prob.H
  pelelmex_prob.cpp
```

各文件职责：

- `GNUmakefile`：指定维度、编译器、MPI、EOS、transport、chemistry model，并包含 PeleLMeX 的 makefile。
- `inputs.axisym_h2o2`：运行时参数，包括网格、边界、时间步、输出频率、点火温度、点火半径、H2/O2 配比、重力等。
- `pelelmex_prob.H`：定义 problem parameters，并实现每个网格单元的初始化逻辑。
- `pelelmex_prob.cpp`：从 `inputs.axisym_h2o2` 读取 `prob.*` 参数。

PeleLMeX 运行后会生成：

```text
plt00020/
plt00040/
chk00100/
```

其中：

- `pltXXXXX/` 是后处理用的 plotfile。
- `chkXXXXX/` 是 restart/checkpoint 文件，可用于续算。

## 5. 计算 1D Cantera H2/O2 预混火焰

先用较容易收敛的 `phi=0.6` 测试：

```bash
python3 scripts/compute_1d_h2o2_flame.py --phi 0.6 --output-dir profiles/h2o2_phi_060
```

输出：

```text
profiles/h2o2_phi_060/
  h2o2_1d_summary.csv
  h2o2_1d_summary.json
  h2o2_1d_profile.csv
  h2o2_1d_profile.h5
```

`h2o2_1d_summary.json` 中包含：

- `S_L_m_s`
- `delta_T_m`
- `Tb_K`
- `rho_u_kg_m3`
- `rho_b_kg_m3`

`h2o2_1d_profile.csv` 和 `h2o2_1d_profile.h5` 保存温度、速度、密度、质量分数和摩尔分数。后续可用于更精细的 PeleLMeX 初始场映射。

注意：非常稀的 H2/O2，例如：

```bash
python3 scripts/compute_1d_h2o2_flame.py --h2-volume-fraction 0.10 --output-dir profiles/h2o2_h2vol_010
```

可能在 Cantera 自由火焰求解中不收敛。脚本会直接失败，不会生成假 profile。

## 6. 准备 PeleLMeX Case

本仓库提供模板：

```text
spherical_flame_dns_like/examples/pelelmex_axisymmetric_h2o2/
```

复制到 PeleLMeX：

```bash
cd /mnt/e/Projects/CombustionSim/FlameDNS-like/spherical_flame_dns_like
cp -r examples/pelelmex_axisymmetric_h2o2 $PELE_HOME/Exec/RegTests/AxisymmetricH2O2
cd $PELE_HOME/Exec/RegTests/AxisymmetricH2O2
```

检查 `GNUmakefile` 中的路径和模型选择：

```makefile
PELE_HOME ?= $(HOME)/PeleLMeX
DIM             = 2
USE_MPI         = TRUE
Eos_Model       := Fuego
Transport_Model := Simple
Chemistry_Model := LiDryer
```

如果你的 PeleLMeX 安装路径不是 `$(HOME)/PeleLMeX`，请修改 `PELE_HOME` 或使用环境变量：

```bash
export PELE_HOME=/path/to/PeleLMeX
```

## 7. 修改 PeleLMeX 输入参数

编辑：

```bash
vim inputs.axisym_h2o2
```

几何和轴对称设置：

```text
geometry.coord_sys = 1
geometry.prob_lo   = 0.0   -0.006
geometry.prob_hi   = 0.006  0.006
amr.n_cell         = 320 640
```

这里 `geometry.coord_sys = 1` 表示 RZ 轴对称坐标。物理域示例为：

```text
r: 0 到 6 mm
z: -6 mm 到 6 mm
```

时间推进和输出：

```text
amr.max_step = 200
amr.stop_time = 0.0035
amr.cfl = 0.3
amr.plot_int = 20
amr.check_int = 100
```

问题参数：

```text
prob.P_mean = 101325.0
prob.T0 = 298.0
prob.phi = 0.6
prob.h2_volume_fraction = 0.5454545454545454
prob.use_phi = 1
prob.ignition_radius = 0.00075
prob.initial_flame_radius = 0.0015
prob.ignition_temperature = 1800.0
prob.gravity_magnitude = 0.0
prob.gravity_direction = 0.0 -1.0
prob.profile_csv = ../../profiles/h2o2_phi_060/h2o2_1d_profile.csv
```

如果使用 `phi`，令：

```text
prob.use_phi = 1
```

如果直接使用 H2 体积分数，令：

```text
prob.use_phi = 0
prob.h2_volume_fraction = 0.10
```

## 8. 初始化逻辑

`pelelmex_prob.H` 中的初始化逻辑做以下事情：

1. 读取当前 cell center 的 RZ 坐标。
2. 计算到原点的球形半径：

```text
radius = sqrt(r^2 + z^2)
```

3. 如果 `radius <= prob.ignition_radius`，设置：

```text
T = prob.ignition_temperature
```

否则设置：

```text
T = prob.T0
```

4. 根据 `phi` 或 `h2_volume_fraction` 得到 H2/O2 摩尔分数。
5. 调用 PelePhysics EOS 将摩尔分数转为质量分数。
6. 根据 `P_mean`、`T` 和质量分数计算密度。
7. 写入 PeleLMeX state，包括密度、动量、温度、焓和组分。

这是一个最小可运行初始化。后续如果要使用 Cantera 1D profile 初始化有限厚度火焰，需要在 `pelelmex_prob.H/.cpp` 中读取或嵌入 profile 并按 `sqrt(r^2+z^2)-initial_flame_radius` 插值。

## 9. 编译 PeleLMeX Case

在 case 目录中：

```bash
cd $PELE_HOME/Exec/RegTests/AxisymmetricH2O2
make TPL
make -j 8
```

如果修改了化学机制、EOS、transport 或编译选项，建议清理后重编：

```bash
make TPLrealclean
make realclean
make TPL
make -j 8
```

编译成功后会生成类似：

```text
PeleLMeX2d.gnu.MPI.ex
```

## 10. 运行 PeleLMeX

串行或小规模 MPI 测试：

```bash
mpiexec -np 2 ./PeleLMeX2d.gnu.MPI.ex inputs.axisym_h2o2
```

运行过程中重点检查：

- 时间步是否稳定
- `cfl` 是否小于 0.5
- 温度是否发散
- 质量分数是否出现明显负值
- plotfile 是否按 `amr.plot_int` 正常输出

短测试建议先用：

```text
amr.max_step = 20
amr.plot_int = 5
```

确认稳定后再提高步数。

## 11. PeleLMeX Restart 工作流

如果运行生成了 checkpoint，例如：

```text
chk00100/
```

可以从 checkpoint 续算。常见方式是在 inputs 中加入或修改：

```text
amr.restart = chk00100
amr.max_step = 300
```

然后重新运行：

```bash
mpiexec -np 2 ./PeleLMeX2d.gnu.MPI.ex inputs.axisym_h2o2
```

restart 时不要删除对应 `chkXXXXX/` 目录。

## 12. DNS 分辨率检查

使用 Cantera 的 `delta_T_m` 检查网格：

```bash
cd /mnt/e/Projects/CombustionSim/FlameDNS-like/spherical_flame_dns_like
python3 scripts/check_dns_resolution.py \
  --summary profiles/h2o2_phi_060/h2o2_1d_summary.json \
  --dx 1.875e-5 \
  --cfl 0.3
```

硬性要求：

```text
dx <= delta_T / 10
CFL < 0.5
```

推荐：

```text
dx <= delta_T / 20
```

如果检查失败，不应称为 DNS。需要加密网格、缩小物理域或降低时间步。

## 13. 后处理真实 PeleLMeX Plotfile

后处理必须读取真实 plotfile：

```bash
python3 scripts/postprocess_pelelmex_plotfiles.py \
  $PELE_HOME/Exec/RegTests/AxisymmetricH2O2/plt00020 \
  $PELE_HOME/Exec/RegTests/AxisymmetricH2O2/plt00040 \
  $PELE_HOME/Exec/RegTests/AxisymmetricH2O2/plt00060 \
  --summary profiles/h2o2_phi_060/h2o2_1d_summary.json \
  --output-dir postprocess_pelelmex \
  --cfl 0.3
```

脚本会读取真实 temperature field，提取：

```text
T = (Tu + Tb) / 2
```

并输出：

```text
postprocess_pelelmex/flame_radii.csv
postprocess_pelelmex/flame_radii.png
postprocess_pelelmex/flame_contours.png
```

`flame_radii.csv` 包含：

- `time_s`
- `radius_upper_m`
- `radius_lower_m`
- `radius_max_m`
- `contour_point_count`

如果 plotfile 不存在，脚本会直接失败，不会生成占位数据。

## 14. 最小完整命令流程

```bash
# 1. 准备 Python 依赖
cd /mnt/e/Projects/CombustionSim/FlameDNS-like/spherical_flame_dns_like
python3 -m pip install -r requirements.txt
python3 scripts/check_dependencies.py

# 2. 计算真实 1D Cantera profile
python3 scripts/compute_1d_h2o2_flame.py --phi 0.6 --output-dir profiles/h2o2_phi_060

# 3. 复制 PeleLMeX case
cp -r examples/pelelmex_axisymmetric_h2o2 $PELE_HOME/Exec/RegTests/AxisymmetricH2O2
cd $PELE_HOME/Exec/RegTests/AxisymmetricH2O2

# 4. 编译
make TPL
make -j 8

# 5. 运行
mpiexec -np 2 ./PeleLMeX2d.gnu.MPI.ex inputs.axisym_h2o2

# 6. 后处理真实 plotfile
cd /mnt/e/Projects/CombustionSim/FlameDNS-like/spherical_flame_dns_like
python3 scripts/postprocess_pelelmex_plotfiles.py \
  $PELE_HOME/Exec/RegTests/AxisymmetricH2O2/plt00060 \
  $PELE_HOME/Exec/RegTests/AxisymmetricH2O2/plt00070 \
  --summary profiles/h2o2_phi_060/h2o2_1d_summary.json \
  --output-dir postprocess_pelelmex \
  --cfl 0.3
```

## 15. 当前限制

当前模板已经包含 PeleLMeX case 文件和球形高温火核初始化，但仍有两个工程限制：

- 本机当前没有可用的 PeleLMeX/AMReX 编译环境，因此不能在本仓库内自动完成 PeleLMeX 编译验证。
- `pelelmex_prob.H` 目前实现的是最小高温火核初始化；如需用 Cantera 1D profile 生成有限厚度初始球形火焰，需要进一步实现 profile 读取和径向插值。

## 16. PeleLMeX 常见编译问题

### CMake 版本过低

如果 `make TPL` 编译 SUNDIALS 时出现：

```text
CMake 3.18 or higher is required. You are running version 3.16.3
```

说明系统默认 CMake 太旧。Ubuntu 20.04 默认常见版本是 3.16.x，而当前 PelePhysics/SUNDIALS 需要至少 3.18。

检查当前版本：

```bash
cmake --version
which cmake
```

推荐用 Kitware 官方 apt 源安装新版 CMake：

```bash
sudo apt update
sudo apt install -y ca-certificates gpg wget
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null \
  | gpg --dearmor \
  | sudo tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null
echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ focal main" \
  | sudo tee /etc/apt/sources.list.d/kitware.list >/dev/null
sudo apt update
sudo apt install -y cmake
cmake --version
```

如果你不是 Ubuntu 20.04，把上面命令里的 `focal` 改成你的发行版代号，例如：

```bash
lsb_release -cs
```

常见代号：

```text
20.04 -> focal
22.04 -> jammy
24.04 -> noble
```

如果之前已经失败过，建议清理第三方库构建目录后重试：

```bash
cd $PELE_HOME/Exec/RegTests/AxisymmetricH2O2
make TPLrealclean
make realclean
make TPL
make -j 8
```

也可以先在官方示例里验证：

```bash
cd $PELE_HOME/Exec/RegTests/FlameSheet
make TPLrealclean
make realclean
make TPL
make -j 8
```

### g++ 不支持 `-std=c++20`

如果编译时出现：

```text
g++: error: unrecognized command line option '-std=c++20'; did you mean '-std=c++2a'?
```

说明当前 `g++` 太旧。当前 PeleLMeX/AMReX 构建链会使用 C++20，建议使用 GCC/G++ 10 或更新版本。

先检查版本：

```bash
g++ --version
gcc --version
mpicxx --version
```

Ubuntu 20.04 可安装 GCC 10：

```bash
sudo apt update
sudo apt install -y gcc-10 g++-10
```

将系统默认 `gcc/g++` 切到新版：

```bash
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100
sudo update-alternatives --config gcc
sudo update-alternatives --config g++
```

确认：

```bash
g++ --version
g++ -std=c++20 -x c++ - -fsyntax-only <<< "int main(){return 0;}"
```

如果使用 OpenMPI，确认 MPI C++ wrapper 调用的也是新版编译器：

```bash
mpicxx --showme:command
```

如果它仍指向旧编译器，可以临时指定：

```bash
export OMPI_CXX=g++-10
```

然后清理并重编：

```bash
cd $PELE_HOME/Exec/RegTests/AxisymmetricH2O2
make TPLrealclean
make realclean
make TPL
make -j 8
```

如果你使用 MPICH，可检查：

```bash
mpicxx -show
```

必要时重新安装 MPI，或确保 MPI wrapper 使用更新后的系统 `g++`。

### `mpiexec` 无法启动 PeleLMeX 可执行文件

如果运行时出现：

```text
mpiexec was unable to launch the specified application as it could not access
or execute an executable:

Executable: ./PeleLMeX2d.gnu.MPI.ex
```

通常说明当前目录下没有这个文件、文件名不一致、没有执行权限，或 MPI wrapper/运行环境不匹配。

先确认当前目录：

```bash
pwd
ls -lh
ls -lh PeleLMeX*ex
```

如果没有 `PeleLMeX2d.gnu.MPI.ex`，查找真实生成的文件名：

```bash
find . -maxdepth 2 -type f -name "*PeleLMeX*ex" -ls
```

常见文件名可能类似：

```text
PeleLMeX2d.gnu.MPI.ex
PeleLMeX2d.gnu.MPI.DEBUG.ex
PeleLMeX2d.gnu.ex
```

用实际文件名运行：

```bash
mpiexec -np 2 ./实际文件名 inputs.axisym_h2o2
```

如果文件存在但没有执行权限：

```bash
chmod +x PeleLMeX2d.gnu.MPI.ex
```

然后先不通过 MPI，直接测试可执行文件能否启动：

```bash
./PeleLMeX2d.gnu.MPI.ex inputs.axisym_h2o2
```

如果直接运行也报：

```text
Permission denied
```

检查当前路径是否在 Windows 挂载目录或 noexec 文件系统中：

```bash
mount | grep "$(df -P . | tail -1 | awk '{print $1}')"
```

在 WSL 中，如果 case 放在 `/mnt/c` 或 `/mnt/e` 下并遇到执行权限问题，建议把 PeleLMeX 和 case 放在 Linux 文件系统中，例如：

```bash
mkdir -p ~/Projects
cp -r /mnt/e/Projects/CombustionSim/FlameDNS-like/spherical_flame_dns_like/examples/pelelmex_axisymmetric_h2o2 \
  $PELE_HOME/Exec/RegTests/AxisymmetricH2O2
```

如果文件存在、权限正常，但 `mpiexec` 仍失败，检查动态库依赖：

```bash
ldd ./PeleLMeX2d.gnu.MPI.ex | grep "not found"
```

如果有缺失库，需要先修复 `LD_LIBRARY_PATH` 或重新编译 TPL。

最后确认不是编译失败后误运行：

```bash
make -j 8
echo $?
```

`make` 必须返回 `0`，并且目录中必须实际生成 `.ex` 文件。

### 找不到 `Exec/Make/PeleLMeX.mak`

如果编译时出现：

```text
GNUmakefile:22: /path/to/PeleLMeX/Exec/Make/PeleLMeX.mak: No such file or directory
make: *** No rule to make target '/path/to/PeleLMeX/Exec/Make/PeleLMeX.mak'. Stop.
```

说明 case 的 `GNUmakefile` include 路径与当前 PeleLMeX 版本不匹配。当前 PeleLMeX 常见路径是：

```text
$PELE_HOME/Exec/Make.PeleLMeX
```

检查你的 PeleLMeX 目录：

```bash
ls $PELE_HOME/Exec/Make.PeleLMeX
ls $PELE_HOME/Exec/Make
```

如果存在 `Exec/Make.PeleLMeX`，把 case 的 `GNUmakefile` 最后一行改成：

```makefile
include $(PELE_HOME)/Exec/Make.PeleLMeX
```

不要使用：

```makefile
include $(PELE_HOME)/Exec/Make/PeleLMeX.mak
```

修改后重新编译：

```bash
make TPL
make -j 8
```

### 找不到 `prob_parm.H`

如果编译时出现：

```text
./pelelmex_prob.H:10:10: fatal error: prob_parm.H: No such file or directory
```

说明 case 的 `pelelmex_prob.H` 里错误地包含了不存在的旧头文件：

```cpp
#include "prob_parm.H"
```

当前 PeleLMeX case 不需要这一行。`pelelmex_prob.H` 是由 PeleLMeX 框架头文件包含的，`ProbParmDefault` 和 `DefaultProblemSpecificFunctions` 由框架提供。

删除这一行，保留类似：

```cpp
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Math.H>

#include "mechanism.H"
```

修改后清理并重编：

```bash
make realclean
make -j 8
```
