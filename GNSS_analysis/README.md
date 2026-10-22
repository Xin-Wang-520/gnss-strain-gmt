# GNSS 数据、输出和数值定义

运行命令见项目根目录的 [README.md](../README.md)。

## 当前默认数据与数量

读取 `data/raw/` 内作者 AHB v6.4 两份 `_eu.dat` 的纯 GNSS 记录，使用2D和3D的水平分量，参考系为欧亚固定；没有额外删除高不确定度台站。`run1_plot_gnss_itrf14.py` 仍可查看 ITRF2014 数据。

| 数量 | 当前值 |
|---|---:|
| 原始输入，含缓冲区67–113°E、22–45°N | 3600 |
| 原始台站，研究区70–110°E、25–42°N | 2571 |
| 30′ blockmedian后，含缓冲区 | 1685 |
| blockmedian代表点，研究区内 | 1167 |
| gpsgridder总模态 | 3370 |
| 保留25%时实际使用模态 | 842 |
| 20′输出网格 | 120列 × 51行＝6120个位置 |

每个输出位置有VE、VN两个速度值；网格位置不是新增的实测台站。

## run2：blockmedian

只做预处理。使用相同地理分块，对VE、VN、sigmaE、sigmaN分别做未加权中位数计算，坐标采用GMT默认经度和纬度中位数，检查四组坐标对齐后合并六列。每块sigma中位数是代表性误差尺度，不是严格误差传播，不除以台站数平方根。

输出在 `output/run2_blockmedian_eu/bm30m/`：

- `gnss_before_blockmedian_eu.txt`：分块前六列数据。
- `gnss_blockmedian_eu.txt`：代表点，列为lon、lat、VE、VN、sigmaE、sigmaN；速度和sigma单位mm/yr。
- `blockmedian_counts.txt`：每块的台站数。
- `blockmedian_{ve,vn,sigma_e,sigma_n}.txt`：四个分量的中位数。
- `blockmedian_summary.json`、`blockmedian_commands.txt`、`run2_summary.json`：参数、命令、计数与指纹。

## run3：gpsgridder

只读取run2已准备的代表点，不在内部重复分块。默认核心参数为 `-R70/110/25/42 -I20m -fg -W -r -Cn25%`，沿用用户run6.1方法，不显式指定-S或-F。本机GMT6.4默认为0.5和-Ff0.01。默认不加距离掩膜或速度滤波。

输出在 `output/run3_gpsgridder_eu/bm30m_c25_i20m_pixel/`：

- `velocity_eu_u.nc`、`velocity_eu_v.nc`：东向、北向速度，mm/yr，像元注册。
- `velocity_eu.npz`：速度数组。
- `run3_velocity_fit_eu.png/.pdf`：速度及代表点拟合残差图。
- `gpsgridder.log`、`gmt_station_misfit.txt`：GMT日志及直接预测残差。
- `mode_summary.json`：实际保留模态数。
- `run_config.json`、`fit_summary.json`、`block_fit_eu.txt`：参数、指纹与拟合误差。

-I控制输出采样，-Cn控制求解模态数量；后续有限差分求导受输出间隔影响。

## run4：Xu应变率

使用用户提供的 `xu/plot_tensor.csh`、`plot_shear.csh`、`plot_dilatation.csh` 的GMT数值方法。速度由mm/yr转m/yr，以 `grdgradient -A270/-A180 -fg` 求梯度，距离尺度半径6378km；不加球面曲率项、速度滤波或地形掩膜。

```text
EE = du/dx
NN = dv/dy
EN = (du/dy + dv/dx)/2
最大剪切率 = sqrt(EN² + (EE-NN)²/4) = abs(e1-e2)/2
面膨胀率 = EE + NN
second_invariant = sqrt(EE² + NN² + 2EN²)
```

最后一项沿用Xu/RSE命名，数值是张量Frobenius范数；EN是张量剪切，工程剪切为2EN。输出单位nstrain/yr。

输出在 `output/run4_strain_eu/bm30m_c25_i20m_pixel/`：各应变率NetCDF、`strain_eu.npz`、`run4_summary.json`、`run4_strain_eu.png/.pdf`、`run4_shear_contributions_eu.png/.pdf`。`xu_intermediate/`保留GMT命令和中间文件。

`run5_plot_strain_eu.py` 只读取已保存数据，把双线性显示插值图另存至 `display_bilinear/`；色标仍来自原始网格，并校验原始文件没有改写。显示插值不增加观测信息或空间分辨率。

## 数值检查

在项目根目录运行：

```bash
python3 GNSS_analysis/checks/check_blockmedian.py
python3 GNSS_analysis/checks/check_mode_logging.py
python3 GNSS_analysis/checks/check_strain.py
```

验证已知分块中位数、坐标对齐、sigma、单点块及计数；用真实GMT小样本检查模态取整和日志；独立执行三份Xu原始计算代码块对比六个应变率分量，并以解析速度场和独立特征值计算验证导数符号、单位和剪切1/2系数。

原项目整理时重新运行了一次完整三步流程，全部应变率数组与整理前选定结果完全一致。本目录仅收录运行所需的代码、输入数据及参考脚本，不收录原项目的历史计算结果。

参数测试在 `output/tests/`。比较统计排除地图四周0.5°，使用原生网格采样。拟合残差和图面光滑程度不属于独立验证；目前未传播完整插值及应变率不确定度。
