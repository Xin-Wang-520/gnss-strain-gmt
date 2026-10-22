# GNSS 速度场与水平应变率

完整流程：**原始 GNSS 绘图 → blockmedian → gpsgridder → Xu 应变率 → 显示插值绘图**。

当前默认参数为欧亚固定参考系、30′分块、25%保留模态、20′输出网格；研究区70–110°E、25–42°N，输入另含67–113°E、22–45°N缓冲区。

## 环境

需要 Python、NumPy、Matplotlib、packaging，以及可从命令行调用的 GMT。处理方法在 GMT 6.4 上验证过；Python 入口已在本机 Python 3.9 和3.13环境检查。独立验证原始 Xu 脚本时还需要 csh。

在仓库根目录安装 Python 依赖并检查 GMT：

```bash
python3 -m pip install -r requirements.txt
gmt --version
```

若 GMT 不在 PATH，可给 run2、run3、run4 传入 `--gmt /path/to/gmt`。

## 正式运行

进入包含这些脚本的仓库根目录，依次执行：

```bash
python3 run1_plot_gnss.py
python3 run2_blockmedian_eu.py
python3 run3_gpsgridder_eu.py --rerun -v
python3 run4_calculate_strain_eu.py
python3 run5_plot_strain_eu.py
```

| 脚本 | 功能 | 输出位置，相对于 GNSS_analysis/output |
|---|---|---|
| run1_plot_gnss.py | 绘制原始台站和水平速度 | run1/ |
| run2_blockmedian_eu.py | 同一地理分块内对速度和sigma分别取中位数 | run2_blockmedian_eu/bm30m/ |
| run3_gpsgridder_eu.py | 从run2代表点插值速度 | run3_gpsgridder_eu/bm30m_c25_i20m_pixel/ |
| run4_calculate_strain_eu.py | 用Xu的GMT方法计算水平应变率 | run4_strain_eu/bm30m_c25_i20m_pixel/ |
| run5_plot_strain_eu.py | 已有run4网格的双线性显示插值图 | 上述run4目录内的display_bilinear/ |

图件输出为PNG和PDF。run5保留原始应变率数据和色标，只改善图面颜色过渡，不增加实际空间分辨率。

`--rerun` 强制run3重新求解；不加时可自动复用输入和参数完全相同的已有结果。`-v`实时显示GMT日志，同时保存日志文件。

仅查看模态数量：

```bash
python3 run3_gpsgridder_eu.py --show-mode-counts
```

## 文件依赖

不要只复制run1–run5，下面这些辅助文件也需要保留：

| 文件或目录 | 用途 |
|---|---|
| GNSS_analysis/__init__.py | Python公共模块所在的包 |
| GNSS_analysis/gnss_processing.py | blockmedian、GMT调用、速度诊断和绘图公共函数 |
| GNSS_analysis/data/raw/ | 作者v6.4原始GNSS表，正式流程使用两份_eu.dat |
| GNSS_analysis/data/basemap/ | 离线海岸线和行政边界 |
| xu/plot_tensor.csh、plot_shear.csh、plot_dilatation.csh | Xu原始参考代码；run4记录它们的指纹，数值检查独立执行其计算部分 |
| GNSS_analysis/parameter_tests.py | 三个test入口共用的参数测试代码 |

run2复用run1的数据读取函数，run3读取run2的默认输出位置，run5复用run4的绘图函数；因此主脚本也需要一起保留。

这份目录已包含必要原始数据、底图和参考脚本，运行时根据仓库所在位置确定路径，不依赖之前的项目目录。

## 可选参数测试

```bash
python3 test_run1_blockmedian_sizes.py --values 15 20 30
python3 test_run2_gpsgridder_percentages.py --values 100 75 50 25
python3 test_run3_gpsgridder_spacings.py --values 5 10 15 20
```

每个测试都自动执行run2→run3→run4，一次只改变指定参数。结果单独写入 `GNSS_analysis/output/tests/`，不覆盖正式结果。加 `--workers 1` 可串行运行，加 `--compare-only` 只重新生成已有结果的对比图。

run2的 `--spacing` 控制分块，run3的 `--spacing` 控制输出网格；`--retain-percent` 是保留模态数量比例，不是台站比例，也不是解释方差比例。1685个代表点对应3370个总模态，目前25%实际保留842个。

## 数值检查与数据来源

完成正式流程后，可运行 `GNSS_analysis/checks/` 下的三个检查脚本。数据来源和应变率定义见 [GNSS_analysis/README.md](GNSS_analysis/README.md)，参考材料的归属说明见 [SOURCES.md](SOURCES.md)。

输出、缓存和本地安装依赖已由 `.gitignore` 排除；原始输入数据随仓库保留，运行结果在本地生成。
