# 输入与参考代码来源

## GNSS

`GNSS_analysis/data/raw/` 为 AHB GPS 数据库v6.4的原始文件：

- 来源仓库：[earjcr1/AHB_GPS](https://github.com/earjcr1/AHB_GPS)
- 版本目录：`ahbgps_v6pt4/`
- 文件日期：29-Nov-2025；包含2D、3D的欧亚固定和ITRF2014参考系文件。
- 本流程仅选GNSS记录，使用水平速度及其sigma。
- 原始文件的Git blob指纹及下载来源记录在 `data/manifest.json` 和 `data/itrf14_sources.json`。

这些输入数据的来源和归属不因收录到本流程而改变。

## 底图

`GNSS_analysis/data/basemap/` 的海岸线和边界为原项目通过GMT 6.4/GSHHG提取的区域底图文本，仅用于地图显示。

## Xu参考脚本

`xu/plot_tensor.csh`、`xu/plot_shear.csh`、`xu/plot_dilatation.csh` 为用户提供的Xu原始参考代码，随包保留原始内容，用于追踪应变率公式及独立核验。相关参考文章为Xu等2021年《Integrated Sentinel-1 InSAR and GNSS Time-Series Along the San Andreas Fault System》。

本流程保留这些参考脚本的归属，不对第三方代码或数据另行声明新的许可证。
