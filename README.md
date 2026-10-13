# 测绘程序设计

武汉大学 / 中国地质大学（北京）测绘工程 MATLAB 课程作业集合。

## 项目列表

### 卫星轨迹绘制（GNSS Broadcast Ephemeris Processing）

使用 RINEX 广播星历文件（N 文件）和观测文件（O 文件），完成：

1. **读取 N 文件**：解析 GPS 广播星历，输出 33×N 星历矩阵
2. **卫星轨道解算**：按 IS-GPS-200 开普勒轨道模型，计算卫星在 ECEF 地固系下的位置
3. **轨迹可视化**：绘制卫星一天内的三维轨道和星下点轨迹
4. **接收机单点定位（SPP）**：利用伪距观测值，最小二乘解算接收机坐标，输出位置散点图

#### 运行环境

- MATLAB R2020b 或更高版本
- 无需额外工具箱

#### 文件说明

| 文件 | 功能 |
|---|---|
| `main_track.m` | 卫星轨迹主脚本（画三维轨道 + 星下点） |
| `plot_circular_track.m` | ECI 惯性系下的圆形轨道对比图 |
| `spp_main.m` | 接收机单点定位主脚本（输出散点图） |
| `Rinex_get_nav.m` | 读取 RINEX N 文件广播星历 |
| `read_obs_gps.m` | 读取 RINEX O 文件 GPS C1C 伪距 |
| `satellite_orbits.m` | 开普勒轨道解算卫星位置 |
| `find_eph.m` | 按时间选最近的星历列 |
| `topocent.m` | 计算方位角/高度角 |
| `cart2geod.m` | ECEF 直角坐标 → 大地经纬度 |
| `sat_clock_error_correction.m` | 卫星钟差改正 |
| `relativistic_error_correction.m` | 相对论效应改正 |
| `earth_rotation_correction.m` | 地球自转（Sagnac）改正 |
| `iono_error_correction.m` | Klobuchar 电离层延迟模型 |
| `glp_ini.m` | 全局物理常数初始化 |

#### 数据文件

仓库**不包含** RINEX 数据文件（体积过大），请从 IGS 官网下载：
- 广播星历：https://cddis.nasa.gov/archive/gnss/products/
- 观测文件：https://cddis.nasa.gov/archive/gnss/data/daily/

将 `.rnx` 和 `.26o` 文件放在项目根目录即可。

#### 运行方法

```matlab
% 画卫星轨迹
main_track

% 画 ECI 惯性系对比图
plot_circular_track

% 接收机单点定位
spp_main
```

## 作者

Ezra240601
