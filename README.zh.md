# UPD 计算工具

面向 NIPT SNP readlist 的独立 UPD (Uniparental Disomy) 分析工具。

[English](README.md)

## 特性

- **自包含**：仅需 5 个第三方库
- **FF 稳健**：胎儿浓度双轨估计，避免全基因组纯合样本因 FF 高估而漏检
- **批量并行**：支持多进程并行处理整批样本

## 系统要求

### 软件依赖

| 组件 | 版本要求 | 说明 |
|---|---|---|
| Python | 3.9 – 3.11 | 已在 3.9 / 3.10 / 3.11 上测试 |
| numpy | >=1.20 | 数值计算 |
| pandas | >=1.3 | readlist 解析与分组统计 |
| scipy | >=1.7 | `scipy.stats.binom` 似然计算 |
| hmmlearn | >=0.3.0 | 必须提供 `CategoricalHMM`（0.2.x 无此类） |
| matplotlib | >=3.4 | 仅绘图子功能需要 |

版本清单见 `upd_tool/requirements.txt`。

### 操作系统

- Linux（CentOS 7 / Ubuntu 20.04、22.04）：生产环境使用的平台
- Windows 10/11、macOS 12+：可运行本文档的 demo；绘图强制使用
  `matplotlib` 的 `Agg` 后端，无需图形界面

### 硬件要求

无特殊硬件，不需要 GPU。代码内部将 BLAS 线程数固定为 1，`batch` 子命令通过
多进程按样本并行。单样本峰值内存约数百 MB（13833 行探针量级）。

## 安装

```bash
git clone https://github.com/biobiggen/upd
cd upd
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r upd_tool/requirements.txt
```

纯 Python 实现，无需编译。在普通台式机上安装耗时主要来自下载
numpy/pandas/scipy 的 wheel 包，通常几分钟内完成。

## Demo

真实血浆数据无法公开分发，因此提供一个模拟数据生成器，输出格式与生产环境
完全一致。

### 1. 生成模拟演示数据

```bash
python -m upd_tool.simulate_demo_data -o demo_data
```

生成内容（默认胎儿浓度 0.10、平均深度 1200×、每个目标区域 120 个 SNP，
共约 1200 个位点）：

| 文件 | 说明 |
|---|---|
| `demo_data/demo_probe.bed` | 探针文件（`#Chr` / `Pos` / `Type`） |
| `demo_data/demoN01..N04_normal_*.reads.list` | 4 个正常二体样本 |
| `demo_data/demoP01_updm_*.reads.list` | `15q11q13` 母源 UPD |
| `demo_data/demoP02_updpi_*.reads.list` | `11p15` 父源异二体 |
| `demo_data/demoP03_updpii_*.reads.list` | `11p15` 父源同二体 |

readlist 列格式：

```
#Chr  Pos  Ref  Alt  Depth  Ref_Dep  Alt_Dep  pA_Ratio  GC  DepRegion  MapQ
```

### 2. 运行单样本计算

```bash
python -m upd_tool.cli single \
    -i demo_data/demoP01_updm_consensus.mapped.clipped.snp.reads.list \
    -o demo_out/demoP01_upd_results.json \
    -p demo_data/demo_probe.bed \
    --plot demo_out/demoP01_upd_regions_pa_ratio.png
```

预期输出：

- 标准输出打印每个区域的 `final_state` 与占比，以及胎儿浓度诊断信息
- `demo_out/demoP01_upd_results.json`：各区域完整结果（见「输出格式」）
- `demo_out/demoP01_upd_regions_pa_ratio.png`：pA_Ratio 散点图

对该模拟样本，`15q11q13` 的 `final_state` 预期为 `UPDM`，其余区域为
`Normal`，`FF_Method` 为 `homozygous`（比值 ≈ 1）。

> 注：模拟数据的位点数远少于真实探针（13833 行），故不要调用
> `validate_row_count()`；`--probe-version` 仅影响 `ignore_snps` 与阈值，
> demo 使用默认的 `NIPT3V4` 即可。

### 3. 批量处理与汇总报告

```bash
python -m upd_tool.cli batch \
    -i demo_data/ -o demo_out/ -p demo_data/demo_probe.bed \
    --threads 4 --plot-dir demo_out/UPD_image/ -v

python -m upd_tool.cli report --results demo_out/ -o demo_out/upd_report.csv
```

预期输出：`demo_out/` 下每个样本一个 `*_upd_results.json`，
`demo_out/UPD_image/` 下每个样本一张 PNG，以及汇总的
`demo_out/upd_report.csv`（每行一个样本，每个区域两列：状态与占比）。

### 预期运行时间

在普通台式机（4 核 CPU、16 GB 内存）上，整个 demo 预期在数分钟内完成：
数据生成约数秒；单样本计算的耗时主要来自
`predict_fetal_genotype_hypergeom` 对位点的逐行循环，demo 规模下约数秒，
真实 13833 行探针约 1–2 分钟；`batch` 按 `--threads` 线性加速。

## 使用

### 单样本计算

```bash
python -m upd_tool.cli single \
    -i sample_consensus.mapped.clipped.snp.reads.list \
    -o result.json \
    -p /path/to/NIPT3V4_CNV-HAP_targets_hg38.bed \
    --probe-version NIPT3V4
```

同时生成散点图：

```bash
python -m upd_tool.cli single -i sample.readslist -o result.json \
    -p probe.bed --plot sample_upd_regions_pa_ratio.png
```

### 批量处理

```bash
python -m upd_tool.cli batch \
    -i ./readslist_dir/ \
    -o ./results/ \
    -p /path/to/probe.bed \
    --threads 8 \
    --plot-dir ./UPD_image/ \
    -v
```

参数说明：

| 参数 | 说明 |
|---|---|
| `-i, --input` | 包含 `*.snp.reads.list` 文件的目录 |
| `-o, --output` | JSON 结果输出目录 |
| `-p, --probe-file` | 探针文件路径 |
| `--probe-version` | `NIPT3V3` / `NIPT3V4`（默认 `NIPT3V4`） |
| `--threads` | 并行进程数（默认 4） |
| `--plot-dir` | 可选，同时生成散点图的输出目录 |
| `--no-recursive` | 不递归扫描子目录 |
| `-v, --verbose` | 显示每个样本的处理结果 |

### 生成报告

```bash
# 从目录扫描
python -m upd_tool.cli report --results ./results/ -o upd_report.csv

# 或指定文件列表
python -m upd_tool.cli report --upd-jsons a.json b.json -o upd_report.csv
```

## 在自有数据上运行

1. **准备 readlist**：由上游比对流程为每个样本生成一个
   `*.snp.reads.list`（制表符分隔）。必需列为 `#Chr`、`Pos`、`Ref`、`Alt`、
   `Depth`、`pA_Ratio`；若存在 `DepRegionGC` 列会优先用于 chrY 深度比值，
   否则回退到 `Depth`。第一列必须是 `#Chr`。
2. **准备探针文件**：`.bed` / `.xls`（制表符分隔）或 `.csv`，必需列
   `#Chr`、`Pos`，以及类型列（识别 `Type` / `SNP_Tag` / `Probe_Type` 等
   常见名称）。类型列的取值需与 `regions.py` 中的区域名一致
   （`6q24`、`7q32`、`11p15`、`14q32`、`15q11q13`、`20q13`、`chrRef`）。
   未在探针文件中出现的位点会被标为 `chrRef`。
3. **选择探针版本**：`--probe-version NIPT3V3` 或 `NIPT3V4`，决定
   `expected_rows`、`ff_threshold`、`depth_threshold` 与被忽略的探针类型
   （`HBA` / `RHD` / `SMN` / `HAP` / `other`）。若上游探针 panel 不同，
   在 `core.PROBE_VERSIONS` 中新增一项即可。
4. **修改目标区域**：若关注的印记区域与默认不同，编辑 `regions.py` 的
   `UPD_TARGETS_HMM`（HMM 用窄区间）与 `UPD_TARGETS_PLOT`（绘图用宽区间），
   并同步 `REPORT_REGIONS`、`PLOT_ORDER`。
5. **批量运行并汇总**：先 `batch` 再 `report`，见上文命令。建议每批样本放在
   同一目录，`--threads` 不超过物理核心数。
6. **可调参数**：判定阈值集中在 `core.py` 顶部常量
   （`MIN_REGION_SNPS`、`MIN_HOMOZYGOUS_SNPS`、`MIN_OBSERVATIONS`、
   `MIN_REGION_LENGTH`、`DEPTH_FILTER`、`NORMAL_RATIO_THRESHOLD`、
   `SIGNIFICANT_UPD_RATIO`），HMM 参数为 `STARTPROB` / `TRANSMAT` /
   `EMISSIONPROB`。低深度或低胎儿浓度样本可适当下调 `DEPTH_FILTER`。
7. **质控**：胎儿浓度低于 `ff_threshold`（0.03）时 UPD 判定不可靠；JSON 中
   区域的 `status` 字段说明该区域为何未能计算（见「输出格式」）。

## 编程接口

```python
from upd_tool import UPDCalculator

calc = UPDCalculator(probe_version='NIPT3V4', probe_file='probe.bed')
calc.load_readlist('sample.snp.reads.list')

ff = calc.get_fetal_fraction()      # 胎儿浓度（已含双轨校正）
results = calc.predict_upd_hmm()    # UPD 预测
print(calc.ff_info)                 # FF 诊断信息（两条路径估计值与比值）
```

便捷函数：

```python
from upd_tool.core import calculate_upd

results = calculate_upd('sample.readslist', probe_file='probe.bed')

# 需要 FF 诊断信息时
results = calculate_upd('sample.readslist', probe_file='probe.bed',
                        with_ff_info=True)
print(results['_ff_info'])
```

## 模块结构

| 文件 | 说明 |
|---|---|
| `core.py` | `UPDCalculator` 类，UPD 计算核心逻辑 |
| `regions.py` | UPD 目标区域坐标定义 |
| `plotting.py` | pA_Ratio 散点图绘制 |
| `cli.py` | 命令行入口 |
| `__main__.py` | 使 `python -m upd_tool` 等价于 `python -m upd_tool.cli` |
| `simulate_demo_data.py` | 小型模拟演示数据集生成器 |
| `requirements.txt` | 依赖清单 |

## 算法说明

### 计算流程

```
readlist 加载
  → get_probe_type       标注 SNP_Tag
  → get_background_gt    推断母亲基因型
  → get_fetal_fraction   计算胎儿浓度（双轨估计 + 自动切换）
  → predict_fetal_genotype_hypergeom   二项似然预测胎儿基因型
  → predict_upd_hmm      HMM 状态预测
```

### 胎儿浓度的双轨估计

胎儿浓度（FF）用两条互相独立的路径估计，并做交叉校验：

**路径 1 — 母亲纯合位点**（`ffAB`）

母亲 BB 位点上血浆 alt 比例的期望为 `AF = ff * 胎儿alt剂量 / 2`。取
`ff = 2 * median(AF)` 隐含了**「胎儿在该位点为杂合」**的假设（即胎儿从父方
获得一个 alt 等位基因，剂量为 1）。

**路径 2 — 母亲杂合（BA）位点**

母亲 BA 位点上：

```
AF = (1 - ff) * 0.5 + ff * 胎儿alt剂量 / 2
   = 0.5 + (ff / 2) * (胎儿alt剂量 - 1)
```

即 `|AF - 0.5| = (ff / 2) * |胎儿alt剂量 - 1|`。胎儿纯合（剂量 0 或 2）时
偏移量恒为 `ff / 2`，取偏移中位数乘 2 即得 FF。**该关系与胎儿基因组来自
双亲还是全部来自单方无关**，因此对全基因组纯合稳健。

**为何全基因组纯合样本必须用路径 2**

若胎儿全基因组纯合，母亲纯合位点上胎儿也是纯合
（AF = `ff` 而非 `ff/2`），故路径 1 会把 FF 系统性**高估约 2 倍**。FF 翻倍后，
`expected_alt_ratio(BB→BA) = ff_est/2` 恰好等于真实观测 AF，导致胎儿的
**AA 被误判为 BA**，观测 `BBAA` 退化为 `BBBA`。而 `BBAA`/`AABB` 是父源同二体
在发射矩阵中的唯一签名（`Normal` 状态下概率为 0），签名一旦消失，区域即被
判为 `Normal`，父源同二体无法检出。

因此当 `ff_hom / ff_het > 1.5`（`FF_RATIO_HOM_UPPER`）时，提示存在全基因组
纯合，自动改用路径 2 的估计值。普通样本两条路径结果一致（比值 ≈ 1），
不受影响；局部 UPD（仅单个区域纯合）也不会触发切换，因为 FF 是全基因组统计量。

两条路径的估计值与比值一并写入结果的 `_ff_info` 字段，可用于复核。

### HMM 状态

| 状态 | 含义 |
|---|---|
| `Normal` | 正常双亲二体 |
| `UPDM` | 母源单亲二体 |
| `UPDPI` | 父源单亲二体（异二体） |
| `UPDPII` | 父源单亲二体（同二体） |

### 目标区域

工具使用两组区域坐标：

- **HMM 计算**（窄区间）：`6q24`、`7q32`、`11p15`、`14q32`、`15q11q13`、`20q13`、`chrRef`
- **绘图**（宽区间）：前 6 个区域（不含 `chrRef`）

### 判定阈值

| 阈值 | 值 | 说明 |
|---|---|---|
| 区域最少 SNP 数 | 50 | 低于则 `insufficient_data` |
| 母亲纯合位点最少数 | 20 | 低于则 `insufficient_homozygous` |
| 深度过滤 | 400 | 仅使用深度 ≥400 的位点 |
| 分段最小长度 | 20 | 短于此长度的分段被丢弃 |
| Normal 判定阈值 | 0.2 | Normal 比例 >0.2 即判为正常 |
| 显著 UPD 比例 | 0.1 | 非 Normal 分段比例 >0.1 计入 `significant_upds` |
| BA 路径最少位点数 | 30 | 低于则不做 FF 交叉校验 |
| FF 比值切换阈值 | 1.5 | `ff_hom/ff_het` 高于此值改用 BA 路径估计 |
| BA 偏移下限 | 0.02 | `\|AF-0.5\|` 低于此值视为胎儿杂合，不参与 FF 估计 |

## 输出格式

### JSON 结果

```json
{
  "6q24": {
    "status": "success",
    "state_ratios": {"Normal": 0.95, "UPDM": 0.05},
    "regions": [
      {"start": 142448249, "end": 145502506, "state": "Normal",
       "length": 120, "ratio": 0.95}
    ],
    "significant_upds": [],
    "final_state": "Normal",
    "final_ratio": 0.95,
    "total_observations": 126,
    "chromosome": "chr6",
    "start_pos": 142448249,
    "end_pos": 145502506,
    "site_details": [...]
  },
  "_ff_info": {
    "ff_used": 0.1012,
    "ff_homozygous": 0.1008,
    "ff_heterozygous": 0.1015,
    "ff_ratio": 0.99,
    "ff_method": "homozygous",
    "het_sites": 412,
    "het_shifted_sites": 208
  }
}
```

除各区域外，顶层还有一个 `_ff_info` 键（下划线前缀，与区域名区分）记录
胎儿浓度诊断信息：

| 字段 | 说明 |
|---|---|
| `ff_used` | 实际用于下游计算的 FF |
| `ff_homozygous` | 母亲纯合位点路径的估计值 |
| `ff_heterozygous` | 母亲杂合（BA）位点路径的估计值；位点不足时为 `null` |
| `ff_ratio` | `ff_homozygous / ff_heterozygous`，接近 2 提示全基因组纯合 |
| `ff_method` | `homozygous`（默认）或 `heterozygous`（已切换校正） |
| `het_sites` | 通过深度过滤的母亲 BA 位点数 |
| `het_shifted_sites` | 其中偏移量超过下限、参与 FF 估计的位点数 |

`status` 可能取值：

| 值 | 说明 |
|---|---|
| `success` | 计算成功 |
| `insufficient_data` | 区域 SNP 数 <50 |
| `insufficient_homozygous` | 母亲纯合位点 <20 |
| `insufficient_observations` | 有效观测 <20 |
| `no_valid_observations` | 无法映射的观测 |
| `insufficient_observation_types` | 观测类型 <2 |
| `observation_conversion_error` | 观测转换失败 |
| `hmm_error` | HMM 预测失败 |

### CSV 报告

各区域两列（状态、占比）之后，追加三列 FF 诊断信息：

```csv
Sample_Name,6q24_UPD,6q24_Ratio,...,chrRef_Ratio,FF_Used,FF_Method,FF_Ratio
A346_cf1439_2516654P1DE,Normal,0.9500,...,0.9600,0.1012,homozygous,0.99
demoP03,UPDPII,0.9800,...,0.9700,0.1015,homozygous,1.01
```

`FF_Method` 为 `heterozygous` 的样本即触发了 FF 校正，配合 `FF_Ratio`
（接近 2）可快速筛出疑似全基因组纯合的样本。

## 许可

本软件采用 **GNU General Public License v3.0**（OSI 批准的开源协议）发布，
完整协议文本见仓库根目录的 `COPYING` 文件。

Copyright (C) 2024 biobiggen

本程序是自由软件：你可以依据自由软件基金会发布的 GNU 通用公共许可证第 3 版
（或你选择的任何更新版本）的条款，重新分发和/或修改本程序。

分发本程序的目的是希望它有用，但**不提供任何担保**，甚至不包含对适销性或
特定用途适用性的默示担保。详见 GNU 通用公共许可证。

你应当已随本程序收到一份 GNU 通用公共许可证的副本；如果没有，请查阅
<https://www.gnu.org/licenses/>。

### 商业授权

GPL-3.0 属于 copyleft 协议：基于本代码的衍生作品在分发时必须同样以 GPL-3.0
开源。若需在闭源产品中集成本软件，或需要不受 copyleft 约束的授权条款，
请与作者联系商谈单独的商业许可。

算法的完整描述（伪代码）见论文 **Methods** 部分。
