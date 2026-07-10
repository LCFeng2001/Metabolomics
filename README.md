# Targeted Metabolomics Rebuild

这是一个用于 **Skyline 导出的广靶代谢组原始 Peak Area 数据** 的完整 R 分析项目。

项目从 Skyline 原始峰面积表出发，自动完成数据清洗、重复代谢物合并、Blank 过滤、缺失值处理、归一化、log2 转换、Pareto scaling、PCA、limma 差异分析、OPLS-DA/VIP、火山图、boxplot、heatmap 以及 KEGG 分析框架。

---

## 1. 项目功能

本项目可以完成：

- 自动读取 Skyline 导出的 Excel / CSV / TSV / TXT 文件
- 自动识别样本列、blank 列和代谢物注释列
- 删除无用列，例如 `序号`、质荷比、`均值` 列
- 将 0 视为缺失值
- 按 `Compound name` 合并重复代谢物或多个 transition
- Blank 过滤
- 缺失值过滤与填补
- 中位数归一化
- log2 转换
- Pareto scaling
- PCA 分析
- 自动生成不同 genotype、不同时间点下的比较
- limma 差异分析
- OPLS-DA / VIP 分析
- ropls 失败时自动生成 fallback VIP-like ranking
- 火山图，VIP 值作为点大小
- 关注代谢物类别标注
- 显著关注代谢物 boxplot
- total heatmap
- differential metabolite heatmap
- heatmap 图与代谢物、样本、分组的对应表
- KEGG 注释模板和 KEGG 富集分析框架

---

## 2. 输入数据格式

原始数据应放在：

```text
data/raw/
```

推荐文件：

```text
data/raw/skyline_peak_area.xlsx
```

然后在：

```text
config/config.R
```

中修改：

```r
raw_file <- file.path(project_dir, "data/raw/skyline_peak_area.xlsx")
raw_sheet <- 1
```

如果是 TSV 文件：

```r
raw_file <- file.path(project_dir, "data/raw/skyline_peak_area.tsv")
```

---

## 3. 原始表格列名要求

典型 Skyline 导出表格应包含：

```text
序号
分子
Compound name
Class
sample1
sample2
...
```

关键列说明：

| 列名 | 说明 |
|---|---|
| `Compound name` | 代谢物名称，默认用于合并重复 transition |
| `分子` | 原始分子 ID 或 transition ID |
| `Class` | 代谢物类别 |

---

## 4. 样本命名规则

默认样本格式为：

```text
Genotype-Treatment Time-Replicate
```

```

---

## 5. 项目结构

```text
Targeted_Metabolomics_Rebuild/
├── README.md
├── config/
│   └── config.R
├── data/
│   ├── raw/
│   └── annotation/
├── functions/
│   ├── io/
│   ├── preprocess/
│   ├── stats/
│   ├── plots/
│   ├── pathway/
│   └── utils/
├── scripts/
│   ├── 00_install_packages.R
│   ├── 00B_Make_KEGG_Annotation_Candidates.R
│   ├── 01_ReadCleanAggregate.R
│   ├── 02_Preprocess.R
│   ├── 03_PCA.R
│   ├── 04_Comparisons.R
│   ├── 05_Limma.R
│   ├── 06_OPLSDA.R
│   ├── 07_Plots.R
│   ├── 08_KEGG.R
│   └── 10_RunAll.R
└── results/
    ├── rds/
    ├── tables/
    ├── pca/
    ├── limma/
    ├── oplsda/
    ├── volcano/
    ├── heatmap/
    ├── boxplot/
    └── kegg/
```

---

## 6. 快速运行

首次运行：

```r
source("scripts/00_install_packages.R")
```

完整运行：

```r
rm(list = ls())
source("scripts/10_RunAll.R")
```

只重新生成图：

```r
rm(list = ls())
source("scripts/07_Plots.R")
```

重新生成 OPLS-DA/VIP 和图：

```r
rm(list = ls())
source("scripts/06_OPLSDA.R")
source("scripts/07_Plots.R")
```

---

## 7. 核心参数

所有主要参数位于：

```text
config/config.R
```

### 7.1 缺失值

Skyline 空值通常为 0，本项目默认将 0 当作缺失值：

```r
zero_as_missing <- TRUE
```

缺失值填补：

```r
impute_method <- "half_min"
```

### 7.2 重复代谢物合并

默认按 `Compound name` 合并：

```r
aggregate_by <- "Compound name"
```

合并后：

- 样本峰面积求和
- `Class` 合并保留
- 原始 `分子` ID 存入 `source_features`

### 7.3 归一化

```r
normalization_method <- "median"
```

可选：

```r
"none"
"TIC"
"median"
```

### 7.4 log2 转换和 scaling

```r
log_base <- 2
pseudo_count <- 1
scale_method <- "pareto"
```

heatmap 当前使用：

```r
heatmap_use_matrix <- "scaled"
```

即：

```r
data02$expr_scaled
```

由于 `scale_method <- "pareto"`，heatmap 使用 Pareto-scaled matrix。

### 7.5 差异分析阈值

```r
logFC_cutoff <- 1
pvalue_cutoff <- 0.05
volcano_stat_col <- "P.Value"
```

火山图默认使用：

```text
P.Value < 0.05
|log2FC| >= 1
```

---

## 8. 脚本说明

### 8.1 `01_ReadCleanAggregate.R`

功能：

- 读取 Skyline 原始数据
- 删除无用列
- 自动识别 blank 和样本列
- 将 0 转为 NA
- 按 `Compound name` 合并重复代谢物
- 生成样本信息表

输出：

```text
results/tables/01_ReadCleanAggregate_output.xlsx
results/rds/01_clean_aggregated.rds
```

---

### 8.2 `02_Preprocess.R`

功能：

- Blank 过滤
- 缺失值过滤
- 缺失值填补
- 归一化
- log2 转换
- Pareto scaling

输出：

```text
results/tables/02_Preprocess_output.xlsx
results/rds/02_preprocessed.rds
```

重要 sheet：

```text
normalized_area
log2_expression
scaled_expression
blank_filter_report
missing_filter_report
```

---

### 8.3 `03_PCA.R`

功能：

- 基于 scaled matrix 做 PCA
- 输出 PCA 图和 PCA 结果表

输出：

```text
results/pca/
results/tables/03_PCA_output.xlsx
```

---

### 8.4 `04_Comparisons.R`

自动生成比较组。

默认生成每个 genotype 在指定时间点下的：

输出：

```text
results/tables/04_Comparisons_output.xlsx
results/rds/04_comparisons.rds
```

---

### 8.5 `05_Limma.R`

功能：

- 使用 limma 做差异分析
- 输出每个 comparison 的差异代谢物结果
- 输出关注类别代谢物结果

输出：

```text
results/tables/05_Limma_output.xlsx
results/rds/05_limma_results.rds
```

---

### 8.6 `06_OPLSDA.R`

功能：

- 优先使用 `ropls::opls()` 做 OPLS-DA
- 如果 ropls 失败，自动使用 fallback VIP-like ranking
- 每个 comparison 输出 VIP 表和 VIP 图

输出：

```text
results/tables/06_OPLSDA_output.xlsx
results/rds/06_oplsda_results.rds
results/oplsda/
```

每个 comparison 的文件夹中包括：

```text
*_OPLSDA_VIP_summary.xlsx
*_OPLSDA_summary.csv
*_VIP.csv
*_VIP_plot.pdf
*_VIP_plot.png
*_VIP_plot.svg
```

如果 ropls 真正建模成功，还会输出：

```text
*_ropls_diagnostic_plots.pdf
```

如果没有 ropls model，会输出：

```text
*_NO_ropls_model_note.txt
```

---

### 8.7 `07_Plots.R`

功能：

- 绘制 total heatmap
- 绘制 differential metabolite heatmap
- 绘制 volcano
- 绘制 boxplot
- 导出 heatmap 对应表格

输出：

```text
results/heatmap/
results/volcano/
results/boxplot/
```

---

### 8.8 `08_KEGG.R`

功能：

- 读取 KEGG 注释表
- 进行 KEGG 富集分析
- 如果没有 KEGG 注释表，则自动生成模板并跳过分析

需要提供：

```text
data/annotation/kegg_annotation.xlsx
```

必须包含列：

```text
Compound name
KEGG
```

示例：

| Compound name | KEGG |
|---|---|
| L-Serine | C00065 |
| L-Proline | C00148 |
| Uracil | C00106 |

---

## 9. 主要输出文件

### 9.1 数据清洗与预处理

```text
results/tables/01_ReadCleanAggregate_output.xlsx
results/tables/02_Preprocess_output.xlsx
```

其中 `02_Preprocess_output.xlsx` 最重要。

重点 sheet：

```text
normalized_area
log2_expression
scaled_expression
```

---

### 9.2 差异分析

```text
results/tables/05_Limma_output.xlsx
```

---

### 9.3 OPLS-DA / VIP

总表：

```text
results/tables/06_OPLSDA_output.xlsx
```

分 comparison 输出：

```text
results/oplsda/<comparison>/
```


---

### 9.4 火山图

```text
results/volcano/
```

颜色含义：

| 类别 | 颜色 |
|---|---|
| Up | 红色 |
| Down | 蓝色 |
| NS | 灰色 |

点大小表示 VIP value。

---

### 9.5 Boxplot

```text
results/boxplot/
```

boxplot 展示显著且属于关注类别的代谢物。

---

## 10. Heatmap 输出说明

heatmap 输出在：

```text
results/heatmap/
```

分为：

```text
results/heatmap/total/
results/heatmap/differential/
```

### 10.1 Total heatmap

位置：

```text
results/heatmap/total/
```

---

### 10.2 Differential heatmap

位置：

```text
results/heatmap/differential/
```

如果差异代谢物数量不足，会生成：

```text
Differential_heatmap_NOT_DRAWN_matrix_and_annotation.xlsx
```

---

### 10.3 Heatmap Excel 关键 sheet

#### `matrix_plot_order_with_info`

最推荐查看。

它按 heatmap 图中顺序排列，并明确每一行是哪一个代谢物。

包含：

```text
heatmap_row_position
row_slice
metabolite_label
metabolite_id
Class
source_features
```

后面的样本列就是 heatmap 使用的矩阵值。

#### `matrix_with_rownames`

带行名版本。

第一列是：

```text
row_name
```

后面是 heatmap 使用的矩阵值。

#### `heatmap_row_order`

图中每一行对应哪个代谢物。

#### `heatmap_column_order`

图中每一列对应哪个样本。

#### `colgrp_*`

按 heatmap 的列分组拆开的矩阵。

例如：

```text
colgrp_M03
colgrp_D3140
colgrp_ZH6218
```

#### `rowgrp_*`

按代谢物 Class 拆开的矩阵。

例如：

```text
rowgrp_BX
rowgrp_Flavone
rowgrp_Phenolamides
```

---

### 10.4 Heatmap 矩阵和配色

当前 heatmap 使用：

```r
heatmap_use_matrix <- "scaled"
```

即：

```r
data02$expr_scaled
```

由于：

```r
scale_method <- "pareto"
```

所以 heatmap 使用 Pareto-scaled matrix。

颜色设置：

```r
heatmap_clip_value <- 3
```

配色：

```r
c("#2166AC", "#F7F7F7", "#B2182B")
```

含义：

| 数值 | 颜色 |
|---|---|
| 低值 | 蓝色 |
| 0 附近 | 白色 |
| 高值 | 红色 |

---

## 11. KEGG 分析

如果需要 KEGG 分析，需要准备：

```text
data/annotation/kegg_annotation.xlsx
```

文件至少包含：

```text
Compound name
KEGG
```

如果没有 KEGG 注释文件，`08_KEGG.R` 会自动跳过 KEGG 分析，并生成模板文件。

---

## 12. 常见问题

### 12.1 为什么 `results/oplsda/` 为空？

旧版本只有 ropls 成功建模时才保存 ropls 诊断图。  
当前版本已修复：即使 ropls 失败，也会输出 VIP 表、summary 和 VIP 图。

---

### 12.2 为什么没有 differential heatmap？

可能是满足阈值的差异代谢物少于 2 个。

默认阈值：

```r
P.Value < 0.05
abs(logFC) >= 1
```

如果不足，会生成：

```text
Differential_heatmap_NOT_DRAWN_matrix_and_annotation.xlsx
```

---

### 12.3 heatmap 表格没生成怎么办？

当前版本已经加入 CSV 兜底输出。

如果 Excel 没生成，请检查同目录下是否有：

```text
*_matrix_plot_order_with_info.csv
*_matrix_with_rownames.csv
*_heatmap_row_order.csv
*_heatmap_column_order.csv
```

如果 CSV 有而 Excel 没有，通常是 `openxlsx` 没安装。

安装：

```r
install.packages("openxlsx")
```

---

### 12.4 如何查某个代谢物在 heatmap 中的位置？

打开：

```text
matrix_plot_order_with_info
```

或：

```text
*_matrix_plot_order_with_info.csv
```

搜索代谢物名称或 ID。

重点看：

```text
heatmap_row_position
metabolite_label
metabolite_id
Class
```

---

### 12.5 如何查某个 heatmap 色块对应哪个样本？

打开：

```text
heatmap_cell_long
```

或：

```text
*_heatmap_cell_long.csv
```

重点列：

```text
heatmap_row_position
heatmap_column_position
metabolite_id
sample
value
```

---

### 12.6 如何修改火山图点大小？

在：

```text
functions/plots/plot_volcano.R
```

查找：

```r
scale_size_continuous(
  range = c(0.35, 3.2)
)
```

如果想让点更小：

```r
range = c(0.25, 2.5)
```

---

### 12.7 如何修改 heatmap 颜色？

在：

```text
functions/plots/plot_heatmap.R
```

查找：

```r
circlize::colorRamp2(
  brks,
  c("#2166AC", "#F7F7F7", "#B2182B")
)
```

修改颜色即可。

---

## 13. 推荐检查顺序

完整运行后，建议依次检查：

```text
results/tables/01_ReadCleanAggregate_output.xlsx
results/tables/02_Preprocess_output.xlsx
results/tables/04_Comparisons_output.xlsx
results/tables/05_Limma_output.xlsx
results/tables/06_OPLSDA_output.xlsx
results/volcano/
results/boxplot/
results/heatmap/total/
results/heatmap/differential/
```

---

## 14. 注意事项

1. 原始数据中的 0 默认会被视为缺失值。
2. 重复代谢物默认按 `Compound name` 合并。
3. heatmap 默认使用 Pareto-scaled matrix，不是原始 Peak Area。
4. 火山图默认使用 nominal `P.Value`，不是 FDR。
5. OPLS-DA 如果失败，会自动使用 fallback VIP-like ranking。
6. KEGG 分析需要用户提供 KEGG Compound ID 注释表。
7. 如果只修改图形参数，不需要从头运行，只需运行 `07_Plots.R`。
