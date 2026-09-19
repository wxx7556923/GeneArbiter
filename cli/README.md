# GeneArbiter

GeneArbiter 是一个证据约束的基因注释集仲裁命令行工具。它以现有注释（current/backbone）和 1–3 套候选 GFF3 为基础，结合可用的 RNA 剪接、长读长和同源蛋白证据，在已有候选结构之间选择，并导出干净 GFF3、精简溯源 GFF3 和旧 ID 对应表。

GeneArbiter 不会自行生成新的 exon 或 CDS 坐标。输出是候选注释的证据约束裁决，不是实验验证、生物学真值或“准确率”。

## 系统要求与安装

- Linux；
- Python 3.10 或更高版本；
- 核心流程不调用 STAR、miniprot、BUSCO 或其他上游软件，用户应先准备好 GFF3 和证据文件。

从源码目录安装：

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install .
genearbiter --version
```

## 直接从文件运行

最简单的 correction 运行不需要手写配置文件：

```bash
genearbiter run-files \
  --current old_annotation.gff3 \
  --candidate ToolA=tool_a.gff3 \
  --candidate ToolB=tool_b.gff3 \
  --evidence splice_junctions=merged_splice_junctions.tsv \
  --genome genome.fa \
  --out-dir results/my_species
```

`--candidate` 可重复 1–3 次。`--evidence` 可重复，裸路径会自动识别类型；正式运行建议显式写出类型。默认 `local_rule` 是无网络、可复现的规则基线。远程模型 profile 需要对应服务的 API key，且结果不应与 `local_rule` 混为同一个运行条件。

如果不提供 `--genome`，流程仍可运行，但会关闭候选 CDS 序列和编码完整性检查；这不建议作为正式结果的默认做法。

## 用户结果

`run-files` 会保存实际运行配置，并将面向用户的结果收敛到一个目录：

```text
results/my_species/
├── genearbiter.config.json
├── results/
│   ├── annotation.clean.gff3
│   ├── annotation.trace.gff3
│   ├── id_mapping.tsv.gz
│   └── run_summary.txt
└── work/                 # 完整中间产物、日志和审计信息
```

- `annotation.clean.gff3`：常规下游分析首选，不含 risk、review、评分或运行参数。
- `annotation.trace.gff3`：结构和 ID 与 clean GFF 完全一致，仅在 gene/transcript 层保留精简溯源字段。
- `id_mapping.tsv.gz`：只记录发生变化的 current lineage 和入选来源；原样保留的旧 ID 不重复列出。
- `work/`：供复现和审计，不建议当作网站下载内容。

## 证据输入

当前公开 CLI 识别以下核心证据类型：

| 类型 | 输入 | 最小识别字段/特征 |
|---|---|---|
| `splice_junctions` | TSV | `chrom, intron_start_1based, intron_end_1based` |
| `model_junction_support` | TSV | `transcript_id, junction_support_fraction` |
| `protein_gff` | GFF3/GTF | miniprot/Spaln 类蛋白比对记录 |
| `long_read_support_dir` | 目录 | GeneArbiter 支持的长读长汇总文件 |

不能识别的 evidence 会立即报错，不会被静默忽略。短读长 TSV 生成示例见 `examples/short_read/README.md`。

## 高级配置模式

需要固定参数、批量运行或远程模型时，使用 JSON/YAML 配置：

```bash
genearbiter validate-config --strict --config configs/example.json
genearbiter plan --config configs/example.json
genearbiter run --config configs/example.json
```

API key 只能通过配置指定的环境变量读取，不应写入配置、日志或代码。

## 科研边界

- 第一个 GFF 是 current/backbone，后续 GFF 是候选来源；同分时依输入顺序确定来源优先级。
- 候选 GFF 默认标准化为 `gene -> mRNA -> exon/CDS`；current/backbone 不改写。
- 整类 evidence 未提供时，评分分母按本次运行的可用证据类型调整。
- 某个位点没有观察到证据支持，不等于该结构被证伪。
- `##gff-version 3` 是 GFF3 规范声明，不是 GeneArbiter 软件版本。

## 测试

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
python3 -m compileall -q genearbiter tests
```

这些测试覆盖代码契约和小型 fixture，不替代真实远程 API、全基因组性能或科研结果验证。

## 发布状态

本项目源码使用 [MIT License](LICENSE)；打包所需的第三方依赖遵循各自许可证。
