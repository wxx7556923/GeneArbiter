# GeneArbiter 0.4.0 Windows 使用说明

GeneArbiter 在固定、可审计的工程流程中，结合现有注释、候选注释和可选生物学证据，对候选基因模型进行仲裁，并输出适合使用的 clean GFF3、用于溯源的 trace GFF3 和精简 ID mapping。

## 启动

1. 将下载的 ZIP 完整解压到本地目录。
2. 保留 `GeneArbiter.exe`、`genearbiter-worker.exe` 和 `_internal` 在同一目录。
3. 双击 `GeneArbiter.exe`。

本成品已经包含 Python 和运行库，用户不需要另行安装 Python。程序适用于64位 Windows 10/11。未签名的科研软件首次启动时可能触发 Windows SmartScreen 提示，请确认下载来源和发布页提供的 SHA256 后再运行。

## 输入

- 当前注释 GFF3：必填，是需要修订的完整注释骨架。
- 候选 GFF3：必填，可提供1–3个；从上到下表示证据接近时的来源优先级。
- 参考基因组 FASTA：可选但推荐，用于候选 CDS 序列和编码质量检查。
- RNA剪接证据：可选，使用 `merged_splice_junctions.tsv`。
- 同源蛋白证据：可选，使用 miniprot/Spaln 类 GFF3/GTF。
- 长读长证据目录：可选，使用 GeneArbiter 支持的长读长汇总目录。
- DeepSeek API Key：必填，仅传给当次运行的 worker 进程。
- 模型：DeepSeek Chat 或 DeepSeek Reasoner。
- 并发请求数：控制 API 并发请求，不是 CPU 线程数。
- 输出目录：必填，必须是尚不存在的新目录。

所有 GFF3、证据和可选 FASTA 必须来自同一基因组组装，并使用一致的染色体名称和坐标系统。候选顺序只是证据接近时的来源优先级，不覆盖明确的生物学证据。

## 证据格式

RNA剪接 TSV 至少包含以下列：

```text
chrom
intron_start_1based
intron_end_1based
strand
junction_read_count
sample_support_count
supporting_samples
```

RNA内含子坐标使用1-based闭区间。同源蛋白 GFF 的主命中特征应为 `mRNA`、`match`、`protein_match` 或 `cDNA_match`，建议包含 `Target` 和 `Identity` 属性。

不提供 FASTA 时程序仍可运行，但会跳过候选 CDS 序列和编码质量检查。不提供 RNA、蛋白或长读长证据时程序也能运行，但生物学证据有限，结果应谨慎审核。

## 运行和结果

填写输入后点击“运行 GeneArbiter”。程序运行期间可以展开日志查看当前阶段。不要在运行期间移动输入文件或关闭网络。

成功后，所选输出目录包含：

```text
results/
  annotation.clean.gff3
  annotation.trace.gff3
  id_mapping.tsv.gz
  run_summary.txt

work/
  完整审计中间结果和分步骤日志
```

- `annotation.clean.gff3`：用于下游分析的整洁注释。
- `annotation.trace.gff3`：与 clean GFF3 保持相同结构，并保留精简溯源信息。
- `id_mapping.tsv.gz`：记录发生变化的旧注释谱系和入选来源，不重复记录原样保留的 ID。
- `run_summary.txt`：记录运行状态、模型配置和结果入口。

## API Key和隐私

API Key不会写入 GeneArbiter 配置、日志或结果。需要模型判断的决策卡会发送给所选 DeepSeek API；GeneArbiter不会向模型发送基因组 FASTA 全序列，也不允许模型生成新坐标。API失败时程序会报错，不会静默改用内部测试规则。

## 使用边界

GeneArbiter输出是证据约束的注释修订结果，不等同于独立真值或准确率证明。正式科研使用应检查运行日志、trace GFF3和mapping，并结合项目自己的质量控制和生物学验证。
