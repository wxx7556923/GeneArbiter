# GeneArbiter

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23113893.svg)](https://doi.org/10.5281/zenodo.23113893)

| 版本 | 内容 | 入口 |
| --- | --- | --- |
| [CLI](cli/README.md) | Linux 命令行源码、示例和测试 | `python3 -m pip install ./cli` |
| [GUI](gui/README.md) | Windows 界面及构建源码，内含相同版本的核心代码 | `gui/build_windows.ps1` |
| [下载文件](downloads/) | CLI 源码包、wheel 和 Windows 可执行 ZIP | 先核对 `downloads/SHA256SUMS` |

Windows 用户可下载 `downloads/GeneArbiter-Windows-x64-0.4.0.zip`，完整解压后按 [使用说明](gui/WINDOWS_USER_GUIDE.md) 启动。CLI 用户可从 `cli/` 安装，或使用 `downloads/` 中的 wheel / 源码包。两个源码目录各自可独立构建；GUI 包含一份与 CLI 相同的 `genearbiter` 核心代码。


本项目原创源码采用 [MIT License](LICENSE)（SPDX：`MIT`）。Windows ZIP 包含第三方运行库，其许可证不因项目使用 MIT 而改变；详见 [第三方说明](gui/THIRD_PARTY_NOTICES.md)。

## Quick Start

GeneArbiter supports both a command-line interface (CLI) and a graphical user interface (GUI).

### 1. CLI (Linux)

First, set your DeepSeek API Key:

```bash
export DEEPSEEK_API_KEY="your_api_key"
```

Run GeneArbiter:

```bash
genearbiter run-files \
  --current current.gff3 \
  --candidate ToolA=annotation_a.gff3 \
  --candidate ToolB=annotation_b.gff3 \
  --genome genome.fa \
  --out-dir results/run1 \
  --profile deepseek_flash
```

- `--current`: Existing reference annotation GFF3.
- `--candidate`: Candidate annotation GFF3 (1–3 files).
- `--genome`: Reference genome FASTA (optional but recommended).
- `--out-dir`: New output directory.
- `--profile deepseek_flash`: Enables DeepSeek API inference instead of the default local rule.

Optional evidence can be provided using `--evidence`, for example:

```bash
--evidence splice_junctions=merged_splice_junctions.tsv
```

### 2. GUI (Windows)

Download `GeneArbiter-Windows-x64-0.4.0.zip` from the [downloads directory](https://github.com/wxx7556923/GeneArbiter/tree/main/downloads).

1. Extract the ZIP file and double-click `GeneArbiter.exe`.
2. Enter your DeepSeek API Key.
3. Select the current GFF3 and add 1–3 candidate GFF3 files.
4. Optionally provide reference genome FASTA and supporting evidence.
5. Select a new output directory and click **运行 GeneArbiter**.

The Windows package includes the required runtime; no additional Python installation is needed.

### 3. GFF3 Input Order

**GFF3 order matters:**

- **Current GFF3**: The existing annotation used as the backbone for correction.
- **Candidate GFF3 files**: Alternative annotations to be evaluated. Their order defines source priority when evidence is comparable.

For example:

```text
Current GFF3:   original.gff3
Candidate 1:    annotation_a.gff3
Candidate 2:    annotation_b.gff3
```

Candidate 1 has higher priority than Candidate 2 when evidence is comparable. However, source priority does not override stronger biological evidence.

All GFF3 files and supporting evidence must correspond to the same genome assembly.

### 4. Output

The main results are saved in `results/run1/results/`:

- `annotation.clean.gff3`: Final annotation for downstream analysis.
- `annotation.trace.gff3`: Final annotation with provenance information.
- `id_mapping.tsv.gz`: Mapping of changed annotation IDs.
- `run_summary.txt`: Run summary.
