# GeneArbiter Windows GUI 构建源码

这是 GeneArbiter 0.4.0 的紧凑型 Windows GUI 构建包。GUI 只负责收集输入并调用同包内的 GeneArbiter 核心；候选构建、证据解释、模型返回校验、ID 对齐和 GFF3 发布均由同一套核心代码完成。

当前是 Windows 构建候选源码，不是已经验证过的 Windows 可执行文件。本目录在 Linux 上完成了源码和命令构造测试；最终 EXE、真实 DeepSeek API 和真实全基因组运行仍需在 Windows 上验证。

## 界面范围

- 单页书页式界面，标题为 GeneArbiter，右上角折页打开帮助。
- 一个当前/骨架 GFF3，1–3 个可排序候选 GFF3。
- 可选参考基因组 FASTA、RNA 剪接 TSV、同源蛋白 GFF、长读长证据目录。
- 第一版仅显示已经接入核心的 DeepSeek Chat 和 DeepSeek Reasoner。
- 固定 Temperature 为 0；固定工程规则不向用户开放。
- API Key 只放入当次 worker 的进程环境，不写入配置、日志或结果。
- 输出 clean GFF3、trace GFF3、精简 ID mapping 和运行摘要。

内部 local_rule 只用于构建 smoke test，不在 GUI 中显示，也不会在 API 失败时静默替代模型服务。

## Windows 构建

建议使用 64 位 Windows 10/11 和 64 位 Python 3.11 或 3.12。首次构建需要联网从 PyPI 安装 PySide6、PyInstaller 和 PyYAML。

1. 解压源码包到较短的英文路径，例如 C:\GeneArbiterBuild。
2. 在该目录打开 PowerShell。
3. 执行：

        Set-ExecutionPolicy -Scope Process Bypass
        .\build_windows.ps1

脚本会：

1. 创建隔离环境 .venv-windows；
2. 安装构建依赖；
3. 运行 GUI 输入模型单元测试和源码编译检查；
4. 生成 GeneArbiter.exe 与内部 genearbiter-worker.exe；
5. 运行 GUI 启动、worker 启动及微型离线完整流程 smoke test；
6. 生成 release\GeneArbiter-Windows-x64-0.4.0.zip 及 SHA256。

最终用户只需解压 Windows 成品 ZIP，保留整个 GeneArbiter 文件夹并双击 GeneArbiter.exe。不要单独移动 EXE，因为它依赖同目录的内部 worker 和运行库。

## 输入约束

- 当前 GFF3 是需要修订的完整注释骨架。
- 候选 GFF3 顺序只在证据接近时作为来源优先级。
- 所有 GFF3、证据和可选 FASTA 必须对应同一基因组组装及坐标系统。
- FASTA 可选；不提供时会关闭候选 CDS 序列和编码质量检查。
- RNA TSV 至少需要：

        chrom  intron_start_1based  intron_end_1based  strand
        junction_read_count  sample_support_count  supporting_samples

- 蛋白证据接受 miniprot/Spaln 类 GFF3/GTF，建议包含 Target 和 Identity。
- 长读长证据应提供汇总目录，而不是单独的中间 TSV。
- 输出目录必须是尚不存在的新目录，防止混合或覆盖旧运行。

## 结果

成功运行后，用户结果位于所选输出目录的 results 子目录：

        annotation.clean.gff3
        annotation.trace.gff3
        id_mapping.tsv.gz
        run_summary.txt

完整审计中间结果和日志保存在所选输出目录的 work 子目录。

## 当前验证边界

已在 Linux/Python 3.12 上验证：

- GUI 的纯 Python 输入检查和 CLI 参数构造；
- 所带 GeneArbiter 核心源码可编译；
- GeneArbiter 核心既有单元测试和微型 CLI 流程。

尚未验证：

- Windows 原生构建及 EXE 启动；
- Windows 中文、空格和超长路径；
- 真实 DeepSeek API 调用；
- 新的真实物种全基因组运行、性能、内存和中断恢复；
- Python 3.10/3.11/3.13 的完整兼容性。

因此第一次 PowerShell 构建完成后，应把构建日志和 smoke test 结果保留下来。本项目源码使用 [MIT License](LICENSE)；Windows 构建包中的 Python、Qt、PySide6 等第三方组件遵循各自许可证，见 [第三方说明](THIRD_PARTY_NOTICES.md)。
