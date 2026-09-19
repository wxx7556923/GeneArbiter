# GeneArbiter 0.4.0 发布候选版

GeneArbiter 用现有注释、候选 GFF3 和可用的转录或蛋白证据进行基因模型仲裁。输出是证据约束的注释修订结果，不代表生物学真值或准确率。

| 版本 | 内容 | 入口 |
| --- | --- | --- |
| [CLI](cli/README.md) | Linux 命令行源码、示例和测试 | `python3 -m pip install ./cli` |
| [GUI](gui/README.md) | Windows 界面及构建源码，内含相同版本的核心代码 | `gui/build_windows.ps1` |
| [下载文件](downloads/) | CLI 源码包、wheel 和 Windows 可执行 ZIP | 先核对 `downloads/SHA256SUMS` |

Windows 用户可下载 `downloads/GeneArbiter-Windows-x64-0.4.0.zip`，完整解压后按 [使用说明](gui/WINDOWS_USER_GUIDE.md) 启动。CLI 用户可从 `cli/` 安装，或使用 `downloads/` 中的 wheel / 源码包。两个源码目录各自可独立构建；GUI 包含一份与 CLI 相同的 `genearbiter` 核心代码。

## 验证范围

源码测试和编译检查覆盖代码契约及小型 fixture。Windows ZIP 是已有构建产物；本仓库不将其标记为经过 Windows 原生、真实 DeepSeek API 或全基因组运行验证的正式科研结果。详情见各目录 README。

本项目原创源码采用 [MIT License](LICENSE)（SPDX：`MIT`）。Windows ZIP 包含第三方运行库，其许可证不因项目使用 MIT 而改变；详见 [第三方说明](gui/THIRD_PARTY_NOTICES.md)。
