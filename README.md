# GeneArbiter

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23113893.svg)](https://doi.org/10.5281/zenodo.23113893)

| 版本 | 内容 | 入口 |
| --- | --- | --- |
| [CLI](cli/README.md) | Linux 命令行源码、示例和测试 | `python3 -m pip install ./cli` |
| [GUI](gui/README.md) | Windows 界面及构建源码，内含相同版本的核心代码 | `gui/build_windows.ps1` |
| [下载文件](downloads/) | CLI 源码包、wheel 和 Windows 可执行 ZIP | 先核对 `downloads/SHA256SUMS` |

Windows 用户可下载 `downloads/GeneArbiter-Windows-x64-0.4.0.zip`，完整解压后按 [使用说明](gui/WINDOWS_USER_GUIDE.md) 启动。CLI 用户可从 `cli/` 安装，或使用 `downloads/` 中的 wheel / 源码包。两个源码目录各自可独立构建；GUI 包含一份与 CLI 相同的 `genearbiter` 核心代码。


本项目原创源码采用 [MIT License](LICENSE)（SPDX：`MIT`）。Windows ZIP 包含第三方运行库，其许可证不因项目使用 MIT 而改变；详见 [第三方说明](gui/THIRD_PARTY_NOTICES.md)。
