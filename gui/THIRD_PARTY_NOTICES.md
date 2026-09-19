# 第三方组件许可证说明

本仓库的 GeneArbiter 原创源码采用 MIT License。Windows 可执行 ZIP 还包含第三方组件，它们不受 GeneArbiter 的 MIT 许可证覆盖。

- Python 3.12：Python Software Foundation License；https://docs.python.org/3/license.html
- PySide6、Shiboken6 和 Qt 6.11.2：Qt for Python 社区版按 LGPLv3/GPLv3 等条款提供。随包附 `LGPL-3.0.txt` 与 `GPL-3.0.txt`。官方说明：https://doc.qt.io/qtforpython-6/ 及 https://doc.qt.io/qt-6/licensing.html
- PyYAML 6.0.3：MIT License；https://github.com/yaml/pyyaml/blob/main/LICENSE
- OpenSSL 3（`libcrypto-3.dll`、`libssl-3.dll`）：Apache License 2.0；https://www.openssl.org/source/license.html
- Qt 和 Python 运行库还可能包含其他上游组件；其具体授权与署名以相应上游发行材料为准。

此说明用于区分项目源码和随附依赖的授权，不改变第三方许可证条款。Windows ZIP 中未使用的 Qt Virtual Keyboard 组件已从公开副本移除；原本地归档保持不变。
