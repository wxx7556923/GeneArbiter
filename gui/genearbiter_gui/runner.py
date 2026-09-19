"""Run the bundled GeneArbiter CLI in an isolated worker process."""

from __future__ import annotations

import os
import sys
from pathlib import Path

from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, QTimer, Signal


class WorkflowRunner(QObject):
    output = Signal(str)
    started = Signal()
    finished = Signal(int)
    failed_to_start = Signal(str)

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read_output)
        self.process.started.connect(self.started)
        self.process.finished.connect(self._finished)
        self.process.errorOccurred.connect(self._error)

    @property
    def running(self) -> bool:
        return self.process.state() != QProcess.ProcessState.NotRunning

    def _command(self, arguments: list[str]) -> tuple[str, list[str]]:
        if getattr(sys, "frozen", False):
            worker = Path(sys.executable).resolve().with_name("genearbiter-worker.exe")
            return str(worker), arguments
        return sys.executable, ["-m", "genearbiter", *arguments]

    def start(self, arguments: list[str], api_key: str, working_directory: str) -> None:
        if self.running:
            raise RuntimeError("GeneArbiter 已在运行。")
        program, process_args = self._command(arguments)
        if getattr(sys, "frozen", False) and not Path(program).is_file():
            self.failed_to_start.emit(f"缺少内部 worker：{program}")
            return
        environment = QProcessEnvironment.systemEnvironment()
        environment.insert("DEEPSEEK_API_KEY", api_key)
        environment.insert("DEEPSEEK_BASE_URL", "https://api.deepseek.com/chat/completions")
        environment.insert("DEEPSEEK_FLASH_MODEL", "deepseek-chat")
        environment.insert("DEEPSEEK_FLASH_THINKING_MODEL", "deepseek-reasoner")
        environment.insert("GENEARBITER_IN_PROCESS_STEPS", "1")
        environment.insert("PYTHONUTF8", "1")
        self.process.setProcessEnvironment(environment)
        self.process.setWorkingDirectory(working_directory or os.getcwd())
        self.process.start(program, process_args)

    def stop(self) -> None:
        if not self.running:
            return
        self.output.emit("正在请求停止…\n")
        self.process.terminate()
        QTimer.singleShot(3000, self._kill_if_running)

    def _kill_if_running(self) -> None:
        if self.running:
            self.output.emit("worker 未及时退出，执行强制停止。\n")
            self.process.kill()

    def _read_output(self) -> None:
        data = bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if data:
            self.output.emit(data)

    def _finished(self, exit_code: int, _status: QProcess.ExitStatus) -> None:
        self._read_output()
        self.finished.emit(exit_code)

    def _error(self, error: QProcess.ProcessError) -> None:
        if error == QProcess.ProcessError.FailedToStart:
            self.failed_to_start.emit(self.process.errorString())
