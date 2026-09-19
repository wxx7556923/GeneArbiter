"""Compact single-page Windows interface for GeneArbiter.

The interface only collects inputs and calls the bundled public CLI. Annotation
scoring, evidence interpretation, ID reconciliation and GFF release stay in the
GeneArbiter core.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, QUrl
from PySide6.QtGui import QColor, QCloseEvent, QDesktopServices, QPainter, QPolygon
from PySide6.QtWidgets import (
    QAbstractButton,
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from . import __version__
from .model import (
    PROFILE_CHOICES,
    CandidateInput,
    RunRequest,
    build_worker_arguments,
    source_name_from_path,
    validate_request,
)
from .runner import WorkflowRunner


GFF_FILTER = "GFF3 annotation (*.gff3 *.gff *.gff3.gz *.gff.gz);;All files (*)"
FASTA_FILTER = "Genome FASTA (*.fa *.fasta *.fna);;All files (*)"
TSV_FILTER = "Tab-separated evidence (*.tsv *.txt);;All files (*)"
PROTEIN_FILTER = "Protein alignment GFF (*.gff3 *.gff *.gtf);;All files (*)"

STYLE = """
QMainWindow { background: #dfe4dc; }
QWidget { color: #26332d; font-family: "Segoe UI", "Microsoft YaHei UI"; font-size: 10pt; }
QFrame#page { background: #fbfaf4; border: 1px solid #c9cec5; border-radius: 7px; }
QFrame#hero { background: #173f36; border: none; border-radius: 6px; }
QLabel#title { color: #fffdf4; font-size: 25pt; font-weight: 650; }
QLabel#subtitle { color: #dbe7df; font-size: 9.5pt; }
QLabel#version { color: #b9ccc2; font-size: 8.5pt; }
QFrame#section { background: #fffef9; border: 1px solid #d9ddd5; border-radius: 6px; }
QLabel#sectionTitle { color: #173f36; font-size: 11pt; font-weight: 650; }
QLabel#hint { color: #69746e; font-size: 8.8pt; }
QLabel#warning { color: #75501d; background: #fff7e4; border: 1px solid #ead7ad; border-radius: 4px; padding: 6px; }
QLineEdit, QComboBox, QSpinBox, QTableWidget, QPlainTextEdit {
  background: white; border: 1px solid #bdc7c0; border-radius: 4px; padding: 5px;
  selection-background-color: #2c7563;
}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QTableWidget:focus { border: 1px solid #2c7563; }
QPushButton { background: #edf1ed; border: 1px solid #b9c4bd; border-radius: 4px; padding: 6px 11px; }
QPushButton:hover { background: #e2eae4; }
QPushButton:disabled { color: #9aa39e; background: #eff1ef; }
QPushButton#primary { color: white; background: #236451; border: 1px solid #1b5141; font-weight: 650; padding: 9px 22px; }
QPushButton#primary:hover { background: #2b735e; }
QPushButton#danger { color: #8c3028; }
QHeaderView::section { background: #eef2ee; color: #425049; border: none; padding: 5px; }
QProgressBar { border: 1px solid #bdc7c0; border-radius: 4px; text-align: center; background: white; }
QProgressBar::chunk { background: #3b816d; }
QScrollArea { border: none; background: transparent; }
"""


class FoldButton(QAbstractButton):
    """Dog-ear help affordance in the page header."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedSize(54, 54)
        self.setToolTip("查看输入格式和软件说明")
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    def paintEvent(self, _event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor("#edf3ee" if self.underMouse() else "#d9e4dc"))
        painter.drawPolygon(QPolygon([self.rect().topRight(), self.rect().bottomRight(), self.rect().topLeft()]))
        painter.setPen(QColor("#315a4d"))
        font = painter.font()
        font.setBold(True)
        font.setPointSize(11)
        painter.setFont(font)
        painter.drawText(self.rect().adjusted(20, 0, 0, -20), Qt.AlignmentFlag.AlignCenter, "?")
        painter.end()


class PathPicker(QWidget):
    def __init__(
        self,
        *,
        file_filter: str = "All files (*)",
        directory: bool = False,
        output_target: bool = False,
        placeholder: str = "",
    ) -> None:
        super().__init__()
        self.file_filter = file_filter
        self.directory = directory
        self.output_target = output_target
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        button = QPushButton("浏览…")
        button.clicked.connect(self.browse)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(self.edit, 1)
        layout.addWidget(button)

    def text(self) -> str:
        return self.edit.text().strip()

    def browse(self) -> None:
        start = self.text() or str(Path.home())
        if self.directory or self.output_target:
            chosen = QFileDialog.getExistingDirectory(self, "选择目录", start)
            if chosen and self.output_target:
                stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                chosen = str(Path(chosen) / f"GeneArbiter_run_{stamp}")
        else:
            chosen, _ = QFileDialog.getOpenFileName(self, "选择文件", start, self.file_filter)
        if chosen:
            self.edit.setText(chosen)


class Section(QFrame):
    def __init__(self, title: str) -> None:
        super().__init__()
        self.setObjectName("section")
        self.content = QVBoxLayout(self)
        self.content.setContentsMargins(14, 11, 14, 12)
        self.content.setSpacing(8)
        heading = QLabel(title)
        heading.setObjectName("sectionTitle")
        self.content.addWidget(heading)


class HelpDialog(QDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("GeneArbiter 输入说明")
        self.resize(670, 570)
        layout = QVBoxLayout(self)
        browser = QTextBrowser()
        browser.setHtml(
            """
            <h2>GeneArbiter 输入说明</h2>
            <h3>GFF3 顺序</h3>
            <p><b>当前注释</b>是需要修订的全基因组骨架；候选 GFF3 可提供 1–3 个。
            候选顺序只在证据接近时作为来源优先级，不覆盖明确的生物学证据。</p>
            <h3>RNA 剪接证据</h3>
            <p>输入 <code>merged_splice_junctions.tsv</code>，至少包含：</p>
            <pre>chrom\tintron_start_1based\tintron_end_1based\tstrand\tjunction_read_count\tsample_support_count\tsupporting_samples</pre>
            <p>坐标为 1-based 闭区间；<code>chrom</code> 必须与 GFF3 一致。</p>
            <h3>同源蛋白证据</h3>
            <p>输入 miniprot/Spaln 类 GFF3/GTF。主命中特征应为
            <code>mRNA</code>、<code>match</code>、<code>protein_match</code> 或
            <code>cDNA_match</code>；建议包含 <code>Target</code> 和
            <code>Identity</code> 属性。</p>
            <h3>长读长证据</h3>
            <p>选择长读长汇总目录。每个来源子目录中使用
            <code>long_read_model_support_by_model.tsv</code> 和
            <code>long_read_set_support_by_locus.tsv</code>。</p>
            <h3>FASTA</h3>
            <p>FASTA 可选但推荐，用于坐标、CDS、内部终止密码子和剪接位点检查；
            必须与 GFF3 对应同一组装版本。</p>
            <h3>结果与密钥</h3>
            <p>结果目录包含 clean GFF3、trace GFF3、ID mapping 和运行摘要。
            API Key 只传入当次 worker 进程，不写入配置、日志或结果。</p>
            """
        )
        layout.addWidget(browser)
        close = QPushButton("关闭")
        close.clicked.connect(self.accept)
        layout.addWidget(close, alignment=Qt.AlignmentFlag.AlignRight)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"GeneArbiter {__version__}")
        self.resize(870, 790)
        self.setMinimumSize(760, 650)
        self.runner = WorkflowRunner(self)
        self.runner.output.connect(self.append_output)
        self.runner.started.connect(self.run_started)
        self.runner.finished.connect(self.run_finished)
        self.runner.failed_to_start.connect(self.run_failed_to_start)
        self.last_output_dir = ""
        self._build_ui()

    def _build_ui(self) -> None:
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        host = QWidget()
        outer = QVBoxLayout(host)
        outer.setContentsMargins(22, 18, 22, 22)
        page = QFrame()
        page.setObjectName("page")
        self.page_layout = QVBoxLayout(page)
        self.page_layout.setContentsMargins(16, 14, 16, 16)
        self.page_layout.setSpacing(10)
        outer.addWidget(page)
        outer.addStretch(1)
        scroll.setWidget(host)
        self.setCentralWidget(scroll)
        self._add_header()
        self._add_annotations()
        self._add_evidence()
        self._add_service()
        self._add_output()
        self._add_actions()

    def _add_header(self) -> None:
        hero = QFrame()
        hero.setObjectName("hero")
        layout = QHBoxLayout(hero)
        layout.setContentsMargins(18, 11, 0, 11)
        text = QVBoxLayout()
        text.setSpacing(1)
        title = QLabel("GeneArbiter")
        title.setObjectName("title")
        subtitle = QLabel("Evidence-guided gene annotation arbitration")
        subtitle.setObjectName("subtitle")
        version = QLabel(f"Windows GUI · Core {__version__}")
        version.setObjectName("version")
        text.addWidget(title)
        text.addWidget(subtitle)
        text.addWidget(version)
        layout.addLayout(text, 1)
        fold = FoldButton()
        fold.clicked.connect(self.show_help)
        layout.addWidget(fold, alignment=Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignRight)
        self.page_layout.addWidget(hero)

    def _add_annotations(self) -> None:
        section = Section("01  注释文件")
        grid = QGridLayout()
        grid.setHorizontalSpacing(9)
        grid.setVerticalSpacing(7)
        self.current = PathPicker(file_filter=GFF_FILTER, placeholder="必填 · 需要修订的当前注释")
        self.fasta = PathPicker(file_filter=FASTA_FILTER, placeholder="可选 · 推荐用于编码质量检查")
        grid.addWidget(QLabel("当前注释 GFF3 *"), 0, 0)
        grid.addWidget(self.current, 0, 1)
        grid.addWidget(QLabel("参考基因组 FASTA"), 1, 0)
        grid.addWidget(self.fasta, 1, 1)
        grid.setColumnStretch(1, 1)
        section.content.addLayout(grid)
        section.content.addWidget(QLabel("候选 GFF3 *（从上到下为证据接近时的来源优先级）"))
        self.candidates = QTableWidget(0, 2)
        self.candidates.setHorizontalHeaderLabels(["来源名", "GFF3 路径"])
        self.candidates.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.candidates.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.candidates.verticalHeader().setVisible(False)
        self.candidates.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.candidates.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.candidates.setMinimumHeight(118)
        self.candidates.setMaximumHeight(150)
        section.content.addWidget(self.candidates)
        buttons = QHBoxLayout()
        for label, callback in (
            ("＋ 添加候选", self.choose_candidates),
            ("↑ 上移", lambda: self.move_candidate(-1)),
            ("↓ 下移", lambda: self.move_candidate(1)),
            ("移除", self.remove_candidate),
        ):
            button = QPushButton(label)
            button.clicked.connect(callback)
            buttons.addWidget(button)
        buttons.addStretch(1)
        section.content.addLayout(buttons)
        self.page_layout.addWidget(section)

    def _add_evidence(self) -> None:
        section = Section("02  支持证据（可选）")
        grid = QGridLayout()
        grid.setHorizontalSpacing(9)
        grid.setVerticalSpacing(7)
        self.rna = PathPicker(file_filter=TSV_FILTER, placeholder="merged_splice_junctions.tsv")
        self.protein = PathPicker(file_filter=PROTEIN_FILTER, placeholder="miniprot/Spaln 类 GFF3")
        self.long_read = PathPicker(directory=True, placeholder="GeneArbiter 长读长汇总目录")
        for row, (label, picker) in enumerate(
            (("RNA 剪接证据", self.rna), ("同源蛋白证据", self.protein), ("长读长证据目录", self.long_read))
        ):
            grid.addWidget(QLabel(label), row, 0)
            grid.addWidget(picker, row, 1)
        grid.setColumnStretch(1, 1)
        section.content.addLayout(grid)
        hint = QLabel("证据坐标和染色体名必须与 GFF3 一致。右上角折页中有格式和必需字段说明。")
        hint.setObjectName("hint")
        section.content.addWidget(hint)
        self.page_layout.addWidget(section)

    def _add_service(self) -> None:
        section = Section("03  模型服务")
        grid = QGridLayout()
        grid.setHorizontalSpacing(9)
        grid.setVerticalSpacing(7)
        self.service = QComboBox()
        self.service.addItem("DeepSeek API")
        self.model = QComboBox()
        for choice in PROFILE_CHOICES:
            self.model.addItem(choice.label, choice.profile)
        self.api_key = QLineEdit()
        self.api_key.setEchoMode(QLineEdit.EchoMode.Password)
        self.api_key.setPlaceholderText("仅传入当次运行，不写入文件")
        show_key = QCheckBox("显示")
        show_key.toggled.connect(
            lambda checked: self.api_key.setEchoMode(
                QLineEdit.EchoMode.Normal if checked else QLineEdit.EchoMode.Password
            )
        )
        key_row = QHBoxLayout()
        key_row.addWidget(self.api_key, 1)
        key_row.addWidget(show_key)
        self.workers = QSpinBox()
        self.workers.setRange(1, 64)
        self.workers.setValue(4)
        self.workers.setToolTip("API 并发请求数，不是 CPU 线程数。")
        grid.addWidget(QLabel("API 服务"), 0, 0)
        grid.addWidget(self.service, 0, 1)
        grid.addWidget(QLabel("模型"), 0, 2)
        grid.addWidget(self.model, 0, 3)
        grid.addWidget(QLabel("API Key *"), 1, 0)
        grid.addLayout(key_row, 1, 1, 1, 3)
        grid.addWidget(QLabel("并发请求数"), 2, 0)
        grid.addWidget(self.workers, 2, 1)
        fixed = QLabel("Temperature 固定为 0；工程规则由 GeneArbiter 固定执行。")
        fixed.setObjectName("hint")
        grid.addWidget(fixed, 2, 2, 1, 2)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 2)
        section.content.addLayout(grid)
        self.page_layout.addWidget(section)

    def _add_output(self) -> None:
        section = Section("04  输出")
        grid = QGridLayout()
        self.output_dir = PathPicker(output_target=True, placeholder="选择父目录后自动生成新的运行目录")
        grid.addWidget(QLabel("新输出目录 *"), 0, 0)
        grid.addWidget(self.output_dir, 0, 1)
        grid.setColumnStretch(1, 1)
        section.content.addLayout(grid)
        self.page_layout.addWidget(section)
        tips = QLabel(
            "① 当前 GFF3 是完整注释骨架　　"
            "② 候选顺序只在证据接近时作为来源优先级　　"
            "③ 所有输入必须对应同一组装版本"
        )
        tips.setWordWrap(True)
        tips.setObjectName("warning")
        self.page_layout.addWidget(tips)

    def _add_actions(self) -> None:
        row = QHBoxLayout()
        self.status = QLabel("准备就绪")
        self.status.setObjectName("hint")
        row.addWidget(self.status, 1)
        self.run_button = QPushButton("运行 GeneArbiter")
        self.run_button.setObjectName("primary")
        self.run_button.clicked.connect(self.start_run)
        self.stop_button = QPushButton("停止")
        self.stop_button.setObjectName("danger")
        self.stop_button.clicked.connect(self.runner.stop)
        self.stop_button.hide()
        row.addWidget(self.stop_button)
        row.addWidget(self.run_button)
        self.page_layout.addLayout(row)
        self.runtime = QFrame()
        self.runtime.setObjectName("section")
        layout = QVBoxLayout(self.runtime)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        layout.addWidget(self.progress)
        controls = QHBoxLayout()
        self.toggle_log = QPushButton("查看运行日志")
        self.toggle_log.setCheckable(True)
        self.toggle_log.toggled.connect(self.toggle_runtime_log)
        self.open_results = QPushButton("打开结果目录")
        self.open_results.clicked.connect(self.open_result_directory)
        self.open_results.hide()
        controls.addWidget(self.toggle_log)
        controls.addStretch(1)
        controls.addWidget(self.open_results)
        layout.addLayout(controls)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(5000)
        self.log.setMinimumHeight(150)
        self.log.hide()
        layout.addWidget(self.log)
        self.runtime.hide()
        self.page_layout.addWidget(self.runtime)

    def choose_candidates(self) -> None:
        remaining = 3 - self.candidates.rowCount()
        if remaining <= 0:
            QMessageBox.information(self, "候选数量", "最多允许 3 个候选 GFF3。")
            return
        paths, _ = QFileDialog.getOpenFileNames(self, "选择候选 GFF3", str(Path.home()), GFF_FILTER)
        for path in paths[:remaining]:
            row = self.candidates.rowCount()
            self.candidates.insertRow(row)
            self.candidates.setItem(row, 0, QTableWidgetItem(source_name_from_path(path, row + 1)))
            self.candidates.setItem(row, 1, QTableWidgetItem(path))

    def move_candidate(self, offset: int) -> None:
        row = self.candidates.currentRow()
        target = row + offset
        if row < 0 or target < 0 or target >= self.candidates.rowCount():
            return
        values = [[self.candidates.item(r, c).text() if self.candidates.item(r, c) else "" for c in range(2)] for r in (row, target)]
        for column in range(2):
            self.candidates.setItem(row, column, QTableWidgetItem(values[1][column]))
            self.candidates.setItem(target, column, QTableWidgetItem(values[0][column]))
        self.candidates.selectRow(target)

    def remove_candidate(self) -> None:
        row = self.candidates.currentRow()
        if row >= 0:
            self.candidates.removeRow(row)

    def request_from_ui(self) -> RunRequest:
        candidates = []
        for row in range(self.candidates.rowCount()):
            name = self.candidates.item(row, 0).text().strip() if self.candidates.item(row, 0) else ""
            path = self.candidates.item(row, 1).text().strip() if self.candidates.item(row, 1) else ""
            candidates.append(CandidateInput(name=name, path=path))
        return RunRequest(
            current_gff=self.current.text(),
            candidates=candidates,
            genome_fasta=self.fasta.text(),
            splice_junctions=self.rna.text(),
            protein_gff=self.protein.text(),
            long_read_support_dir=self.long_read.text(),
            output_dir=self.output_dir.text(),
            api_key=self.api_key.text().strip(),
            profile=str(self.model.currentData()),
            workers=self.workers.value(),
        )

    def start_run(self) -> None:
        request = self.request_from_ui()
        validation = validate_request(request)
        if validation.errors:
            QMessageBox.critical(self, "输入检查未通过", "\n".join(f"• {item}" for item in validation.errors))
            return
        if validation.warnings:
            answer = QMessageBox.warning(
                self,
                "运行前提示",
                "\n".join(f"• {item}" for item in validation.warnings) + "\n\n是否继续？",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel,
            )
            if answer != QMessageBox.StandardButton.Yes:
                return
        self.last_output_dir = str(Path(request.output_dir).expanduser().resolve())
        self.log.clear()
        self.runtime.show()
        self.open_results.hide()
        self.progress.setRange(0, 0)
        self.status.setText("正在启动…")
        self.set_inputs_enabled(False)
        self.runner.start(
            build_worker_arguments(request),
            request.api_key,
            str(Path.home()),
        )
        self.api_key.clear()

    def set_inputs_enabled(self, enabled: bool) -> None:
        for widget in (
            self.current, self.fasta, self.candidates, self.rna, self.protein, self.long_read,
            self.service, self.model, self.api_key, self.workers, self.output_dir,
        ):
            widget.setEnabled(enabled)
        self.run_button.setEnabled(enabled)
        self.stop_button.setVisible(not enabled)

    def run_started(self) -> None:
        self.status.setText("正在运行 GeneArbiter…")

    def append_output(self, text: str) -> None:
        self.log.appendPlainText(text.rstrip("\n"))
        match = re.search(r"\] run ([a-z_]+)", text)
        if match:
            labels = {
                "normalize_gffs": "标准化 GFF3", "validate_candidates": "检查候选模型",
                "audit_eligibility": "审核候选资格", "cards": "构建仲裁卡",
                "ai_cards": "整理模型输入", "api": "请求模型服务",
                "final_calls": "汇总最终决定", "export_gff": "生成候选结果",
                "reconcile_ids": "对齐注释 ID", "release_gff_qc": "检查并发布 GFF3",
            }
            self.status.setText(labels.get(match.group(1), match.group(1)))
        match = re.search(r"\]\s+(\d+)/(\d+)\s+\(", text)
        if match:
            done, total = int(match.group(1)), int(match.group(2))
            self.progress.setRange(0, max(total, 1))
            self.progress.setValue(done)
            self.status.setText(f"模型判断：{done} / {total}")

    def run_finished(self, exit_code: int) -> None:
        self.set_inputs_enabled(True)
        self.progress.setRange(0, 1)
        if exit_code == 0:
            self.progress.setValue(1)
            self.status.setText("运行完成")
            results = Path(self.last_output_dir) / "results"
            expected = (
                results / "annotation.clean.gff3", results / "annotation.trace.gff3",
                results / "id_mapping.tsv.gz", results / "run_summary.txt",
            )
            if all(path.is_file() for path in expected):
                self.open_results.show()
            else:
                self.status.setText("进程结束，但公开结果文件不完整")
                self.toggle_runtime_log(True)
        else:
            self.progress.setValue(0)
            self.status.setText(f"运行失败（退出码 {exit_code}）")
            self.toggle_runtime_log(True)

    def run_failed_to_start(self, message: str) -> None:
        self.set_inputs_enabled(True)
        self.status.setText("无法启动 GeneArbiter worker")
        QMessageBox.critical(self, "启动失败", message)

    def toggle_runtime_log(self, visible: bool) -> None:
        if self.toggle_log.isChecked() != visible:
            self.toggle_log.blockSignals(True)
            self.toggle_log.setChecked(visible)
            self.toggle_log.blockSignals(False)
        self.log.setVisible(visible)
        self.toggle_log.setText("收起运行日志" if visible else "查看运行日志")

    def open_result_directory(self) -> None:
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(self.last_output_dir) / "results")))

    def show_help(self) -> None:
        HelpDialog(self).exec()

    def closeEvent(self, event: QCloseEvent) -> None:  # noqa: N802
        if not self.runner.running:
            event.accept()
            return
        answer = QMessageBox.question(
            self, "停止运行？", "关闭窗口将停止当前 GeneArbiter worker。",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel,
        )
        if answer == QMessageBox.StandardButton.Yes:
            self.runner.stop()
            event.accept()
        else:
            event.ignore()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="GeneArbiter Windows GUI")
    parser.add_argument("--smoke-test", action="store_true", help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.smoke_test:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QApplication(sys.argv[:1])
    app.setApplicationName("GeneArbiter")
    app.setApplicationVersion(__version__)
    app.setStyleSheet(STYLE)
    window = MainWindow()
    if args.smoke_test:
        window.show()
        app.processEvents()
        QTimer.singleShot(0, app.quit)
        return app.exec()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
