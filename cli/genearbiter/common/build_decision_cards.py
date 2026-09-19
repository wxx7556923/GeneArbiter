#!/usr/bin/env python3
# script_id_md5: 79f27873b1228710c2174333b2e969f3
# created: 2026-06-26
# modified: 2026-06-26
# owner: project
# status: project_code
# purpose: 当前入口：调用共享 core 构建 AI annotation cards。
# inputs: 命令行参数指定的 full model_arbitration_cards.jsonl。
# outputs: 命令行参数指定的 AI annotation card JSONL、summary TSV、set trace JSONL 和 auto decisions JSONL。
# notes: 输出不暴露 correction/deletion/mode 语义。

"""Build AI annotation cards from full model-arbitration cards.

This wrapper delegates to decision_card_core.py. The generated card is the
AI-facing annotation-set decision layer. Gate A/B auto decisions and set traces
are written separately from API-bound AI cards.
"""

from decision_card_core import main


if __name__ == "__main__":
    main()
