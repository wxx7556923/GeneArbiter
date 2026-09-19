# Short-read Evidence Demo

本示例说明只有 short-read RNA-seq 时，如何先生成 GeneArbiter evidence TSV，再接入 `genearbiter run`。

## 1. 提取 short-read evidence

```bash
cd /path/to/GeneArbiter-Linux

genearbiter extract-short-read \
  --sample-manifest tests/fixtures/sample_manifest.tsv \
  --gff current=tests/fixtures/tiny_current.gff3 \
  --sj-out sample1=tests/fixtures/sample1.SJ.out.tab \
  --out-dir /tmp/genearbiter_short_read_demo
```

输出：

```text
/tmp/genearbiter_short_read_demo/sample_manifest.tsv
/tmp/genearbiter_short_read_demo/merged_splice_junctions.tsv
/tmp/genearbiter_short_read_demo/junction_support_by_model.tsv
```

## 2. 校验 evidence

```bash
genearbiter validate-evidence \
  --strict \
  --evidence-dir /tmp/genearbiter_short_read_demo
```

## 3. 在 config 中引用

```yaml
evidence:
  - name: short_read_splice_junctions
    type: splice_junctions
    path: /tmp/genearbiter_short_read_demo/merged_splice_junctions.tsv
  - name: rna_junction_support_by_model
    type: model_junction_support
    path: /tmp/genearbiter_short_read_demo/junction_support_by_model.tsv
```
