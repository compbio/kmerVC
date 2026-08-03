# AGENTS.md — kmerVC

Ahmed co-authored **kmerVC** ([NAR Cancer 2020](https://doi.org/10.1093/narcan/zcaa013)): a k-mer **validator** for somatic variants (not a de-novo caller). Validates candidate variants using jellyfish k-mer counts from tumor/normal (+ optional reference).

Upstream/public history: Stanford lab software page + original scripts in-repo. **Revival** is a showcase/usability project **in this repo** (not a side tree).

## North star (priority order)

1. Run on modern systems (Python 3.10+)
2. Public data end-to-end
3. Head-to-head benchmark vs standard tools
4. Simple standalone web app on his domain

Accuracy extras / RNA / germline = later. **Bioconda packaging is a rewrite constraint.**

## Repo map

| Path | Role |
|------|------|
| `kmervc.py`, `kmervclib.py` | Original CLI + library (Python-2-era; partially broken on 3.10+) |
| `revival/baseline/` | Patched runnable baseline |
| `revival/data/`, `slice/`, `runs/` | Scratch / golden fixtures (large data gitignored) |
| `revival/bench/` | Pipeline / slice scripts |
| `revival/research/` | Refactor design notes (hexagonal rewrite plan) |
| `examples/`, `TEST/` | Original tutorial paths |

## Locks

- **Benchmark core:** SEQC2 HCC1395 vs HCC1395BL; score with som.py / hap.py / vcfeval. WES reads vs WGS truth → many "Insufficient" is expected; judge **analyzable** variants.
- **Machine:** bound jellyfish/RAM; no two heavy jobs at once. Scratch under `revival/` (gitignored large data).
- **Proven path (chr22):** patched baseline runs; golden under `revival/slice/run_chr22/`. In rewrite: missing k-mers → count **0** (not KeyError); jellyfish `-s` must exceed distinct k-mer count.
- External deps: **Jellyfish**, **bedtools**, pandas/numpy/scipy.

## Smoke (baseline)

From `revival/` (after data/jellyfish setup — see `bench/` scripts):

```bash
# chr22 golden path is the regression fixture
# Invalidated / Insufficient / Validated counts in out_variant_summary_table.txt
```

Original CLI shape (for reference):

```bash
python kmervc.py compare -k 31 -j1 tumor.jf -j2 normal.jf -b variants.bed -o out -fi ref.fa
```

## Style for rewrite

Prefer pure core (window math, k-mer enum, stats) behind adapters for jellyfish/bedtools/IO. Thin CLI over pipeline. Golden-file tests against the chr22 slice. Don't thrash the laptop on full WGS until bounded jobs are solid.
