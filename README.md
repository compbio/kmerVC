# kmerVC

k-mer **validator** for somatic variant calls (not a de-novo caller). Paper: *Feasibility of k-mers counts for somatic mutation calls* (NAR Cancer 2020). Lab page: [dna-discovery.stanford.edu](https://dna-discovery.stanford.edu/research/software/kmervc/).

Agent / revival context: see [AGENTS.md](AGENTS.md).

## Quick smoke (original CLI)

Needs Jellyfish, bedtools, and a Python env (see `python3_requirements.txt`).

```bash
python kmervc.py compare -k 31 -j1 tumor.jf -j2 normal.jf -b variants.bed -o out -fi ref.fa
```

Or from FASTQs (`-t1/-t2` tumor, `-c1/-c2` normal). Examples under `examples/` (download companion data from the lab software page).
