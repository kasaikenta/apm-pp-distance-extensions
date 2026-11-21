# Increasing Quantum LDPC Code Distance through Explicit Extensions: data

Concrete code data accompanying Kenta Kasai's research draft
*Increasing Quantum LDPC Code Distance through Explicit Extensions*.
The original v1.0.0 data remain unchanged. The current revision includes five
default-rule CPM–PP examples; cite the specific commit used.

## Complete exact-distance evidence

Release v1.3.0 adds the archived lower-bound search records and physical upper
witnesses for **15 exact distances**, with an auditable claim-to-matrix index.
See [the evidence table and verification instructions](exact-distance/README.md)
and [download the raw archives](https://github.com/kasaikenta/apm-pp-distance-extensions/releases/tag/v1.3.0).
It also adds the fixed-rule CPM (4,10) P37 and P73 matrices with exact distances
22 and 24. Historical data files and release tags retain their original content.

## Included codes

| Family | P | n | k | Distance | Check weight | Column weight per CSS sector |
|---|---:|---:|---:|---|---:|---:|
| Minimal-twist APM seed | 192 | 1536 | 388 | 26 (exact) | 8 | 3 |
| Fixed-data APM extension | 384 | 3072 | 772 | 32 (exact) | 8 | 3 |
| Selected power-9 series | 16 | 128 | 36 | 6 (exact) | 8 | 3 |
| Selected power-9 series | 32 | 256 | 68 | 12 (exact) | 8 | 3 |
| Selected power-9 series | 64 | 512 | 132 | 16 (exact) | 8 | 3 |
| Selected power-9 series | 128 | 1024 | 260 | 18 (exact) | 8 | 3 |
| Selected power-9 series | 256 | 2048 | 516 | 20 (exact) | 8 | 3 |
| Power-9 next-size failure | 512 | 4096 | 1088 | 4 (exact) | 8 | 3 |
| CPM extension | 127 | 1270 | 260 | 32 ≤ d ≤ 70 | 10 | 4 |

The main P192/P384 pair uses one fixed M, S and theta, with
P=192·2^t and D=96·2^floor(t/2). The small power-9 series used
distance-guided selection of integer representatives preserving preceding
matrices. Its next-size failure is included: monotonic distance growth is not
a theorem. The P127 interval is frozen; ongoing computations do not change
this release, and its distance must not be cited as exactly 32.

## Additional default-rule CPM series

| P | n | k | Exact distance | Column weight | Check weight |
|---:|---:|---:|---:|---:|---:|
| 23 | 276 | 98 | 14 | 4 | 12 |
| 37 | 444 | 154 | 16 | 4 | 12 |
| 43 | 516 | 178 | 16 | 4 | 12 |
| 47 | 564 | 194 | 14 | 4 | 12 |
| 53 | 636 | 218 | 20 | 4 | 12 |

The P43 and P47 cases test shorter lengths against the distance 20 P53 reference.
Neither retains distance 20; both results are included without reselection.

See [the self-contained data and verifier](default-cpm12/README.md). The original
seed is the P23 example of [Design Principles, Table 8](https://arxiv.org/html/2609.30069v1#A9.T8).
Run `python3 default-cpm12/verify.py` in addition to the original verifier below.
The unchanged general algorithm fixes parameters at the seed and performs no
candidate screening at the extension sizes. The original nine instances below
and their frozen tag are preserved; not every archived instance is in the paper.

## Files and coordinates

- `data/index.json`: code parameters and per-code file locations.
- `data/codes/*.json`: full sparse HX and HZ, affine multipliers/shifts,
  ranks, distance status, and nontrivial logical-operator supports.
- `data/rules/*.json`: matching arrays, integer data and fixed theta.
- `scripts/zsnf.py`: pinned deterministic Smith-coordinate implementation.
- `evidence/`: archived audit summaries and physical verification records.
- `SHA256SUMS.json`: hashes of the released data, evidence and scripts.

All indices are zero-based. Every row of `hx` or `hz` lists the positions of
its nonzero entries. For APM codes q=8u+l; for the CPM example q=lP+u.
Block (f,l) maps a qubit lift index u to check lift index a[f,l]u+b[f,l]
modulo P. The first J block rows of a,b describe HX and the next J describe HZ.
An X witness belongs to ker(HZ) but not row(HX); a Z witness belongs to
ker(HX) but not row(HZ). Different coordinate conventions are explicit in
each code file.

## Reproduce the matrices and physical checks

Python 3.10 or newer, standard library only:

```sh
python3 scripts/verify.py
```

This checks file hashes, reconstructs all nine matrix pairs from the fixed
rules, compares every matrix row and affine label, computes binary ranks,
checks CSS orthogonality, row/column weights and even kernels, and verifies
zero syndrome and stabilizer nonmembership of the supplied logical operators.

Exact distance and lower-bound claims also depend on the previously completed
exhaustive exclusion audits. The short script does **not** rerun those searches.
The Git tree includes compact summaries. The complete archived search records
for the exact-distance cases are available as release assets, with receipt and
coverage auditors in [exact-distance/](exact-distance/README.md). These execution
records are not a formal per-node proof trace; independently rerunning the
exclusions requires the supplied search inputs and implementation.

## Attribution and citation

The P192 minimal-twist seed was publicly listed as **Okada's construction** in
the [PP construction catalogue](https://kasai.ict.eng.isct.ac.jp/pair_partition_cpm_css_codes_20260714.html).
Its original defining JSON is preserved under `evidence/` (SHA-256
`5a079d93fb3e9b5daa6f73665ac59b44da05be0a35464e07501ff64e1c2bec8e`).
The general APM–PP cancellation condition is prior work: Okada and Kasai,
*Pair-Partition Constructions for CPM-Based Quantum LDPC Codes*,
[arXiv:2607.14091v4](https://arxiv.org/abs/2607.14091v4), Section VII.
The CPM template derives from the P37 [[370,80,22]] example in
[Design Principles, Table 8](https://arxiv.org/abs/2609.30069v1).
Distance searches used [Jong Yeon Lee's QEDist](https://github.com/yeon-lee/QEDist)
at commit `aeb7006147b777c121fd3ab8c6841da99238a84e`; the matching engine and its
license are included in the evidence assets. The earlier CPM P37 certificates
come from the original PP verifier, as documented in the evidence index.

Use `CITATION.cff` and cite a specific commit when reproducing these data.
Repository: https://github.com/kasaikenta/apm-pp-distance-extensions

## 日本語

当初収録した9個の符号について、検査行列・構成パラメータ・論理演算子を公開します。
主例は [[1536,388,26]] → [[3072,772,32]] で、チェック重み8・各側列重み3を維持します。
CPM例 [[1270,260,d]] は32≤d≤70の区間であり、距離32が確定した符号ではありません。
再構成スクリプトは物理行列、階数、直交性、重み、上界論理を検証します。
下界の全探索は再実行せず、保存済み監査に依存します。公開種の岡田氏への帰属を保持します。
引用時はリポジトリ名に加えて使用したコミットを指定してください。
