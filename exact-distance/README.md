# Exact-distance evidence

This release publishes the archived evidence for **15 exact code distances**,
including every exact example in the manuscript and the previously published
small-series and shorter-length examples. The non-exact P127 instance is excluded
from this list. Its later bound is 34 ≤ d ≤ 70; its frozen matrix-data record
still retains the earlier bound 32 ≤ d ≤ 70.

| Code | Exact parameters | Lower-bound evidence |
|---|---|---|
| APM P192 | [[1536,388,26]] | Both sectors excluded through 24; even kernels |
| APM P384 | [[3072,772,32]] | Both sectors excluded through 30; even kernels |
| CPM (4,10) P37 | [[370,80,22]] | Original PP-verifier exclusions through 20; parity and X/Z exchange; verified seed equivalence |
| CPM (4,10) P73 | [[730,152,24]] | Exclusion through 22; parity and verified X/Z exchange |
| CPM (4,12) P23 | [[276,98,14]] | Complete both-sector exclusion through 13 |
| CPM (4,12) P37 | [[444,154,16]] | Complete both-sector exclusion through 15 |
| CPM (4,12) P43 | [[516,178,16]] | Complete both-sector exclusion through 15 |
| CPM (4,12) P47 | [[564,194,14]] | Complete both-sector exclusion through 13 |
| CPM (4,12) P53 | [[636,218,20]] | Complete both-sector exclusion through 19 |
| Power-9 P16 | [[128,36,6]] | Complete both-sector exclusion through 5 |
| Power-9 P32 | [[256,68,12]] | Complete both-sector exclusion through 11 |
| Power-9 P64 | [[512,132,16]] | Complete both-sector exclusion through 15 |
| Power-9 P128 | [[1024,260,18]] | Complete both-sector exclusion through 17 |
| Power-9 P256 | [[2048,516,20]] | Complete both-sector exclusion through 19 |
| Power-9 P512 (failed extension) | [[4096,1088,4]] | Complete both-sector exclusion through 3 |

Each row also has a physical nontrivial logical operator of the stated weight.
`index.json` binds each claim to its exact matrix file and SHA-256. Equal n and k
do not identify the same code. The original construction provenance remains in
the repository README.

## Download and verify

The large raw records are **GitHub Release assets**, not ordinary Git blobs:
https://github.com/kasaikenta/apm-pp-distance-extensions/releases/tag/v1.3.0

Clone this repository at tag `v1.3.0`. With Python 3.10+ and GitHub CLI:

```sh
mkdir -p downloaded-evidence
gh release download v1.3.0 --repo kasaikenta/apm-pp-distance-extensions \
  --pattern '*.tar.gz' --dir downloaded-evidence
python3 exact-distance/verify.py --assets downloaded-evidence
for archive in downloaded-evidence/*.tar.gz; do
  tar -xzf "$archive" -C downloaded-evidence
done
python3 exact-distance/audit_release.py \
  --evidence downloaded-evidence/exact-distance-evidence --case all
```

The first verifier checks the public matrices, ranks, commuting checks, physical
logical witnesses, and all downloaded asset hashes. The second audits the
archived native receipts, input/source/executable identities where recorded,
partition coverage, coordinate permutations, parity and X/Z exchange. It does
**not** start new distance searches, execute archived binaries, connect to a
cluster, or purchase compute time. P384 has 286,445 receipts and its audit takes
longer than the others. Case selection: `standard`, `p192`, `p384`, `cpm10-P37`,
`cpm10-P73`, or `all`.

All four archives unpack into the same `exact-distance-evidence/` directory.
`common.tar.gz` contains the pinned audit software and physical matrices;
`apm-P192.tar.gz` and `apm-P384.tar.gz` contain adaptive partition trees and
SQLite receipt archives; `other-exact-codes.tar.gz` contains the remaining
campaigns and the original P37 shard certificates. `files.json` records the
SHA-256 of every uncompressed file. Archive hard links only deduplicate identical
bytes. Original historical records are preserved, including incomplete attempts.

These are **exhaustive-execution records**, not a formally checked per-node proof
trace. Hashes alone do not establish a distance. The auditors check the recorded
coverage and computation identities; independently rerunning the enumerations
uses the bundled inputs and search sources. The P37 lower bound is inherited
from the original PP-verifier computation: its 5,005 shard certificates, parity,
duality and coordinate-equivalence evidence are supplied; it is not described as
a newly executed QEDist exclusion. Its legacy certificates do not bind a native
source/executable digest; the release does not claim to recover that historical
binary. An independent new exclusion can instead use the bundled general search
engine on the supplied P37 physical matrices. The earlier frozen release tags
are unchanged.

## Algorithm provenance

Support-growing searches for sparse codes are described in Dumer, Kovalev and
Pryadko, *IEEE Transactions on Information Theory* 63, 4675–4686 (2017),
https://doi.org/10.1109/TIT.2017.2690381 . PP-specific distance verification and
symmetry reduction are described in
https://arxiv.org/abs/2607.14091v4 and https://arxiv.org/abs/2609.30069v1 .
The recorded later computations use the QEDist implementation and its pruning
rules, https://arxiv.org/abs/2611.21709v1 , at commit
`aeb7006147b777c121fd3ab8c6841da99238a84e` of
https://github.com/yeon-lee/QEDist . Its MIT license and notices are included.
Adaptive coverage auditing and the release verification scripts are supplied
with these data. Proposed two-generator and T-join enhancements in the draft
were not used to establish the archived distances.

## 日本語

最小距離が確定した15符号について、物理検査行列、同じ重みの非自明論理演算子、
下界を与える全探索の記録、分割範囲の監査手順を公開します。距離が未確定の
P127はこの15符号に含めません。監査は保存記録の整合性と全範囲の網羅を検査し、
距離探索を新しく実行するものではありません。記録は形式検証済みの逐次証明では
なく、独立した探索の再実行には同梱した入力とソースを用います。
