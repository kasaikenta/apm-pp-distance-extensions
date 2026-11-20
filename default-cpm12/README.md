# Default-rule CPM–PP extensions, column weight 4 and check weight 12

Starting from the published [[276,98,14]] code, the unchanged general extraction
and construction algorithms give [[444,154,16]] and [[636,218,20]]. Parameters
are frozen once at P=23; no descendant-based screening or reselection is used.
The P=37 code here differs from the published [[444,154,18]] code.

Run `python3 default-cpm12/verify.py` from the repository root. It checks hashes,
re-extracts the frozen rule from the original seed, reconstructs every physical
matrix, and verifies ranks, orthogonality, row/column weights and nontrivial
logical operators. All indices are zero-based; q=lP+u.

The exact distances additionally rely on independent audits of the complete
both-sector QEDist exclusions through weights 13, 15 and 19. The compact audit
summary is included; full native search records are archived separately.
This script checks the summary's consistency, not the underlying search again.
The examples demonstrate improvement at these three sizes, not monotonic growth
for arbitrary inputs or sizes. Source attribution is in provenance.json.
