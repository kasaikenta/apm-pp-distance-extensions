# Default-rule CPM–PP extensions, column weight 4 and check weight 12

Starting from the published [[276,98,14]] code, the unchanged general extraction
and construction algorithms give [[444,154,16]] and [[636,218,20]]. Parameters
are frozen once at P=23; no descendant-based screening or reselection is used.
The P=37 code here differs from the published [[444,154,18]] code.

Two requested shorter sizes were also tested with exactly the same frozen data:
P=43 gives [[516,178,16]], and P=47 gives [[564,194,14]]. Neither retains the
distance 20 of the P=53 reference. Both tested sizes are included, with no
parameter reselection. The P=43 code here differs from the published
[[516,178,20]] code; equal n and k do not identify a code.

Run `python3 default-cpm12/verify.py` from the repository root. It checks hashes,
re-extracts the frozen rule from the original seed, reconstructs every physical
matrix, and verifies ranks, orthogonality, row/column weights and nontrivial
logical operators. All indices are zero-based; q=lP+u.

The exact distances additionally rely on independent audits of the complete
both-sector QEDist exclusions through weight d-1 for every listed exact distance d.
The compact audit summary is included; full native search records are archived separately.
This script checks the summary's consistency, not the underlying search again.
The original sizes P=23,37,53 demonstrate distance improvement. The additional
P=43,47 results explicitly show that distance is not monotone in P under this
same fixed rule. Source attribution is in provenance.json.
