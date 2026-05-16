# Partial Kostka Collapse Probe

## Idea

The GT-chain LR DP stores the current partition plus the prefix sums of the
previous strip.  A partial Kostka collapse keeps only selected Yamanouchi prefix
rows.  If a small row mask gives the same coefficient as the all-row LR mask,
then the omitted prefix rows are unnecessary for that coefficient.

## Current Hook

- `src/lr_gt.rs`: `lrcoef_gt_yamanouchi_mask_stats`
- `src/lr_gt.rs`: `lrcoef_gt_yamanouchi_mask_certified_stats`
- `src/lr_gt.rs`: `lrcoef_gt_partial_collapse_rows`
- `src/lr_gt.rs`: `lrcoef_gt_partial_collapse_stats`
- `src/bin/partial_collapse_probe.rs`: small diagnostic suite

The all-row mask agrees with the usual LR count on all small triples tested.
Smaller masks are relaxations and can overcount.
The certified variant gives a sufficient condition: if every omitted
Yamanouchi inequality is forced on the relaxed reachable transitions, the
masked count is exact.

## Probe Result

Command:

```bash
timeout 60s nice -n 10 cargo run --release --bin partial_collapse_probe
```

Observed on 2026-05-16:

| case | rows | LR value | all peak | candidate | candidate value | certified | safe mode | safe peak |
|---|---:|---:|---:|---|---:|---|---|---:|
| exact tiny | 3 | 2 | 7 | `-` | 15 | false | gt-fallback | 3 |
| exact sparse | 4 | 4 | 11 | `-` | 25 | false | gt-fallback | 9 |
| exact medium x2 | 5 | 106 | 2225 | `-` | 109115 | false | gt-fallback | 1644 |
| one tail defect | 2 | 1 | 2 | `-` | 2 | false | gt-fallback | 1 |
| left extension defect | 4 | 4 | 8 | `0,1` | 4 | true | certified-mask | 6 |
| right gap defect | 4 | 2 | 8 | `-` | 19 | false | gt-fallback | 3 |
| irregular mixed | 7 | 56 | 84 | `-` | 929 | false | gt-fallback | 84 |

Small survey through outer size `8`:

- triples: `4136`
- nonempty candidates: `1400`
- exact candidates: `2724`
- improved exact candidates: `5`
- certified candidates: `2662`
- certified inexact candidates: `0`
- certified improved candidates: `5`

## Takeaway

Exact Kostka translations should keep using the direct packed Kostka dispatch
from the hybrid path.  The masked DP is not as strong there, and the certificate
correctly rejects the empty-mask relaxation.

For near-Kostka shapes, the useful pattern is different: low row masks can
recover the LR coefficient while dropping some prefix coordinates.  The current
candidate only fires on underfull last defects whose mask is at most half of
the row set.  This keeps the useful left-extension cases while avoiding
exact-but-slower masks and expensive failed certificates for tail, excess-gap,
and large irregular defects.

## Next Step

The full-count selector is now:

```text
exact Kostka translation -> packed Kostka DP
certified partial collapse -> row-masked LR/GT DP
otherwise -> ordinary GT-chain LR DP
```

Next:

1. Benchmark the certified candidate on larger near-Kostka families.
2. Extend the same idea to paired full/interior counts only after the full-count
   mask behavior is stable.

Current mixed full-count benchmark:

```bash
timeout 60s nice -n 10 cargo run --release --bin lr_hybrid_bench -- 5 3
```

The production tableau selector is about `382x` faster than raw GT-chain full
counts and about `1.08x` faster than Buch full counts on the exact/near-Kostka
benchmark suite.
