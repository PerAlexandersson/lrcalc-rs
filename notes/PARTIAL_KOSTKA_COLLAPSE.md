# Partial Kostka Collapse Probe

## Idea

The GT-chain LR DP stores the current partition plus the prefix sums of the
previous strip.  A partial Kostka collapse keeps only selected Yamanouchi prefix
rows.  If a small row mask gives the same coefficient as the all-row LR mask,
then the omitted prefix rows are unnecessary for that coefficient.

## Current Hook

- `src/lr_gt.rs`: `lrcoef_gt_yamanouchi_mask_stats`
- `src/bin/partial_collapse_probe.rs`: small diagnostic suite

The all-row mask agrees with the usual LR count on all small triples tested.
Smaller masks are relaxations and can overcount.

## Probe Result

Command:

```bash
timeout 60s nice -n 10 cargo run --release --bin partial_collapse_probe
```

Observed on 2026-05-16:

| case | rows | LR value | all peak | empty value | empty peak | best mask | best peak |
|---|---:|---:|---:|---:|---:|---|---:|
| exact tiny | 3 | 2 | 7 | 15 | 7 | `0,1` | 7 |
| exact sparse | 4 | 4 | 11 | 25 | 11 | `0,1` | 11 |
| exact medium x2 | 5 | 106 | 2225 | 109115 | 411 | `0,1,2,3` | 1429 |
| one tail defect | 2 | 1 | 2 | 2 | 2 | `0` | 2 |
| left extension defect | 4 | 4 | 8 | 14 | 6 | `0,1` | 6 |
| right gap defect | 4 | 2 | 8 | 19 | 8 | `0,1` | 8 |
| irregular mixed | 7 | 56 | 84 | 929 | 54 | `0,1,2,3,4` | 76 |

## Takeaway

Exact Kostka translations should keep using the direct packed Kostka dispatch.
The masked DP is not as strong there.

For near-Kostka shapes, the useful pattern is different: low row masks often
recover the LR coefficient while dropping some prefix coordinates.  The simple
defect-row mask is not sufficient by itself; the closure seems to include rows
above the first defect, especially rows `0` and `1`.

## Next Step

Turn the probe into a safe heuristic:

1. Build a candidate mask from low rows plus the defect closure.
2. Run the masked DP and a cheap certificate on small/medium cases.
3. Use all-row GT fallback when the mask is not certified.
4. Extend the same idea to paired full/interior counts only after the full-count
   mask behavior is stable.
