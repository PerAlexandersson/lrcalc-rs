use lrcalc::kostka::kostka_lr_triple;
use lrcalc::lr_gt::{lrcoef_gt_yamanouchi_mask_stats, LrGtYamanouchiMaskStats};

#[derive(Clone)]
struct Case {
    label: String,
    outer: Vec<i32>,
    inner: Vec<i32>,
    content: Vec<i32>,
}

fn main() {
    let cases = cases();
    println!("suite: partial_collapse_probe");
    println!(
        "case\trows\tvalue\tall_peak\tempty_value\tempty_peak\tdefects\tdefect_value\tdefect_peak\tbest_mask\tbest_peak"
    );

    for case in &cases {
        print_case(case);
    }
}

fn print_case(case: &Case) {
    let rows = partition_length(&case.outer);
    let all_rows = (0..rows).collect::<Vec<_>>();
    let full = mask_stats(case, &all_rows);
    let empty = mask_stats(case, &[]);
    let defects = defect_rows(&case.outer, &case.inner);
    let defect_stats = mask_stats(case, &defects);
    let best = best_matching_mask(case, full.value, rows);
    let (best_mask, best_peak) = best.as_ref().map_or_else(
        || ("skip".to_string(), 0),
        |(mask, stats)| (format_mask(mask), stats.peak_states),
    );

    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        case.label,
        rows,
        full.value,
        full.peak_states,
        empty.value,
        empty.peak_states,
        format_mask(&defects),
        defect_stats.value,
        defect_stats.peak_states,
        best_mask,
        best_peak
    );
}

fn best_matching_mask(
    case: &Case,
    full_value: u128,
    rows: usize,
) -> Option<(Vec<usize>, LrGtYamanouchiMaskStats)> {
    if rows > 10 {
        return None;
    }

    let mut best = None::<(Vec<usize>, LrGtYamanouchiMaskStats)>;
    for bits in 0u64..(1u64 << rows) {
        let mask = mask_from_bits(bits, rows);
        let stats = mask_stats(case, &mask);
        if stats.value != full_value {
            continue;
        }
        let replace = best.as_ref().is_none_or(|(best_mask, best_stats)| {
            mask.len() < best_mask.len()
                || (mask.len() == best_mask.len() && stats.peak_states < best_stats.peak_states)
        });
        if replace {
            best = Some((mask, stats));
        }
    }
    best
}

fn mask_stats(case: &Case, mask: &[usize]) -> LrGtYamanouchiMaskStats {
    lrcoef_gt_yamanouchi_mask_stats(&case.outer, &case.inner, &case.content, mask)
        .unwrap_or_else(|_| panic!("masked LR failed for {}", case.label))
}

fn cases() -> Vec<Case> {
    let mut cases = Vec::new();
    push_kostka_case(&mut cases, "exact tiny", &[3, 2, 1], &[2, 2, 2], 1);
    push_kostka_case(&mut cases, "exact sparse", &[4, 2, 1], &[3, 2, 1, 1], 1);
    push_kostka_case(
        &mut cases,
        "exact medium x2",
        &[5, 4, 2, 1],
        &[4, 3, 2, 2, 1],
        2,
    );

    cases.push(Case {
        label: "one tail defect".to_string(),
        outer: vec![4, 2],
        inner: vec![2, 1],
        content: vec![2, 1],
    });
    cases.push(Case {
        label: "left extension defect".to_string(),
        outer: vec![7, 4, 2, 1],
        inner: vec![4, 2],
        content: vec![5, 2, 1],
    });
    cases.push(Case {
        label: "right gap defect".to_string(),
        outer: vec![7, 4, 2, 1],
        inner: vec![4, 3, 1],
        content: vec![3, 2, 1],
    });
    cases.push(Case {
        label: "irregular mixed".to_string(),
        outer: vec![7, 6, 5, 4, 3, 2, 1],
        inner: vec![4, 4, 3, 2, 1],
        content: vec![5, 4, 3, 2],
    });

    cases
}

fn push_kostka_case(cases: &mut Vec<Case>, label: &str, shape: &[i32], weight: &[i32], scale: i32) {
    let (outer, inner, content) = kostka_lr_triple(shape, weight)
        .unwrap_or_else(|_| panic!("Kostka-to-LR conversion failed for {label}"));
    cases.push(Case {
        label: label.to_string(),
        outer: scale_parts(&outer, scale),
        inner: scale_parts(&inner, scale),
        content: scale_parts(&content, scale),
    });
}

fn defect_rows(outer: &[i32], inner: &[i32]) -> Vec<usize> {
    let rows = partition_length(outer);
    (0..rows)
        .filter(|&row| part(inner, row) != part(outer, row + 1))
        .collect()
}

fn partition_length(parts: &[i32]) -> usize {
    parts
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1)
}

fn part(parts: &[i32], row: usize) -> i32 {
    parts.get(row).copied().unwrap_or(0)
}

fn mask_from_bits(bits: u64, rows: usize) -> Vec<usize> {
    (0..rows).filter(|&row| bits & (1u64 << row) != 0).collect()
}

fn format_mask(mask: &[usize]) -> String {
    if mask.is_empty() {
        return "-".to_string();
    }
    mask.iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn scale_parts(parts: &[i32], scale: i32) -> Vec<i32> {
    parts.iter().map(|part| part * scale).collect()
}
