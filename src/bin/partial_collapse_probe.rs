use lrcalc::kostka::kostka_lr_triple;
use lrcalc::lr_gt::{
    lrcoef_gt_partial_collapse_rows, lrcoef_gt_partial_collapse_stats,
    lrcoef_gt_yamanouchi_mask_certified_stats, lrcoef_gt_yamanouchi_mask_stats,
    LrGtYamanouchiMaskCertifiedStats, LrGtYamanouchiMaskStats,
};

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
        "case\trows\tvalue\tall_peak\tempty_value\tempty_peak\tcandidate\tcandidate_value\tcandidate_peak\tcandidate_cert\tsafe_mode\tsafe_peak\tbest_mask\tbest_peak"
    );

    for case in &cases {
        print_case(case);
    }

    print_small_survey(8);
}

fn print_case(case: &Case) {
    let rows = partition_length(&case.outer);
    let all_rows = (0..rows).collect::<Vec<_>>();
    let full = mask_stats(case, &all_rows);
    let empty = mask_stats(case, &[]);
    let candidate = lrcoef_gt_partial_collapse_rows(&case.outer, &case.inner)
        .unwrap_or_else(|_| panic!("candidate rows failed for {}", case.label));
    let candidate_certified = certified_mask_stats(case, &candidate);
    let safe = lrcoef_gt_partial_collapse_stats(&case.outer, &case.inner, &case.content)
        .unwrap_or_else(|_| panic!("safe partial collapse failed for {}", case.label));
    let best = best_matching_mask(case, full.value, rows);
    let (best_mask, best_peak) = best.as_ref().map_or_else(
        || ("skip".to_string(), 0),
        |(mask, stats)| (format_mask(mask), stats.peak_states),
    );

    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        case.label,
        rows,
        full.value,
        full.peak_states,
        empty.value,
        empty.peak_states,
        format_mask(&candidate),
        candidate_certified.stats.value,
        candidate_certified.stats.peak_states,
        candidate_certified.omitted_constraints_forced,
        safe.mode.label(),
        safe.stats.peak_states,
        best_mask,
        best_peak
    );
}

fn print_small_survey(max_outer_size: i32) {
    let mut triples = 0usize;
    let mut nonempty_candidates = 0usize;
    let mut exact_candidates = 0usize;
    let mut improved_exact_candidates = 0usize;
    let mut certified_candidates = 0usize;
    let mut certified_inexact_candidates = 0usize;
    let mut certified_improved_candidates = 0usize;
    let mut peak_saved = 0usize;
    let mut first_certified_inexact = None::<Case>;

    for outer_size in 0..=max_outer_size {
        for outer in partitions_of(outer_size) {
            let rows = partition_length(&outer);
            let all_rows = (0..rows).collect::<Vec<_>>();
            for inner_size in 0..=outer_size {
                for inner in partitions_of(inner_size) {
                    if !partition_less_equal(&inner, &outer) {
                        continue;
                    }
                    let content_size = outer_size - inner_size;
                    for content in partitions_of(content_size) {
                        let case = Case {
                            label: "survey".to_string(),
                            outer: outer.clone(),
                            inner: inner.clone(),
                            content,
                        };
                        let full = mask_stats(&case, &all_rows);
                        let candidate =
                            lrcoef_gt_partial_collapse_rows(&case.outer, &case.inner).unwrap();
                        let candidate_stats = certified_mask_stats(&case, &candidate);
                        triples += 1;
                        if !candidate.is_empty() {
                            nonempty_candidates += 1;
                        }
                        if candidate_stats.stats.value == full.value {
                            exact_candidates += 1;
                            if candidate_stats.stats.peak_states < full.peak_states {
                                improved_exact_candidates += 1;
                                peak_saved += full.peak_states - candidate_stats.stats.peak_states;
                            }
                        }
                        if candidate_stats.omitted_constraints_forced {
                            certified_candidates += 1;
                            if candidate_stats.stats.value != full.value {
                                certified_inexact_candidates += 1;
                                if first_certified_inexact.is_none() {
                                    first_certified_inexact = Some(case.clone());
                                }
                            }
                            if candidate_stats.stats.peak_states < full.peak_states {
                                certified_improved_candidates += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    println!("small_survey_max_outer_size: {max_outer_size}");
    println!("small_survey_triples: {triples}");
    println!("small_survey_nonempty_candidates: {nonempty_candidates}");
    println!("small_survey_exact_candidates: {exact_candidates}");
    println!("small_survey_improved_exact_candidates: {improved_exact_candidates}");
    println!("small_survey_certified_candidates: {certified_candidates}");
    println!("small_survey_certified_inexact_candidates: {certified_inexact_candidates}");
    println!("small_survey_certified_improved_candidates: {certified_improved_candidates}");
    println!("small_survey_total_peak_saved: {peak_saved}");
    if let Some(case) = first_certified_inexact {
        println!(
            "small_survey_first_certified_inexact: outer={:?} inner={:?} content={:?}",
            case.outer, case.inner, case.content
        );
    }
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

fn certified_mask_stats(case: &Case, mask: &[usize]) -> LrGtYamanouchiMaskCertifiedStats {
    lrcoef_gt_yamanouchi_mask_certified_stats(&case.outer, &case.inner, &case.content, mask)
        .unwrap_or_else(|_| panic!("certified masked LR failed for {}", case.label))
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

fn partition_length(parts: &[i32]) -> usize {
    parts
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1)
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

fn partitions_of(n: i32) -> Vec<Vec<i32>> {
    fn go(remaining: i32, max_part: i32, current: &mut Vec<i32>, out: &mut Vec<Vec<i32>>) {
        if remaining == 0 {
            out.push(current.clone());
            return;
        }
        for part in (1..=remaining.min(max_part)).rev() {
            current.push(part);
            go(remaining - part, part, current, out);
            current.pop();
        }
    }

    let mut out = Vec::new();
    go(n, n, &mut Vec::new(), &mut out);
    out
}

fn partition_less_equal(inner: &[i32], outer: &[i32]) -> bool {
    let len = inner.len().max(outer.len());
    for index in 0..len {
        if *inner.get(index).unwrap_or(&0) > *outer.get(index).unwrap_or(&0) {
            return false;
        }
    }
    true
}
