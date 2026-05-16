use lrcalc::kostka::kostka_lr_triple;
use lrcalc::lr_gt::{
    lrcoef_gt_hybrid_stats, lrcoef_gt_stats, lrcoef_tableau_hybrid_stats, LrGtHybridStats,
    LrGtStats, LrTableauHybridStats,
};
use lrcalc::lrcoef::lrcoef;
use std::hint::black_box;
use std::time::{Duration, Instant};

#[derive(Clone)]
struct BaseCase {
    label: &'static str,
    outer: Vec<i32>,
    inner: Vec<i32>,
    content: Vec<i32>,
}

#[derive(Clone)]
struct ScaledCase {
    label: String,
    outer: Vec<i32>,
    inner: Vec<i32>,
    content: Vec<i32>,
}

struct CaseStats {
    gt: LrGtStats,
    gt_hybrid: LrGtHybridStats,
    tableau_hybrid: LrTableauHybridStats,
    buch: u128,
}

fn main() {
    let mut repeat = 10usize;
    let mut max_scale = 3i32;
    let mut args = std::env::args().skip(1);
    if let Some(arg) = args.next() {
        repeat = arg
            .parse::<usize>()
            .unwrap_or_else(|_| panic!("expected integer repeat count, got {arg}"));
    }
    if let Some(arg) = args.next() {
        max_scale = arg
            .parse::<i32>()
            .unwrap_or_else(|_| panic!("expected integer max scale, got {arg}"));
    }
    assert!(repeat > 0, "repeat must be positive");
    assert!(max_scale > 0, "max scale must be positive");

    let cases = scaled_cases(max_scale);
    println!("suite: lr_hybrid_bench");
    println!("repeat: {repeat}");
    println!("max_scale: {max_scale}");
    verify_cases(&cases);
    print_case_stats(&cases);

    let (buch_time, buch_sink) = time_loop(repeat, &cases, |case| {
        lrcoef(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch LR failed for {}", case.label))
    });
    let (gt_time, gt_sink) = time_loop(repeat, &cases, |case| {
        lrcoef_gt_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT LR failed for {}", case.label))
            .value
    });
    let (gt_hybrid_time, gt_hybrid_sink) = time_loop(repeat, &cases, |case| {
        lrcoef_gt_hybrid_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT hybrid LR failed for {}", case.label))
            .stats
            .value
    });
    let (tableau_hybrid_time, tableau_hybrid_sink) = time_loop(repeat, &cases, |case| {
        lrcoef_tableau_hybrid_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("tableau hybrid LR failed for {}", case.label))
            .value
    });
    black_box((buch_sink, gt_sink, gt_hybrid_sink, tableau_hybrid_sink));

    println!(
        "Buch full:     {}  ({} evals)",
        format_duration(buch_time),
        repeat * cases.len()
    );
    println!(
        "GT-chain full: {}  ({} evals)",
        format_duration(gt_time),
        repeat * cases.len()
    );
    println!(
        "GT hybrid:     {}  ({} evals)",
        format_duration(gt_hybrid_time),
        repeat * cases.len()
    );
    println!(
        "Tableau hybrid:{}  ({} evals)",
        format_duration(tableau_hybrid_time),
        repeat * cases.len()
    );
    println!(
        "comparison: GT/tableau-hybrid = {:.3}x",
        gt_time.as_secs_f64() / tableau_hybrid_time.as_secs_f64()
    );
    println!(
        "comparison: Buch/tableau-hybrid = {:.3}x",
        buch_time.as_secs_f64() / tableau_hybrid_time.as_secs_f64()
    );
    println!(
        "comparison: GT-hybrid/tableau-hybrid = {:.3}x",
        gt_hybrid_time.as_secs_f64() / tableau_hybrid_time.as_secs_f64()
    );
}

fn scaled_cases(max_scale: i32) -> Vec<ScaledCase> {
    let mut cases = Vec::new();
    for base in base_cases() {
        for scale in 1..=max_scale {
            cases.push(ScaledCase {
                label: format!("{} x{}", base.label, scale),
                outer: scale_parts(&base.outer, scale),
                inner: scale_parts(&base.inner, scale),
                content: scale_parts(&base.content, scale),
            });
        }
    }
    cases
}

fn verify_cases(cases: &[ScaledCase]) {
    for case in cases {
        let stats = case_stats(case);
        assert_eq!(
            stats.gt.value, stats.gt_hybrid.stats.value,
            "GT hybrid mismatch for {}",
            case.label
        );
        assert_eq!(
            stats.buch, stats.tableau_hybrid.value,
            "tableau hybrid mismatch for {}",
            case.label
        );
    }
    println!("correctness: ok ({} scaled cases)", cases.len());
}

fn print_case_stats(cases: &[ScaledCase]) {
    println!(
        "case\tvalue\tgt_peak\tgt_hybrid_mode\tgt_hybrid_peak\ttableau_mode\ttableau_peak\tenforced_rows"
    );
    for case in cases {
        let stats = case_stats(case);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}",
            case.label,
            stats.gt.value,
            stats.gt.peak_states,
            stats.gt_hybrid.mode.label(),
            stats.gt_hybrid.stats.peak_states,
            stats.tableau_hybrid.mode.label(),
            format_optional_usize(stats.tableau_hybrid.peak_states),
            stats.tableau_hybrid.enforced_rows
        );
    }
}

fn case_stats(case: &ScaledCase) -> CaseStats {
    CaseStats {
        gt: lrcoef_gt_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT LR failed for {}", case.label)),
        gt_hybrid: lrcoef_gt_hybrid_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT hybrid LR failed for {}", case.label)),
        tableau_hybrid: lrcoef_tableau_hybrid_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("tableau hybrid LR failed for {}", case.label)),
        buch: lrcoef(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch LR failed for {}", case.label)),
    }
}

fn time_loop<F>(repeat: usize, cases: &[ScaledCase], mut f: F) -> (Duration, u128)
where
    F: FnMut(&ScaledCase) -> u128,
{
    let start = Instant::now();
    let mut sink = 0u128;
    for _ in 0..repeat {
        for case in cases {
            sink ^= black_box(f(case));
        }
    }
    (start.elapsed(), sink)
}

fn base_cases() -> Vec<BaseCase> {
    let mut cases = Vec::new();
    push_kostka_case(&mut cases, "exact medium", &[5, 4, 2, 1], &[4, 3, 2, 2, 1]);
    push_kostka_case(&mut cases, "exact large", &[8, 6, 4, 2], &[6, 5, 4, 3, 2]);
    cases.push(BaseCase {
        label: "one tail defect",
        outer: vec![4, 2],
        inner: vec![2, 1],
        content: vec![2, 1],
    });
    cases.push(BaseCase {
        label: "left extension defect",
        outer: vec![7, 4, 2, 1],
        inner: vec![4, 2],
        content: vec![5, 2, 1],
    });
    cases.push(BaseCase {
        label: "right gap defect",
        outer: vec![7, 4, 2, 1],
        inner: vec![4, 3, 1],
        content: vec![3, 2, 1],
    });
    cases.push(BaseCase {
        label: "irregular mixed",
        outer: vec![7, 6, 5, 4, 3, 2, 1],
        inner: vec![4, 4, 3, 2, 1],
        content: vec![5, 4, 3, 2],
    });
    cases
}

fn push_kostka_case(cases: &mut Vec<BaseCase>, label: &'static str, shape: &[i32], weight: &[i32]) {
    let (outer, inner, content) = kostka_lr_triple(shape, weight)
        .unwrap_or_else(|_| panic!("Kostka-to-LR conversion failed for {label}"));
    cases.push(BaseCase {
        label,
        outer,
        inner,
        content,
    });
}

fn scale_parts(parts: &[i32], scale: i32) -> Vec<i32> {
    parts.iter().map(|part| part * scale).collect()
}

fn format_duration(duration: Duration) -> String {
    format!("{:.6}s", duration.as_secs_f64())
}

fn format_optional_usize(value: Option<usize>) -> String {
    value.map_or_else(|| "-".to_string(), |value| value.to_string())
}
