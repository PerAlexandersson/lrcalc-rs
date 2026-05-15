use lrcalc::kostka::kostka_lr_triple;
use lrcalc::kostka_fast::{kostka_counts_stats, KostkaCountsStats};
use lrcalc::lr_gt::{lrcoef_gt_counts_stats, LrGtCountsStats};
use std::hint::black_box;
use std::time::{Duration, Instant};

#[derive(Clone)]
struct BaseCase {
    label: &'static str,
    shape: &'static [i32],
    weight: &'static [i32],
}

#[derive(Clone)]
struct ScaledCase {
    label: String,
    shape: Vec<i32>,
    weight: Vec<i32>,
    outer: Vec<i32>,
    inner: Vec<i32>,
    content: Vec<i32>,
}

fn main() {
    let mut repeat = 5usize;
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
    println!("suite: stretched_dp_bench");
    println!("repeat: {repeat}");
    println!("max_scale: {max_scale}");
    verify_cases(&cases);
    print_case_stats(&cases);

    let (kostka_time, kostka_sink) = time_loop(repeat, &cases, |case| {
        let stats = kostka_counts_stats(&case.shape, &case.weight)
            .unwrap_or_else(|_| panic!("Kostka counts failed for {}", case.label));
        stats.full ^ stats.interior
    });
    let (lr_time, lr_sink) = time_loop(repeat, &cases, |case| {
        let stats = lrcoef_gt_counts_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("LR GT counts failed for {}", case.label));
        stats.full ^ stats.interior
    });
    black_box((kostka_sink, lr_sink));

    println!(
        "Kostka DP counts: {}  ({} evals)",
        format_duration(kostka_time),
        repeat * cases.len()
    );
    println!(
        "LR GT DP counts:  {}  ({} evals)",
        format_duration(lr_time),
        repeat * cases.len()
    );
    println!(
        "comparison: LR/Kostka = {:.3}x",
        lr_time.as_secs_f64() / kostka_time.as_secs_f64()
    );
}

fn scaled_cases(max_scale: i32) -> Vec<ScaledCase> {
    let mut cases = Vec::new();
    for base in BASE_CASES {
        let (outer, inner, content) = kostka_lr_triple(base.shape, base.weight)
            .unwrap_or_else(|_| panic!("Kostka-to-LR conversion failed for {}", base.label));
        for scale in 1..=max_scale {
            cases.push(ScaledCase {
                label: format!("{} x{}", base.label, scale),
                shape: scale_parts(base.shape, scale),
                weight: scale_parts(base.weight, scale),
                outer: scale_parts(&outer, scale),
                inner: scale_parts(&inner, scale),
                content: scale_parts(&content, scale),
            });
        }
    }
    cases
}

fn verify_cases(cases: &[ScaledCase]) {
    for case in cases {
        let kostka = kostka_counts_stats(&case.shape, &case.weight)
            .unwrap_or_else(|_| panic!("Kostka counts failed for {}", case.label));
        let lr = lrcoef_gt_counts_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("LR GT counts failed for {}", case.label));
        assert_eq!(kostka.full, lr.full, "full mismatch for {}", case.label);
        assert_eq!(
            kostka.interior, lr.interior,
            "interior mismatch for {}",
            case.label
        );
    }
    println!("correctness: ok ({} scaled cases)", cases.len());
}

fn print_case_stats(cases: &[ScaledCase]) {
    println!("case\tfull\tinterior\tk_full_peak\tk_int_peak\tlr_full_peak\tlr_int_peak");
    for case in cases {
        let kostka = kostka_counts_stats(&case.shape, &case.weight)
            .unwrap_or_else(|_| panic!("Kostka counts failed for {}", case.label));
        let lr = lrcoef_gt_counts_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("LR GT counts failed for {}", case.label));
        print_stats(case, &kostka, &lr);
    }
}

fn print_stats(case: &ScaledCase, kostka: &KostkaCountsStats, lr: &LrGtCountsStats) {
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        case.label,
        kostka.full,
        kostka.interior,
        kostka.full_peak_states,
        kostka.interior_peak_states,
        lr.full_peak_states,
        lr.interior_peak_states
    );
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

fn scale_parts(parts: &[i32], scale: i32) -> Vec<i32> {
    parts.iter().map(|part| part * scale).collect()
}

fn format_duration(duration: Duration) -> String {
    format!("{:.6}s", duration.as_secs_f64())
}

const BASE_CASES: &[BaseCase] = &[
    BaseCase {
        label: "tiny 321 / 222",
        shape: &[3, 2, 1],
        weight: &[2, 2, 2],
    },
    BaseCase {
        label: "small sparse",
        shape: &[4, 2, 1],
        weight: &[3, 2, 1, 1],
    },
    BaseCase {
        label: "balanced 642 / 444",
        shape: &[6, 4, 2],
        weight: &[4, 4, 4],
    },
    BaseCase {
        label: "five-part medium",
        shape: &[5, 4, 2, 1],
        weight: &[4, 3, 2, 2, 1],
    },
    BaseCase {
        label: "five-part large",
        shape: &[8, 6, 4, 2],
        weight: &[6, 5, 4, 3, 2],
    },
];
