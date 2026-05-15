use lrcalc::kostka::kostka_lr_triple;
use lrcalc::lr_ehrhart::lr_stretch_h_vector;
use lrcalc::lr_gt::{lrcoef_gt_interior_dfs_u128, lrcoef_gt_interior_stats};
use lrcalc::lrcoef::{
    lrcoef, lrcoef_buch_interior_memo_u128, lrcoef_buch_interior_stats, lrcoef_buch_interior_u128,
};
use std::hint::black_box;
use std::time::{Duration, Instant};

struct KostkaCase {
    label: &'static str,
    shape: &'static [i32],
    weight: &'static [i32],
}

struct LrCase {
    label: &'static str,
    outer: Vec<i32>,
    inner: Vec<i32>,
    content: Vec<i32>,
}

fn main() {
    let mut repeat = 100usize;
    for arg in std::env::args().skip(1) {
        repeat = arg
            .parse::<usize>()
            .unwrap_or_else(|_| panic!("expected integer repeat count, got {arg}"));
    }

    let cases = LR_HSTAR_CASES
        .iter()
        .map(|case| {
            let (outer, inner, content) = kostka_lr_triple(case.shape, case.weight)
                .unwrap_or_else(|_| panic!("Kostka-to-LR conversion failed for {}", case.label));
            LrCase {
                label: case.label,
                outer,
                inner,
                content,
            }
        })
        .collect::<Vec<_>>();

    println!("suite: lr_hstar_bench");
    println!("repeat: {repeat}");
    verify_cases(&cases);
    print_case_stats(&cases);

    let (full_time, full_sink) = time_loop(repeat, &cases, |case| {
        lrcoef(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("full LR failed for {}", case.label))
    });
    let (buch_interior_time, buch_sink) = time_loop(repeat, &cases, |case| {
        lrcoef_buch_interior_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch interior failed for {}", case.label))
    });
    let (gt_interior_time, gt_sink) = time_loop(repeat, &cases, |case| {
        lrcoef_gt_interior_dfs_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT interior failed for {}", case.label))
    });
    let (memo_interior_time, memo_sink) = time_loop(repeat, &cases, |case| {
        lrcoef_buch_interior_memo_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch memo interior failed for {}", case.label))
    });
    let (hstar_time, hstar_sink) = time_loop(repeat, &cases, |case| {
        lr_stretch_h_vector(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("h* failed for {}", case.label))
            .h_vector
            .len() as u128
    });
    black_box((full_sink, buch_sink, gt_sink, memo_sink, hstar_sink));

    println!(
        "full Buch LR:       {}  ({} evals)",
        format_duration(full_time),
        repeat * cases.len()
    );
    println!(
        "Buch strict LR:     {}  ({} evals)",
        format_duration(buch_interior_time),
        repeat * cases.len()
    );
    println!(
        "GT strict LR:       {}  ({} evals)",
        format_duration(gt_interior_time),
        repeat * cases.len()
    );
    println!(
        "Buch memo strict:   {}  ({} evals)",
        format_duration(memo_interior_time),
        repeat * cases.len()
    );
    println!(
        "LR h*:              {}  ({} evals)",
        format_duration(hstar_time),
        repeat * cases.len()
    );
    println!(
        "strict comparison: GT/Buch = {:.3}x",
        gt_interior_time.as_secs_f64() / buch_interior_time.as_secs_f64()
    );
    println!(
        "memo comparison: memo/Buch = {:.3}x",
        memo_interior_time.as_secs_f64() / buch_interior_time.as_secs_f64()
    );
}

fn verify_cases(cases: &[LrCase]) {
    for case in cases {
        let buch = lrcoef_buch_interior_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch interior failed for {}", case.label));
        let gt = lrcoef_gt_interior_dfs_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT interior failed for {}", case.label));
        assert_eq!(buch, gt, "strict count mismatch for {}", case.label);
        let memo = lrcoef_buch_interior_memo_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch memo interior failed for {}", case.label));
        assert_eq!(buch, memo, "memo strict count mismatch for {}", case.label);
    }
    println!("correctness: ok ({} cases)", cases.len());
}

fn print_case_stats(cases: &[LrCase]) {
    for case in cases {
        let hstar = lr_stretch_h_vector(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("h* failed for {}", case.label));
        let buch_stats = lrcoef_buch_interior_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch interior stats failed for {}", case.label));
        let gt_stats = lrcoef_gt_interior_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT interior stats failed for {}", case.label));
        println!(
            "{}: dim {}, h* len {}, Buch weak {}, Buch strict {}, GT peak {}",
            case.label,
            hstar.dimension,
            hstar.h_vector.len(),
            buch_stats.weak_tableaux,
            buch_stats.strict_tableaux,
            gt_stats.peak_states
        );
    }
}

fn time_loop<F>(repeat: usize, cases: &[LrCase], mut f: F) -> (Duration, u128)
where
    F: FnMut(&LrCase) -> u128,
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

fn format_duration(duration: Duration) -> String {
    format!("{:.6}s", duration.as_secs_f64())
}

const LR_HSTAR_CASES: &[KostkaCase] = &[
    KostkaCase {
        label: "tiny 321 / 222",
        shape: &[3, 2, 1],
        weight: &[2, 2, 2],
    },
    KostkaCase {
        label: "small sparse",
        shape: &[4, 2, 1],
        weight: &[3, 2, 1, 1],
    },
    KostkaCase {
        label: "stretched 642 / 444",
        shape: &[6, 4, 2],
        weight: &[4, 4, 4],
    },
    KostkaCase {
        label: "five-part medium",
        shape: &[5, 4, 2, 1],
        weight: &[4, 3, 2, 2, 1],
    },
    KostkaCase {
        label: "five-part large",
        shape: &[8, 6, 4, 2],
        weight: &[6, 5, 4, 3, 2],
    },
];
