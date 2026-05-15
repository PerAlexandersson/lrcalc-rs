use lrcalc::lr_gt::{lrcoef_gt_stats, lrcoef_gt_u128};
use lrcalc::lrcoef::lrcoef;
use std::hint::black_box;
use std::time::{Duration, Instant};

struct Case {
    label: &'static str,
    outer: &'static [i32],
    inner: &'static [i32],
    content: &'static [i32],
}

fn main() {
    let mut repeat = 10_000usize;
    let mut suite = "mixed";
    for arg in std::env::args().skip(1) {
        if let Ok(value) = arg.parse::<usize>() {
            repeat = value;
        } else {
            suite = Box::leak(arg.into_boxed_str());
        }
    }

    let cases = cases(suite);
    for case in cases {
        let buch = lrcoef(case.outer, case.inner, case.content)
            .unwrap_or_else(|_| panic!("Buch LR failed for {}", case.label));
        let gt = lrcoef_gt_u128(case.outer, case.inner, case.content)
            .unwrap_or_else(|_| panic!("GT LR failed for {}", case.label));
        if buch != gt {
            eprintln!("mismatch: {}", case.label);
            eprintln!("outer:   {:?}", case.outer);
            eprintln!("inner:   {:?}", case.inner);
            eprintln!("content: {:?}", case.content);
            eprintln!("Buch: {buch}");
            eprintln!("GT:   {gt}");
            std::process::exit(1);
        }
    }

    println!("correctness: ok ({} cases)", cases.len());
    println!("suite: {suite}");
    println!("repeat: {repeat}");
    print_case_stats(cases);

    let (buch_time, buch_sink) = time_loop(repeat, cases, |case| {
        lrcoef(case.outer, case.inner, case.content)
            .unwrap_or_else(|_| panic!("Buch LR failed for {}", case.label))
    });
    let (gt_time, gt_sink) = time_loop(repeat, cases, |case| {
        lrcoef_gt_u128(case.outer, case.inner, case.content)
            .unwrap_or_else(|_| panic!("GT LR failed for {}", case.label))
    });
    black_box((buch_sink, gt_sink));

    println!(
        "Buch-port LR: {}  ({} evals)",
        format_duration(buch_time),
        repeat * cases.len()
    );
    println!(
        "GT-chain LR:  {}  ({} evals)",
        format_duration(gt_time),
        repeat * cases.len()
    );
    println!(
        "comparison: GT/Buch = {:.3}x",
        gt_time.as_secs_f64() / buch_time.as_secs_f64()
    );
}

fn time_loop<F>(repeat: usize, cases: &[Case], mut f: F) -> (Duration, u128)
where
    F: FnMut(&Case) -> u128,
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

fn print_case_stats(cases: &[Case]) {
    println!("GT peak states:");
    for case in cases {
        let stats = lrcoef_gt_stats(case.outer, case.inner, case.content)
            .unwrap_or_else(|_| panic!("GT stats failed for {}", case.label));
        println!(
            "  {}: value {}, peak {}, levels {:?}",
            case.label, stats.value, stats.peak_states, stats.levels
        );
    }
}

fn format_duration(duration: Duration) -> String {
    format!("{:.6}s", duration.as_secs_f64())
}

fn cases(suite: &str) -> &'static [Case] {
    match suite {
        "small" => SMALL_CASES,
        "mixed" => MIXED_CASES,
        "high-coeff" => HIGH_COEFF_CASES,
        "large-few-parts" | "large" => LARGE_FEW_PARTS,
        "three-part" => THREE_PART_CASES,
        other => {
            panic!(
                "unknown suite '{other}', expected small, mixed, high-coeff, large-few-parts, or three-part"
            )
        }
    }
}

const SMALL_CASES: &[Case] = &[
    Case {
        label: "s21 s21 coefficient 2",
        outer: &[3, 2, 1],
        inner: &[2, 1],
        content: &[2, 1],
    },
    Case {
        label: "s21 s21 coefficient 1 A",
        outer: &[4, 2],
        inner: &[2, 1],
        content: &[2, 1],
    },
    Case {
        label: "vertical strip short",
        outer: &[4, 3, 2, 1],
        inner: &[3, 2, 1],
        content: &[1, 1, 1],
    },
    Case {
        label: "medium",
        outer: &[5, 4, 3, 2, 1],
        inner: &[3, 2, 1],
        content: &[4, 3, 1, 1],
    },
];

const MIXED_CASES: &[Case] = &[
    Case {
        label: "s21 s21 coefficient 2",
        outer: &[3, 2, 1],
        inner: &[2, 1],
        content: &[2, 1],
    },
    Case {
        label: "upstream oom larger",
        outer: &[7, 6, 5, 4, 3, 2, 1],
        inner: &[4, 4, 3, 2, 1],
        content: &[5, 4, 3, 2],
    },
    Case {
        label: "stretched 321 scale 3",
        outer: &[9, 6, 3],
        inner: &[6, 3],
        content: &[6, 3],
    },
    Case {
        label: "wide compact medium",
        outer: &[30, 20, 10],
        inner: &[20, 10],
        content: &[20, 10],
    },
    Case {
        label: "irregular D",
        outer: &[12, 10, 8, 5, 3, 1],
        inner: &[7, 6, 4, 2],
        content: &[9, 8, 5, 3, 1],
    },
];

const LARGE_FEW_PARTS: &[Case] = &[
    Case {
        label: "2-row wide A",
        outer: &[1000, 700],
        inner: &[600],
        content: &[700, 400],
    },
    Case {
        label: "2-row wide B",
        outer: &[2000, 1000],
        inner: &[1500],
        content: &[1000, 500],
    },
    Case {
        label: "2-row wide C",
        outer: &[5000, 3000],
        inner: &[3500],
        content: &[3000, 1500],
    },
    Case {
        label: "3-row balanced A",
        outer: &[300, 200, 100],
        inner: &[200, 100],
        content: &[200, 100],
    },
    Case {
        label: "3-row balanced B",
        outer: &[450, 300, 150],
        inner: &[300, 150],
        content: &[300, 150],
    },
];

const HIGH_COEFF_CASES: &[Case] = &[
    Case {
        label: "(20,15,10)^2 top coeff",
        outer: &[30, 24, 18, 12, 6],
        inner: &[20, 15, 10],
        content: &[20, 15, 10],
    },
    Case {
        label: "(40,30,20)^2 top coeff",
        outer: &[61, 48, 37, 23, 11],
        inner: &[40, 30, 20],
        content: &[40, 30, 20],
    },
];

const THREE_PART_CASES: &[Case] = &[
    Case {
        label: "3-part stretched scale 40",
        outer: &[120, 80, 40],
        inner: &[60, 40, 20],
        content: &[60, 40, 20],
    },
    Case {
        label: "3-part stretched scale 60",
        outer: &[180, 120, 60],
        inner: &[90, 60, 30],
        content: &[90, 60, 30],
    },
    Case {
        label: "3-part stretched scale 80",
        outer: &[240, 160, 80],
        inner: &[120, 80, 40],
        content: &[120, 80, 40],
    },
];
