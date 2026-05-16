use lrcalc::lr_shifted::{lrcoef_shifted_interval_stats, lrcoef_shifted_interval_u128};
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
    let mut repeat = 100usize;
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
        let shifted = lrcoef_shifted_interval_u128(case.outer, case.inner, case.content)
            .unwrap_or_else(|_| panic!("shifted interval LR failed for {}", case.label));
        if buch != shifted {
            eprintln!("mismatch: {}", case.label);
            eprintln!("outer:   {:?}", case.outer);
            eprintln!("inner:   {:?}", case.inner);
            eprintln!("content: {:?}", case.content);
            eprintln!("Buch:    {buch}");
            eprintln!("Shifted: {shifted}");
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
    let (shifted_time, shifted_sink) = time_loop(repeat, cases, |case| {
        lrcoef_shifted_interval_u128(case.outer, case.inner, case.content)
            .unwrap_or_else(|_| panic!("shifted interval LR failed for {}", case.label))
    });
    black_box((buch_sink, shifted_sink));

    println!(
        "Buch-port LR:       {}  ({} evals)",
        format_duration(buch_time),
        repeat * cases.len()
    );
    println!(
        "Shifted interval:   {}  ({} evals)",
        format_duration(shifted_time),
        repeat * cases.len()
    );
    println!(
        "comparison: shifted/Buch = {:.3}x",
        shifted_time.as_secs_f64() / buch_time.as_secs_f64()
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
    println!("shifted interval states:");
    for case in cases {
        let stats = lrcoef_shifted_interval_stats(case.outer, case.inner, case.content)
            .unwrap_or_else(|_| panic!("shifted stats failed for {}", case.label));
        println!(
            "  {}: value {}, partitions {}, peak {}, total {}, layers {:?}, fixed {:?}",
            case.label,
            stats.value,
            stats.interval_partitions,
            stats.peak_layer_states,
            stats.total_pair_states,
            stats.layer_states,
            stats.fixed
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
        other => panic!("unknown suite '{other}', expected small or mixed"),
    }
}

const SMALL_CASES: &[Case] = &[
    Case {
        label: "s1 s1 to s2",
        outer: &[2],
        inner: &[1],
        content: &[1],
    },
    Case {
        label: "s1 s1 to s11",
        outer: &[1, 1],
        inner: &[1],
        content: &[1],
    },
    Case {
        label: "s21 s21 coefficient 2",
        outer: &[3, 2, 1],
        inner: &[2, 1],
        content: &[2, 1],
    },
    Case {
        label: "s21 s21 coefficient 1",
        outer: &[4, 2],
        inner: &[2, 1],
        content: &[2, 1],
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
        label: "irregular D",
        outer: &[12, 10, 8, 5, 3, 1],
        inner: &[7, 6, 4, 2],
        content: &[9, 8, 5, 3, 1],
    },
];
