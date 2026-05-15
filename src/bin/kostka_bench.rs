use lrcalc::kostka::kostka_via_lr;
use lrcalc::kostka_fast::kostka_fast_u128;
use std::hint::black_box;
use std::time::{Duration, Instant};

#[derive(Clone)]
struct Case {
    label: &'static str,
    shape: Vec<i32>,
    weight: Vec<i32>,
}

fn main() {
    let repeat = std::env::args()
        .nth(1)
        .map(|arg| arg.parse::<usize>().expect("repeat must be an integer"))
        .unwrap_or(20_000);
    let cases = cases();

    for case in &cases {
        let fast = kostka_fast_u128(&case.shape, &case.weight).expect(case.label);
        let lr = kostka_via_lr(&case.shape, &case.weight).expect(case.label);
        if fast != lr {
            eprintln!("mismatch: {}", case.label);
            eprintln!("shape:  {:?}", case.shape);
            eprintln!("weight: {:?}", case.weight);
            eprintln!("fast: {fast}");
            eprintln!("lr:   {lr}");
            std::process::exit(1);
        }
    }

    println!("correctness: ok ({} cases)", cases.len());
    println!("repeat: {repeat}");

    let (fast_time, fast_sink) = time_loop(repeat, &cases, |case| {
        kostka_fast_u128(&case.shape, &case.weight).expect(case.label)
    });
    let (lr_time, lr_sink) = time_loop(repeat, &cases, |case| {
        kostka_via_lr(&case.shape, &case.weight).expect(case.label)
    });
    black_box((fast_sink, lr_sink));

    println!(
        "fast u128 DP: {}  ({} evals)",
        format_duration(fast_time),
        repeat * cases.len()
    );
    println!(
        "LR/lrcalc path: {}  ({} evals)",
        format_duration(lr_time),
        repeat * cases.len()
    );
    println!(
        "comparison: fast/LR = {:.3}x",
        fast_time.as_secs_f64() / lr_time.as_secs_f64()
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

fn format_duration(duration: Duration) -> String {
    format!("{:.3}s", duration.as_secs_f64())
}

fn cases() -> Vec<Case> {
    vec![
        case("row shape small", &[3], &[2, 1]),
        case("standard 21", &[2, 1], &[1, 1, 1]),
        case("dominance zero", &[1, 1, 1], &[2, 1]),
        case("single row all ones", &[6], &[1, 1, 1, 1, 1, 1]),
        case("single column heavy zero", &[1, 1, 1, 1], &[3, 1]),
        case("balanced 321", &[3, 2, 1], &[2, 2, 2]),
        case("example 421", &[4, 2, 1], &[3, 2, 1, 1]),
        case("standard 321", &[3, 2, 1], &[1, 1, 1, 1, 1, 1]),
        case(
            "standard 4321",
            &[4, 3, 2, 1],
            &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        ),
        case("stretched 321 x2", &[6, 4, 2], &[4, 4, 4]),
        case("stretched 321 x3", &[9, 6, 3], &[6, 6, 6]),
        case("stretched 421 x2", &[8, 4, 2], &[6, 4, 2, 2]),
        case("stretched 421 x3", &[12, 6, 3], &[9, 6, 3]),
        case("irregular A", &[5, 4, 2, 1], &[4, 3, 2, 2, 1]),
        case("irregular B", &[6, 4, 3, 1], &[5, 3, 3, 2, 1]),
        case("irregular C", &[7, 5, 3, 2], &[6, 4, 3, 2, 2]),
        case("irregular D", &[8, 6, 4, 2], &[6, 5, 4, 3, 2]),
        case("irregular E", &[9, 7, 4, 2, 1], &[7, 5, 4, 3, 2, 2]),
        case("many labels 42", &[4, 2], &[1, 1, 1, 1, 1, 1]),
        case(
            "many labels 5321",
            &[5, 3, 2, 1],
            &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        ),
        case(
            "many labels 7431",
            &[7, 4, 3, 1],
            &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        ),
    ]
}

fn case(label: &'static str, shape: &[i32], weight: &[i32]) -> Case {
    Case {
        label,
        shape: shape.to_vec(),
        weight: weight.to_vec(),
    }
}
