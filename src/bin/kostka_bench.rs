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

#[derive(Clone, Copy, Eq, PartialEq)]
enum OutputFormat {
    Plain,
    Markdown,
}

struct BenchResult {
    repeat: usize,
    fast_time: Duration,
    lr_time: Duration,
    diagnostics: Vec<CaseDiagnostic>,
}

struct CaseDiagnostic {
    label: &'static str,
    value: u128,
    shape_size: i32,
    labels: usize,
    fast_time: Duration,
    lr_time: Duration,
}

fn main() {
    let (repeat, format) = parse_options();
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

    let result = run_benchmark(repeat, &cases);
    match format {
        OutputFormat::Plain => print_plain_result(&result),
        OutputFormat::Markdown => print_markdown_result(&result),
    }
}

fn parse_options() -> (usize, OutputFormat) {
    let mut repeat = 20_000usize;
    let mut format = OutputFormat::Plain;
    for arg in std::env::args().skip(1) {
        match arg.as_str() {
            "--markdown" | "-m" => format = OutputFormat::Markdown,
            "--plain" => format = OutputFormat::Plain,
            "--help" | "-h" => {
                println!("usage: kostka_bench [--markdown] [repeat]");
                std::process::exit(0);
            }
            _ => repeat = arg.parse::<usize>().expect("repeat must be an integer"),
        }
    }
    (repeat, format)
}

fn run_benchmark(repeat: usize, cases: &[Case]) -> BenchResult {
    let (fast_time, fast_sink) = time_loop(repeat, cases, |case| {
        kostka_fast_u128(&case.shape, &case.weight).expect(case.label)
    });
    let (lr_time, lr_sink) = time_loop(repeat, cases, |case| {
        kostka_via_lr(&case.shape, &case.weight).expect(case.label)
    });
    black_box((fast_sink, lr_sink));

    let diagnostics = case_diagnostics(repeat, cases);
    BenchResult {
        repeat,
        fast_time,
        lr_time,
        diagnostics,
    }
}

fn print_plain_result(result: &BenchResult) {
    println!("correctness: ok ({} cases)", result.diagnostics.len());
    println!("repeat: {}", result.repeat);
    println!(
        "fast u128 DP: {}  ({} evals)",
        format_duration(result.fast_time),
        result.evals()
    );
    println!(
        "LR/lrcalc path: {}  ({} evals)",
        format_duration(result.lr_time),
        result.evals()
    );
    println!("comparison: fast/LR = {:.3}x", result.ratio());
}

fn print_markdown_result(result: &BenchResult) {
    println!("## `kostka_bench`: direct Kostka DP");
    println!();
    println!("| Field | Value |");
    println!("|---|---:|");
    println!("| Repeat | `{}` |", result.repeat);
    println!("| Cases | `{}` |", result.diagnostics.len());
    println!("| Total evaluations | `{}` |", result.evals());
    println!("| Correctness | `ok` |");
    println!();
    println!("| Metric | Direct Kostka DP | LR-translation path | Direct/LR |");
    println!("|---|---:|---:|---:|");
    println!(
        "| Total wall time | `{}` | `{}` | `{:.3}x` |",
        format_duration(result.fast_time),
        format_duration(result.lr_time),
        result.ratio()
    );
    println!();
    println!("| Per-case statistic | Value |");
    println!("|---|---:|");
    println!(
        "| Direct faster cases | `{}/{}` |",
        result.fast_faster_cases(),
        result.diagnostics.len()
    );
    println!("| Median direct/LR | `{:.3}x` |", result.median_ratio());
    println!(
        "| Geometric mean direct/LR | `{:.3}x` |",
        result.geometric_mean_ratio()
    );
    println!();
    println!("| Case | Value | Shape size | Labels | Direct DP | LR path | Direct/LR |");
    println!("|---|---:|---:|---:|---:|---:|---:|");
    for diagnostic in &result.diagnostics {
        println!(
            "| {} | `{}` | `{}` | `{}` | `{}` | `{}` | `{:.3}x` |",
            markdown_cell(diagnostic.label),
            diagnostic.value,
            diagnostic.shape_size,
            diagnostic.labels,
            format_duration(diagnostic.fast_time),
            format_duration(diagnostic.lr_time),
            diagnostic.ratio()
        );
    }
}

fn case_diagnostics(repeat: usize, cases: &[Case]) -> Vec<CaseDiagnostic> {
    cases
        .iter()
        .map(|case| {
            let value = kostka_fast_u128(&case.shape, &case.weight).expect(case.label);
            let (fast_time, fast_sink) = time_one(repeat, case, |case| {
                kostka_fast_u128(&case.shape, &case.weight).expect(case.label)
            });
            let (lr_time, lr_sink) = time_one(repeat, case, |case| {
                kostka_via_lr(&case.shape, &case.weight).expect(case.label)
            });
            black_box((fast_sink, lr_sink));
            CaseDiagnostic {
                label: case.label,
                value,
                shape_size: case.shape.iter().sum(),
                labels: case.weight.iter().take_while(|&&part| part != 0).count(),
                fast_time,
                lr_time,
            }
        })
        .collect()
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

fn time_one<F>(repeat: usize, case: &Case, mut f: F) -> (Duration, u128)
where
    F: FnMut(&Case) -> u128,
{
    let start = Instant::now();
    let mut sink = 0u128;
    for _ in 0..repeat {
        sink ^= black_box(f(case));
    }
    (start.elapsed(), sink)
}

fn format_duration(duration: Duration) -> String {
    let seconds = duration.as_secs_f64();
    if seconds >= 1.0 {
        format!("{seconds:.3}s")
    } else if seconds >= 0.001 {
        format!("{:.3}ms", seconds * 1_000.0)
    } else {
        format!("{:.3}us", seconds * 1_000_000.0)
    }
}

impl BenchResult {
    fn evals(&self) -> usize {
        self.repeat * self.diagnostics.len()
    }

    fn ratio(&self) -> f64 {
        self.fast_time.as_secs_f64() / self.lr_time.as_secs_f64()
    }

    fn fast_faster_cases(&self) -> usize {
        self.diagnostics
            .iter()
            .filter(|diagnostic| diagnostic.ratio() < 1.0)
            .count()
    }

    fn median_ratio(&self) -> f64 {
        let mut ratios = self
            .diagnostics
            .iter()
            .map(CaseDiagnostic::ratio)
            .collect::<Vec<_>>();
        ratios.sort_by(|left, right| left.total_cmp(right));
        let mid = ratios.len() / 2;
        if ratios.len() % 2 == 0 {
            (ratios[mid - 1] + ratios[mid]) / 2.0
        } else {
            ratios[mid]
        }
    }

    fn geometric_mean_ratio(&self) -> f64 {
        let log_sum = self
            .diagnostics
            .iter()
            .map(|diagnostic| diagnostic.ratio().ln())
            .sum::<f64>();
        (log_sum / self.diagnostics.len() as f64).exp()
    }
}

impl CaseDiagnostic {
    fn ratio(&self) -> f64 {
        self.fast_time.as_secs_f64() / self.lr_time.as_secs_f64()
    }
}

fn markdown_cell(value: &str) -> String {
    value.replace('|', "\\|").replace('\n', " ")
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
