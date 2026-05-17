use lrcalc::lr_ehrhart::{evaluate_h_vector, lr_stretch_h_vector, LrStretchPolynomial};
use lrcalc::lr_gt::{lrcoef_gt_interior_dfs_u128, lrcoef_gt_interior_stats};
use lrcalc::lrcoef::{
    lrcoef, lrcoef_buch_interior_memo_u128, lrcoef_buch_interior_stats, lrcoef_buch_interior_u128,
};
use num_bigint::BigInt;
use num_traits::ToPrimitive;
use std::hint::black_box;
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Eq, PartialEq)]
enum OutputFormat {
    Plain,
    Markdown,
}

#[derive(Clone)]
struct LrCase {
    label: &'static str,
    outer: Vec<i32>,
    inner: Vec<i32>,
    content: Vec<i32>,
}

struct PreparedCase {
    case: LrCase,
    scaled_outer: Vec<i32>,
    scaled_inner: Vec<i32>,
    scaled_content: Vec<i32>,
    polynomial: LrStretchPolynomial,
    direct_scaled_value: u128,
    evaluated_scaled_value: BigInt,
}

struct CaseReport {
    label: &'static str,
    dimension: usize,
    h_len: usize,
    samples: usize,
    base_value: u128,
    scaled_value: u128,
    buch_weak: u128,
    buch_strict: u128,
    gt_peak: usize,
}

struct BenchResult {
    repeat: usize,
    eval_stretch: u64,
    cases: Vec<CaseReport>,
    full_time: Duration,
    buch_interior_time: Duration,
    gt_interior_time: Duration,
    memo_interior_time: Duration,
    hstar_build_time: Duration,
    hstar_eval_time: Duration,
    direct_scaled_time: Duration,
}

fn main() {
    let (repeat, eval_stretch, format) = parse_options();
    assert!(repeat > 0, "repeat must be positive");
    assert!(eval_stretch > 0, "evaluation stretch must be positive");

    let cases = cases();
    let prepared = prepare_cases(&cases, eval_stretch);
    verify_cases(&prepared);
    let result = run_benchmark(repeat, eval_stretch, &cases, &prepared);

    match format {
        OutputFormat::Plain => print_plain_result(&result),
        OutputFormat::Markdown => print_markdown_result(&result),
    }
}

fn parse_options() -> (usize, u64, OutputFormat) {
    let mut repeat = 10usize;
    let mut eval_stretch = 5u64;
    let mut format = OutputFormat::Plain;
    let mut positional = Vec::new();

    for arg in std::env::args().skip(1) {
        match arg.as_str() {
            "--markdown" | "-m" => format = OutputFormat::Markdown,
            "--plain" => format = OutputFormat::Plain,
            "--help" | "-h" => {
                println!("usage: lr_hstar_bench [--markdown] [repeat] [eval-stretch]");
                std::process::exit(0);
            }
            _ => positional.push(arg),
        }
    }
    if let Some(arg) = positional.first() {
        repeat = arg
            .parse::<usize>()
            .unwrap_or_else(|_| panic!("expected integer repeat count, got {arg}"));
    }
    if let Some(arg) = positional.get(1) {
        eval_stretch = arg
            .parse::<u64>()
            .unwrap_or_else(|_| panic!("expected integer evaluation stretch, got {arg}"));
    }

    (repeat, eval_stretch, format)
}

fn prepare_cases(cases: &[LrCase], eval_stretch: u64) -> Vec<PreparedCase> {
    cases
        .iter()
        .map(|case| {
            let polynomial = lr_stretch_h_vector(&case.outer, &case.inner, &case.content)
                .unwrap_or_else(|_| panic!("h* failed for {}", case.label));
            let scaled_outer = scale_partition(&case.outer, eval_stretch);
            let scaled_inner = scale_partition(&case.inner, eval_stretch);
            let scaled_content = scale_partition(&case.content, eval_stretch);
            let direct_scaled_value = lrcoef(&scaled_outer, &scaled_inner, &scaled_content)
                .unwrap_or_else(|_| panic!("scaled direct LR failed for {}", case.label));
            let evaluated_scaled_value =
                evaluate_h_vector(&polynomial.h_vector, polynomial.dimension, eval_stretch);
            PreparedCase {
                case: case.clone(),
                scaled_outer,
                scaled_inner,
                scaled_content,
                polynomial,
                direct_scaled_value,
                evaluated_scaled_value,
            }
        })
        .collect()
}

fn verify_cases(prepared: &[PreparedCase]) {
    for prepared_case in prepared {
        let case = &prepared_case.case;
        let buch = lrcoef_buch_interior_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch interior failed for {}", case.label));
        let gt = lrcoef_gt_interior_dfs_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT interior failed for {}", case.label));
        assert_eq!(buch, gt, "strict count mismatch for {}", case.label);
        let memo = lrcoef_buch_interior_memo_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch memo interior failed for {}", case.label));
        assert_eq!(buch, memo, "memo strict count mismatch for {}", case.label);
        assert_eq!(
            BigInt::from(prepared_case.direct_scaled_value),
            prepared_case.evaluated_scaled_value,
            "h* evaluation mismatch for {}",
            case.label
        );
    }
}

fn run_benchmark(
    repeat: usize,
    eval_stretch: u64,
    cases: &[LrCase],
    prepared: &[PreparedCase],
) -> BenchResult {
    let case_reports = case_reports(prepared);

    let (full_time, full_sink) = time_case_loop(repeat, cases, |case| {
        lrcoef(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("full LR failed for {}", case.label))
    });
    let (buch_interior_time, buch_sink) = time_case_loop(repeat, cases, |case| {
        lrcoef_buch_interior_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch interior failed for {}", case.label))
    });
    let (gt_interior_time, gt_sink) = time_case_loop(repeat, cases, |case| {
        lrcoef_gt_interior_dfs_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("GT interior failed for {}", case.label))
    });
    let (memo_interior_time, memo_sink) = time_case_loop(repeat, cases, |case| {
        lrcoef_buch_interior_memo_u128(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("Buch memo interior failed for {}", case.label))
    });
    let (hstar_build_time, hstar_sink) = time_case_loop(repeat, cases, |case| {
        let polynomial = lr_stretch_h_vector(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("h* failed for {}", case.label));
        (polynomial.dimension ^ polynomial.h_vector.len()) as u128
    });
    let (hstar_eval_time, eval_sink) = time_prepared_loop(repeat, prepared, |case| {
        let value = evaluate_h_vector(
            &case.polynomial.h_vector,
            case.polynomial.dimension,
            eval_stretch,
        );
        bigint_sink(&value)
    });
    let (direct_scaled_time, direct_scaled_sink) = time_prepared_loop(repeat, prepared, |case| {
        lrcoef(&case.scaled_outer, &case.scaled_inner, &case.scaled_content)
            .unwrap_or_else(|_| panic!("scaled direct LR failed for {}", case.case.label))
    });
    black_box((
        full_sink,
        buch_sink,
        gt_sink,
        memo_sink,
        hstar_sink,
        eval_sink,
        direct_scaled_sink,
    ));

    BenchResult {
        repeat,
        eval_stretch,
        cases: case_reports,
        full_time,
        buch_interior_time,
        gt_interior_time,
        memo_interior_time,
        hstar_build_time,
        hstar_eval_time,
        direct_scaled_time,
    }
}

fn print_plain_result(result: &BenchResult) {
    println!("suite: lr_hstar_bench");
    println!("repeat: {}", result.repeat);
    println!("eval_stretch: {}", result.eval_stretch);
    println!("correctness: ok ({} cases)", result.cases.len());
    for case in &result.cases {
        println!(
            "{}: dim {}, h* len {}, base {}, scale {} value {}, Buch weak {}, Buch strict {}, GT peak {}",
            case.label,
            case.dimension,
            case.h_len,
            case.base_value,
            result.eval_stretch,
            case.scaled_value,
            case.buch_weak,
            case.buch_strict,
            case.gt_peak
        );
    }
    println!(
        "full Buch LR:       {}  ({} evals)",
        format_duration(result.full_time),
        result.evals()
    );
    println!(
        "Buch strict LR:     {}  ({} evals)",
        format_duration(result.buch_interior_time),
        result.evals()
    );
    println!(
        "GT strict LR:       {}  ({} evals)",
        format_duration(result.gt_interior_time),
        result.evals()
    );
    println!(
        "Buch memo strict:   {}  ({} evals)",
        format_duration(result.memo_interior_time),
        result.evals()
    );
    println!(
        "LR h* build:        {}  ({} evals)",
        format_duration(result.hstar_build_time),
        result.evals()
    );
    println!(
        "LR h* eval-only:    {}  ({} evals)",
        format_duration(result.hstar_eval_time),
        result.evals()
    );
    println!(
        "direct scaled LR:   {}  ({} evals)",
        format_duration(result.direct_scaled_time),
        result.evals()
    );
    println!(
        "strict comparison: GT/Buch = {:.3}x",
        result.gt_vs_buch_strict()
    );
    println!(
        "stretched comparison: direct/eval = {:.3}x",
        result.direct_scaled_vs_hstar_eval()
    );
}

fn print_markdown_result(result: &BenchResult) {
    println!("## `lr_hstar_bench`: explicit stretched LR h*-evaluation");
    println!();
    println!("| Field | Value |");
    println!("|---|---:|");
    println!("| Repeat | `{}` |", result.repeat);
    println!("| Evaluation stretch | `{}` |", result.eval_stretch);
    println!("| Cases | `{}` |", result.cases.len());
    println!("| Total evaluations per metric | `{}` |", result.evals());
    println!("| Correctness | `ok` |");
    println!();
    println!("| Metric | Wall time |");
    println!("|---|---:|");
    println!(
        "| Base full Buch LR | `{}` |",
        format_duration(result.full_time)
    );
    println!(
        "| Base Buch strict LR | `{}` |",
        format_duration(result.buch_interior_time)
    );
    println!(
        "| Base GT strict LR | `{}` |",
        format_duration(result.gt_interior_time)
    );
    println!(
        "| Base Buch memo strict LR | `{}` |",
        format_duration(result.memo_interior_time)
    );
    println!(
        "| Build LR h* polynomial | `{}` |",
        format_duration(result.hstar_build_time)
    );
    println!(
        "| Evaluate cached h* at stretch | `{}` |",
        format_duration(result.hstar_eval_time)
    );
    println!(
        "| Direct LR at stretch | `{}` |",
        format_duration(result.direct_scaled_time)
    );
    println!();
    println!("| Ratio | Value |");
    println!("|---|---:|");
    println!(
        "| Base GT strict / Buch strict | `{:.3}x` |",
        result.gt_vs_buch_strict()
    );
    println!(
        "| Direct stretched LR / cached h* eval | `{:.3}x` |",
        result.direct_scaled_vs_hstar_eval()
    );
    println!(
        "| h* build / direct stretched LR | `{:.3}x` |",
        result.hstar_build_vs_direct_scaled()
    );
    println!();
    println!("| Case | Dimension | h* len | Samples | Base value | Value at stretch | Buch weak | Buch strict | GT peak |");
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|");
    for case in &result.cases {
        println!(
            "| {} | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` |",
            markdown_cell(case.label),
            case.dimension,
            case.h_len,
            case.samples,
            case.base_value,
            case.scaled_value,
            case.buch_weak,
            case.buch_strict,
            case.gt_peak
        );
    }
}

fn case_reports(prepared: &[PreparedCase]) -> Vec<CaseReport> {
    prepared
        .iter()
        .map(|prepared_case| {
            let case = &prepared_case.case;
            let buch_stats = lrcoef_buch_interior_stats(&case.outer, &case.inner, &case.content)
                .unwrap_or_else(|_| panic!("Buch interior stats failed for {}", case.label));
            let gt_stats = lrcoef_gt_interior_stats(&case.outer, &case.inner, &case.content)
                .unwrap_or_else(|_| panic!("GT interior stats failed for {}", case.label));
            CaseReport {
                label: case.label,
                dimension: prepared_case.polynomial.dimension,
                h_len: prepared_case.polynomial.h_vector.len(),
                samples: prepared_case.polynomial.sample_points.len(),
                base_value: lrcoef(&case.outer, &case.inner, &case.content)
                    .unwrap_or_else(|_| panic!("full LR failed for {}", case.label)),
                scaled_value: prepared_case.direct_scaled_value,
                buch_weak: buch_stats.weak_tableaux,
                buch_strict: buch_stats.strict_tableaux,
                gt_peak: gt_stats.peak_states,
            }
        })
        .collect()
}

fn time_case_loop<F>(repeat: usize, cases: &[LrCase], mut f: F) -> (Duration, u128)
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

fn time_prepared_loop<F>(repeat: usize, cases: &[PreparedCase], mut f: F) -> (Duration, u128)
where
    F: FnMut(&PreparedCase) -> u128,
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
    let seconds = duration.as_secs_f64();
    if seconds >= 1.0 {
        format!("{seconds:.3}s")
    } else if seconds >= 0.001 {
        format!("{:.3}ms", seconds * 1_000.0)
    } else {
        format!("{:.3}us", seconds * 1_000_000.0)
    }
}

fn scale_partition(parts: &[i32], stretch: u64) -> Vec<i32> {
    parts
        .iter()
        .map(|&part| {
            i32::try_from(i64::from(part) * stretch as i64)
                .expect("benchmark stretch should fit i32")
        })
        .collect()
}

fn bigint_sink(value: &BigInt) -> u128 {
    value
        .to_u128()
        .unwrap_or_else(|| value.to_string().len() as u128)
}

fn markdown_cell(value: &str) -> String {
    value.replace('|', "\\|").replace('\n', " ")
}

impl BenchResult {
    fn evals(&self) -> usize {
        self.repeat * self.cases.len()
    }

    fn gt_vs_buch_strict(&self) -> f64 {
        self.gt_interior_time.as_secs_f64() / self.buch_interior_time.as_secs_f64()
    }

    fn direct_scaled_vs_hstar_eval(&self) -> f64 {
        self.direct_scaled_time.as_secs_f64() / self.hstar_eval_time.as_secs_f64()
    }

    fn hstar_build_vs_direct_scaled(&self) -> f64 {
        self.hstar_build_time.as_secs_f64() / self.direct_scaled_time.as_secs_f64()
    }
}

fn cases() -> Vec<LrCase> {
    vec![
        case("s21 square", &[3, 2, 1], &[2, 1], &[2, 1]),
        case(
            "medium generic",
            &[5, 4, 3, 2, 1],
            &[3, 2, 1],
            &[4, 3, 1, 1],
        ),
        case(
            "larger generic",
            &[7, 6, 5, 4, 3, 2, 1],
            &[4, 4, 3, 2, 1],
            &[5, 4, 3, 2],
        ),
        case(
            "kostka medium translated",
            &[12, 8, 5, 3, 1],
            &[8, 5, 3, 1],
            &[5, 4, 2, 1],
        ),
        case(
            "kostka large translated",
            &[20, 14, 9, 5, 2],
            &[14, 9, 5, 2],
            &[8, 6, 4, 2],
        ),
    ]
}

fn case(label: &'static str, outer: &[i32], inner: &[i32], content: &[i32]) -> LrCase {
    LrCase {
        label,
        outer: outer.to_vec(),
        inner: inner.to_vec(),
        content: content.to_vec(),
    }
}
