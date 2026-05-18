use lrcalc::kostka::kostka_lr_triple;
use lrcalc::kostka_fast::{kostka_counts_stats, KostkaCountsStats};
use lrcalc::lr_gt::{
    lrcoef_gt_counts_stats, lrcoef_hybrid_counts_stats, LrGtCountsStats, LrHybridCountsMode,
    LrHybridCountsStats,
};
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

struct CaseStats {
    kostka: KostkaCountsStats,
    lr: LrGtCountsStats,
    hybrid: LrHybridCountsStats,
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum OutputFormat {
    Plain,
    Markdown,
}

struct BenchResult {
    repeat: usize,
    max_scale: i32,
    kostka_time: Duration,
    lr_time: Duration,
    hybrid_time: Duration,
    cases: Vec<CaseReport>,
}

struct CaseReport {
    label: String,
    coefficient: String,
    full: u128,
    interior: u128,
    k_full_peak: usize,
    k_interior_peak: usize,
    lr_full_peak: usize,
    lr_interior_peak: usize,
    hybrid_mode: &'static str,
    hybrid_full_peak: usize,
    hybrid_interior_peak: usize,
}

fn main() {
    let (repeat, max_scale, format) = parse_options();
    assert!(repeat > 0, "repeat must be positive");
    assert!(max_scale > 0, "max scale must be positive");

    let cases = scaled_cases(max_scale);
    verify_cases(&cases);
    let result = run_benchmark(repeat, max_scale, &cases);
    match format {
        OutputFormat::Plain => print_plain_result(&result),
        OutputFormat::Markdown => print_markdown_result(&result),
    }
}

fn parse_options() -> (usize, i32, OutputFormat) {
    let mut repeat = 5usize;
    let mut max_scale = 3i32;
    let mut format = OutputFormat::Plain;
    let mut positional = Vec::new();
    for arg in std::env::args().skip(1) {
        match arg.as_str() {
            "--markdown" | "-m" => format = OutputFormat::Markdown,
            "--plain" => format = OutputFormat::Plain,
            "--help" | "-h" => {
                println!("usage: stretched_dp_bench [--markdown] [repeat] [max-scale]");
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
        max_scale = arg
            .parse::<i32>()
            .unwrap_or_else(|_| panic!("expected integer max scale, got {arg}"));
    }
    (repeat, max_scale, format)
}

fn run_benchmark(repeat: usize, max_scale: i32, cases: &[ScaledCase]) -> BenchResult {
    let reports = case_reports(cases);

    let (kostka_time, kostka_sink) = time_loop(repeat, cases, |case| {
        let stats = kostka_counts_stats(&case.shape, &case.weight)
            .unwrap_or_else(|_| panic!("Kostka counts failed for {}", case.label));
        stats.full ^ stats.interior
    });
    let (lr_time, lr_sink) = time_loop(repeat, cases, |case| {
        let stats = lrcoef_gt_counts_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("LR GT counts failed for {}", case.label));
        stats.full ^ stats.interior
    });
    let (hybrid_time, hybrid_sink) = time_loop(repeat, cases, |case| {
        let stats = lrcoef_hybrid_counts_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("hybrid LR counts failed for {}", case.label));
        stats.counts.full ^ stats.counts.interior
    });
    black_box((kostka_sink, lr_sink, hybrid_sink));

    BenchResult {
        repeat,
        max_scale,
        kostka_time,
        lr_time,
        hybrid_time,
        cases: reports,
    }
}

fn print_plain_result(result: &BenchResult) {
    println!("suite: stretched_dp_bench");
    println!("repeat: {}", result.repeat);
    println!("max_scale: {}", result.max_scale);
    println!("correctness: ok ({} scaled cases)", result.cases.len());
    print_case_report_table(&result.cases);

    println!(
        "Kostka DP counts: {}  ({} evals)",
        format_duration(result.kostka_time),
        result.evals()
    );
    println!(
        "LR GT DP counts:  {}  ({} evals)",
        format_duration(result.lr_time),
        result.evals()
    );
    println!(
        "Hybrid LR counts: {}  ({} evals)",
        format_duration(result.hybrid_time),
        result.evals()
    );
    println!("comparison: LR/Kostka = {:.3}x", result.lr_vs_kostka());
    println!(
        "comparison: hybrid/Kostka = {:.3}x",
        result.hybrid_vs_kostka()
    );
    println!("comparison: LR/hybrid = {:.3}x", result.lr_vs_hybrid());
}

fn print_markdown_result(result: &BenchResult) {
    println!("## `stretched_dp_bench`: stretched Kostka-as-LR counts");
    println!();
    println!("| Field | Value |");
    println!("|---|---:|");
    println!("| Repeat | `{}` |", result.repeat);
    println!("| Max scale | `{}` |", result.max_scale);
    println!("| Scaled cases | `{}` |", result.cases.len());
    println!("| Total evaluations per method | `{}` |", result.evals());
    println!("| Correctness | `ok` |");
    println!();
    println!("| Metric | Kostka DP | LR GT-chain DP | Hybrid LR |");
    println!("|---|---:|---:|---:|");
    println!(
        "| Total wall time | `{}` | `{}` | `{}` |",
        format_duration(result.kostka_time),
        format_duration(result.lr_time),
        format_duration(result.hybrid_time)
    );
    println!();
    println!("| Ratio | Value |");
    println!("|---|---:|");
    println!(
        "| LR GT-chain / Kostka DP | `{:.3}x` |",
        result.lr_vs_kostka()
    );
    println!(
        "| Hybrid LR / Kostka DP | `{:.3}x` |",
        result.hybrid_vs_kostka()
    );
    println!(
        "| LR GT-chain / Hybrid LR | `{:.3}x` |",
        result.lr_vs_hybrid()
    );
    println!();
    println!("| Case | Coefficient | Full | Interior | Kostka full peak | Kostka interior peak | LR full peak | LR interior peak | Hybrid mode | Hybrid full peak | Hybrid interior peak |");
    println!("|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|");
    for case in &result.cases {
        println!(
            "| {} | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` |",
            markdown_cell(&case.label),
            markdown_cell(&case.coefficient),
            case.full,
            case.interior,
            case.k_full_peak,
            case.k_interior_peak,
            case.lr_full_peak,
            case.lr_interior_peak,
            case.hybrid_mode,
            case.hybrid_full_peak,
            case.hybrid_interior_peak
        );
    }
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
        let stats = case_stats(case);
        assert_eq!(
            stats.kostka.full, stats.lr.full,
            "full mismatch for {}",
            case.label
        );
        assert_eq!(
            stats.kostka.interior, stats.lr.interior,
            "interior mismatch for {}",
            case.label
        );
        assert_eq!(
            stats.kostka.full, stats.hybrid.counts.full,
            "hybrid full mismatch for {}",
            case.label
        );
        assert_eq!(
            stats.kostka.interior, stats.hybrid.counts.interior,
            "hybrid interior mismatch for {}",
            case.label
        );
        assert_eq!(
            stats.hybrid.mode,
            LrHybridCountsMode::KostkaTranslation,
            "hybrid missed Kostka translation for {}",
            case.label
        );
    }
}

fn print_case_report_table(cases: &[CaseReport]) {
    println!(
        "case\tfull\tinterior\tk_full_peak\tk_int_peak\tlr_full_peak\tlr_int_peak\thybrid_mode\thybrid_full_peak\thybrid_int_peak"
    );
    for case in cases {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            case.label,
            case.full,
            case.interior,
            case.k_full_peak,
            case.k_interior_peak,
            case.lr_full_peak,
            case.lr_interior_peak,
            case.hybrid_mode,
            case.hybrid_full_peak,
            case.hybrid_interior_peak
        );
    }
}

fn case_reports(cases: &[ScaledCase]) -> Vec<CaseReport> {
    cases
        .iter()
        .map(|case| {
            let stats = case_stats(case);
            CaseReport {
                label: case.label.clone(),
                coefficient: stretched_kostka_lr_expr(case),
                full: stats.kostka.full,
                interior: stats.kostka.interior,
                k_full_peak: stats.kostka.full_peak_states,
                k_interior_peak: stats.kostka.interior_peak_states,
                lr_full_peak: stats.lr.full_peak_states,
                lr_interior_peak: stats.lr.interior_peak_states,
                hybrid_mode: stats.hybrid.mode.label(),
                hybrid_full_peak: stats.hybrid.counts.full_peak_states,
                hybrid_interior_peak: stats.hybrid.counts.interior_peak_states,
            }
        })
        .collect()
}

fn case_stats(case: &ScaledCase) -> CaseStats {
    CaseStats {
        kostka: kostka_counts_stats(&case.shape, &case.weight)
            .unwrap_or_else(|_| panic!("Kostka counts failed for {}", case.label)),
        lr: lrcoef_gt_counts_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("LR GT counts failed for {}", case.label)),
        hybrid: lrcoef_hybrid_counts_stats(&case.outer, &case.inner, &case.content)
            .unwrap_or_else(|_| panic!("hybrid LR counts failed for {}", case.label)),
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

fn scale_parts(parts: &[i32], scale: i32) -> Vec<i32> {
    parts.iter().map(|part| part * scale).collect()
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
        self.repeat * self.cases.len()
    }

    fn lr_vs_kostka(&self) -> f64 {
        self.lr_time.as_secs_f64() / self.kostka_time.as_secs_f64()
    }

    fn hybrid_vs_kostka(&self) -> f64 {
        self.hybrid_time.as_secs_f64() / self.kostka_time.as_secs_f64()
    }

    fn lr_vs_hybrid(&self) -> f64 {
        self.lr_time.as_secs_f64() / self.hybrid_time.as_secs_f64()
    }
}

fn markdown_cell(value: &str) -> String {
    value.replace('|', "\\|").replace('\n', " ")
}

fn stretched_kostka_lr_expr(case: &ScaledCase) -> String {
    format!(
        "K_{},{} = c^{}_{},{}",
        format_partition(&case.shape),
        format_partition(&case.weight),
        format_partition(&case.outer),
        format_partition(&case.inner),
        format_partition(&case.content)
    )
}

fn format_partition(parts: &[i32]) -> String {
    let trimmed = parts
        .iter()
        .copied()
        .take_while(|&part| part != 0)
        .collect::<Vec<_>>();
    if trimmed.is_empty() {
        "()".to_string()
    } else {
        let body = trimmed
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(",");
        format!("({body})")
    }
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
