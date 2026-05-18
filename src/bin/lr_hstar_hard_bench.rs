use lrcalc::lr_ehrhart::{evaluate_h_vector, lr_stretch_h_vector, LrStretchPolynomial};
use num_bigint::BigInt;
use std::hint::black_box;
use std::process::Command;
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Eq, PartialEq)]
enum OutputFormat {
    Plain,
    Markdown,
}

struct Options {
    repeat: usize,
    timeout_secs: u64,
    format: OutputFormat,
    upstream_bin: String,
}

struct Family {
    label: &'static str,
    outer: &'static [i32],
    inner: &'static [i32],
    content: &'static [i32],
    stretches: &'static [u64],
}

struct FamilyReport {
    label: &'static str,
    coefficient: String,
    dimension: usize,
    h_len: usize,
    samples: usize,
    hstar_build_time: Duration,
    stretches: Vec<StretchReport>,
}

struct StretchReport {
    stretch: u64,
    value: BigInt,
    hstar_eval_time: Duration,
    upstream: UpstreamReport,
}

struct UpstreamReport {
    elapsed: Duration,
    outcome: UpstreamOutcome,
}

enum UpstreamOutcome {
    Match,
    Timeout,
    Mismatch { output: String },
    Failed { status: String, stderr: String },
}

fn main() {
    let options = parse_options();
    assert!(options.repeat > 0, "repeat must be positive");
    assert!(options.timeout_secs > 0, "timeout must be positive");

    let reports = families()
        .into_iter()
        .map(|family| run_family(&family, &options))
        .collect::<Vec<_>>();

    match options.format {
        OutputFormat::Plain => print_plain(&reports, &options),
        OutputFormat::Markdown => print_markdown(&reports, &options),
    }
}

fn parse_options() -> Options {
    let mut repeat = 1_000usize;
    let mut timeout_secs = 5u64;
    let mut format = OutputFormat::Plain;
    let mut positional = Vec::new();
    let upstream_bin = std::env::var("UPSTREAM_BIN")
        .unwrap_or_else(|_| "/tmp/lrcalc-upstream/src/lrcalc".to_string());

    for arg in std::env::args().skip(1) {
        match arg.as_str() {
            "--markdown" | "-m" => format = OutputFormat::Markdown,
            "--plain" => format = OutputFormat::Plain,
            "--help" | "-h" => {
                println!("usage: lr_hstar_hard_bench [--markdown] [repeat] [timeout-secs]");
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
        timeout_secs = arg
            .parse::<u64>()
            .unwrap_or_else(|_| panic!("expected integer timeout seconds, got {arg}"));
    }

    Options {
        repeat,
        timeout_secs,
        format,
        upstream_bin,
    }
}

fn run_family(family: &Family, options: &Options) -> FamilyReport {
    let start = Instant::now();
    let polynomial = lr_stretch_h_vector(family.outer, family.inner, family.content)
        .unwrap_or_else(|_| panic!("h* failed for {}", family.label));
    let hstar_build_time = start.elapsed();

    let stretches = family
        .stretches
        .iter()
        .copied()
        .map(|stretch| run_stretch(family, &polynomial, stretch, options))
        .collect();

    FamilyReport {
        label: family.label,
        coefficient: lrcoef_expr(family.outer, family.inner, family.content),
        dimension: polynomial.dimension,
        h_len: polynomial.h_vector.len(),
        samples: polynomial.sample_points.len(),
        hstar_build_time,
        stretches,
    }
}

fn run_stretch(
    family: &Family,
    polynomial: &LrStretchPolynomial,
    stretch: u64,
    options: &Options,
) -> StretchReport {
    let (hstar_eval_time, value) = time_hstar_eval(polynomial, stretch, options.repeat);
    let upstream = run_upstream_direct(
        &options.upstream_bin,
        options.timeout_secs,
        family.outer,
        family.inner,
        family.content,
        stretch,
        &value,
    );

    StretchReport {
        stretch,
        value,
        hstar_eval_time,
        upstream,
    }
}

fn time_hstar_eval(
    polynomial: &LrStretchPolynomial,
    stretch: u64,
    repeat: usize,
) -> (Duration, BigInt) {
    let start = Instant::now();
    let mut value = BigInt::from(0u8);
    for _ in 0..repeat {
        value = black_box(evaluate_h_vector(
            &polynomial.h_vector,
            polynomial.dimension,
            stretch,
        ));
    }
    (start.elapsed(), value)
}

fn run_upstream_direct(
    upstream_bin: &str,
    timeout_secs: u64,
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    stretch: u64,
    expected: &BigInt,
) -> UpstreamReport {
    let mut args = Vec::<String>::new();
    args.push(format!("{}s", timeout_secs));
    args.push("nice".to_string());
    args.push("-n".to_string());
    args.push("10".to_string());
    args.push(upstream_bin.to_string());
    args.push("coef".to_string());
    push_scaled_partition_args(&mut args, outer, stretch);
    args.push("-".to_string());
    push_scaled_partition_args(&mut args, inner, stretch);
    args.push("-".to_string());
    push_scaled_partition_args(&mut args, content, stretch);

    let start = Instant::now();
    let output = Command::new("timeout")
        .args(&args)
        .output()
        .expect("failed to run timeout/upstream lrcalc");
    let elapsed = start.elapsed();

    let outcome = if output.status.code() == Some(124) {
        UpstreamOutcome::Timeout
    } else if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
        if stdout == expected.to_string() {
            UpstreamOutcome::Match
        } else {
            UpstreamOutcome::Mismatch { output: stdout }
        }
    } else {
        UpstreamOutcome::Failed {
            status: output
                .status
                .code()
                .map_or_else(|| "signal".to_string(), |code| code.to_string()),
            stderr: String::from_utf8_lossy(&output.stderr).trim().to_string(),
        }
    };

    UpstreamReport { elapsed, outcome }
}

fn push_scaled_partition_args(args: &mut Vec<String>, parts: &[i32], stretch: u64) {
    for &part in parts {
        let scaled = i64::from(part)
            .checked_mul(stretch as i64)
            .expect("benchmark scale should not overflow");
        args.push(scaled.to_string());
    }
}

fn print_plain(reports: &[FamilyReport], options: &Options) {
    println!("suite: lr_hstar_hard_bench");
    println!("repeat: {}", options.repeat);
    println!("upstream_timeout_secs: {}", options.timeout_secs);
    println!("upstream_bin: {}", options.upstream_bin);
    for family in reports {
        println!(
            "{}: dim {}, h* len {}, samples {}, build {}",
            family.label,
            family.dimension,
            family.h_len,
            family.samples,
            format_duration(family.hstar_build_time)
        );
        for stretch in &family.stretches {
            println!(
                "  x{}: value {}, h* eval {}, upstream {}",
                stretch.stretch,
                stretch.value,
                format_duration(stretch.hstar_eval_time),
                upstream_plain(&stretch.upstream, options.timeout_secs)
            );
        }
    }
}

fn print_markdown(reports: &[FamilyReport], options: &Options) {
    println!("## `lr_hstar_hard_bench`: hard stretched LR examples");
    println!();
    println!("| Field | Value |");
    println!("|---|---:|");
    println!("| h* evaluation repeat | `{}` |", options.repeat);
    println!("| Upstream direct cap | `{}s` |", options.timeout_secs);
    println!("| Families | `{}` |", reports.len());
    println!(
        "| Upstream executable | `{}` |",
        markdown_cell(&options.upstream_bin)
    );
    println!();
    println!("| Family | Base coefficient | Stretch | Dimension | h* len | Samples | h* build | h* eval total | Value | Upstream direct |");
    println!("|---|---|---:|---:|---:|---:|---:|---:|---:|---|");
    for family in reports {
        for stretch in &family.stretches {
            println!(
                "| {} | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | `{}` | {} |",
                markdown_cell(family.label),
                markdown_cell(&family.coefficient),
                stretch.stretch,
                family.dimension,
                family.h_len,
                family.samples,
                format_duration(family.hstar_build_time),
                format_duration(stretch.hstar_eval_time),
                stretch.value,
                upstream_markdown(&stretch.upstream, options.timeout_secs)
            );
        }
    }
}

fn upstream_plain(report: &UpstreamReport, timeout_secs: u64) -> String {
    match &report.outcome {
        UpstreamOutcome::Match => format!("match in {}", format_duration(report.elapsed)),
        UpstreamOutcome::Timeout => format!("timeout >{timeout_secs}s"),
        UpstreamOutcome::Mismatch { output } => format!("mismatch output {output}"),
        UpstreamOutcome::Failed { status, stderr } => {
            format!("failed status {status}: {stderr}")
        }
    }
}

fn upstream_markdown(report: &UpstreamReport, timeout_secs: u64) -> String {
    match &report.outcome {
        UpstreamOutcome::Match => format!("`match in {}`", format_duration(report.elapsed)),
        UpstreamOutcome::Timeout => format!("`timeout >{}s`", timeout_secs),
        UpstreamOutcome::Mismatch { output } => {
            format!("`mismatch: {}`", markdown_cell(output))
        }
        UpstreamOutcome::Failed { status, stderr } => {
            format!(
                "`failed {}` {}",
                markdown_cell(status),
                markdown_cell(stderr)
            )
        }
    }
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

fn markdown_cell(value: &str) -> String {
    value.replace('|', "\\|").replace('\n', " ")
}

fn lrcoef_expr(outer: &[i32], inner: &[i32], content: &[i32]) -> String {
    format!(
        "c^{}_{},{}",
        format_partition(outer),
        format_partition(inner),
        format_partition(content)
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

fn families() -> Vec<Family> {
    vec![
        Family {
            label: "seven-row generic",
            outer: &[7, 6, 5, 4, 3, 2, 1],
            inner: &[4, 4, 3, 2, 1],
            content: &[5, 4, 3, 2],
            stretches: &[7, 10, 20],
        },
        Family {
            label: "six-row generic",
            outer: &[6, 5, 4, 3, 2, 1],
            inner: &[4, 3, 2, 1],
            content: &[4, 3, 2, 2],
            stretches: &[10, 20],
        },
    ]
}
