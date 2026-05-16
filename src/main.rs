use lrcalc::kostka::{kostka_lr_triple, kostka_via_lr};
use lrcalc::kostka_fast::{
    kostka_fast_stats, kostka_fast_u128, kostka_interior_stats, kostka_interior_u128,
    KostkaFastError,
};
use lrcalc::lr_ehrhart::{
    format_bigint_vector, format_error as format_lr_stretch_error, lr_stretch_coefficient,
    lr_stretch_h_vector,
};
use lrcalc::lr_gt::{
    lrcoef_gt_hybrid_stats, lrcoef_gt_interior_dfs_stats, lrcoef_gt_interior_dfs_u128,
    lrcoef_gt_interior_stats, lrcoef_gt_interior_u128, lrcoef_gt_stats, lrcoef_gt_u128, LrGtError,
};
use lrcalc::lr_signed::{lrcoef_signed_kostka, lrcoef_signed_kostka_stats, SignedLrError};
use lrcalc::lrcoef::{
    lrcoef, lrcoef_buch_counts_u128, lrcoef_buch_dimension, lrcoef_buch_interior_memo_stats,
    lrcoef_buch_interior_memo_u128, lrcoef_buch_interior_stats, lrcoef_buch_interior_u128,
    LrCoefError,
};

fn main() {
    let mut args = std::env::args();
    let program = args.next().unwrap_or_else(|| "lrcalc".to_string());
    match args.next().as_deref() {
        Some("--version") | Some("-V") => {
            println!("lrcalc-new 0.1.0");
        }
        Some("coef" | "lrcoef") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef(&parts[0], &parts[1], &parts[2]).map_err(format_lrcoef_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("kostka-lr") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_pair(&rest)
                .and_then(|parts| kostka_via_lr(&parts[0], &parts[1]).map_err(format_lrcoef_error))
            {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("kostka-lr-triple") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_pair(&rest).and_then(|parts| {
                kostka_lr_triple(&parts[0], &parts[1]).map_err(format_lrcoef_error)
            }) {
                Ok((outer, inner, content)) => {
                    println!(
                        "{} - {} - {}",
                        format_partition(&outer),
                        format_partition(&inner),
                        format_partition(&content)
                    );
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("kostka-fast") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_pair(&rest).and_then(|parts| {
                kostka_fast_u128(&parts[0], &parts[1]).map_err(format_kostka_fast_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("kostka-fast-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_pair(&rest).and_then(|parts| {
                kostka_fast_stats(&parts[0], &parts[1]).map_err(format_kostka_fast_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("peak_states: {}", stats.peak_states);
                    println!("cached_transitions: {}", stats.cached_transitions);
                    println!("levels: {:?}", stats.levels);
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("kostka-fast-interior") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_pair(&rest).and_then(|parts| {
                kostka_interior_u128(&parts[0], &parts[1]).map_err(format_kostka_fast_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("kostka-fast-interior-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_pair(&rest).and_then(|parts| {
                kostka_interior_stats(&parts[0], &parts[1]).map_err(format_kostka_fast_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("peak_states: {}", stats.peak_states);
                    println!("levels: {:?}", stats.levels);
                    println!("full_reachable_levels: {:?}", stats.full_reachable_levels);
                    println!(
                        "strict_lower_constraints: {}",
                        stats.strict_lower_constraints
                    );
                    println!(
                        "strict_diagonal_constraints: {}",
                        stats.strict_diagonal_constraints
                    );
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-signed-kostka") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_signed_kostka(&parts[0], &parts[1], &parts[2])
                    .map_err(format_signed_lr_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-signed-kostka-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_signed_kostka_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_signed_lr_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("mode: {:?}", stats.mode);
                    println!("determinant_len: {}", stats.determinant_len);
                    println!("valid_permutations: {}", stats.valid_permutations);
                    println!("aggregated_terms: {}", stats.aggregated_terms);
                    println!("cached_transitions: {}", stats.cached_transitions);
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-gt") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_gt_u128(&parts[0], &parts[1], &parts[2]).map_err(format_lr_gt_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-gt-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_gt_stats(&parts[0], &parts[1], &parts[2]).map_err(format_lr_gt_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("peak_states: {}", stats.peak_states);
                    println!("levels: {:?}", stats.levels);
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-gt-hybrid") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_gt_hybrid_stats(&parts[0], &parts[1], &parts[2]).map_err(format_lr_gt_error)
            }) {
                Ok(stats) => println!("{}", stats.stats.value),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-gt-hybrid-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_gt_hybrid_stats(&parts[0], &parts[1], &parts[2]).map_err(format_lr_gt_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.stats.value);
                    println!("mode: {}", stats.mode.label());
                    println!("enforced_rows: {:?}", stats.enforced_rows);
                    println!("peak_states: {}", stats.stats.peak_states);
                    println!("levels: {:?}", stats.stats.levels);
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-gt-interior") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_gt_interior_u128(&parts[0], &parts[1], &parts[2]).map_err(format_lr_gt_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-gt-interior-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_gt_interior_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_gt_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("peak_states: {}", stats.peak_states);
                    println!("levels: {:?}", stats.levels);
                    println!("full_reachable_levels: {:?}", stats.full_reachable_levels);
                    println!(
                        "strict_lower_constraints: {}",
                        stats.strict_lower_constraints
                    );
                    println!(
                        "strict_diagonal_constraints: {}",
                        stats.strict_diagonal_constraints
                    );
                    println!(
                        "strict_yamanouchi_constraints: {}",
                        stats.strict_yamanouchi_constraints
                    );
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-gt-interior-dfs") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_gt_interior_dfs_u128(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_gt_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-gt-interior-dfs-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_gt_interior_dfs_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_gt_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("weak_nodes: {}", stats.weak_nodes);
                    println!("strict_nodes: {}", stats.strict_nodes);
                    println!(
                        "strict_lower_constraints: {}",
                        stats.strict_lower_constraints
                    );
                    println!(
                        "strict_diagonal_constraints: {}",
                        stats.strict_diagonal_constraints
                    );
                    println!(
                        "strict_yamanouchi_constraints: {}",
                        stats.strict_yamanouchi_constraints
                    );
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-buch-interior") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_buch_interior_u128(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lrcoef_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-buch-counts") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_buch_counts_u128(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lrcoef_error)
            }) {
                Ok(counts) => {
                    println!("full: {}", counts.full);
                    println!("interior: {}", counts.interior);
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-buch-interior-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_buch_interior_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lrcoef_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("weak_tableaux: {}", stats.weak_tableaux);
                    println!("strict_tableaux: {}", stats.strict_tableaux);
                    println!(
                        "strict_lower_constraints: {}",
                        stats.strict_lower_constraints
                    );
                    println!(
                        "strict_diagonal_constraints: {}",
                        stats.strict_diagonal_constraints
                    );
                    println!(
                        "strict_yamanouchi_constraints: {}",
                        stats.strict_yamanouchi_constraints
                    );
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-buch-interior-memo") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_buch_interior_memo_u128(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lrcoef_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-buch-interior-memo-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_buch_interior_memo_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lrcoef_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("weak_tableaux: {}", stats.weak_tableaux);
                    println!("memo_states: {}", stats.memo_states);
                    println!("cache_hits: {}", stats.cache_hits);
                    println!("generated_transitions: {}", stats.generated_transitions);
                    println!(
                        "strict_lower_constraints: {}",
                        stats.strict_lower_constraints
                    );
                    println!(
                        "strict_diagonal_constraints: {}",
                        stats.strict_diagonal_constraints
                    );
                    println!(
                        "strict_yamanouchi_constraints: {}",
                        stats.strict_yamanouchi_constraints
                    );
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-buch-dimension") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_buch_dimension(&parts[0], &parts[1], &parts[2]).map_err(format_lrcoef_error)
            }) {
                Ok(Some(dimension)) => println!("{dimension}"),
                Ok(None) => println!("empty"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-stretch-hvector") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lr_stretch_h_vector(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_stretch_error)
            }) {
                Ok(polynomial) => {
                    println!("dimension: {}", polynomial.dimension);
                    println!("h_vector: {}", format_bigint_vector(&polynomial.h_vector));
                    println!("sample_points: {:?}", polynomial.sample_points);
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-stretch-eval") => {
            let rest: Vec<String> = args.collect();
            match parse_stretch_and_triple(&rest).and_then(|(stretch, parts)| {
                lr_stretch_coefficient(&parts[0], &parts[1], &parts[2], stretch)
                    .map_err(format_lr_stretch_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some(command) => {
            eprintln!("{program}: command '{command}' is not implemented yet");
            std::process::exit(2);
        }
        None => {
            eprintln!("Usage: {program} <coef|mult|skew|coprod|tab> [arguments]");
            std::process::exit(2);
        }
    }
}

fn parse_partition_triple(args: &[String]) -> Result<[Vec<i32>; 3], String> {
    let mut parts = [Vec::new(), Vec::new(), Vec::new()];
    let mut section = 0usize;

    for token in args {
        if token == "-" {
            section += 1;
            if section >= parts.len() {
                return Err("expected exactly two '-' separators".to_string());
            }
            continue;
        }

        for piece in token.split(',') {
            let piece = piece.trim_matches(|ch| matches!(ch, '(' | ')' | '[' | ']'));
            if piece.is_empty() {
                continue;
            }
            let value = piece
                .parse::<i32>()
                .map_err(|_| format!("invalid integer '{piece}'"))?;
            parts[section].push(value);
        }
    }

    if section != 2 {
        return Err("usage: coef OUTER - INNER1 - INNER2".to_string());
    }

    Ok(parts)
}

fn parse_stretch_and_triple(args: &[String]) -> Result<(u64, [Vec<i32>; 3]), String> {
    let Some((stretch, rest)) = args.split_first() else {
        return Err("usage: lr-stretch-eval STRETCH OUTER - INNER - CONTENT".to_string());
    };
    let stretch = stretch
        .parse::<u64>()
        .map_err(|_| format!("invalid stretch factor '{stretch}'"))?;
    Ok((stretch, parse_partition_triple(rest)?))
}

fn parse_partition_pair(args: &[String]) -> Result<[Vec<i32>; 2], String> {
    let mut parts = [Vec::new(), Vec::new()];
    let mut section = 0usize;

    for token in args {
        if token == "-" {
            section += 1;
            if section >= parts.len() {
                return Err("expected exactly one '-' separator".to_string());
            }
            continue;
        }

        for piece in token.split(',') {
            let piece = piece.trim_matches(|ch| matches!(ch, '(' | ')' | '[' | ']'));
            if piece.is_empty() {
                continue;
            }
            let value = piece
                .parse::<i32>()
                .map_err(|_| format!("invalid integer '{piece}'"))?;
            parts[section].push(value);
        }
    }

    if section != 1 {
        return Err("usage: kostka-lr SHAPE - WEIGHT".to_string());
    }

    Ok(parts)
}

fn format_partition(partition: &[i32]) -> String {
    if partition.is_empty() {
        return "0".to_string();
    }
    partition
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(" ")
}

fn format_lrcoef_error(error: LrCoefError) -> String {
    match error {
        LrCoefError::InvalidPartition => "invalid partition".to_string(),
        LrCoefError::ArithmeticOverflow => "arithmetic overflow".to_string(),
    }
}

fn format_kostka_fast_error(error: KostkaFastError) -> String {
    match error {
        KostkaFastError::InvalidInput => "invalid input".to_string(),
        KostkaFastError::ArithmeticOverflow => "arithmetic overflow".to_string(),
        KostkaFastError::StateTooWide => "state does not fit in u128".to_string(),
    }
}

fn format_signed_lr_error(error: SignedLrError) -> String {
    match error {
        SignedLrError::InvalidInput => "invalid input".to_string(),
        SignedLrError::ArithmeticOverflow => "arithmetic overflow".to_string(),
        SignedLrError::StateTooWide => "state does not fit in u128".to_string(),
    }
}

fn format_lr_gt_error(error: LrGtError) -> String {
    match error {
        LrGtError::InvalidInput => "invalid input".to_string(),
        LrGtError::ArithmeticOverflow => "arithmetic overflow".to_string(),
        LrGtError::StateTooWide => "state does not fit in u128".to_string(),
    }
}
