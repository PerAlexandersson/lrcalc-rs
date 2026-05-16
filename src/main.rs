use lrcalc::kostka::{kostka_lr_triple, kostka_via_lr};
use lrcalc::kostka_fast::{
    kostka_fast_stats, kostka_fast_u128, kostka_interior_stats, kostka_interior_u128,
    skew_kostka_fast_stats, skew_kostka_fast_u128, KostkaFastError,
};
use lrcalc::lr_ehrhart::{
    beta_lr_stretch_coefficient, beta_lr_stretch_dimension, beta_lr_stretch_h_vector,
    format_bigint_vector, format_error as format_lr_stretch_error, lr_stretch_coefficient,
    lr_stretch_dimension, lr_stretch_h_vector,
};
use lrcalc::lr_gt::{
    lrcoef_gt_hybrid_stats, lrcoef_gt_interior_dfs_stats, lrcoef_gt_interior_dfs_u128,
    lrcoef_gt_interior_stats, lrcoef_gt_interior_u128, lrcoef_gt_stats, lrcoef_gt_u128,
    lrcoef_tableau_hybrid_counts_stats, lrcoef_tableau_hybrid_stats, LrGtError,
};
use lrcalc::lr_signed::{lrcoef_signed_kostka, lrcoef_signed_kostka_stats, SignedLrError};
use lrcalc::lrcoef::{
    beta_lrcoef, beta_lrcoef_buch_counts_u128, beta_lrcoef_buch_dimension,
    beta_lrcoef_buch_interior_stats, beta_lrcoef_buch_interior_u128, lrcoef,
    lrcoef_buch_counts_u128, lrcoef_buch_dimension, lrcoef_buch_interior_memo_stats,
    lrcoef_buch_interior_memo_u128, lrcoef_buch_interior_stats, lrcoef_buch_interior_u128,
    LrBuchInteriorStats, LrCoefError,
};
use lrcalc::schur::{
    schur_product_expansion, schur_skew_expansion, SchurExpansionError, SchurTerm,
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
        Some("mult") => {
            let rest: Vec<String> = args.collect();
            match parse_mult_args(&rest).and_then(|parsed| {
                schur_product_expansion(&parsed.left, &parsed.right, parsed.rows, parsed.cols)
                    .map(|terms| (terms, parsed.maple))
                    .map_err(format_schur_error)
            }) {
                Ok((terms, maple)) => print_schur_terms(&terms, maple),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("skew") => {
            let rest: Vec<String> = args.collect();
            match parse_skew_args(&rest).and_then(|parsed| {
                schur_skew_expansion(&parsed.outer, &parsed.inner, parsed.rows)
                    .map(|terms| (terms, parsed.maple))
                    .map_err(format_schur_error)
            }) {
                Ok((terms, maple)) => print_schur_terms(&terms, maple),
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
        Some("skew-kostka-fast") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                skew_kostka_fast_u128(&parts[0], &parts[1], &parts[2])
                    .map_err(format_kostka_fast_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("skew-kostka-fast-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                skew_kostka_fast_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_kostka_fast_error)
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
        Some("lr-tableau-hybrid") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_tableau_hybrid_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_gt_error)
            }) {
                Ok(stats) => println!("{}", stats.value),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-tableau-hybrid-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_tableau_hybrid_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_gt_error)
            }) {
                Ok(stats) => {
                    println!("value: {}", stats.value);
                    println!("mode: {}", stats.mode.label());
                    println!("enforced_rows: {:?}", stats.enforced_rows);
                    println!("peak_states: {}", format_optional_usize(stats.peak_states));
                    println!(
                        "levels: {}",
                        format_optional_usize_vec(stats.levels.as_deref())
                    );
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-tableau-hybrid-counts") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_tableau_hybrid_counts_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_gt_error)
            }) {
                Ok(stats) => {
                    println!("full: {}", stats.full);
                    println!("interior: {}", stats.interior);
                }
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("lr-tableau-hybrid-counts-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lrcoef_tableau_hybrid_counts_stats(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_gt_error)
            }) {
                Ok(stats) => {
                    println!("full: {}", stats.full);
                    println!("interior: {}", stats.interior);
                    println!("mode: {}", stats.mode.label());
                    println!(
                        "full_peak_states: {}",
                        format_optional_usize(stats.full_peak_states)
                    );
                    println!(
                        "interior_peak_states: {}",
                        format_optional_usize(stats.interior_peak_states)
                    );
                    println!(
                        "full_levels: {}",
                        format_optional_usize_vec(stats.full_levels.as_deref())
                    );
                    println!(
                        "interior_levels: {}",
                        format_optional_usize_vec(stats.interior_levels.as_deref())
                    );
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
        Some("beta-lr") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_quad(&rest).and_then(|parts| {
                beta_lrcoef(&parts[0], &parts[1], &parts[2], &parts[3]).map_err(format_lrcoef_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("beta-lr-buch-interior") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_quad(&rest).and_then(|parts| {
                beta_lrcoef_buch_interior_u128(&parts[0], &parts[1], &parts[2], &parts[3])
                    .map_err(format_lrcoef_error)
            }) {
                Ok(coef) => println!("{coef}"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("beta-lr-buch-counts") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_quad(&rest).and_then(|parts| {
                beta_lrcoef_buch_counts_u128(&parts[0], &parts[1], &parts[2], &parts[3])
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
        Some("beta-lr-buch-interior-stats") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_quad(&rest).and_then(|parts| {
                beta_lrcoef_buch_interior_stats(&parts[0], &parts[1], &parts[2], &parts[3])
                    .map_err(format_lrcoef_error)
            }) {
                Ok(stats) => print_buch_interior_stats(&stats),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("beta-lr-buch-dimension") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_quad(&rest).and_then(|parts| {
                beta_lrcoef_buch_dimension(&parts[0], &parts[1], &parts[2], &parts[3])
                    .map_err(format_lrcoef_error)
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
        Some("lr-stretch-dimension") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_triple(&rest).and_then(|parts| {
                lr_stretch_dimension(&parts[0], &parts[1], &parts[2])
                    .map_err(format_lr_stretch_error)
            }) {
                Ok(Some(dimension)) => println!("{dimension}"),
                Ok(None) => println!("empty"),
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
        Some("beta-lr-stretch-hvector") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_quad(&rest).and_then(|parts| {
                beta_lr_stretch_h_vector(&parts[0], &parts[1], &parts[2], &parts[3])
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
        Some("beta-lr-stretch-dimension") => {
            let rest: Vec<String> = args.collect();
            match parse_partition_quad(&rest).and_then(|parts| {
                beta_lr_stretch_dimension(&parts[0], &parts[1], &parts[2], &parts[3])
                    .map_err(format_lr_stretch_error)
            }) {
                Ok(Some(dimension)) => println!("{dimension}"),
                Ok(None) => println!("empty"),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("beta-lr-stretch-eval") => {
            let rest: Vec<String> = args.collect();
            match parse_stretch_and_quad(&rest).and_then(|(stretch, parts)| {
                beta_lr_stretch_coefficient(&parts[0], &parts[1], &parts[2], &parts[3], stretch)
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
            eprintln!("Usage: {program} <coef|lr-buch-counts|beta-lr|...> [arguments]");
            std::process::exit(2);
        }
    }
}

struct MultArgs {
    left: Vec<i32>,
    right: Vec<i32>,
    rows: i32,
    cols: i32,
    maple: bool,
}

struct SkewArgs {
    outer: Vec<i32>,
    inner: Vec<i32>,
    rows: i32,
    maple: bool,
}

fn parse_mult_args(args: &[String]) -> Result<MultArgs, String> {
    let mut rows = -1;
    let mut cols = -1;
    let mut maple = false;
    let mut parts = Vec::new();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "-m" => {
                maple = true;
                index += 1;
            }
            "-r" => {
                let value = args
                    .get(index + 1)
                    .ok_or_else(|| "missing value after -r".to_string())?;
                rows = parse_i32_option(value, "rows")?;
                index += 2;
            }
            "-c" => {
                let value = args
                    .get(index + 1)
                    .ok_or_else(|| "missing value after -c".to_string())?;
                cols = parse_i32_option(value, "cols")?;
                index += 2;
            }
            "-q" | "-f" => {
                return Err(
                    "quantum/fusion Schur multiplication is not implemented yet".to_string()
                );
            }
            token => {
                parts.push(token.to_string());
                index += 1;
            }
        }
    }
    let [left, right] = parse_partition_pair_with_separator(
        &parts,
        "-",
        "usage: mult [-m] [-r rows] [-c cols] PART1 - PART2",
    )?;
    Ok(MultArgs {
        left,
        right,
        rows,
        cols,
        maple,
    })
}

fn parse_skew_args(args: &[String]) -> Result<SkewArgs, String> {
    let mut rows = -1;
    let mut maple = false;
    let mut parts = Vec::new();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "-m" => {
                maple = true;
                index += 1;
            }
            "-r" => {
                let value = args
                    .get(index + 1)
                    .ok_or_else(|| "missing value after -r".to_string())?;
                rows = parse_i32_option(value, "rows")?;
                index += 2;
            }
            token => {
                parts.push(token.to_string());
                index += 1;
            }
        }
    }
    let [outer, inner] = parse_partition_pair_with_separator(
        &parts,
        "/",
        "usage: skew [-m] [-r rows] OUTER / INNER",
    )?;
    Ok(SkewArgs {
        outer,
        inner,
        rows,
        maple,
    })
}

fn parse_i32_option(value: &str, name: &str) -> Result<i32, String> {
    value
        .parse::<i32>()
        .map_err(|_| format!("invalid {name} value '{value}'"))
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

fn parse_partition_quad(args: &[String]) -> Result<[Vec<i32>; 4], String> {
    let mut parts = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    let mut section = 0usize;

    for token in args {
        if token == "-" {
            section += 1;
            if section >= parts.len() {
                return Err("expected exactly three '-' separators".to_string());
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

    if section != 3 {
        return Err("usage: beta-lr OUTER - INNER - CONTENT - BETA".to_string());
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

fn parse_stretch_and_quad(args: &[String]) -> Result<(u64, [Vec<i32>; 4]), String> {
    let Some((stretch, rest)) = args.split_first() else {
        return Err(
            "usage: beta-lr-stretch-eval STRETCH OUTER - INNER - CONTENT - BETA".to_string(),
        );
    };
    let stretch = stretch
        .parse::<u64>()
        .map_err(|_| format!("invalid stretch factor '{stretch}'"))?;
    Ok((stretch, parse_partition_quad(rest)?))
}

fn parse_partition_pair(args: &[String]) -> Result<[Vec<i32>; 2], String> {
    parse_partition_pair_with_separator(args, "-", "usage: kostka-lr SHAPE - WEIGHT")
}

fn parse_partition_pair_with_separator(
    args: &[String],
    separator: &str,
    usage: &str,
) -> Result<[Vec<i32>; 2], String> {
    let mut parts = [Vec::new(), Vec::new()];
    let mut section = 0usize;

    for token in args {
        if token == separator {
            section += 1;
            if section >= parts.len() {
                return Err(format!("expected exactly one '{separator}' separator"));
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
        return Err(usage.to_string());
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

fn print_schur_terms(terms: &[SchurTerm], maple: bool) {
    if maple {
        print!("0");
        for term in terms {
            print!(
                "+{}*s[{}]",
                term.coefficient,
                format_comma_partition(&term.partition)
            );
        }
        println!();
        return;
    }

    for term in terms {
        println!(
            "{}  ({})",
            term.coefficient,
            format_comma_partition(&term.partition)
        );
    }
}

fn format_comma_partition(partition: &[i32]) -> String {
    partition
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn format_optional_usize(value: Option<usize>) -> String {
    value.map_or_else(|| "-".to_string(), |value| value.to_string())
}

fn format_optional_usize_vec(value: Option<&[usize]>) -> String {
    value.map_or_else(|| "-".to_string(), |value| format!("{value:?}"))
}

fn print_buch_interior_stats(stats: &LrBuchInteriorStats) {
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

fn format_lrcoef_error(error: LrCoefError) -> String {
    match error {
        LrCoefError::InvalidPartition => "invalid partition".to_string(),
        LrCoefError::ArithmeticOverflow => "arithmetic overflow".to_string(),
    }
}

fn format_schur_error(error: SchurExpansionError) -> String {
    match error {
        SchurExpansionError::InvalidPartition => "invalid partition".to_string(),
        SchurExpansionError::ArithmeticOverflow => "arithmetic overflow".to_string(),
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
