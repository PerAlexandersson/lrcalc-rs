#![allow(clippy::items_after_test_module)]

use lrcalc::abi::{
    iv_free, iv_new_zero, lrit_free, lrit_good, lrit_new, lrit_next, IVector, LrTabIter,
};
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
    schur_coproduct_expansion, schur_product_expansion, schur_product_fusion_expansion,
    schur_skew_expansion, SchurExpansionError, SchurTerm, SignedSchurTerm,
};
use std::{ptr, slice};

fn main() {
    let mut args = std::env::args();
    let program = args.next().unwrap_or_else(|| "lrcalc".to_string());
    match args.next().as_deref() {
        Some("--version") | Some("-V") => {
            println!("lrcalc-rs 0.1.0");
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
                match parsed.mode {
                    MultMode::Ordinary { rows, cols } => {
                        let terms =
                            schur_product_expansion(&parsed.left, &parsed.right, rows, cols)
                                .map_err(format_schur_error)?;
                        print_schur_terms(&terms, parsed.maple);
                    }
                    MultMode::Fusion { rows, level } => {
                        let terms = schur_product_fusion_expansion(
                            &parsed.left,
                            &parsed.right,
                            rows,
                            level,
                        )
                        .map_err(format_schur_error)?;
                        print_signed_schur_terms(&terms, parsed.maple);
                    }
                    MultMode::Quantum { rows, level } => {
                        let terms = schur_product_fusion_expansion(
                            &parsed.left,
                            &parsed.right,
                            rows,
                            level,
                        )
                        .map_err(format_schur_error)?;
                        print_quantum_schur_terms(&terms, parsed.maple, rows, level);
                    }
                }
                Ok(())
            }) {
                Ok(()) => {}
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
        Some("coprod") => {
            let rest: Vec<String> = args.collect();
            match parse_coprod_args(&rest).and_then(|parsed| {
                schur_coproduct_expansion(&parsed.shape, parsed.rows, parsed.cols, parsed.all)
                    .map(|terms| (terms, parsed.rows, parsed.cols))
                    .map_err(format_schur_error)
            }) {
                Ok((terms, rows, cols)) => print_coproduct_terms(&terms, rows, cols),
                Err(message) => {
                    eprintln!("{program}: {message}");
                    std::process::exit(2);
                }
            }
        }
        Some("tab") => {
            let rest: Vec<String> = args.collect();
            match parse_tab_args(&rest).and_then(run_tab_command) {
                Ok(()) => {}
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
    mode: MultMode,
    maple: bool,
}

enum MultMode {
    Ordinary { rows: i32, cols: i32 },
    Fusion { rows: i32, level: i32 },
    Quantum { rows: i32, level: i32 },
}

struct SkewArgs {
    outer: Vec<i32>,
    inner: Vec<i32>,
    rows: i32,
    maple: bool,
}

struct CoprodArgs {
    shape: Vec<i32>,
    rows: i32,
    cols: i32,
    all: bool,
}

struct TabArgs {
    outer: Vec<i32>,
    inner: Vec<i32>,
    weight: Option<Vec<i32>>,
    rows: i32,
}

fn expand_short_options(args: &[String], value_options: &[char]) -> Vec<String> {
    let mut expanded = Vec::with_capacity(args.len());
    for arg in args {
        if arg == "-" || !arg.starts_with('-') || arg.starts_with("--") || arg.len() <= 2 {
            expanded.push(arg.clone());
            continue;
        }

        let mut split_value = None;
        for (offset, option) in arg[1..].char_indices() {
            if !option.is_ascii_alphabetic() {
                expanded.push(arg.clone());
                split_value = Some(arg.len());
                break;
            }
            expanded.push(format!("-{option}"));
            if value_options.contains(&option) {
                let value_start = 1 + offset + option.len_utf8();
                if value_start < arg.len() {
                    expanded.push(arg[value_start..].to_string());
                }
                split_value = Some(arg.len());
                break;
            }
        }
        if split_value.is_none() && arg.len() == 2 {
            expanded.push(arg.clone());
        }
    }
    expanded
}

fn parse_mult_args(args: &[String]) -> Result<MultArgs, String> {
    let args = expand_short_options(args, &['r', 'c', 'f', 'q']);
    let mut rows = -1;
    let mut cols = -1;
    let mut fusion = false;
    let mut quantum = false;
    let mut maple = false;
    let mut parts = Vec::new();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--" => {
                parts.extend(args[index + 1..].iter().cloned());
                break;
            }
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
            "-f" => {
                let value = args
                    .get(index + 1)
                    .ok_or_else(|| "missing value after -f".to_string())?;
                let (fusion_rows, level) = parse_i32_pair_option(value, "rows", "level")?;
                if fusion_rows < 0 || level < 0 {
                    return Err("fusion rows and level must be nonnegative".to_string());
                }
                fusion = true;
                rows = fusion_rows;
                cols = level;
                index += 2;
            }
            "-q" => {
                let value = args
                    .get(index + 1)
                    .ok_or_else(|| "missing value after -q".to_string())?;
                let (quantum_rows, level) = parse_i32_pair_option(value, "rows", "level")?;
                if quantum_rows < 0 || level < 0 {
                    return Err("quantum rows and level must be nonnegative".to_string());
                }
                quantum = true;
                rows = quantum_rows;
                cols = level;
                index += 2;
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
        "usage: mult [-m] [-r rows] [-c cols] [-q rows,level] [-f rows,level] PART1 - PART2",
    )?;
    let mode = if quantum {
        if rows < 0 || cols < 0 {
            return Err("quantum rows and level must be nonnegative".to_string());
        }
        MultMode::Quantum { rows, level: cols }
    } else if fusion {
        if rows < 0 || cols < 0 {
            return Err("fusion rows and level must be nonnegative".to_string());
        }
        MultMode::Fusion { rows, level: cols }
    } else {
        MultMode::Ordinary { rows, cols }
    };
    Ok(MultArgs {
        left,
        right,
        mode,
        maple,
    })
}

fn parse_skew_args(args: &[String]) -> Result<SkewArgs, String> {
    let args = expand_short_options(args, &['r']);
    let mut rows = -1;
    let mut maple = false;
    let mut parts = Vec::new();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--" => {
                parts.extend(args[index + 1..].iter().cloned());
                break;
            }
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

fn parse_coprod_args(args: &[String]) -> Result<CoprodArgs, String> {
    let args = expand_short_options(args, &[]);
    let mut all = false;
    let mut parts = Vec::new();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--" => {
                parts.extend(args[index + 1..].iter().cloned());
                break;
            }
            "-a" => {
                all = true;
                index += 1;
            }
            token => {
                parts.push(token.to_string());
                index += 1;
            }
        }
    }
    let shape = parse_partition_args(&parts, "usage: coprod [-a] PART")?;
    let rows = i32::try_from(partition_length(&shape))
        .map_err(|_| "partition length overflow".to_string())?;
    let cols = shape.first().copied().unwrap_or(0);
    Ok(CoprodArgs {
        shape,
        rows,
        cols,
        all,
    })
}

fn parse_tab_args(args: &[String]) -> Result<TabArgs, String> {
    let args = expand_short_options(args, &['r']);
    let mut rows = -1;
    let mut parts = Vec::new();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--" => {
                parts.extend(args[index + 1..].iter().cloned());
                break;
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
    let parsed = parse_tab_partitions(&parts)?;
    Ok(TabArgs {
        outer: parsed[0].clone(),
        inner: parsed[1].clone(),
        weight: if parsed[2].is_empty() {
            None
        } else {
            Some(parsed[2].clone())
        },
        rows,
    })
}

fn parse_i32_option(value: &str, name: &str) -> Result<i32, String> {
    value
        .parse::<i32>()
        .map_err(|_| format!("invalid {name} value '{value}'"))
}

fn parse_i32_pair_option(
    value: &str,
    first_name: &str,
    second_name: &str,
) -> Result<(i32, i32), String> {
    let Some((left, right)) = value.split_once(',') else {
        return Err(format!(
            "expected {first_name},{second_name} value, got '{value}'"
        ));
    };
    let left = parse_i32_option(left, first_name)?;
    let right = parse_i32_option(right, second_name)?;
    Ok((left, right))
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

fn parse_tab_partitions(args: &[String]) -> Result<[Vec<i32>; 3], String> {
    let usage = "usage: tab [-r rows] OUTER / INNER [- WEIGHT]";
    let mut parts = [Vec::new(), Vec::new(), Vec::new()];
    let mut section = 0usize;
    let mut saw_slash = false;
    let mut saw_weight = false;

    for token in args {
        if token == "/" {
            if saw_slash || saw_weight {
                return Err(usage.to_string());
            }
            saw_slash = true;
            section = 1;
            continue;
        }
        if token == "-" {
            if !saw_slash || saw_weight {
                return Err(usage.to_string());
            }
            saw_weight = true;
            section = 2;
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

    if !saw_slash {
        return Err(usage.to_string());
    }
    Ok(parts)
}

fn parse_partition_args(args: &[String], usage: &str) -> Result<Vec<i32>, String> {
    if args.is_empty() {
        return Err(usage.to_string());
    }
    let mut part = Vec::new();
    for token in args {
        for piece in token.split(',') {
            let piece = piece.trim_matches(|ch| matches!(ch, '(' | ')' | '[' | ']'));
            if piece.is_empty() {
                continue;
            }
            let value = piece
                .parse::<i32>()
                .map_err(|_| format!("invalid integer '{piece}'"))?;
            part.push(value);
        }
    }
    Ok(part)
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

fn partition_length(partition: &[i32]) -> usize {
    partition
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1)
}

fn valid_cli_partition(partition: &[i32]) -> bool {
    let mut previous = 0;
    for &part in partition.iter().rev() {
        if part < previous {
            return false;
        }
        previous = part;
    }
    true
}

fn run_tab_command(parsed: TabArgs) -> Result<(), String> {
    if !valid_cli_partition(&parsed.outer)
        || !valid_cli_partition(&parsed.inner)
        || parsed
            .weight
            .as_deref()
            .is_some_and(|weight| !valid_cli_partition(weight))
    {
        return Err("invalid partition".to_string());
    }

    unsafe {
        let outer = abi_vector_from_partition(&parsed.outer)?;
        let inner = abi_vector_from_partition(&parsed.inner)?;
        let lrit = lrit_new(outer, inner, ptr::null(), parsed.rows, -1, -1);
        if lrit.is_null() {
            iv_free(inner);
            iv_free(outer);
            return Err("out of memory".to_string());
        }

        while lrit_good(lrit) != 0 {
            if parsed
                .weight
                .as_deref()
                .is_none_or(|weight| tab_weight_matches(lrit, weight))
            {
                print_lrit_skewtab(lrit, &parsed.outer, &parsed.inner);
                println!();
            }
            lrit_next(lrit);
        }

        lrit_free(lrit);
        iv_free(inner);
        iv_free(outer);
    }
    Ok(())
}

unsafe fn abi_vector_from_partition(partition: &[i32]) -> Result<*mut IVector, String> {
    let length =
        u32::try_from(partition.len()).map_err(|_| "partition length overflow".to_string())?;
    let vector = iv_new_zero(length);
    if vector.is_null() {
        return Err("out of memory".to_string());
    }
    unsafe {
        abi_vector_values_mut(vector).copy_from_slice(partition);
    }
    Ok(vector)
}

unsafe fn abi_vector_values_mut<'a>(vector: *mut IVector) -> &'a mut [i32] {
    let length = unsafe { (*vector).length as usize };
    let data = unsafe { ptr::addr_of_mut!((*vector).array).cast::<i32>() };
    unsafe { slice::from_raw_parts_mut(data, length) }
}

unsafe fn abi_vector_values<'a>(vector: *const IVector) -> &'a [i32] {
    let length = unsafe { (*vector).length as usize };
    let data = unsafe { ptr::addr_of!((*vector).array).cast::<i32>() };
    unsafe { slice::from_raw_parts(data, length) }
}

unsafe fn lrit_array<'a>(lrit: *const LrTabIter) -> &'a [lrcalc::abi::LritBox] {
    let length = unsafe { (*lrit).array_len as usize };
    let data = unsafe { ptr::addr_of!((*lrit).array).cast::<lrcalc::abi::LritBox>() };
    unsafe { slice::from_raw_parts(data, length) }
}

unsafe fn tab_weight_matches(lrit: *const LrTabIter, weight: &[i32]) -> bool {
    let content = unsafe { abi_vector_values((*lrit).cont) };
    let weight_len = partition_length(weight);
    if partition_length(content) != weight_len {
        return false;
    }
    content
        .iter()
        .take(weight_len)
        .eq(weight.iter().take(weight_len))
}

unsafe fn print_lrit_skewtab(lrit: *const LrTabIter, outer: &[i32], inner: &[i32]) {
    let array = unsafe { lrit_array(lrit) };
    let mut size = unsafe { (*lrit).size };
    let ilen = inner.len();
    let mut len = partition_length(outer);
    if len <= ilen {
        while len > 0 && inner.get(len - 1).copied().unwrap_or(0) == outer[len - 1] {
            len -= 1;
        }
    }
    if len == 0 {
        return;
    }

    let col_first = if ilen < len {
        0
    } else {
        inner.get(len - 1).copied().unwrap_or(0)
    };
    let mut row = 0usize;
    while row < ilen && inner[row] == outer[row] {
        row += 1;
    }
    while row < len {
        let inn_r = inner.get(row).copied().unwrap_or(0);
        let out_r = outer[row];
        let row_size = out_r - inn_r;
        size -= row_size;
        for _ in col_first..inn_r {
            print!("  ");
        }
        for col in 0..row_size {
            let index = usize::try_from(size + col).unwrap_or(0);
            print!("{:2}", array[index].value);
        }
        println!();
        row += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
    }

    #[test]
    fn mult_parser_accepts_clustered_short_options() {
        let parsed = parse_mult_args(&args(&["-mr", "3", "2", "1", "-", "1"])).unwrap();
        assert!(parsed.maple);
        assert_eq!(parsed.left, vec![2, 1]);
        assert_eq!(parsed.right, vec![1]);
        match parsed.mode {
            MultMode::Ordinary { rows, cols } => {
                assert_eq!(rows, 3);
                assert_eq!(cols, -1);
            }
            _ => panic!("expected ordinary mode"),
        }

        let parsed = parse_mult_args(&args(&["-q0,2", "1", "-", "1"])).unwrap();
        match parsed.mode {
            MultMode::Quantum { rows, level } => {
                assert_eq!(rows, 0);
                assert_eq!(level, 2);
            }
            _ => panic!("expected quantum mode"),
        }

        let parsed = parse_mult_args(&args(&["-q3,2", "-r2", "2", "1", "-", "2", "1"])).unwrap();
        match parsed.mode {
            MultMode::Quantum { rows, level } => {
                assert_eq!(rows, 2);
                assert_eq!(level, 2);
            }
            _ => panic!("expected quantum mode"),
        }

        let parsed =
            parse_mult_args(&args(&["-m", "-q3,2", "-f3,2", "2", "1", "-", "2", "1"])).unwrap();
        assert!(parsed.maple);
        match parsed.mode {
            MultMode::Quantum { rows, level } => {
                assert_eq!(rows, 3);
                assert_eq!(level, 2);
            }
            _ => panic!("expected quantum mode"),
        }

        let parsed = parse_mult_args(&args(&["--", "1", "-", "1"])).unwrap();
        assert_eq!(parsed.left, vec![1]);
        assert_eq!(parsed.right, vec![1]);
    }

    #[test]
    fn skew_and_tab_parsers_accept_clustered_row_option() {
        let skew = parse_skew_args(&args(&["-mr", "2", "3", "2", "/", "1"])).unwrap();
        assert!(skew.maple);
        assert_eq!(skew.rows, 2);

        let tab = parse_tab_args(&args(&["-r2", "3", "2", "/", "1"])).unwrap();
        assert_eq!(tab.rows, 2);
    }
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

fn print_signed_schur_terms(terms: &[SignedSchurTerm], maple: bool) {
    if maple {
        print!("0");
        for term in terms {
            print_signed_maple_prefix(term.coefficient);
            print!("*s[{}]", format_comma_partition(&term.partition));
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

fn print_quantum_schur_terms(terms: &[SignedSchurTerm], maple: bool, rows: i32, level: i32) {
    let rows = usize::try_from(rows).unwrap_or(0);
    if maple {
        print!("0");
        for term in terms {
            let (degree, partition) = quantum_partition_and_degree(&term.partition, rows, level);
            print_signed_maple_prefix(term.coefficient);
            print!("*q^{degree}*s[{}]", format_comma_i64_partition(&partition));
        }
        println!();
        return;
    }

    for term in terms {
        let (_, partition) = quantum_partition_and_degree(&term.partition, rows, level);
        println!(
            "{}  ({})",
            term.coefficient,
            format_comma_i64_partition(&partition)
        );
    }
}

fn print_signed_maple_prefix(coefficient: i128) {
    if coefficient < 0 {
        print!("-{}", coefficient.unsigned_abs());
    } else {
        print!("+{coefficient}");
    }
}

fn print_coproduct_terms(terms: &[SchurTerm], rows: i32, cols: i32) {
    let rows = usize::try_from(rows).unwrap_or(0);
    for term in terms {
        let (left, right) = coproduct_pair(&term.partition, rows, cols);
        println!(
            "{}  ({})  ({})",
            term.coefficient,
            format_comma_partition(&left),
            format_comma_partition(&right)
        );
    }
}

fn coproduct_pair(partition: &[i32], rows: usize, cols: i32) -> (Vec<i32>, Vec<i32>) {
    let mut left = Vec::new();
    for row in 0..rows {
        let part = partition.get(row).copied().unwrap_or(0);
        if part <= cols {
            break;
        }
        left.push(part - cols);
    }

    let mut right = Vec::new();
    for &part in partition.iter().skip(rows) {
        if part == 0 {
            break;
        }
        right.push(part);
    }
    (left, right)
}

fn format_comma_partition(partition: &[i32]) -> String {
    partition
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn format_comma_i64_partition(partition: &[i64]) -> String {
    partition
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn quantum_partition_and_degree(partition: &[i32], rows: usize, level: i32) -> (i64, Vec<i64>) {
    let degree = quantum_degree(partition, rows, level);
    let rows_i64 = i64::try_from(rows).unwrap_or(1).max(1);
    let level_i64 = i64::from(level);
    let mut out = Vec::new();
    for index in 0..rows {
        let shifted = i64::try_from(index).unwrap_or(0) + degree;
        let source = usize::try_from(shifted.rem_euclid(rows_i64)).unwrap_or(0);
        let entry = i64::from(partition.get(source).copied().unwrap_or(0))
            - (shifted / rows_i64) * level_i64
            - degree;
        if entry == 0 {
            break;
        }
        out.push(entry);
    }
    (degree, out)
}

fn quantum_degree(partition: &[i32], rows: usize, level: i32) -> i64 {
    let rows_i64 = i64::try_from(rows).unwrap_or(1).max(1);
    let n = rows_i64 + i64::from(level);
    let mut degree = 0i64;
    for index in 0..rows {
        let index_i64 = i64::try_from(index).unwrap_or(0);
        let a = i64::from(partition.get(index).copied().unwrap_or(0)) + rows_i64 - index_i64 - 1;
        degree += a.div_euclid(n);
    }
    degree
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
