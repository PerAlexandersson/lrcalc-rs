use lrcalc::lr_ehrhart::{compute_beta_lr_stretch_polynomial, evaluate_h_vector};
use num_bigint::BigInt;

#[test]
fn dominant_beta_shortcut_preserves_inputs_wider_than_the_packed_engine() {
    // A six-box one-dimensional fiber, but 9 * 17 = 153 state bits.
    let mut outer = vec![65_536; 9];
    outer[0] += 3;
    outer[1] += 2;
    outer[2] += 1;
    let inner = vec![65_536; 9];
    let polynomial =
        compute_beta_lr_stretch_polynomial(&outer, &inner, &[1, 1, 3, 1], &[9, 7, 3]).unwrap();
    assert_eq!(polynomial.dimension, 1);
    for n in [0, 1, 2, 7] {
        assert_eq!(
            evaluate_h_vector(&polynomial.h_vector, polynomial.dimension, n),
            BigInt::from(n + 1)
        );
    }
}

#[test]
fn cli_reports_the_manifest_version() {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_lrcalc"))
        .arg("--version")
        .output()
        .unwrap();
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap().trim(),
        format!("lrcalc-rs {}", env!("CARGO_PKG_VERSION"))
    );
}
