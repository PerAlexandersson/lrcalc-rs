use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=src/native/lrcoef_count.c");

    let out_dir = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR is set by Cargo"));
    let object = out_dir.join("lrcoef_count.o");
    let library = out_dir.join("liblrcalc_native.a");
    let cc = env::var_os("CC").unwrap_or_else(|| "cc".into());
    let ar = env::var_os("AR").unwrap_or_else(|| "ar".into());

    let cc_status = Command::new(&cc)
        .args(["-O3", "-fPIC", "-std=c99", "-Wall", "-Wextra", "-c"])
        .arg("src/native/lrcoef_count.c")
        .arg("-o")
        .arg(&object)
        .status()
        .expect("failed to launch C compiler");
    assert!(cc_status.success(), "C compiler failed");

    let ar_status = Command::new(&ar)
        .arg("crs")
        .arg(&library)
        .arg(&object)
        .status()
        .expect("failed to launch ar");
    assert!(ar_status.success(), "ar failed");

    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=lrcalc_native");
}
