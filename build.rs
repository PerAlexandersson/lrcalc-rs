use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=src/native/lrcoef_count.c");
    println!("cargo:rerun-if-changed=src/native/abi_variadic.c");

    if env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("linux") {
        println!("cargo:rustc-cdylib-link-arg=-Wl,-soname,liblrcalc.so.2");
    }

    let out_dir = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR is set by Cargo"));
    let library = out_dir.join("liblrcalc_native.a");
    let cc = env::var_os("CC").unwrap_or_else(|| "cc".into());
    let ar = env::var_os("AR").unwrap_or_else(|| "ar".into());

    let sources = ["src/native/lrcoef_count.c", "src/native/abi_variadic.c"];
    let mut objects = Vec::with_capacity(sources.len());
    for source in sources {
        let stem = PathBuf::from(source)
            .file_stem()
            .expect("C source has a file stem")
            .to_owned();
        let object = out_dir.join(stem).with_extension("o");
        let cc_status = Command::new(&cc)
            .args(["-O3", "-fPIC", "-std=c99", "-Wall", "-Wextra", "-c"])
            .arg(source)
            .arg("-o")
            .arg(&object)
            .status()
            .expect("failed to launch C compiler");
        assert!(cc_status.success(), "C compiler failed on {source}");
        objects.push(object);
    }

    let ar_status = Command::new(&ar)
        .arg("crs")
        .arg(&library)
        .args(&objects)
        .status()
        .expect("failed to launch ar");
    assert!(ar_status.success(), "ar failed");

    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static:+whole-archive=lrcalc_native");
}
