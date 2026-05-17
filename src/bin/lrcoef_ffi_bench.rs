use lrcalc::lrcoef::lrcoef_i64;
use std::ffi::{CStr, CString};
use std::hint::black_box;
use std::ptr;
use std::time::{Duration, Instant};

#[repr(C)]
struct CVector {
    length: u32,
    array: [i32; 1],
}

type IvNew = unsafe extern "C" fn(u32) -> *mut CVector;
type IvFree = unsafe extern "C" fn(*mut CVector);
type SchurLrcoef = unsafe extern "C" fn(*mut CVector, *mut CVector, *mut CVector) -> i64;

struct Upstream {
    handle: *mut libc::c_void,
    iv_new: IvNew,
    iv_free: IvFree,
    schur_lrcoef: SchurLrcoef,
}

impl Upstream {
    fn load() -> Self {
        let path = std::env::var("UPSTREAM_LIB")
            .unwrap_or_else(|_| "/tmp/lrcalc-upstream/src/.libs/liblrcalc.so".to_string());
        let path = CString::new(path).expect("UPSTREAM_LIB contained an interior NUL");
        let flags = libc::RTLD_NOW | libc::RTLD_LOCAL | rtld_deepbind();
        let handle = unsafe { libc::dlopen(path.as_ptr(), flags) };
        if handle.is_null() {
            panic!("dlopen failed: {}", dl_error());
        }
        unsafe {
            Self {
                handle,
                iv_new: symbol(handle, b"iv_new\0"),
                iv_free: symbol(handle, b"iv_free\0"),
                schur_lrcoef: symbol(handle, b"schur_lrcoef\0"),
            }
        }
    }
}

impl Drop for Upstream {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe {
                libc::dlclose(self.handle);
            }
        }
    }
}

struct OwnedCVector {
    ptr: *mut CVector,
    iv_free: IvFree,
}

impl OwnedCVector {
    fn new(upstream: &Upstream, values: &[i32]) -> Self {
        let values = if values.is_empty() { &[0][..] } else { values };
        let length = u32::try_from(values.len()).expect("input vector too long for upstream");
        let ptr = unsafe { (upstream.iv_new)(length) };
        if ptr.is_null() {
            panic!("upstream iv_new returned null");
        }
        for (index, &value) in values.iter().enumerate() {
            unsafe {
                *ptr::addr_of_mut!((*ptr).array).cast::<i32>().add(index) = value;
            }
        }
        Self {
            ptr,
            iv_free: upstream.iv_free,
        }
    }

    fn as_mut_ptr(&self) -> *mut CVector {
        self.ptr
    }
}

impl Drop for OwnedCVector {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { (self.iv_free)(self.ptr) };
        }
    }
}

struct Case {
    label: &'static str,
    outer: Vec<i32>,
    inner: Vec<i32>,
    content: Vec<i32>,
    c_outer: OwnedCVector,
    c_inner: OwnedCVector,
    c_content: OwnedCVector,
}

fn main() {
    let mut repeat = 5_000usize;
    let mut suite = "mixed".to_string();
    for arg in std::env::args().skip(1) {
        if let Ok(value) = arg.parse::<usize>() {
            repeat = value;
        } else {
            suite = arg;
        }
    }

    let upstream = Upstream::load();
    let cases = build_cases(&upstream, &suite);

    for case in &cases {
        let rust_value = rust_lrcoef(case);
        let c_value = c_lrcoef(&upstream, case);
        if rust_value != c_value {
            eprintln!("mismatch: {}", case.label);
            eprintln!("outer:   {:?}", case.outer);
            eprintln!("inner:   {:?}", case.inner);
            eprintln!("content: {:?}", case.content);
            eprintln!("Rust: {rust_value}");
            eprintln!("C:    {c_value}");
            std::process::exit(1);
        }
    }

    println!("suite: lrcoef_ffi_bench");
    println!("case-suite: {suite}");
    println!("repeat: {repeat}");
    println!(
        "upstream: {}",
        std::env::var("UPSTREAM_LIB")
            .unwrap_or_else(|_| "/tmp/lrcalc-upstream/src/.libs/liblrcalc.so".to_string())
    );

    let (rust_time, rust_sink) = time_loop(repeat, &cases, rust_lrcoef);
    let (c_time, c_sink) = time_loop(repeat, &cases, |case| c_lrcoef(&upstream, case));
    black_box((rust_sink, c_sink));

    println!("correctness: ok ({} cases)", cases.len());
    println!(
        "Rust lrcoef_i64: {}  ({} evals)",
        format_duration(rust_time),
        repeat * cases.len()
    );
    println!(
        "C schur_lrcoef: {}  ({} evals)",
        format_duration(c_time),
        repeat * cases.len()
    );
    println!(
        "ratio Rust/C: {:.3}x",
        rust_time.as_secs_f64() / c_time.as_secs_f64()
    );
    print_case_diagnostics(repeat, &cases, &upstream);
}

fn rust_lrcoef(case: &Case) -> i64 {
    lrcoef_i64(&case.outer, &case.inner, &case.content)
        .unwrap_or_else(|| panic!("Rust LR failed for {}", case.label))
}

fn c_lrcoef(upstream: &Upstream, case: &Case) -> i64 {
    let value = unsafe {
        (upstream.schur_lrcoef)(
            case.c_outer.as_mut_ptr(),
            case.c_inner.as_mut_ptr(),
            case.c_content.as_mut_ptr(),
        )
    };
    if value < 0 {
        panic!("C schur_lrcoef returned {value} for {}", case.label);
    }
    value
}

fn time_loop<F>(repeat: usize, cases: &[Case], mut f: F) -> (Duration, i64)
where
    F: FnMut(&Case) -> i64,
{
    let start = Instant::now();
    let mut sink = 0i64;
    for _ in 0..repeat {
        for case in cases {
            sink ^= black_box(f(case));
        }
    }
    (start.elapsed(), sink)
}

fn time_one<F>(repeat: usize, case: &Case, mut f: F) -> (Duration, i64)
where
    F: FnMut(&Case) -> i64,
{
    let start = Instant::now();
    let mut sink = 0i64;
    for _ in 0..repeat {
        sink ^= black_box(f(case));
    }
    (start.elapsed(), sink)
}

fn print_case_diagnostics(repeat: usize, cases: &[Case], upstream: &Upstream) {
    println!("per-case diagnostics:");
    println!(
        "  {:<28} {:>7} {:>6} {:>6} {:>9} {:>9} {:>8}",
        "case", "value", "skew", "labels", "Rust", "C", "Rust/C"
    );
    for case in cases {
        let value = rust_lrcoef(case);
        let skew_size = case.outer.iter().sum::<i32>() - case.inner.iter().sum::<i32>();
        let labels = case.content.iter().take_while(|&&part| part != 0).count();
        let (rust_time, rust_sink) = time_one(repeat, case, rust_lrcoef);
        let (c_time, c_sink) = time_one(repeat, case, |case| c_lrcoef(upstream, case));
        black_box((rust_sink, c_sink));
        println!(
            "  {:<28} {:>7} {:>6} {:>6} {:>9} {:>9} {:>8.3}x",
            truncate_label(case.label, 28),
            value,
            skew_size,
            labels,
            format_duration(rust_time),
            format_duration(c_time),
            rust_time.as_secs_f64() / c_time.as_secs_f64()
        );
    }
}

fn build_cases(upstream: &Upstream, suite: &str) -> Vec<Case> {
    raw_cases(suite)
        .into_iter()
        .map(|(label, outer, inner, content)| Case {
            label,
            c_outer: OwnedCVector::new(upstream, &outer),
            c_inner: OwnedCVector::new(upstream, &inner),
            c_content: OwnedCVector::new(upstream, &content),
            outer,
            inner,
            content,
        })
        .collect()
}

fn raw_cases(suite: &str) -> Vec<(&'static str, Vec<i32>, Vec<i32>, Vec<i32>)> {
    match suite {
        "small" => raw_small_cases(),
        "mixed" => raw_mixed_cases(),
        "large-few-parts" | "large-short" | "large" => raw_large_short_cases(),
        "large-three-part" | "large-3part" => raw_large_three_part_cases(),
        other => {
            panic!(
                "unknown suite '{other}', expected small, mixed, large-few-parts, or large-three-part"
            )
        }
    }
}

fn raw_small_cases() -> Vec<(&'static str, Vec<i32>, Vec<i32>, Vec<i32>)> {
    vec![
        case("empty coefficient", &[0], &[0], &[0]),
        case("empty outer nonempty", &[0], &[0], &[1]),
        case("size mismatch zero", &[1], &[0], &[0]),
        case("single box Pieri", &[2, 1], &[2], &[1]),
        case("s21 s21 coefficient 2", &[3, 2, 1], &[2, 1], &[2, 1]),
        case("s21 s21 coefficient 1 A", &[4, 2], &[2, 1], &[2, 1]),
        case("s21 absent shape", &[5, 1], &[2, 1], &[2, 1]),
        case(
            "vertical strip short",
            &[4, 3, 2, 1],
            &[3, 2, 1],
            &[1, 1, 1],
        ),
    ]
}

fn raw_mixed_cases() -> Vec<(&'static str, Vec<i32>, Vec<i32>, Vec<i32>)> {
    vec![
        case("empty coefficient", &[0], &[0], &[0]),
        case("empty outer nonempty", &[0], &[0], &[1]),
        case("size mismatch zero", &[1], &[0], &[0]),
        case("inner not contained", &[2, 2, 1], &[3], &[2]),
        case("single box Pieri", &[2, 1], &[2], &[1]),
        case("s21 s21 coefficient 2", &[3, 2, 1], &[2, 1], &[2, 1]),
        case("s21 coefficient 1 A", &[4, 2], &[2, 1], &[2, 1]),
        case("s21 coefficient 1 B", &[3, 3], &[2, 1], &[2, 1]),
        case("s21 absent shape", &[5, 1], &[2, 1], &[2, 1]),
        case("horizontal strip short", &[5, 3, 1], &[3, 2, 1], &[2, 1]),
        case(
            "vertical strip short",
            &[4, 3, 2, 1],
            &[3, 2, 1],
            &[1, 1, 1],
        ),
        case("row times row", &[9], &[5], &[4]),
        case(
            "column times column",
            &[2, 2, 2, 1, 1],
            &[1, 1, 1, 1],
            &[1, 1, 1],
        ),
        case("medium", &[5, 4, 3, 2, 1], &[3, 2, 1], &[4, 3, 1, 1]),
        case(
            "larger",
            &[7, 6, 5, 4, 3, 2, 1],
            &[4, 4, 3, 2, 1],
            &[5, 4, 3, 2],
        ),
        case("stretched 321 x2", &[6, 4, 2], &[4, 2], &[4, 2]),
        case("stretched 321 x3", &[9, 6, 3], &[6, 3], &[6, 3]),
        case(
            "wide compact small",
            &[20, 14, 8, 2],
            &[12, 8, 4],
            &[14, 6, 2],
        ),
        case("wide compact medium", &[30, 20, 10], &[20, 10], &[20, 10]),
        case(
            "wide compact skewed",
            &[28, 21, 15, 6],
            &[18, 12, 6],
            &[20, 14],
        ),
        case(
            "irregular A",
            &[9, 7, 5, 3, 1],
            &[5, 3, 2, 1],
            &[6, 5, 3, 2],
        ),
        case(
            "irregular B",
            &[10, 8, 6, 4, 2],
            &[6, 4, 3, 1],
            &[7, 6, 4, 3],
        ),
        case(
            "irregular C",
            &[11, 9, 7, 4, 2, 1],
            &[7, 5, 3, 2],
            &[8, 7, 4, 2, 1],
        ),
        case(
            "irregular D",
            &[12, 10, 8, 5, 3, 1],
            &[7, 6, 4, 2],
            &[9, 8, 5, 3, 1],
        ),
    ]
}

fn raw_large_short_cases() -> Vec<(&'static str, Vec<i32>, Vec<i32>, Vec<i32>)> {
    vec![
        case("2-row wide A", &[1000, 700], &[600], &[700, 400]),
        case("2-row wide B", &[2000, 1000], &[1500], &[1000, 500]),
        case("2-row wide C", &[5000, 3000], &[3500], &[3000, 1500]),
        case("2-row wide D", &[8000, 5000], &[7000], &[4000, 2000]),
        case(
            "3-row balanced A",
            &[300, 200, 100],
            &[200, 100],
            &[200, 100],
        ),
        case(
            "3-row balanced B",
            &[450, 300, 150],
            &[300, 150],
            &[300, 150],
        ),
    ]
}

fn raw_large_three_part_cases() -> Vec<(&'static str, Vec<i32>, Vec<i32>, Vec<i32>)> {
    vec![
        case(
            "3-part stretched x40",
            &[120, 80, 40],
            &[60, 40, 20],
            &[60, 40, 20],
        ),
        case(
            "3-part stretched x60",
            &[180, 120, 60],
            &[90, 60, 30],
            &[90, 60, 30],
        ),
        case(
            "3-part stretched x80",
            &[240, 160, 80],
            &[120, 80, 40],
            &[120, 80, 40],
        ),
    ]
}

fn case(
    label: &'static str,
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> (&'static str, Vec<i32>, Vec<i32>, Vec<i32>) {
    (label, outer.to_vec(), inner.to_vec(), content.to_vec())
}

fn truncate_label(label: &str, width: usize) -> String {
    if label.len() <= width {
        label.to_string()
    } else {
        let mut out = label[..width.saturating_sub(1)].to_string();
        out.push('~');
        out
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

unsafe fn symbol<T: Copy>(handle: *mut libc::c_void, name: &[u8]) -> T {
    let ptr = unsafe { libc::dlsym(handle, name.as_ptr().cast()) };
    if ptr.is_null() {
        let name = CStr::from_bytes_with_nul(name)
            .expect("symbol name should be NUL-terminated")
            .to_string_lossy();
        panic!("dlsym failed for {name}: {}", dl_error());
    }
    unsafe { std::mem::transmute_copy(&ptr) }
}

fn dl_error() -> String {
    let error = unsafe { libc::dlerror() };
    if error.is_null() {
        "unknown dynamic loader error".to_string()
    } else {
        unsafe { CStr::from_ptr(error) }
            .to_string_lossy()
            .into_owned()
    }
}

fn rtld_deepbind() -> libc::c_int {
    #[cfg(target_os = "linux")]
    {
        0x00008
    }
    #[cfg(not(target_os = "linux"))]
    {
        0
    }
}
