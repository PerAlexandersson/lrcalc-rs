use lrcalc::schubert::{multiply_schubert, multiply_schubert_strings};
use lrcalc::schur::{
    schur_coproduct_expansion, schur_product_expansion, schur_product_fusion_expansion,
    schur_skew_expansion,
};
use std::ffi::{CStr, CString};
use std::hint::black_box;
use std::os::raw::c_int;
use std::ptr;
use std::time::{Duration, Instant};

type Terms = Vec<(Vec<i32>, i128)>;

#[derive(Clone, Copy)]
enum KeyMode {
    Exact,
    TrimPartition,
}

#[repr(C)]
struct CVector {
    length: u32,
    array: [i32; 1],
}

#[repr(C)]
struct CIvLinComb {
    _private: [u8; 0],
}

#[repr(C)]
#[derive(Default)]
struct CIvLcIter {
    ht: *mut CIvLinComb,
    index: usize,
    i: usize,
}

type IvNew = unsafe extern "C" fn(u32) -> *mut CVector;
type IvFree = unsafe extern "C" fn(*mut CVector);
type IvLcFreeAll = unsafe extern "C" fn(*mut CIvLinComb);
type IvLcFirst = unsafe extern "C" fn(*mut CIvLinComb, *mut CIvLcIter);
type IvLcGood = unsafe extern "C" fn(*mut CIvLcIter) -> c_int;
type IvLcNext = unsafe extern "C" fn(*mut CIvLcIter);
type IvLcKey = unsafe extern "C" fn(*mut CIvLcIter) -> *mut CVector;
type IvLcValue = unsafe extern "C" fn(*mut CIvLcIter) -> i32;
type SchurMult =
    unsafe extern "C" fn(*mut CVector, *mut CVector, c_int, c_int, c_int) -> *mut CIvLinComb;
type SchurSkew = unsafe extern "C" fn(*mut CVector, *mut CVector, c_int, c_int) -> *mut CIvLinComb;
type SchurCoprod =
    unsafe extern "C" fn(*mut CVector, c_int, c_int, c_int, c_int) -> *mut CIvLinComb;
type SchurMultFusion =
    unsafe extern "C" fn(*mut CVector, *mut CVector, c_int, c_int) -> *mut CIvLinComb;
type MultSchubert = unsafe extern "C" fn(*mut CVector, *mut CVector, c_int) -> *mut CIvLinComb;
type MultSchubertStr = unsafe extern "C" fn(*mut CVector, *mut CVector) -> *mut CIvLinComb;

struct Upstream {
    handle: *mut libc::c_void,
    iv_new: IvNew,
    iv_free: IvFree,
    ivlc_free_all: IvLcFreeAll,
    ivlc_first: IvLcFirst,
    ivlc_good: IvLcGood,
    ivlc_next: IvLcNext,
    ivlc_key: IvLcKey,
    ivlc_value: IvLcValue,
    schur_mult: SchurMult,
    schur_skew: SchurSkew,
    schur_coprod: SchurCoprod,
    schur_mult_fusion: SchurMultFusion,
    mult_schubert: MultSchubert,
    mult_schubert_str: MultSchubertStr,
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
                ivlc_free_all: symbol(handle, b"ivlc_free_all\0"),
                ivlc_first: symbol(handle, b"ivlc_first\0"),
                ivlc_good: symbol(handle, b"ivlc_good\0"),
                ivlc_next: symbol(handle, b"ivlc_next\0"),
                ivlc_key: symbol(handle, b"ivlc_key\0"),
                ivlc_value: symbol(handle, b"ivlc_value\0"),
                schur_mult: symbol(handle, b"schur_mult\0"),
                schur_skew: symbol(handle, b"schur_skew\0"),
                schur_coprod: symbol(handle, b"schur_coprod\0"),
                schur_mult_fusion: symbol(handle, b"schur_mult_fusion\0"),
                mult_schubert: symbol(handle, b"mult_schubert\0"),
                mult_schubert_str: symbol(handle, b"mult_schubert_str\0"),
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

struct ProductCase {
    label: &'static str,
    sh1: Vec<i32>,
    sh2: Vec<i32>,
    rows: i32,
    cols: i32,
    c_sh1: OwnedCVector,
    c_sh2: OwnedCVector,
}

struct SkewCase {
    label: &'static str,
    outer: Vec<i32>,
    inner: Vec<i32>,
    rows: i32,
    c_outer: OwnedCVector,
    c_inner: OwnedCVector,
}

struct CoprodCase {
    label: &'static str,
    shape: Vec<i32>,
    rows: i32,
    cols: i32,
    all: bool,
    c_shape: OwnedCVector,
}

struct FusionCase {
    label: &'static str,
    sh1: Vec<i32>,
    sh2: Vec<i32>,
    rows: i32,
    level: i32,
    c_sh1: OwnedCVector,
    c_sh2: OwnedCVector,
}

struct SchubertCase {
    label: &'static str,
    left: Vec<i32>,
    right: Vec<i32>,
    rank: i32,
    c_left: OwnedCVector,
    c_right: OwnedCVector,
}

struct SchubertStringCase {
    label: &'static str,
    left: Vec<i32>,
    right: Vec<i32>,
    c_left: OwnedCVector,
    c_right: OwnedCVector,
}

trait BenchCase {
    fn label(&self) -> &str;
}

impl BenchCase for ProductCase {
    fn label(&self) -> &str {
        self.label
    }
}

impl BenchCase for SkewCase {
    fn label(&self) -> &str {
        self.label
    }
}

impl BenchCase for CoprodCase {
    fn label(&self) -> &str {
        self.label
    }
}

impl BenchCase for FusionCase {
    fn label(&self) -> &str {
        self.label
    }
}

impl BenchCase for SchubertCase {
    fn label(&self) -> &str {
        self.label
    }
}

impl BenchCase for SchubertStringCase {
    fn label(&self) -> &str {
        self.label
    }
}

fn main() {
    let repeat = std::env::args()
        .nth(1)
        .map(|arg| arg.parse::<usize>().expect("repeat must be an integer"))
        .unwrap_or(500);
    let upstream = Upstream::load();

    let product_cases = product_cases(&upstream);
    let skew_cases = skew_cases(&upstream);
    let coprod_cases = coprod_cases(&upstream);
    let fusion_cases = fusion_cases(&upstream);
    let schubert_cases = schubert_cases(&upstream);
    let schubert_string_cases = schubert_string_cases(&upstream);

    println!("suite: schur_schubert_ffi_bench");
    println!("repeat: {repeat}");
    println!(
        "upstream: {}",
        std::env::var("UPSTREAM_LIB")
            .unwrap_or_else(|_| { "/tmp/lrcalc-upstream/src/.libs/liblrcalc.so".to_string() })
    );

    run_suite(
        "Schur product",
        repeat,
        &product_cases,
        KeyMode::TrimPartition,
        |case| {
            schur_terms_to_terms(
                schur_product_expansion(&case.sh1, &case.sh2, case.rows, case.cols)
                    .expect(case.label),
            )
        },
        |case| unsafe {
            collect_lc(
                &upstream,
                (upstream.schur_mult)(
                    case.c_sh1.as_mut_ptr(),
                    case.c_sh2.as_mut_ptr(),
                    case.rows,
                    case.cols,
                    -1,
                ),
                KeyMode::TrimPartition,
            )
        },
        false,
    );

    run_suite(
        "Skew Schur",
        repeat,
        &skew_cases,
        KeyMode::TrimPartition,
        |case| {
            schur_terms_to_terms(
                schur_skew_expansion(&case.outer, &case.inner, case.rows).expect(case.label),
            )
        },
        |case| unsafe {
            collect_lc(
                &upstream,
                (upstream.schur_skew)(
                    case.c_outer.as_mut_ptr(),
                    case.c_inner.as_mut_ptr(),
                    case.rows,
                    -1,
                ),
                KeyMode::TrimPartition,
            )
        },
        false,
    );

    run_suite(
        "Schur coproduct",
        repeat,
        &coprod_cases,
        KeyMode::TrimPartition,
        |case| {
            schur_terms_to_terms(
                schur_coproduct_expansion(&case.shape, case.rows, case.cols, case.all)
                    .expect(case.label),
            )
        },
        |case| unsafe {
            collect_lc(
                &upstream,
                (upstream.schur_coprod)(
                    case.c_shape.as_mut_ptr(),
                    case.rows,
                    case.cols,
                    -1,
                    case.all as c_int,
                ),
                KeyMode::TrimPartition,
            )
        },
        false,
    );

    run_suite(
        "Fusion product",
        repeat,
        &fusion_cases,
        KeyMode::TrimPartition,
        |case| {
            signed_schur_terms_to_terms(
                schur_product_fusion_expansion(&case.sh1, &case.sh2, case.rows, case.level)
                    .expect(case.label),
            )
        },
        |case| unsafe {
            collect_lc(
                &upstream,
                (upstream.schur_mult_fusion)(
                    case.c_sh1.as_mut_ptr(),
                    case.c_sh2.as_mut_ptr(),
                    case.rows,
                    case.level,
                ),
                KeyMode::TrimPartition,
            )
        },
        false,
    );

    run_suite(
        "Schubert permutations",
        repeat,
        &schubert_cases,
        KeyMode::Exact,
        |case| {
            schubert_terms_to_terms(
                multiply_schubert(&case.left, &case.right, case.rank).expect(case.label),
            )
        },
        |case| unsafe {
            collect_lc(
                &upstream,
                (upstream.mult_schubert)(
                    case.c_left.as_mut_ptr(),
                    case.c_right.as_mut_ptr(),
                    case.rank,
                ),
                KeyMode::Exact,
            )
        },
        true,
    );

    run_suite(
        "Schubert strings",
        repeat,
        &schubert_string_cases,
        KeyMode::Exact,
        |case| {
            schubert_terms_to_terms(
                multiply_schubert_strings(&case.left, &case.right).expect(case.label),
            )
        },
        |case| unsafe {
            collect_lc(
                &upstream,
                (upstream.mult_schubert_str)(case.c_left.as_mut_ptr(), case.c_right.as_mut_ptr()),
                KeyMode::Exact,
            )
        },
        true,
    );
}

fn run_suite<C, FRust, FC>(
    name: &str,
    repeat: usize,
    cases: &[C],
    mode: KeyMode,
    mut rust_eval: FRust,
    mut c_eval: FC,
    detail: bool,
) where
    C: BenchCase,
    FRust: FnMut(&C) -> Terms,
    FC: FnMut(&C) -> Terms,
{
    for case in cases {
        let rust_terms = normalize_terms(rust_eval(case), mode);
        let c_terms = normalize_terms(c_eval(case), mode);
        if rust_terms != c_terms {
            eprintln!("{name}: mismatch in {}", case.label());
            eprintln!("rust: {rust_terms:?}");
            eprintln!("C:    {c_terms:?}");
            std::process::exit(1);
        }
    }

    let (rust_time, rust_sink) = time_loop(repeat, cases, |case| {
        checksum(&normalize_terms(rust_eval(case), mode))
    });
    let (c_time, c_sink) = time_loop(repeat, cases, |case| {
        checksum(&normalize_terms(c_eval(case), mode))
    });
    black_box((rust_sink, c_sink));

    println!("{name}: correctness ok ({} cases)", cases.len());
    println!(
        "  Rust:      {}  ({} evals)",
        format_duration(rust_time),
        repeat * cases.len()
    );
    println!(
        "  C lrcalc:  {}  ({} evals)",
        format_duration(c_time),
        repeat * cases.len()
    );
    println!(
        "  ratio Rust/C: {:.3}x",
        rust_time.as_secs_f64() / c_time.as_secs_f64()
    );

    if detail {
        print_case_diagnostics(name, repeat, cases, mode, rust_eval, c_eval);
    }
}

fn print_case_diagnostics<C, FRust, FC>(
    name: &str,
    repeat: usize,
    cases: &[C],
    mode: KeyMode,
    mut rust_eval: FRust,
    mut c_eval: FC,
) where
    C: BenchCase,
    FRust: FnMut(&C) -> Terms,
    FC: FnMut(&C) -> Terms,
{
    println!("  {name} per-case diagnostics:");
    println!(
        "    {:<24} {:>6} {:>7} {:>9} {:>9} {:>9} {:>8}",
        "case", "terms", "keylen", "abscoef", "Rust", "C", "Rust/C"
    );
    for case in cases {
        let terms = normalize_terms(rust_eval(case), mode);
        let term_count = terms.len();
        let max_key_len = terms.iter().map(|(key, _)| key.len()).max().unwrap_or(0);
        let abs_coeff_sum: i128 = terms.iter().map(|(_, coefficient)| coefficient.abs()).sum();
        let (rust_time, rust_sink) = time_one(repeat, case, |case| {
            checksum(&normalize_terms(rust_eval(case), mode))
        });
        let (c_time, c_sink) = time_one(repeat, case, |case| {
            checksum(&normalize_terms(c_eval(case), mode))
        });
        black_box((rust_sink, c_sink));
        println!(
            "    {:<24} {:>6} {:>7} {:>9} {:>9} {:>9} {:>8.3}x",
            truncate_label(case.label(), 24),
            term_count,
            max_key_len,
            abs_coeff_sum,
            format_duration(rust_time),
            format_duration(c_time),
            rust_time.as_secs_f64() / c_time.as_secs_f64()
        );
    }
}

fn time_loop<C, F>(repeat: usize, cases: &[C], mut f: F) -> (Duration, u64)
where
    F: FnMut(&C) -> u64,
{
    let start = Instant::now();
    let mut sink = 0u64;
    for _ in 0..repeat {
        for case in cases {
            sink ^= black_box(f(case));
        }
    }
    (start.elapsed(), sink)
}

fn time_one<C, F>(repeat: usize, case: &C, mut f: F) -> (Duration, u64)
where
    F: FnMut(&C) -> u64,
{
    let start = Instant::now();
    let mut sink = 0u64;
    for _ in 0..repeat {
        sink ^= black_box(f(case));
    }
    (start.elapsed(), sink)
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

unsafe fn collect_lc(upstream: &Upstream, lc: *mut CIvLinComb, mode: KeyMode) -> Terms {
    if lc.is_null() {
        panic!("upstream returned null linear combination");
    }
    let mut terms = Vec::new();
    let mut iter = CIvLcIter::default();
    unsafe {
        (upstream.ivlc_first)(lc, &mut iter);
        while (upstream.ivlc_good)(&mut iter) != 0 {
            let coefficient = (upstream.ivlc_value)(&mut iter);
            if coefficient != 0 {
                let mut key = vector_to_vec((upstream.ivlc_key)(&mut iter));
                if let KeyMode::TrimPartition = mode {
                    trim_trailing_zeroes(&mut key);
                }
                terms.push((key, i128::from(coefficient)));
            }
            (upstream.ivlc_next)(&mut iter);
        }
        (upstream.ivlc_free_all)(lc);
    }
    terms
}

unsafe fn vector_to_vec(vector: *mut CVector) -> Vec<i32> {
    let length = unsafe { (*vector).length as usize };
    let entries = unsafe { ptr::addr_of!((*vector).array).cast::<i32>() };
    (0..length)
        .map(|index| unsafe { *entries.add(index) })
        .collect()
}

fn schur_terms_to_terms(terms: Vec<lrcalc::schur::SchurTerm>) -> Terms {
    terms
        .into_iter()
        .map(|term| {
            (
                term.partition,
                i128::try_from(term.coefficient).expect("coefficient exceeds i128"),
            )
        })
        .collect()
}

fn signed_schur_terms_to_terms(terms: Vec<lrcalc::schur::SignedSchurTerm>) -> Terms {
    terms
        .into_iter()
        .map(|term| (term.partition, term.coefficient))
        .collect()
}

fn schubert_terms_to_terms(terms: lrcalc::schubert::LinearCombination) -> Terms {
    terms
        .into_iter()
        .map(|(key, coefficient)| (key, i128::from(coefficient)))
        .collect()
}

fn normalize_terms(mut terms: Terms, mode: KeyMode) -> Terms {
    if let KeyMode::TrimPartition = mode {
        for (key, _) in &mut terms {
            trim_trailing_zeroes(key);
        }
    }
    terms.retain(|(_, coefficient)| *coefficient != 0);
    terms.sort();
    terms
}

fn trim_trailing_zeroes(values: &mut Vec<i32>) {
    while values.last() == Some(&0) {
        values.pop();
    }
}

fn checksum(terms: &Terms) -> u64 {
    let mut hash = terms.len() as u64;
    for (key, coefficient) in terms {
        hash = hash.wrapping_mul(1_099_511_628_211);
        hash ^= *coefficient as u64;
        for &entry in key {
            hash = hash.wrapping_mul(16_777_619) ^ entry as u64;
        }
    }
    hash
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

fn rtld_deepbind() -> c_int {
    #[cfg(target_os = "linux")]
    {
        0x00008
    }
    #[cfg(not(target_os = "linux"))]
    {
        0
    }
}

fn product_cases(upstream: &Upstream) -> Vec<ProductCase> {
    vec![
        product_case(upstream, "single boxes", &[1], &[1], -1, -1),
        product_case(upstream, "two hooks", &[2, 1], &[2, 1], -1, -1),
        product_case(upstream, "row bounded", &[2, 1], &[2, 1], 3, -1),
        product_case(upstream, "column bounded", &[2, 1], &[2, 1], -1, 2),
        product_case(upstream, "balanced medium", &[4, 2], &[3, 1], -1, -1),
        product_case(upstream, "dense medium", &[5, 3, 1], &[4, 2], -1, -1),
        product_case(upstream, "dense larger", &[6, 4, 2], &[5, 3, 1], -1, -1),
    ]
}

fn skew_cases(upstream: &Upstream) -> Vec<SkewCase> {
    vec![
        skew_case(upstream, "two cells", &[2, 1], &[1], -1),
        skew_case(upstream, "larger example", &[3, 2, 1], &[2, 1], -1),
        skew_case(upstream, "row bounded", &[3, 2, 1], &[2, 1], 2),
        skew_case(upstream, "balanced skew", &[6, 4, 2], &[3, 2, 1], -1),
        skew_case(
            upstream,
            "optimized A",
            &[20, 18, 16, 14, 12],
            &[10, 8, 6, 4, 2],
            -1,
        ),
        skew_case(upstream, "optimized B", &[24, 22, 20, 18], &[12, 10, 8], -1),
    ]
}

fn coprod_cases(upstream: &Upstream) -> Vec<CoprodCase> {
    vec![
        coprod_case(upstream, "single box", &[1], 1, 1, false),
        coprod_case(upstream, "single box all", &[1], 1, 1, true),
        coprod_case(upstream, "two rows", &[2, 1], 1, 2, false),
        coprod_case(upstream, "rectangle filter", &[3, 2], 2, 2, false),
        coprod_case(upstream, "medium all", &[4, 2, 1], 2, 3, true),
    ]
}

fn fusion_cases(upstream: &Upstream) -> Vec<FusionCase> {
    vec![
        fusion_case(upstream, "small level", &[1], &[1], 2, 1),
        fusion_case(upstream, "example 321", &[2, 1], &[2, 1], 3, 2),
        fusion_case(upstream, "example 42", &[3, 1], &[2, 1], 3, 2),
        fusion_case(upstream, "level three", &[3, 2], &[2, 1], 3, 3),
        fusion_case(upstream, "four rows", &[3, 2, 1], &[2, 2], 4, 3),
    ]
}

fn schubert_cases(upstream: &Upstream) -> Vec<SchubertCase> {
    vec![
        schubert_case(upstream, "s3 identity/id", &[1, 2, 3], &[1, 2, 3], 0),
        schubert_case(upstream, "s3 simple/id", &[2, 1, 3], &[1, 2, 3], 0),
        schubert_case(upstream, "s1 squared", &[2, 1], &[2, 1], 0),
        schubert_case(upstream, "rank trunc", &[2, 1], &[2, 1], 2),
        schubert_case(upstream, "rank three simple", &[2, 1], &[2, 1], 3),
        schubert_case(upstream, "s3 adjacent", &[2, 1, 3], &[1, 3, 2], 0),
        schubert_case(upstream, "s3 pair", &[3, 1, 2], &[2, 3, 1], 0),
        schubert_case(upstream, "s3 longest/simple", &[3, 2, 1], &[2, 1, 3], 0),
        schubert_case(upstream, "s4 simple/id", &[1, 3, 2, 4], &[1, 2, 3, 4], 0),
        schubert_case(upstream, "s4 pair", &[4, 2, 3, 1], &[3, 4, 1, 2], 0),
        schubert_case(upstream, "s4 grassmannian", &[2, 4, 1, 3], &[3, 1, 4, 2], 0),
        schubert_case(upstream, "s4 bounded", &[4, 2, 3, 1], &[3, 4, 1, 2], 4),
        schubert_case(upstream, "s4 long/simple", &[4, 3, 2, 1], &[2, 1, 4, 3], 0),
        schubert_case(upstream, "s5 pair", &[3, 5, 1, 4, 2], &[2, 4, 5, 1, 3], 0),
        schubert_case(upstream, "s5 mixed", &[5, 1, 3, 2, 4], &[2, 5, 4, 1, 3], 0),
        schubert_case(
            upstream,
            "s6 sparse",
            &[2, 1, 4, 3, 6, 5],
            &[3, 2, 1, 6, 5, 4],
            0,
        ),
        schubert_case(
            upstream,
            "s6 medium",
            &[4, 1, 6, 2, 5, 3],
            &[3, 6, 1, 5, 2, 4],
            0,
        ),
        schubert_case(
            upstream,
            "s7 medium",
            &[3, 7, 1, 5, 2, 6, 4],
            &[4, 1, 6, 2, 7, 3, 5],
            0,
        ),
    ]
}

fn schubert_string_cases(upstream: &Upstream) -> Vec<SchubertStringCase> {
    vec![
        schubert_string_case(upstream, "binary two", &[0, 1], &[1, 0]),
        schubert_string_case(upstream, "binary four", &[0, 1, 0, 1], &[1, 0, 1, 0]),
        schubert_string_case(
            upstream,
            "binary six",
            &[0, 1, 0, 1, 0, 1],
            &[1, 0, 1, 0, 1, 0],
        ),
        schubert_string_case(
            upstream,
            "binary eight",
            &[0, 1, 0, 1, 0, 1, 0, 1],
            &[1, 0, 1, 0, 1, 0, 1, 0],
        ),
        schubert_string_case(
            upstream,
            "ternary six",
            &[0, 1, 2, 0, 1, 2],
            &[2, 1, 0, 2, 1, 0],
        ),
        schubert_string_case(
            upstream,
            "ternary mixed",
            &[0, 2, 1, 0, 2, 1],
            &[1, 0, 2, 1, 0, 2],
        ),
        schubert_string_case(
            upstream,
            "ternary nine",
            &[0, 1, 2, 0, 1, 2, 0, 1, 2],
            &[2, 1, 0, 2, 1, 0, 2, 1, 0],
        ),
        schubert_string_case(
            upstream,
            "quaternary eight",
            &[0, 1, 2, 3, 0, 1, 2, 3],
            &[3, 2, 1, 0, 3, 2, 1, 0],
        ),
    ]
}

fn product_case(
    upstream: &Upstream,
    label: &'static str,
    sh1: &[i32],
    sh2: &[i32],
    rows: i32,
    cols: i32,
) -> ProductCase {
    ProductCase {
        label,
        sh1: sh1.to_vec(),
        sh2: sh2.to_vec(),
        rows,
        cols,
        c_sh1: OwnedCVector::new(upstream, sh1),
        c_sh2: OwnedCVector::new(upstream, sh2),
    }
}

fn skew_case(
    upstream: &Upstream,
    label: &'static str,
    outer: &[i32],
    inner: &[i32],
    rows: i32,
) -> SkewCase {
    SkewCase {
        label,
        outer: outer.to_vec(),
        inner: inner.to_vec(),
        rows,
        c_outer: OwnedCVector::new(upstream, outer),
        c_inner: OwnedCVector::new(upstream, inner),
    }
}

fn coprod_case(
    upstream: &Upstream,
    label: &'static str,
    shape: &[i32],
    rows: i32,
    cols: i32,
    all: bool,
) -> CoprodCase {
    CoprodCase {
        label,
        shape: shape.to_vec(),
        rows,
        cols,
        all,
        c_shape: OwnedCVector::new(upstream, shape),
    }
}

fn fusion_case(
    upstream: &Upstream,
    label: &'static str,
    sh1: &[i32],
    sh2: &[i32],
    rows: i32,
    level: i32,
) -> FusionCase {
    FusionCase {
        label,
        sh1: sh1.to_vec(),
        sh2: sh2.to_vec(),
        rows,
        level,
        c_sh1: OwnedCVector::new(upstream, sh1),
        c_sh2: OwnedCVector::new(upstream, sh2),
    }
}

fn schubert_case(
    upstream: &Upstream,
    label: &'static str,
    left: &[i32],
    right: &[i32],
    rank: i32,
) -> SchubertCase {
    SchubertCase {
        label,
        left: left.to_vec(),
        right: right.to_vec(),
        rank,
        c_left: OwnedCVector::new(upstream, left),
        c_right: OwnedCVector::new(upstream, right),
    }
}

fn schubert_string_case(
    upstream: &Upstream,
    label: &'static str,
    left: &[i32],
    right: &[i32],
) -> SchubertStringCase {
    SchubertStringCase {
        label,
        left: left.to_vec(),
        right: right.to_vec(),
        c_left: OwnedCVector::new(upstream, left),
        c_right: OwnedCVector::new(upstream, right),
    }
}
