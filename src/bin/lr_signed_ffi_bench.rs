#![allow(unexpected_cfgs)]

#[cfg(not(upstream_lrcalc_ffi))]
fn main() {
    eprintln!(
        "lr_signed_ffi_bench requires cfg upstream_lrcalc_ffi and a link path to upstream liblrcalc"
    );
    eprintln!("use scripts/lr_signed_ffi_bench.sh");
    std::process::exit(2);
}

#[cfg(upstream_lrcalc_ffi)]
fn main() {
    enabled::main();
}

#[cfg(upstream_lrcalc_ffi)]
mod enabled {
    use lrcalc::lr_signed::lrcoef_signed_kostka;
    use std::hint::black_box;
    use std::ptr;
    use std::time::{Duration, Instant};

    #[repr(C)]
    struct CVector {
        length: u32,
        array: [i32; 1],
    }

    #[link(name = "lrcalc", kind = "static")]
    extern "C" {
        fn iv_new(length: u32) -> *mut CVector;
        fn iv_free(v: *mut CVector);
        fn schur_lrcoef(outer: *mut CVector, inner1: *mut CVector, inner2: *mut CVector) -> i64;
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

    struct OwnedCVector {
        ptr: *mut CVector,
    }

    impl OwnedCVector {
        fn new(values: &[i32]) -> Self {
            let values = if values.is_empty() { &[0][..] } else { values };
            let length = u32::try_from(values.len()).expect("input vector too long for upstream");
            let ptr = unsafe { iv_new(length) };
            if ptr.is_null() {
                panic!("upstream iv_new returned null");
            }
            for (index, &value) in values.iter().enumerate() {
                unsafe {
                    *ptr::addr_of_mut!((*ptr).array).cast::<i32>().add(index) = value;
                }
            }
            Self { ptr }
        }

        fn as_mut_ptr(&self) -> *mut CVector {
            self.ptr
        }
    }

    impl Drop for OwnedCVector {
        fn drop(&mut self) {
            if !self.ptr.is_null() {
                unsafe { iv_free(self.ptr) };
            }
        }
    }

    pub fn main() {
        let mut repeat = 5_000usize;
        let mut suite = "mixed".to_string();
        for arg in std::env::args().skip(1) {
            if let Ok(value) = arg.parse::<usize>() {
                repeat = value;
            } else {
                suite = arg;
            }
        }
        let cases = build_cases(&suite);

        for case in &cases {
            let signed = lrcoef_signed_kostka(&case.outer, &case.inner, &case.content)
                .unwrap_or_else(|_| panic!("signed LR failed for {}", case.label));
            let c_value = c_lrcoef(case);
            if c_value < 0 || signed != c_value as u128 {
                eprintln!("mismatch: {}", case.label);
                eprintln!("outer:   {:?}", case.outer);
                eprintln!("inner:   {:?}", case.inner);
                eprintln!("content: {:?}", case.content);
                eprintln!("signed: {signed}");
                eprintln!("C:      {c_value}");
                std::process::exit(1);
            }
        }

        println!("correctness: ok ({} cases)", cases.len());
        println!("suite: {suite}");
        println!("repeat: {repeat}");

        let (signed_time, signed_sink) = time_loop(repeat, &cases, |case| {
            lrcoef_signed_kostka(&case.outer, &case.inner, &case.content)
                .unwrap_or_else(|_| panic!("signed LR failed for {}", case.label))
        });
        let (c_time, c_sink) = time_loop(repeat, &cases, |case| {
            let value = c_lrcoef(case);
            if value < 0 {
                panic!("C schur_lrcoef returned {value}");
            }
            value as u128
        });
        black_box((signed_sink, c_sink));

        println!(
            "signed Kostka LR: {}  ({} evals)",
            format_duration(signed_time),
            repeat * cases.len()
        );
        println!(
            "C lrcalc FFI:    {}  ({} evals)",
            format_duration(c_time),
            repeat * cases.len()
        );
        println!(
            "comparison: signed/C-lrcalc = {:.3}x",
            signed_time.as_secs_f64() / c_time.as_secs_f64()
        );
    }

    fn c_lrcoef(case: &Case) -> i64 {
        unsafe {
            schur_lrcoef(
                case.c_outer.as_mut_ptr(),
                case.c_inner.as_mut_ptr(),
                case.c_content.as_mut_ptr(),
            )
        }
    }

    fn time_loop<F>(repeat: usize, cases: &[Case], mut f: F) -> (Duration, u128)
    where
        F: FnMut(&Case) -> u128,
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

    fn format_duration(duration: Duration) -> String {
        format!("{:.6}s", duration.as_secs_f64())
    }

    fn build_cases(suite: &str) -> Vec<Case> {
        raw_cases(suite)
            .into_iter()
            .map(|(label, outer, inner, content)| Case {
                label,
                c_outer: OwnedCVector::new(&outer),
                c_inner: OwnedCVector::new(&inner),
                c_content: OwnedCVector::new(&content),
                outer,
                inner,
                content,
            })
            .collect()
    }

    fn raw_cases(suite: &str) -> Vec<(&'static str, Vec<i32>, Vec<i32>, Vec<i32>)> {
        match suite {
            "mixed" => raw_mixed_cases(),
            "large-few-parts" | "large-short" | "large" => raw_large_short_cases(),
            "large-three-part" | "large-3part" => raw_large_three_part_cases(),
            other => {
                panic!(
                    "unknown suite '{other}', expected mixed, large-few-parts, or large-three-part"
                )
            }
        }
    }

    fn raw_mixed_cases() -> Vec<(&'static str, Vec<i32>, Vec<i32>, Vec<i32>)> {
        vec![
            case("empty coefficient", &[0], &[0], &[0]),
            case("empty outer nonempty content", &[0], &[0], &[1]),
            case("size mismatch zero", &[1], &[0], &[0]),
            case("inner not contained zero", &[2, 2, 1], &[3], &[2]),
            case("single box Pieri", &[2, 1], &[2], &[1]),
            case("s21 s21 coefficient 2", &[3, 2, 1], &[2, 1], &[2, 1]),
            case("s21 s21 coefficient 1 A", &[4, 2], &[2, 1], &[2, 1]),
            case("s21 s21 coefficient 1 B", &[3, 3], &[2, 1], &[2, 1]),
            case("s21 s21 absent shape", &[5, 1], &[2, 1], &[2, 1]),
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
            case(
                "upstream oom medium",
                &[5, 4, 3, 2, 1],
                &[3, 2, 1],
                &[4, 3, 1, 1],
            ),
            case(
                "upstream oom larger",
                &[7, 6, 5, 4, 3, 2, 1],
                &[4, 4, 3, 2, 1],
                &[5, 4, 3, 2],
            ),
            case("stretched 321 scale 2", &[6, 4, 2], &[4, 2], &[4, 2]),
            case("stretched 321 scale 3", &[9, 6, 3], &[6, 3], &[6, 3]),
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
                "3-part stretched scale 40",
                &[120, 80, 40],
                &[60, 40, 20],
                &[60, 40, 20],
            ),
            case(
                "3-part stretched scale 60",
                &[180, 120, 60],
                &[90, 60, 30],
                &[90, 60, 30],
            ),
            case(
                "3-part stretched scale 80",
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
}
