#![allow(unexpected_cfgs)]

#[cfg(not(upstream_lrcalc_ffi))]
fn main() {
    eprintln!(
        "kostka_ffi_bench requires cfg upstream_lrcalc_ffi and a link path to upstream liblrcalc"
    );
    eprintln!(
        "use scripts/kostka_ffi_bench.sh, or set RUSTFLAGS='--cfg upstream_lrcalc_ffi -L native=/tmp/lrcalc-upstream/src/.libs'"
    );
    std::process::exit(2);
}

#[cfg(upstream_lrcalc_ffi)]
fn main() {
    enabled::main();
}

#[cfg(upstream_lrcalc_ffi)]
mod enabled {
    use lrcalc::kostka::kostka_lr_triple;
    use lrcalc::kostka_fast::kostka_fast_u128;
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
        shape: Vec<i32>,
        weight: Vec<i32>,
        c_outer: OwnedCVector,
        c_inner: OwnedCVector,
        c_content: OwnedCVector,
    }

    struct OwnedCVector {
        ptr: *mut CVector,
    }

    impl OwnedCVector {
        fn new(values: &[i32]) -> Self {
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
        let repeat = std::env::args()
            .nth(1)
            .map(|arg| arg.parse::<usize>().expect("repeat must be an integer"))
            .unwrap_or(5_000);
        let cases = build_cases();

        for case in &cases {
            let fast = kostka_fast_u128(&case.shape, &case.weight).expect(case.label);
            let c_value = c_lrcoef(case);
            if c_value < 0 || fast != c_value as u128 {
                eprintln!("mismatch: {}", case.label);
                eprintln!("shape:  {:?}", case.shape);
                eprintln!("weight: {:?}", case.weight);
                eprintln!("fast: {fast}");
                eprintln!("C:    {c_value}");
                std::process::exit(1);
            }
        }

        println!("correctness: ok ({} cases)", cases.len());
        println!("repeat: {repeat}");

        let (fast_time, fast_sink) = time_loop(repeat, &cases, |case| {
            kostka_fast_u128(&case.shape, &case.weight).expect(case.label)
        });
        let (c_time, c_sink) = time_loop(repeat, &cases, |case| {
            let value = c_lrcoef(case);
            if value < 0 {
                panic!("C schur_lrcoef returned {value}");
            }
            value as u128
        });
        black_box((fast_sink, c_sink));

        println!(
            "fast u128 DP: {}  ({} evals)",
            format_duration(fast_time),
            repeat * cases.len()
        );
        println!(
            "C lrcalc FFI: {}  ({} evals)",
            format_duration(c_time),
            repeat * cases.len()
        );
        println!(
            "comparison: fast/C-lrcalc = {:.3}x",
            fast_time.as_secs_f64() / c_time.as_secs_f64()
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
        format!("{:.3}s", duration.as_secs_f64())
    }

    fn build_cases() -> Vec<Case> {
        raw_cases()
            .into_iter()
            .map(|(label, shape, weight)| {
                let (outer, inner, content) =
                    kostka_lr_triple(&shape, &weight).expect("Kostka-to-LR triple");
                Case {
                    label,
                    shape,
                    weight,
                    c_outer: OwnedCVector::new(&outer),
                    c_inner: OwnedCVector::new(&inner),
                    c_content: OwnedCVector::new(&content),
                }
            })
            .collect()
    }

    fn raw_cases() -> Vec<(&'static str, Vec<i32>, Vec<i32>)> {
        vec![
            case("row shape small", &[3], &[2, 1]),
            case("standard 21", &[2, 1], &[1, 1, 1]),
            case("dominance zero", &[1, 1, 1], &[2, 1]),
            case("single row all ones", &[6], &[1, 1, 1, 1, 1, 1]),
            case("single column heavy zero", &[1, 1, 1, 1], &[3, 1]),
            case("balanced 321", &[3, 2, 1], &[2, 2, 2]),
            case("example 421", &[4, 2, 1], &[3, 2, 1, 1]),
            case("standard 321", &[3, 2, 1], &[1, 1, 1, 1, 1, 1]),
            case(
                "standard 4321",
                &[4, 3, 2, 1],
                &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            ),
            case("stretched 321 x2", &[6, 4, 2], &[4, 4, 4]),
            case("stretched 321 x3", &[9, 6, 3], &[6, 6, 6]),
            case("stretched 421 x2", &[8, 4, 2], &[6, 4, 2, 2]),
            case("stretched 421 x3", &[12, 6, 3], &[9, 6, 3]),
            case("irregular A", &[5, 4, 2, 1], &[4, 3, 2, 2, 1]),
            case("irregular B", &[6, 4, 3, 1], &[5, 3, 3, 2, 1]),
            case("irregular C", &[7, 5, 3, 2], &[6, 4, 3, 2, 2]),
            case("irregular D", &[8, 6, 4, 2], &[6, 5, 4, 3, 2]),
            case("irregular E", &[9, 7, 4, 2, 1], &[7, 5, 4, 3, 2, 2]),
            case("many labels 42", &[4, 2], &[1, 1, 1, 1, 1, 1]),
            case(
                "many labels 5321",
                &[5, 3, 2, 1],
                &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            ),
            case(
                "many labels 7431",
                &[7, 4, 3, 1],
                &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            ),
        ]
    }

    fn case(
        label: &'static str,
        shape: &[i32],
        weight: &[i32],
    ) -> (&'static str, Vec<i32>, Vec<i32>) {
        (label, shape.to_vec(), weight.to_vec())
    }
}
