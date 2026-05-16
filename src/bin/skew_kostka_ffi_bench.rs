fn main() {
    enabled::main();
}

mod enabled {
    use lrcalc::kostka_fast::skew_kostka_fast_u128;
    use std::ffi::{CStr, CString};
    use std::hint::black_box;
    use std::os::raw::c_int;
    use std::ptr;
    use std::time::{Duration, Instant};

    const LC_COPY_KEY: c_int = 1;
    const LC_FREE_ZERO: c_int = 2;
    type RawCase = (&'static str, Vec<i32>, Vec<i32>, Vec<i32>);

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
    type IvHash = unsafe extern "C" fn(*mut CVector) -> u32;
    type IvLcNew = unsafe extern "C" fn(u32, u32) -> *mut CIvLinComb;
    type IvLcFreeAll = unsafe extern "C" fn(*mut CIvLinComb);
    type IvLcFirst = unsafe extern "C" fn(*mut CIvLinComb, *mut CIvLcIter);
    type IvLcGood = unsafe extern "C" fn(*mut CIvLcIter) -> c_int;
    type IvLcNext = unsafe extern "C" fn(*mut CIvLcIter);
    type IvLcKey = unsafe extern "C" fn(*mut CIvLcIter) -> *mut CVector;
    type IvLcValue = unsafe extern "C" fn(*mut CIvLcIter) -> i32;
    type IvLcAddElement =
        unsafe extern "C" fn(*mut CIvLinComb, i32, *mut CVector, u32, c_int) -> c_int;
    type SchurMult =
        unsafe extern "C" fn(*mut CVector, *mut CVector, c_int, c_int, c_int) -> *mut CIvLinComb;

    struct Upstream {
        handle: *mut libc::c_void,
        iv_new: IvNew,
        iv_free: IvFree,
        iv_hash: IvHash,
        ivlc_new: IvLcNew,
        ivlc_free_all: IvLcFreeAll,
        ivlc_first: IvLcFirst,
        ivlc_good: IvLcGood,
        ivlc_next: IvLcNext,
        ivlc_key: IvLcKey,
        ivlc_value: IvLcValue,
        ivlc_add_element: IvLcAddElement,
        schur_mult: SchurMult,
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
                    iv_hash: symbol(handle, b"iv_hash\0"),
                    ivlc_new: symbol(handle, b"ivlc_new\0"),
                    ivlc_free_all: symbol(handle, b"ivlc_free_all\0"),
                    ivlc_first: symbol(handle, b"ivlc_first\0"),
                    ivlc_good: symbol(handle, b"ivlc_good\0"),
                    ivlc_next: symbol(handle, b"ivlc_next\0"),
                    ivlc_key: symbol(handle, b"ivlc_key\0"),
                    ivlc_value: symbol(handle, b"ivlc_value\0"),
                    ivlc_add_element: symbol(handle, b"ivlc_add_element\0"),
                    schur_mult: symbol(handle, b"schur_mult\0"),
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

    unsafe fn symbol<T: Copy>(handle: *mut libc::c_void, name: &[u8]) -> T {
        let ptr = libc::dlsym(handle, name.as_ptr().cast());
        if ptr.is_null() {
            let name = CStr::from_bytes_with_nul(name)
                .expect("symbol name should be NUL-terminated")
                .to_string_lossy();
            panic!("dlsym failed for {name}: {}", dl_error());
        }
        std::mem::transmute_copy(&ptr)
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

    struct Case {
        label: &'static str,
        outer: Vec<i32>,
        inner: Vec<i32>,
        weight: Vec<i32>,
        c_inner: OwnedCVector,
        c_rows: Vec<OwnedCVector>,
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

    struct OwnedLinearCombination {
        ptr: *mut CIvLinComb,
        ivlc_free_all: IvLcFreeAll,
    }

    impl OwnedLinearCombination {
        fn new(upstream: &Upstream) -> Self {
            let ptr = unsafe { (upstream.ivlc_new)(2003, 100) };
            if ptr.is_null() {
                panic!("upstream ivlc_new returned null");
            }
            Self {
                ptr,
                ivlc_free_all: upstream.ivlc_free_all,
            }
        }

        fn from_partition(upstream: &Upstream, partition: &OwnedCVector) -> Self {
            let result = Self::new(upstream);
            let status = unsafe {
                (upstream.ivlc_add_element)(
                    result.ptr,
                    1,
                    partition.as_mut_ptr(),
                    (upstream.iv_hash)(partition.as_mut_ptr()),
                    LC_COPY_KEY,
                )
            };
            if status != 0 {
                panic!("upstream ivlc_add_element failed");
            }
            result
        }

        fn as_mut_ptr(&self) -> *mut CIvLinComb {
            self.ptr
        }
    }

    impl Drop for OwnedLinearCombination {
        fn drop(&mut self) {
            if !self.ptr.is_null() {
                unsafe { (self.ivlc_free_all)(self.ptr) };
            }
        }
    }

    pub fn main() {
        let repeat = std::env::args()
            .nth(1)
            .map(|arg| arg.parse::<usize>().expect("repeat must be an integer"))
            .unwrap_or(2_000);
        let upstream = Upstream::load();
        let cases = build_cases(&upstream);

        for case in &cases {
            let fast =
                skew_kostka_fast_u128(&case.outer, &case.inner, &case.weight).expect(case.label);
            let c_value = c_skew_kostka(&upstream, case);
            if fast != c_value as u128 {
                eprintln!("mismatch: {}", case.label);
                eprintln!("outer:  {:?}", case.outer);
                eprintln!("inner:  {:?}", case.inner);
                eprintln!("weight: {:?}", case.weight);
                eprintln!("fast: {fast}");
                eprintln!("C:    {c_value}");
                std::process::exit(1);
            }
        }

        println!("correctness: ok ({} cases)", cases.len());
        println!("repeat: {repeat}");

        let (fast_time, fast_sink) = time_loop(repeat, &cases, |case| {
            skew_kostka_fast_u128(&case.outer, &case.inner, &case.weight).expect(case.label)
        });
        let (c_time, c_sink) = time_loop(repeat, &cases, |case| {
            c_skew_kostka(&upstream, case) as u128
        });
        black_box((fast_sink, c_sink));

        println!(
            "fast skew DP:        {}  ({} evals)",
            format_duration(fast_time),
            repeat * cases.len()
        );
        println!(
            "C repeated Schur:    {}  ({} evals)",
            format_duration(c_time),
            repeat * cases.len()
        );
        println!(
            "comparison: fast/C-lrcalc = {:.3}x",
            fast_time.as_secs_f64() / c_time.as_secs_f64()
        );
    }

    fn c_skew_kostka(upstream: &Upstream, case: &Case) -> i32 {
        let mut current = OwnedLinearCombination::from_partition(upstream, &case.c_inner);
        for row in &case.c_rows {
            current = multiply_by_row(upstream, current, row);
        }

        let mut value = 0i32;
        let mut iter = CIvLcIter::default();
        unsafe {
            (upstream.ivlc_first)(current.as_mut_ptr(), &mut iter);
            while (upstream.ivlc_good)(&mut iter) != 0 {
                if vector_matches_partition((upstream.ivlc_key)(&mut iter), &case.outer) {
                    value = value
                        .checked_add((upstream.ivlc_value)(&mut iter))
                        .expect("upstream i32 coefficient overflow");
                }
                (upstream.ivlc_next)(&mut iter);
            }
        }
        value
    }

    fn multiply_by_row(
        upstream: &Upstream,
        current: OwnedLinearCombination,
        row: &OwnedCVector,
    ) -> OwnedLinearCombination {
        let result = OwnedLinearCombination::new(upstream);
        let mut iter = CIvLcIter::default();
        unsafe {
            (upstream.ivlc_first)(current.as_mut_ptr(), &mut iter);
            while (upstream.ivlc_good)(&mut iter) != 0 {
                let source_coef = (upstream.ivlc_value)(&mut iter);
                let product = (upstream.schur_mult)(
                    (upstream.ivlc_key)(&mut iter),
                    row.as_mut_ptr(),
                    -1,
                    -1,
                    -1,
                );
                if product.is_null() {
                    panic!("upstream schur_mult returned null");
                }
                add_scaled_product(upstream, result.as_mut_ptr(), product, source_coef);
                (upstream.ivlc_free_all)(product);
                (upstream.ivlc_next)(&mut iter);
            }
        }
        drop(current);
        result
    }

    fn add_scaled_product(
        upstream: &Upstream,
        dst: *mut CIvLinComb,
        product: *mut CIvLinComb,
        scale: i32,
    ) {
        let mut iter = CIvLcIter::default();
        unsafe {
            (upstream.ivlc_first)(product, &mut iter);
            while (upstream.ivlc_good)(&mut iter) != 0 {
                let coef = scale
                    .checked_mul((upstream.ivlc_value)(&mut iter))
                    .expect("upstream i32 coefficient overflow");
                let key = (upstream.ivlc_key)(&mut iter);
                let status = (upstream.ivlc_add_element)(
                    dst,
                    coef,
                    key,
                    (upstream.iv_hash)(key),
                    LC_COPY_KEY | LC_FREE_ZERO,
                );
                if status != 0 {
                    panic!("upstream ivlc_add_element failed");
                }
                (upstream.ivlc_next)(&mut iter);
            }
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

    fn build_cases(upstream: &Upstream) -> Vec<Case> {
        raw_cases()
            .into_iter()
            .map(|(label, outer, inner, weight)| {
                let rows = weight
                    .iter()
                    .copied()
                    .filter(|&part| part > 0)
                    .map(|part| OwnedCVector::new(upstream, &[part]))
                    .collect::<Vec<_>>();
                Case {
                    label,
                    c_inner: OwnedCVector::new(upstream, &inner),
                    c_rows: rows,
                    outer,
                    inner,
                    weight,
                }
            })
            .collect()
    }

    unsafe fn vector_matches_partition(vector: *mut CVector, partition: &[i32]) -> bool {
        let length = (*vector).length as usize;
        let entries = ptr::addr_of!((*vector).array).cast::<i32>();
        let max_len = length.max(partition.len());
        for index in 0..max_len {
            let left = if index < length {
                *entries.add(index)
            } else {
                0
            };
            let right = partition.get(index).copied().unwrap_or(0);
            if left != right {
                return false;
            }
        }
        true
    }

    fn raw_cases() -> Vec<RawCase> {
        vec![
            case("horizontal strip", &[5, 3, 1], &[3, 2, 1], &[2, 1]),
            case("small skew", &[4, 2, 1], &[2, 1], &[2, 1, 1]),
            case("two row skew", &[5, 4, 2], &[3, 1], &[3, 3, 1]),
            case("balanced skew", &[6, 4, 2], &[3, 2, 1], &[3, 2, 1, 1]),
            case("stretched skew x2", &[8, 6, 4], &[4, 2], &[4, 4, 4]),
            case(
                "standard skew",
                &[5, 3, 2, 1],
                &[2, 1],
                &[1, 1, 1, 1, 1, 1, 1, 1],
            ),
            case("irregular A", &[7, 5, 3, 1], &[4, 2, 1], &[4, 3, 2]),
            case("irregular B", &[8, 6, 4, 2], &[5, 3, 1], &[5, 4, 3, 1]),
            case("irregular C", &[9, 7, 4, 2], &[5, 4, 2], &[5, 4, 3, 1]),
            case(
                "many labels",
                &[7, 5, 3, 2],
                &[3, 2],
                &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            ),
        ]
    }

    fn case(label: &'static str, outer: &[i32], inner: &[i32], weight: &[i32]) -> RawCase {
        (label, outer.to_vec(), inner.to_vec(), weight.to_vec())
    }
}
