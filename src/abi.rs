//! C ABI types and functions matching the original `lrcalc` headers.

use crate::lrcoef::{lrcoef_i64, optim_coef as native_optim_coef, OptimizedCoef};
use crate::schubert::{
    monk_product, multiply_poly_schubert, multiply_schubert, multiply_schubert_strings,
    trans_polynomial, LinearCombination,
};
use crate::schur::{
    fusion_reduce_values, schur_coproduct_expansion, schur_product_expansion,
    schur_product_fusion_expansion, visit_schur_skew_expansion_with_len, SchurTerm,
    SignedSchurTerm,
};
use libc::{c_char, c_int, c_longlong, c_void};
use std::cell::Cell;
use std::ffi::CStr;
use std::mem;
use std::ptr;
use std::slice;

unsafe extern "C" {
    static mut optind: c_int;
}

#[repr(C)]
pub struct IVector {
    pub length: u32,
    pub array: [i32; 1],
}

#[repr(C)]
pub struct IvlcKeyVal {
    pub key: *mut IVector,
    pub value: i32,
    pub hash: u32,
    pub next: u32,
}

#[repr(C)]
pub struct IvLinComb {
    pub table: *mut u32,
    pub elts: *mut IvlcKeyVal,
    pub card: u32,
    pub free_elts: u32,
    pub elts_len: u32,
    pub elts_sz: u32,
    pub table_sz: u32,
}

#[repr(C)]
pub struct IList {
    pub array: *mut c_int,
    pub allocated: usize,
    pub length: usize,
}

#[repr(C)]
pub struct IvList {
    pub array: *mut *mut IVector,
    pub allocated: usize,
    pub length: usize,
}

#[repr(C)]
pub struct PartIter {
    pub part: *mut IVector,
    pub outer: *mut IVector,
    pub inner: *mut IVector,
    pub length: c_int,
    pub rows: c_int,
    pub opt: c_int,
}

#[repr(C)]
pub struct SkewShapeAbi {
    pub outer: *mut IVector,
    pub inner: *mut IVector,
    pub cont: *mut IVector,
    pub sign: c_int,
}

#[repr(C)]
pub struct IvlcIter {
    pub ht: *mut IvLinComb,
    pub index: u32,
    pub i: u32,
}

const IVLC_HASHTABLE_SZ: u32 = 2003;
const IVLC_ARRAY_SZ: u32 = 100;
const IVLC_USE_FACTOR: u32 = 2;
const LC_COPY_KEY: c_int = 1;
const LC_FREE_ZERO: c_int = 2;
const PITR_USE_OUTER: c_int = 1;
const PITR_USE_INNER: c_int = 2;
const PITR_USE_SIZE: c_int = 4;

#[repr(C)]
pub struct LritBox {
    pub value: c_int,
    pub max: c_int,
    pub above: c_int,
    pub right: c_int,
}

#[repr(C)]
pub struct LrTabIter {
    pub cont: *mut IVector,
    pub size: c_int,
    pub array_len: c_int,
    pub array: [LritBox; 1],
}

fn lrtab_iter_alloc_size(array_len: usize) -> Option<usize> {
    let entries = array_len.max(1);
    mem::size_of::<LrTabIter>().checked_add(
        entries
            .checked_sub(1)?
            .checked_mul(mem::size_of::<LritBox>())?,
    )
}

unsafe fn lrit_array_mut<'a>(lrit: *mut LrTabIter) -> &'a mut [LritBox] {
    let length = unsafe { (*lrit).array_len as usize };
    let data = unsafe { ptr::addr_of_mut!((*lrit).array).cast::<LritBox>() };
    unsafe { slice::from_raw_parts_mut(data, length) }
}

fn ivector_alloc_size(length: u32) -> Option<usize> {
    let entries = usize::try_from(length).ok()?.max(1);
    mem::size_of::<u32>().checked_add(entries.checked_mul(mem::size_of::<i32>())?)
}

unsafe fn ivector_values_mut<'a>(v: *mut IVector) -> &'a mut [i32] {
    let length = unsafe { (*v).length as usize };
    let data = unsafe { ptr::addr_of_mut!((*v).array).cast::<i32>() };
    unsafe { slice::from_raw_parts_mut(data, length) }
}

unsafe fn ivector_values_mut_len<'a>(v: *mut IVector, length: usize) -> &'a mut [i32] {
    let data = unsafe { ptr::addr_of_mut!((*v).array).cast::<i32>() };
    unsafe { slice::from_raw_parts_mut(data, length) }
}

unsafe fn ivector_values<'a>(v: *const IVector) -> &'a [i32] {
    let length = unsafe { (*v).length as usize };
    let data = unsafe { ptr::addr_of!((*v).array).cast::<i32>() };
    unsafe { slice::from_raw_parts(data, length) }
}

unsafe fn ivector_data(v: *const IVector) -> *const i32 {
    unsafe { ptr::addr_of!((*v).array).cast::<i32>() }
}

unsafe fn ivector_data_mut(v: *mut IVector) -> *mut i32 {
    unsafe { ptr::addr_of_mut!((*v).array).cast::<i32>() }
}

unsafe fn ivlist_values_mut<'a>(lst: *mut IvList) -> &'a mut [*mut IVector] {
    let length = unsafe { (*lst).length };
    unsafe { slice::from_raw_parts_mut((*lst).array, length) }
}

fn abi_part_length(partition: &[i32]) -> usize {
    partition
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1)
}

fn abi_valid_partition(partition: &[i32]) -> bool {
    let mut previous = 0;
    for &part in partition.iter().rev() {
        if part < previous {
            return false;
        }
        previous = part;
    }
    true
}

fn abi_partition_leq(left: &[i32], right: &[i32]) -> bool {
    let len = abi_part_length(left);
    if len > abi_part_length(right) {
        return false;
    }
    (0..len).all(|index| left[index] <= right[index])
}

fn ivlc_table_alloc_size(length: u32) -> Option<usize> {
    usize::try_from(length)
        .ok()?
        .checked_mul(mem::size_of::<u32>())
}

fn ivlc_elts_alloc_size(length: u32) -> Option<usize> {
    usize::try_from(length)
        .ok()?
        .checked_mul(mem::size_of::<IvlcKeyVal>())
}

fn ivlc_table_index(hash: u32, table_sz: u32) -> u32 {
    hash % table_sz
}

fn ivlc_new_table_size(sz: u32) -> Option<u32> {
    let mut new_sz = IVLC_USE_FACTOR
        .checked_mul(2)?
        .checked_mul(sz)?
        .checked_add(1)?;
    if new_sz % 3 == 0 {
        new_sz = new_sz.checked_add(2)?;
    }
    if new_sz % 5 == 0 {
        new_sz = new_sz.checked_add(6)?;
    }
    if new_sz % 7 == 0 {
        new_sz = new_sz.checked_add(30)?;
    }
    Some(new_sz)
}

unsafe fn ivlc_lookup_index(ht: *mut IvLinComb, key: *const IVector, hash: u32) -> u32 {
    if ht.is_null() || key.is_null() || unsafe { (*ht).table_sz } == 0 {
        return 0;
    }
    let table_index = ivlc_table_index(hash, unsafe { (*ht).table_sz });
    let mut i = unsafe { *(*ht).table.add(table_index as usize) };
    while i != 0 {
        let elt = unsafe { (*ht).elts.add(i as usize) };
        if unsafe { iv_cmp(key, (*elt).key) } == 0 {
            return i;
        }
        i = unsafe { (*elt).next };
    }
    0
}

unsafe fn ivlc_makeroom_inner(ht: *mut IvLinComb, sz: u32) -> c_int {
    if ht.is_null() {
        return -1;
    }
    let needs_table = match IVLC_USE_FACTOR.checked_mul(sz) {
        Some(needed) => needed > unsafe { (*ht).table_sz },
        None => true,
    };
    if needs_table && unsafe { ivlc__grow_table(ht, sz) } != 0 {
        return -1;
    }
    let Some(needed_elts) = sz.checked_add(1) else {
        return -1;
    };
    if needed_elts > unsafe { (*ht).elts_sz } && unsafe { ivlc__grow_elts(ht, needed_elts) } != 0 {
        return -1;
    }
    0
}

unsafe fn ivector_from_partition(partition: &[i32], length: usize) -> *mut IVector {
    let Ok(length) = u32::try_from(length.max(partition.len())) else {
        return ptr::null_mut();
    };
    let key = iv_new_zero(length);
    if key.is_null() {
        return ptr::null_mut();
    }
    let dst = unsafe { ivector_values_mut(key) };
    for (index, &part) in partition.iter().enumerate() {
        dst[index] = part;
    }
    key
}

unsafe fn ivlc_from_terms(terms: &[SchurTerm], key_len: usize) -> *mut IvLinComb {
    let initial_elts = u32::try_from(terms.len().saturating_add(1))
        .unwrap_or(u32::MAX)
        .max(IVLC_ARRAY_SZ);
    let lc = ivlc_new(IVLC_HASHTABLE_SZ, initial_elts);
    if lc.is_null() {
        return ptr::null_mut();
    }
    for term in terms {
        let Ok(value) = i32::try_from(term.coefficient) else {
            unsafe { ivlc_free_all(lc) };
            return ptr::null_mut();
        };
        let key = unsafe { ivector_from_partition(&term.partition, key_len) };
        if key.is_null() {
            unsafe { ivlc_free_all(lc) };
            return ptr::null_mut();
        }
        let hash = unsafe { iv_hash(key) } as u32;
        if unsafe { ivlc_add_element(lc, value, key, hash, LC_FREE_ZERO) } != 0 {
            unsafe { ivlc_free_all(lc) };
            return ptr::null_mut();
        }
    }
    lc
}

unsafe fn ivlc_from_skew_expansion_direct(
    outer: &[i32],
    inner: &[i32],
    rows: c_int,
    key_len: usize,
) -> *mut IvLinComb {
    let lc = Cell::new(ptr::null_mut());
    let failed = Cell::new(false);
    let result = visit_schur_skew_expansion_with_len(
        outer,
        inner,
        rows,
        |term_count| {
            let initial_elts = u32::try_from(term_count.saturating_add(1))
                .unwrap_or(u32::MAX)
                .max(IVLC_ARRAY_SZ);
            let new_lc = ivlc_new(IVLC_HASHTABLE_SZ, initial_elts);
            lc.set(new_lc);
            if new_lc.is_null() {
                failed.set(true);
            }
        },
        |partition, coefficient| {
            if failed.get() {
                return;
            }
            let Ok(value) = i32::try_from(coefficient) else {
                failed.set(true);
                return;
            };
            let key = unsafe { ivector_from_partition(&partition, key_len) };
            if key.is_null() {
                failed.set(true);
                return;
            }
            let hash = unsafe { iv_hash(key) } as u32;
            if unsafe { ivlc_add_element(lc.get(), value, key, hash, LC_FREE_ZERO) } != 0 {
                failed.set(true);
            }
        },
    );
    if result.is_err() || failed.get() {
        if !lc.get().is_null() {
            unsafe { ivlc_free_all(lc.get()) };
        }
        return ptr::null_mut();
    }
    lc.get()
}

unsafe fn ivlc_from_signed_terms(terms: &[SignedSchurTerm], key_len: usize) -> *mut IvLinComb {
    let initial_elts = u32::try_from(terms.len().saturating_add(1))
        .unwrap_or(u32::MAX)
        .max(IVLC_ARRAY_SZ);
    let lc = ivlc_new(IVLC_HASHTABLE_SZ, initial_elts);
    if lc.is_null() {
        return ptr::null_mut();
    }
    for term in terms {
        let Ok(value) = i32::try_from(term.coefficient) else {
            unsafe { ivlc_free_all(lc) };
            return ptr::null_mut();
        };
        let key = unsafe { ivector_from_partition(&term.partition, key_len) };
        if key.is_null() {
            unsafe { ivlc_free_all(lc) };
            return ptr::null_mut();
        }
        let hash = unsafe { iv_hash(key) } as u32;
        if unsafe { ivlc_add_element(lc, value, key, hash, LC_FREE_ZERO) } != 0 {
            unsafe { ivlc_free_all(lc) };
            return ptr::null_mut();
        }
    }
    lc
}

unsafe fn ivlc_from_i32_terms(terms: &LinearCombination) -> *mut IvLinComb {
    let initial_elts = u32::try_from(terms.len().saturating_add(1))
        .unwrap_or(u32::MAX)
        .max(IVLC_ARRAY_SZ);
    let lc = ivlc_new(IVLC_HASHTABLE_SZ, initial_elts);
    if lc.is_null() {
        return ptr::null_mut();
    }
    if unsafe { ivlc_fill_i32_terms(lc, terms) } != 0 {
        unsafe { ivlc_free_all(lc) };
        return ptr::null_mut();
    }
    lc
}

unsafe fn ivlc_fill_i32_terms(lc: *mut IvLinComb, terms: &LinearCombination) -> c_int {
    for (key_values, &value) in terms {
        let key = unsafe { ivector_from_partition(key_values, key_values.len()) };
        if key.is_null() {
            return -1;
        }
        let hash = unsafe { iv_hash(key) } as u32;
        if unsafe { ivlc_add_element(lc, value, key, hash, LC_FREE_ZERO) } != 0 {
            return -1;
        }
    }
    0
}

unsafe fn ivlc_collect_i32_terms(lc: *mut IvLinComb) -> LinearCombination {
    let mut terms = LinearCombination::new();
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(lc, &mut itr);
        while ivlc_good(&itr) != 0 {
            let key = ivlc_key(&itr);
            terms.insert(ivector_values(key).to_vec(), ivlc_value(&itr));
            ivlc_next(&mut itr);
        }
    }
    terms
}

unsafe fn ivlc_collect_owned_keys(lc: *mut IvLinComb) -> Vec<*mut IVector> {
    let mut keys = Vec::new();
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(lc, &mut itr);
        while ivlc_good(&itr) != 0 {
            keys.push(ivlc_key(&itr));
            ivlc_next(&mut itr);
        }
    }
    keys
}

unsafe fn lrit_alloc(array_len: usize) -> *mut LrTabIter {
    let Ok(array_len_i32) = c_int::try_from(array_len) else {
        return ptr::null_mut();
    };
    let Some(size) = lrtab_iter_alloc_size(array_len) else {
        return ptr::null_mut();
    };
    let lrit = unsafe { libc::malloc(size) }.cast::<LrTabIter>();
    if lrit.is_null() {
        return ptr::null_mut();
    }
    unsafe {
        (*lrit).cont = ptr::null_mut();
        (*lrit).size = -1;
        (*lrit).array_len = array_len_i32;
    }
    lrit
}

unsafe fn lrit_empty() -> *mut LrTabIter {
    let lrit = unsafe { lrit_alloc(1) };
    if lrit.is_null() {
        return ptr::null_mut();
    }
    let cont = iv_new(1);
    if cont.is_null() {
        unsafe { libc::free(lrit.cast::<c_void>()) };
        return ptr::null_mut();
    }
    unsafe {
        (*lrit).cont = cont;
    }
    lrit
}

fn default_product_key_len(sh1: &[i32], sh2: &[i32], rows: c_int, partsz: c_int) -> usize {
    if partsz >= 0 {
        partsz as usize
    } else if rows >= 0 {
        rows as usize
    } else {
        let len1 = sh1.iter().rposition(|&part| part != 0).map_or(0, |i| i + 1);
        let len2 = sh2.iter().rposition(|&part| part != 0).map_or(0, |i| i + 1);
        len1.saturating_add(len2)
    }
}

fn default_skew_key_len(outer: &[i32], inner: &[i32], rows: c_int, partsz: c_int) -> usize {
    if partsz >= 0 {
        partsz as usize
    } else if rows >= 0 {
        rows as usize
    } else {
        let outer_sum: i64 = outer.iter().map(|&part| i64::from(part)).sum();
        let inner_sum: i64 = inner.iter().map(|&part| i64::from(part)).sum();
        usize::try_from(outer_sum.saturating_sub(inner_sum).max(0)).unwrap_or(0)
    }
}

fn default_coproduct_key_len(sh: &[i32], rows: c_int, partsz: c_int) -> usize {
    if partsz >= 0 {
        partsz as usize
    } else if rows >= 0 {
        let len = sh.iter().rposition(|&part| part != 0).map_or(0, |i| i + 1);
        (rows as usize).saturating_add(len)
    } else {
        sh.iter().map(|&part| i64::from(part)).sum::<i64>().max(0) as usize
    }
}

#[no_mangle]
pub extern "C" fn iv_new(length: u32) -> *mut IVector {
    let Some(size) = ivector_alloc_size(length) else {
        return ptr::null_mut();
    };
    let raw = unsafe { libc::malloc(size) }.cast::<IVector>();
    if raw.is_null() {
        return ptr::null_mut();
    }
    unsafe {
        (*raw).length = length;
    }
    raw
}

#[no_mangle]
pub extern "C" fn iv_new_zero(length: u32) -> *mut IVector {
    let Some(size) = ivector_alloc_size(length) else {
        return ptr::null_mut();
    };
    let raw = unsafe { libc::calloc(1, size) }.cast::<IVector>();
    if raw.is_null() {
        return ptr::null_mut();
    }
    unsafe {
        (*raw).length = length;
    }
    raw
}

/// # Safety
///
/// `v` must be null or a pointer returned by an `iv_*` allocation function
/// that has not already been freed.
#[no_mangle]
pub unsafe extern "C" fn iv_free(v: *mut IVector) {
    if !v.is_null() {
        unsafe { libc::free(v.cast::<c_void>()) };
    }
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn iv_new_copy(v: *const IVector) -> *mut IVector {
    if v.is_null() {
        return ptr::null_mut();
    }
    let out = iv_new(unsafe { (*v).length });
    if out.is_null() {
        return ptr::null_mut();
    }
    let src = unsafe { ivector_values(v) };
    let dst = unsafe { ivector_values_mut(out) };
    dst.copy_from_slice(src);
    out
}

#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub extern "C" fn iv_new_init(
    length: u32,
    x0: i32,
    x1: i32,
    x2: i32,
    x3: i32,
    x4: i32,
    x5: i32,
    x6: i32,
    x7: i32,
) -> *mut IVector {
    let v = iv_new(length);
    if v.is_null() {
        return ptr::null_mut();
    }
    unsafe {
        let data = ivector_data_mut(v);
        if length > 0 {
            *data.add(0) = x0;
        }
        if length > 1 {
            *data.add(1) = x1;
        }
        if length > 2 {
            *data.add(2) = x2;
        }
        if length > 3 {
            *data.add(3) = x3;
        }
        if length > 4 {
            *data.add(4) = x4;
        }
        if length > 5 {
            *data.add(5) = x5;
        }
        if length > 6 {
            *data.add(6) = x6;
        }
        if length > 7 {
            *data.add(7) = x7;
        }
        if length > 8 {
            for index in 8..length as usize {
                *data.add(index) = 0;
            }
        }
    }
    v
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid mutable `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn iv_set_zero(v: *mut IVector) {
    if !v.is_null() {
        unsafe { ivector_values_mut(v) }.fill(0);
    }
}

/// # Safety
///
/// Non-null pointers must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn iv_cmp(v1: *const IVector, v2: *const IVector) -> c_int {
    if v1.is_null() || v2.is_null() {
        return (v1 as isize).cmp(&(v2 as isize)) as c_int;
    }
    let len1 = unsafe { (*v1).length };
    let len2 = unsafe { (*v2).length };
    if len1 != len2 {
        return len1 as c_int - len2 as c_int;
    }
    for (&x, &y) in unsafe { ivector_values(v1) }
        .iter()
        .zip(unsafe { ivector_values(v2) })
    {
        if x != y {
            return x - y;
        }
    }
    0
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn iv_hash(v: *const IVector) -> i32 {
    if v.is_null() {
        return 0;
    }
    let mut h = unsafe { (*v).length };
    for &x in unsafe { ivector_values(v) } {
        h = ((h << 5) ^ (h >> 27)).wrapping_add(x as u32);
    }
    h as i32
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn iv_sum(v: *const IVector) -> i32 {
    if v.is_null() {
        return 0;
    }
    unsafe { ivector_values(v) }.iter().copied().sum()
}

/// # Safety
///
/// `dst` and `src` must point to valid vectors of the same length.
#[no_mangle]
pub unsafe extern "C" fn iv_copy(dst: *mut IVector, src: *const IVector) {
    if dst.is_null() || src.is_null() || unsafe { (*dst).length != (*src).length } {
        return;
    }
    let src_values = unsafe { ivector_values(src) };
    let dst_values = unsafe { ivector_values_mut(dst) };
    dst_values.copy_from_slice(src_values);
}

/// # Safety
///
/// Non-null pointers must point to valid vectors of the same length.
#[no_mangle]
pub unsafe extern "C" fn iv_lesseq(v1: *const IVector, v2: *const IVector) -> c_int {
    if v1.is_null() || v2.is_null() || unsafe { (*v1).length != (*v2).length } {
        return 0;
    }
    let left = unsafe { ivector_values(v1) };
    let right = unsafe { ivector_values(v2) };
    left.iter().zip(right).all(|(&x, &y)| x <= y) as c_int
}

/// # Safety
///
/// `dst` and `src` must point to valid vectors of the same length.
#[no_mangle]
pub unsafe extern "C" fn iv_mult(dst: *mut IVector, c: i32, src: *const IVector) {
    if dst.is_null() || src.is_null() || unsafe { (*dst).length != (*src).length } {
        return;
    }
    let src_values = unsafe { ivector_values(src) };
    let dst_values = unsafe { ivector_values_mut(dst) };
    for (out, &value) in dst_values.iter_mut().zip(src_values) {
        *out = c.wrapping_mul(value);
    }
}

/// # Safety
///
/// `dst` and `src` must point to valid vectors of the same length.
#[no_mangle]
pub unsafe extern "C" fn iv_div(dst: *mut IVector, src: *const IVector, c: i32) {
    if dst.is_null() || src.is_null() || c == 0 || unsafe { (*dst).length != (*src).length } {
        return;
    }
    let src_values = unsafe { ivector_values(src) };
    let dst_values = unsafe { ivector_values_mut(dst) };
    for (out, &value) in dst_values.iter_mut().zip(src_values) {
        *out = value / c;
    }
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn iv_max(v: *const IVector) -> i32 {
    if v.is_null() {
        return i32::MIN;
    }
    unsafe { ivector_values(v) }
        .iter()
        .copied()
        .max()
        .unwrap_or(i32::MIN)
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn iv_min(v: *const IVector) -> i32 {
    if v.is_null() {
        return i32::MAX;
    }
    unsafe { ivector_values(v) }
        .iter()
        .copied()
        .min()
        .unwrap_or(i32::MAX)
}

/// # Safety
///
/// `dst` and `src` must point to valid vectors of the same length.
#[no_mangle]
pub unsafe extern "C" fn iv_reverse(dst: *mut IVector, src: *const IVector) {
    if dst.is_null() || src.is_null() || unsafe { (*dst).length != (*src).length } {
        return;
    }
    let len = unsafe { (*dst).length as usize };
    for index in 0..len / 2 {
        let left = unsafe { *ivector_data(src).add(index) };
        let right = unsafe { *ivector_data(src).add(len - 1 - index) };
        unsafe {
            *ivector_data_mut(dst).add(index) = right;
            *ivector_data_mut(dst).add(len - 1 - index) = left;
        }
    }
    if len % 2 == 1 {
        unsafe {
            *ivector_data_mut(dst).add(len / 2) = *ivector_data(src).add(len / 2);
        }
    }
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn iv_gcd(v: *const IVector) -> i32 {
    if v.is_null() {
        return 0;
    }
    let mut gcd = 0i32;
    for &value in unsafe { ivector_values(v) } {
        let mut x = value;
        let mut y = gcd;
        while y != 0 {
            let z = x % y;
            x = y;
            y = z;
        }
        gcd = x;
    }
    gcd.wrapping_abs()
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn iv_print(v: *const IVector) {
    if v.is_null() {
        print!("()");
        return;
    }
    print!("(");
    for (index, value) in unsafe { ivector_values(v) }.iter().enumerate() {
        if index != 0 {
            print!(",");
        }
        print!("{value}");
    }
    print!(")");
}

/// # Safety
///
/// If `v` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn iv_printnl(v: *const IVector) {
    unsafe { iv_print(v) };
    println!();
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn part_valid(p: *const IVector) -> c_int {
    if p.is_null() {
        return 0;
    }
    abi_valid_partition(unsafe { ivector_values(p) }) as c_int
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn part_decr(p: *const IVector) -> c_int {
    if p.is_null() {
        return 0;
    }
    let values = unsafe { ivector_values(p) };
    values.windows(2).all(|pair| pair[0] >= pair[1]) as c_int
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn part_length(p: *const IVector) -> c_int {
    if p.is_null() {
        return 0;
    }
    c_int::try_from(abi_part_length(unsafe { ivector_values(p) })).unwrap_or(c_int::MAX)
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn part_entry(p: *const IVector, i: c_int) -> c_int {
    if p.is_null() || i < 0 {
        return 0;
    }
    let values = unsafe { ivector_values(p) };
    values.get(i as usize).copied().unwrap_or(0)
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid mutable vector.
#[no_mangle]
pub unsafe extern "C" fn part_chop(p: *mut IVector) {
    if p.is_null() {
        return;
    }
    let len = unsafe { part_length(p) };
    unsafe {
        (*p).length = u32::try_from(len.max(0)).unwrap_or(0);
    }
}

/// # Safety
///
/// `p` must have been allocated with enough capacity for `len` entries.
#[no_mangle]
pub unsafe extern "C" fn part_unchop(p: *mut IVector, len: c_int) {
    if p.is_null() || len < 0 {
        return;
    }
    let old_len = unsafe { (*p).length as usize };
    let new_len = len as usize;
    if new_len < old_len {
        unsafe {
            (*p).length = len as u32;
        }
        return;
    }
    unsafe {
        let values = ivector_values_mut_len(p, new_len);
        values[old_len..new_len].fill(0);
        (*p).length = len as u32;
    }
}

/// # Safety
///
/// Non-null pointers must point to valid vectors.
#[no_mangle]
pub unsafe extern "C" fn part_leq(p1: *const IVector, p2: *const IVector) -> c_int {
    if p1.is_null() || p2.is_null() {
        return 0;
    }
    abi_partition_leq(unsafe { ivector_values(p1) }, unsafe { ivector_values(p2) }) as c_int
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid partition vector.
#[no_mangle]
pub unsafe extern "C" fn part_conj(p: *const IVector) -> *mut IVector {
    if p.is_null() || unsafe { part_valid(p) } == 0 {
        return ptr::null_mut();
    }
    let values = unsafe { ivector_values(p) };
    let rows = abi_part_length(values);
    let cols = if rows == 0 {
        0
    } else {
        values[0].max(0) as usize
    };
    let conj = iv_new(cols as u32);
    if conj.is_null() {
        return ptr::null_mut();
    }
    let out = unsafe { ivector_values_mut(conj) };
    for (col, entry) in out.iter_mut().enumerate() {
        let threshold = (col + 1) as i32;
        *entry = values
            .iter()
            .take(rows)
            .filter(|&&part| part >= threshold)
            .count() as i32;
    }
    conj
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn part_print(p: *const IVector) {
    if p.is_null() {
        print!("()");
        return;
    }
    print!("(");
    for (index, value) in unsafe { ivector_values(p) }
        .iter()
        .take_while(|&&value| value != 0)
        .enumerate()
    {
        if index != 0 {
            print!(",");
        }
        print!("{value}");
    }
    print!(")");
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn part_printnl(p: *const IVector) {
    unsafe { part_print(p) };
    println!();
}

/// # Safety
///
/// `p` must point to a valid `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn part_qdegree(p: *const IVector, level: c_int) -> c_int {
    if p.is_null() {
        return 0;
    }
    let values = unsafe { ivector_values(p) };
    let rows = values.len();
    let Ok(rows_i64) = i64::try_from(rows) else {
        return 0;
    };
    let n = rows_i64 + i64::from(level);
    if rows == 0 || n <= 0 {
        return 0;
    }

    let mut degree = 0i64;
    for (index, &value) in values.iter().enumerate() {
        let Ok(index_i64) = i64::try_from(index) else {
            return 0;
        };
        let Some(a) = i64::from(value)
            .checked_add(rows_i64)
            .and_then(|a| a.checked_sub(index_i64))
            .and_then(|a| a.checked_sub(1))
        else {
            return 0;
        };
        degree += a.div_euclid(n);
    }
    i32::try_from(degree).unwrap_or(0)
}

/// # Safety
///
/// `p` must point to a valid `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn part_qentry(p: *const IVector, i: c_int, d: c_int, level: c_int) -> c_int {
    if p.is_null() {
        return 0;
    }
    let values = unsafe { ivector_values(p) };
    let rows = values.len();
    let Ok(rows_i64) = i64::try_from(rows) else {
        return 0;
    };
    if rows == 0 {
        return 0;
    }

    let shifted = i64::from(i) + i64::from(d);
    let Ok(source) = usize::try_from(shifted.rem_euclid(rows_i64)) else {
        return 0;
    };
    let value = i64::from(values[source]) - (shifted / rows_i64) * i64::from(level) - i64::from(d);
    i32::try_from(value).unwrap_or(0)
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn part_qprint(p: *const IVector, level: c_int) {
    if p.is_null() {
        print!("()");
        return;
    }
    let d = unsafe { part_qdegree(p, level) };
    print!("(");
    let values = unsafe { ivector_values(p) };
    for index in 0..values.len() {
        let x = unsafe { part_qentry(p, index as c_int, d, level) };
        if x == 0 {
            break;
        }
        if index != 0 {
            print!(",");
        }
        print!("{x}");
    }
    print!(")");
}

/// # Safety
///
/// If `p` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn part_qprintnl(p: *const IVector, level: c_int) {
    unsafe { part_qprint(p, level) };
    println!();
}

unsafe fn il_realloc_array_inner(lst: *mut IList, sz: usize) -> c_int {
    if lst.is_null() {
        return -1;
    }
    let Some(new_allocated) = sz.checked_mul(2) else {
        return -1;
    };
    let Some(bytes) = new_allocated.checked_mul(mem::size_of::<c_int>()) else {
        return -1;
    };
    let raw = unsafe { libc::realloc((*lst).array.cast::<c_void>(), bytes) }.cast::<c_int>();
    if raw.is_null() {
        return -1;
    }
    unsafe {
        (*lst).array = raw;
        (*lst).allocated = new_allocated;
    }
    0
}

unsafe fn ivl_realloc_array_inner(lst: *mut IvList, sz: usize) -> c_int {
    if lst.is_null() {
        return -1;
    }
    let Some(new_allocated) = sz.checked_mul(2) else {
        return -1;
    };
    let Some(bytes) = new_allocated.checked_mul(mem::size_of::<*mut IVector>()) else {
        return -1;
    };
    let raw = unsafe { libc::realloc((*lst).array.cast::<c_void>(), bytes) }.cast::<*mut IVector>();
    if raw.is_null() {
        return -1;
    }
    unsafe {
        (*lst).array = raw;
        (*lst).allocated = new_allocated;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn il_init(lst: *mut IList, sz: usize) -> c_int {
    if lst.is_null() {
        return -1;
    }
    let Some(bytes) = sz.checked_mul(mem::size_of::<c_int>()) else {
        return -1;
    };
    let array = unsafe { libc::malloc(bytes) }.cast::<c_int>();
    if array.is_null() && bytes != 0 {
        return -1;
    }
    unsafe {
        (*lst).array = array;
        (*lst).allocated = sz;
        (*lst).length = 0;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn il_new(sz: usize) -> *mut IList {
    let lst = unsafe { libc::malloc(mem::size_of::<IList>()) }.cast::<IList>();
    if lst.is_null() {
        return ptr::null_mut();
    }
    if unsafe { il_init(lst, sz) } != 0 {
        unsafe { libc::free(lst.cast::<c_void>()) };
        return ptr::null_mut();
    }
    lst
}

#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn il_new_init(
    sz: usize,
    count: usize,
    x0: c_int,
    x1: c_int,
    x2: c_int,
    x3: c_int,
    x4: c_int,
    x5: c_int,
    x6: c_int,
    x7: c_int,
) -> *mut IList {
    let list = unsafe { il_new(sz) };
    if list.is_null() {
        return ptr::null_mut();
    }
    if count > 0 && unsafe { il_append(list, x0) } != 0
        || count > 1 && unsafe { il_append(list, x1) } != 0
        || count > 2 && unsafe { il_append(list, x2) } != 0
        || count > 3 && unsafe { il_append(list, x3) } != 0
        || count > 4 && unsafe { il_append(list, x4) } != 0
        || count > 5 && unsafe { il_append(list, x5) } != 0
        || count > 6 && unsafe { il_append(list, x6) } != 0
        || count > 7 && unsafe { il_append(list, x7) } != 0
    {
        unsafe { il_free(list) };
        return ptr::null_mut();
    }
    list
}

#[no_mangle]
pub unsafe extern "C" fn il_dealloc(lst: *mut IList) {
    if !lst.is_null() {
        unsafe {
            libc::free((*lst).array.cast::<c_void>());
            (*lst).array = ptr::null_mut();
            (*lst).allocated = 0;
            (*lst).length = 0;
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn il_free(lst: *mut IList) {
    if !lst.is_null() {
        unsafe {
            il_dealloc(lst);
            libc::free(lst.cast::<c_void>());
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn il_reset(lst: *mut IList) {
    if !lst.is_null() {
        unsafe {
            (*lst).length = 0;
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn il__realloc_array(lst: *mut IList, sz: usize) -> c_int {
    unsafe { il_realloc_array_inner(lst, sz) }
}

#[no_mangle]
pub unsafe extern "C" fn il_makeroom(lst: *mut IList, sz: usize) -> c_int {
    if lst.is_null() {
        return -1;
    }
    if sz <= unsafe { (*lst).allocated } {
        0
    } else {
        unsafe { il_realloc_array_inner(lst, sz) }
    }
}

#[no_mangle]
pub unsafe extern "C" fn il_append(lst: *mut IList, x: c_int) -> c_int {
    if lst.is_null() {
        return -1;
    }
    let new_len = unsafe { (*lst).length }.saturating_add(1);
    if unsafe { il_makeroom(lst, new_len) } != 0 {
        return -1;
    }
    unsafe {
        *(*lst).array.add((*lst).length) = x;
        (*lst).length = new_len;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn il_poplast(lst: *mut IList) -> c_int {
    if lst.is_null() || unsafe { (*lst).length == 0 } {
        return 0;
    }
    unsafe {
        (*lst).length -= 1;
        *(*lst).array.add((*lst).length)
    }
}

#[no_mangle]
pub unsafe extern "C" fn il_insert(lst: *mut IList, i: usize, x: c_int) -> c_int {
    if lst.is_null() || i > unsafe { (*lst).length } {
        return -1;
    }
    let len = unsafe { (*lst).length };
    if unsafe { il_makeroom(lst, len.saturating_add(1)) } != 0 {
        return -1;
    }
    unsafe {
        ptr::copy((*lst).array.add(i), (*lst).array.add(i + 1), len - i);
        *(*lst).array.add(i) = x;
        (*lst).length = len + 1;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn il_delete(lst: *mut IList, i: usize) -> c_int {
    if lst.is_null() || i >= unsafe { (*lst).length } {
        return 0;
    }
    unsafe {
        let value = *(*lst).array.add(i);
        (*lst).length -= 1;
        ptr::copy(
            (*lst).array.add(i + 1),
            (*lst).array.add(i),
            (*lst).length - i,
        );
        value
    }
}

#[no_mangle]
pub unsafe extern "C" fn il_fastdelete(lst: *mut IList, i: usize) -> c_int {
    if lst.is_null() || i >= unsafe { (*lst).length } {
        return 0;
    }
    unsafe {
        let value = *(*lst).array.add(i);
        (*lst).length -= 1;
        *(*lst).array.add(i) = *(*lst).array.add((*lst).length);
        value
    }
}

#[no_mangle]
pub unsafe extern "C" fn il_extend(dst: *mut IList, src: *mut IList) -> c_int {
    if dst.is_null() || src.is_null() {
        return -1;
    }
    let dlen = unsafe { (*dst).length };
    let slen = unsafe { (*src).length };
    if unsafe { il_makeroom(dst, dlen.saturating_add(slen)) } != 0 {
        return -1;
    }
    unsafe {
        ptr::copy((*src).array, (*dst).array.add(dlen), slen);
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn il_copy(dst: *mut IList, src: *mut IList) -> c_int {
    if dst.is_null() || src.is_null() {
        return -1;
    }
    let slen = unsafe { (*src).length };
    if unsafe { il_makeroom(dst, slen) } != 0 {
        return -1;
    }
    unsafe {
        (*dst).length = slen;
        ptr::copy_nonoverlapping((*src).array, (*dst).array, slen);
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn il_new_copy(lst: *mut IList) -> *mut IList {
    if lst.is_null() {
        return ptr::null_mut();
    }
    let out = unsafe { il_new((*lst).length) };
    if out.is_null() {
        return ptr::null_mut();
    }
    if unsafe { il_copy(out, lst) } != 0 {
        unsafe { il_free(out) };
        return ptr::null_mut();
    }
    out
}

#[no_mangle]
pub unsafe extern "C" fn il_reverse(dst: *mut IList, src: *mut IList) -> c_int {
    if dst.is_null() || src.is_null() {
        return -1;
    }
    let n = unsafe { (*src).length };
    if dst != src && unsafe { il_makeroom(dst, n) } != 0 {
        return -1;
    }
    for i in 0..n / 2 {
        let left = unsafe { *(*src).array.add(i) };
        let right = unsafe { *(*src).array.add(n - 1 - i) };
        unsafe {
            *(*dst).array.add(i) = right;
            *(*dst).array.add(n - 1 - i) = left;
        }
    }
    if n % 2 == 1 {
        unsafe {
            *(*dst).array.add(n / 2) = *(*src).array.add(n / 2);
        }
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn ivl_init(lst: *mut IvList, sz: usize) -> c_int {
    if lst.is_null() {
        return -1;
    }
    let Some(bytes) = sz.checked_mul(mem::size_of::<*mut IVector>()) else {
        return -1;
    };
    let array = unsafe { libc::malloc(bytes) }.cast::<*mut IVector>();
    if array.is_null() && bytes != 0 {
        return -1;
    }
    unsafe {
        (*lst).array = array;
        (*lst).allocated = sz;
        (*lst).length = 0;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn ivl_new(sz: usize) -> *mut IvList {
    let lst = unsafe { libc::malloc(mem::size_of::<IvList>()) }.cast::<IvList>();
    if lst.is_null() {
        return ptr::null_mut();
    }
    if unsafe { ivl_init(lst, sz) } != 0 {
        unsafe { libc::free(lst.cast::<c_void>()) };
        return ptr::null_mut();
    }
    lst
}

#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn ivl_new_init(
    sz: usize,
    count: usize,
    x0: *mut IVector,
    x1: *mut IVector,
    x2: *mut IVector,
    x3: *mut IVector,
    x4: *mut IVector,
    x5: *mut IVector,
    x6: *mut IVector,
    x7: *mut IVector,
) -> *mut IvList {
    let list = unsafe { ivl_new(sz) };
    if list.is_null() {
        return ptr::null_mut();
    }
    if count > 0 && unsafe { ivl_append(list, x0) } != 0
        || count > 1 && unsafe { ivl_append(list, x1) } != 0
        || count > 2 && unsafe { ivl_append(list, x2) } != 0
        || count > 3 && unsafe { ivl_append(list, x3) } != 0
        || count > 4 && unsafe { ivl_append(list, x4) } != 0
        || count > 5 && unsafe { ivl_append(list, x5) } != 0
        || count > 6 && unsafe { ivl_append(list, x6) } != 0
        || count > 7 && unsafe { ivl_append(list, x7) } != 0
    {
        unsafe { ivl_free(list) };
        return ptr::null_mut();
    }
    list
}

#[no_mangle]
pub unsafe extern "C" fn ivl_dealloc(lst: *mut IvList) {
    if !lst.is_null() {
        unsafe {
            libc::free((*lst).array.cast::<c_void>());
            (*lst).array = ptr::null_mut();
            (*lst).allocated = 0;
            (*lst).length = 0;
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn ivl_free(lst: *mut IvList) {
    if !lst.is_null() {
        unsafe {
            ivl_dealloc(lst);
            libc::free(lst.cast::<c_void>());
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn ivl_free_all(lst: *mut IvList) {
    if lst.is_null() {
        return;
    }
    unsafe {
        for &item in ivlist_values_mut(lst).iter() {
            iv_free(item);
        }
        ivl_free(lst);
    }
}

#[no_mangle]
pub unsafe extern "C" fn ivl_reset(lst: *mut IvList) {
    if !lst.is_null() {
        unsafe {
            (*lst).length = 0;
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn ivl__realloc_array(lst: *mut IvList, sz: usize) -> c_int {
    unsafe { ivl_realloc_array_inner(lst, sz) }
}

#[no_mangle]
pub unsafe extern "C" fn ivl_makeroom(lst: *mut IvList, sz: usize) -> c_int {
    if lst.is_null() {
        return -1;
    }
    if sz <= unsafe { (*lst).allocated } {
        0
    } else {
        unsafe { ivl_realloc_array_inner(lst, sz) }
    }
}

#[no_mangle]
pub unsafe extern "C" fn ivl_append(lst: *mut IvList, x: *mut IVector) -> c_int {
    if lst.is_null() {
        return -1;
    }
    let new_len = unsafe { (*lst).length }.saturating_add(1);
    if unsafe { ivl_makeroom(lst, new_len) } != 0 {
        return -1;
    }
    unsafe {
        *(*lst).array.add((*lst).length) = x;
        (*lst).length = new_len;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn ivl_poplast(lst: *mut IvList) -> *mut IVector {
    if lst.is_null() || unsafe { (*lst).length == 0 } {
        return ptr::null_mut();
    }
    unsafe {
        (*lst).length -= 1;
        *(*lst).array.add((*lst).length)
    }
}

#[no_mangle]
pub unsafe extern "C" fn ivl_insert(lst: *mut IvList, i: usize, x: *mut IVector) -> c_int {
    if lst.is_null() || i > unsafe { (*lst).length } {
        return -1;
    }
    let len = unsafe { (*lst).length };
    if unsafe { ivl_makeroom(lst, len.saturating_add(1)) } != 0 {
        return -1;
    }
    unsafe {
        ptr::copy((*lst).array.add(i), (*lst).array.add(i + 1), len - i);
        *(*lst).array.add(i) = x;
        (*lst).length = len + 1;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn ivl_delete(lst: *mut IvList, i: usize) -> *mut IVector {
    if lst.is_null() || i >= unsafe { (*lst).length } {
        return ptr::null_mut();
    }
    unsafe {
        let value = *(*lst).array.add(i);
        (*lst).length -= 1;
        ptr::copy(
            (*lst).array.add(i + 1),
            (*lst).array.add(i),
            (*lst).length - i,
        );
        value
    }
}

#[no_mangle]
pub unsafe extern "C" fn ivl_fastdelete(lst: *mut IvList, i: usize) -> *mut IVector {
    if lst.is_null() || i >= unsafe { (*lst).length } {
        return ptr::null_mut();
    }
    unsafe {
        let value = *(*lst).array.add(i);
        (*lst).length -= 1;
        *(*lst).array.add(i) = *(*lst).array.add((*lst).length);
        value
    }
}

#[no_mangle]
pub unsafe extern "C" fn ivl_extend(dst: *mut IvList, src: *mut IvList) -> c_int {
    if dst.is_null() || src.is_null() {
        return -1;
    }
    let dlen = unsafe { (*dst).length };
    let slen = unsafe { (*src).length };
    if unsafe { ivl_makeroom(dst, dlen.saturating_add(slen)) } != 0 {
        return -1;
    }
    unsafe {
        ptr::copy((*src).array, (*dst).array.add(dlen), slen);
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn ivl_copy(dst: *mut IvList, src: *mut IvList) -> c_int {
    if dst.is_null() || src.is_null() {
        return -1;
    }
    let slen = unsafe { (*src).length };
    if unsafe { ivl_makeroom(dst, slen) } != 0 {
        return -1;
    }
    unsafe {
        (*dst).length = slen;
        ptr::copy_nonoverlapping((*src).array, (*dst).array, slen);
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn ivl_new_copy(lst: *mut IvList) -> *mut IvList {
    if lst.is_null() {
        return ptr::null_mut();
    }
    let out = unsafe { ivl_new((*lst).length) };
    if out.is_null() {
        return ptr::null_mut();
    }
    if unsafe { ivl_copy(out, lst) } != 0 {
        unsafe { ivl_free(out) };
        return ptr::null_mut();
    }
    out
}

#[no_mangle]
pub unsafe extern "C" fn ivl_reverse(dst: *mut IvList, src: *mut IvList) -> c_int {
    if dst.is_null() || src.is_null() {
        return -1;
    }
    let n = unsafe { (*src).length };
    if dst != src && unsafe { ivl_makeroom(dst, n) } != 0 {
        return -1;
    }
    for i in 0..n / 2 {
        let left = unsafe { *(*src).array.add(i) };
        let right = unsafe { *(*src).array.add(n - 1 - i) };
        unsafe {
            *(*dst).array.add(i) = right;
            *(*dst).array.add(n - 1 - i) = left;
        }
    }
    if n % 2 == 1 {
        unsafe {
            *(*dst).array.add(n / 2) = *(*src).array.add(n / 2);
        }
    }
    0
}

fn abi_perm_group(values: &[i32]) -> usize {
    let mut len = values.len();
    while len > 0 && values[len - 1] == len as i32 {
        len -= 1;
    }
    len
}

fn abi_perm_length(values: &[i32]) -> c_int {
    let mut inversions = 0i32;
    for i in 0..values.len().saturating_sub(1) {
        for j in i + 1..values.len() {
            if values[i] > values[j] {
                inversions = inversions.saturating_add(1);
            }
        }
    }
    inversions
}

fn abi_dimvec(values: &[i32]) -> Option<Vec<i32>> {
    let mut classes = 0usize;
    for &value in values {
        if value < 0 {
            return None;
        }
        classes = classes.max(usize::try_from(value).ok()?.saturating_add(1));
    }
    let mut out = vec![0; classes];
    for &value in values {
        out[usize::try_from(value).ok()?] += 1;
    }
    for index in 1..out.len() {
        out[index] += out[index - 1];
    }
    Some(out)
}

/// # Safety
///
/// If `w` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn perm_valid(w: *const IVector) -> c_int {
    if w.is_null() {
        return 0;
    }
    let values = unsafe { ivector_values(w) };
    let n = values.len();
    let mut seen = vec![false; n];
    for &value in values {
        let Ok(index) = usize::try_from(value - 1) else {
            return 0;
        };
        if index >= n || seen[index] {
            return 0;
        }
        seen[index] = true;
    }
    1
}

/// # Safety
///
/// If `w` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn perm_length(w: *const IVector) -> c_int {
    if w.is_null() {
        return 0;
    }
    abi_perm_length(unsafe { ivector_values(w) })
}

/// # Safety
///
/// If `w` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn perm_group(w: *const IVector) -> c_int {
    if w.is_null() {
        return 0;
    }
    c_int::try_from(abi_perm_group(unsafe { ivector_values(w) })).unwrap_or(c_int::MAX)
}

/// # Safety
///
/// If `dv` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn dimvec_valid(dv: *const IVector) -> c_int {
    if dv.is_null() {
        return 0;
    }
    let values = unsafe { ivector_values(dv) };
    if values.is_empty() || values[0] < 0 {
        return 0;
    }
    values.windows(2).all(|pair| pair[0] <= pair[1]) as c_int
}

/// # Safety
///
/// `w1` and `w2` must point to valid vectors.
#[no_mangle]
pub unsafe extern "C" fn bruhat_leq(w1: *const IVector, w2: *const IVector) -> c_int {
    if w1.is_null() || w2.is_null() {
        return 0;
    }
    let left = unsafe { ivector_values(w1) };
    let right = unsafe { ivector_values(w2) };
    let n = abi_perm_group(left);
    if n > abi_perm_group(right) {
        return 0;
    }
    for q in 1..n {
        let mut r1 = 0;
        let mut r2 = 0;
        for p in 0..n.saturating_sub(1) {
            if left[p] <= q as i32 {
                r1 += 1;
            }
            if right[p] <= q as i32 {
                r2 += 1;
            }
            if r1 < r2 {
                return 0;
            }
        }
    }
    1
}

/// # Safety
///
/// `w1` and `w2` must point to valid vectors.
#[no_mangle]
pub unsafe extern "C" fn bruhat_zero(w1: *const IVector, w2: *const IVector, rank: c_int) -> c_int {
    if w1.is_null() || w2.is_null() || rank < 0 {
        return 1;
    }
    let mut left = unsafe { ivector_values(w1) };
    let mut right = unsafe { ivector_values(w2) };
    let mut n1 = abi_perm_group(left);
    let n2 = abi_perm_group(right);
    if n1 > rank as usize || n2 > rank as usize {
        return 1;
    }
    if n1 > n2 {
        std::mem::swap(&mut left, &mut right);
        n1 = n2;
    }
    for q in 1..n1 {
        let q2 = rank - q as i32;
        let mut r1 = 0;
        let mut r2 = 0;
        for p in 0..n1.saturating_sub(1) {
            if left[p] <= q as i32 {
                r1 += 1;
            }
            if right[p] > q2 {
                r2 += 1;
            }
            if r1 < r2 {
                return 1;
            }
        }
    }
    0
}

/// # Safety
///
/// If `str` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn str2dimvec(str: *const IVector) -> *mut IVector {
    if str.is_null() {
        return ptr::null_mut();
    }
    let Some(dimvec) = abi_dimvec(unsafe { ivector_values(str) }) else {
        return ptr::null_mut();
    };
    unsafe { ivector_from_partition(&dimvec, dimvec.len()) }
}

/// # Safety
///
/// `str1` and `str2` must point to valid vectors.
#[no_mangle]
pub unsafe extern "C" fn str_iscompat(str1: *const IVector, str2: *const IVector) -> c_int {
    if str1.is_null() || str2.is_null() || unsafe { (*str1).length != (*str2).length } {
        return 0;
    }
    let left = abi_dimvec(unsafe { ivector_values(str1) });
    let right = abi_dimvec(unsafe { ivector_values(str2) });
    (left.is_some() && left == right) as c_int
}

/// # Safety
///
/// If `str` is non-null, it must point to a valid vector.
#[no_mangle]
pub unsafe extern "C" fn string2perm(str: *const IVector) -> *mut IVector {
    if str.is_null() {
        return ptr::null_mut();
    }
    let string = unsafe { ivector_values(str) };
    let Some(mut dimvec) = abi_dimvec(string) else {
        return ptr::null_mut();
    };
    let perm = iv_new(string.len() as u32);
    if perm.is_null() {
        return ptr::null_mut();
    }
    let out = unsafe { ivector_values_mut(perm) };
    for index in (0..string.len()).rev() {
        let Ok(class) = usize::try_from(string[index]) else {
            unsafe { iv_free(perm) };
            return ptr::null_mut();
        };
        dimvec[class] -= 1;
        let Ok(target) = usize::try_from(dimvec[class]) else {
            unsafe { iv_free(perm) };
            return ptr::null_mut();
        };
        out[target] = (index + 1) as i32;
    }
    perm
}

/// # Safety
///
/// `perm` and `dimvec` must point to valid vectors.
#[no_mangle]
pub unsafe extern "C" fn perm2string(perm: *const IVector, dimvec: *const IVector) -> *mut IVector {
    if perm.is_null() || dimvec.is_null() {
        return ptr::null_mut();
    }
    let perm_values = unsafe { ivector_values(perm) };
    let dim_values = unsafe { ivector_values(dimvec) };
    let n = dim_values.last().copied().unwrap_or(0);
    if n < 0 {
        return ptr::null_mut();
    }
    let out = iv_new(n as u32);
    if out.is_null() {
        return ptr::null_mut();
    }
    let out_values = unsafe { ivector_values_mut(out) };
    let mut j = 0usize;
    for (class, &limit) in dim_values.iter().enumerate() {
        let Ok(limit) = usize::try_from(limit) else {
            unsafe { iv_free(out) };
            return ptr::null_mut();
        };
        while j < limit {
            let wj = perm_values.get(j).copied().unwrap_or((j + 1) as i32);
            let Ok(target) = usize::try_from(wj - 1) else {
                unsafe { iv_free(out) };
                return ptr::null_mut();
            };
            if target >= out_values.len() {
                unsafe { iv_free(out) };
                return ptr::null_mut();
            }
            out_values[target] = class as i32;
            j += 1;
        }
    }
    out
}

/// # Safety
///
/// `dimvec` must point to a valid dimension vector.
#[no_mangle]
pub unsafe extern "C" fn all_strings(dimvec: *const IVector) -> *mut IvList {
    if dimvec.is_null() || unsafe { dimvec_valid(dimvec) } == 0 {
        return ptr::null_mut();
    }
    let dim_values = unsafe { ivector_values(dimvec) };
    let n = dim_values.last().copied().unwrap_or(0);
    if n < 0 {
        return ptr::null_mut();
    }
    let ld = dim_values.len();
    let mut counts = vec![0i32; ld];
    let mut current = Vec::with_capacity(n as usize);
    let mut j = 0i32;
    for (class, &limit) in dim_values.iter().enumerate() {
        while j < limit {
            current.push(class as i32);
            j += 1;
        }
    }

    let res = unsafe { ivl_new(200) };
    if res.is_null() {
        return ptr::null_mut();
    }
    if n == 0 {
        let str_vec = unsafe { ivector_from_partition(&current, current.len()) };
        if str_vec.is_null() || unsafe { ivl_append(res, str_vec) } != 0 {
            unsafe {
                iv_free(str_vec);
                ivl_free_all(res);
            }
            return ptr::null_mut();
        }
        return res;
    }

    loop {
        let str_vec = unsafe { ivector_from_partition(&current, current.len()) };
        if str_vec.is_null() || unsafe { ivl_append(res, str_vec) } != 0 {
            unsafe {
                iv_free(str_vec);
                ivl_free_all(res);
            }
            return ptr::null_mut();
        }

        let mut pos = current.len() - 1;
        counts[current[pos] as usize] += 1;
        while pos > 0 && current[pos - 1] >= current[pos] {
            pos -= 1;
            counts[current[pos] as usize] += 1;
        }
        if pos == 0 {
            break;
        }

        let mut class = current[pos - 1] as usize;
        counts[class] += 1;
        class += 1;
        while class < counts.len() && counts[class] == 0 {
            class += 1;
        }
        if class == counts.len() {
            break;
        }
        current[pos - 1] = class as i32;
        counts[class] -= 1;

        for (class, count) in counts.iter_mut().enumerate() {
            for _ in 0..*count {
                current[pos] = class as i32;
                pos += 1;
            }
            *count = 0;
        }
    }

    res
}

#[no_mangle]
pub extern "C" fn all_perms(n: c_int) -> *mut IvList {
    if n < 0 {
        return ptr::null_mut();
    }
    let dimvec = iv_new((n + 1) as u32);
    if dimvec.is_null() {
        return ptr::null_mut();
    }
    unsafe {
        for (index, entry) in ivector_values_mut(dimvec).iter_mut().enumerate() {
            *entry = index as i32;
        }
    }
    let res = unsafe { all_strings(dimvec) };
    unsafe { iv_free(dimvec) };
    res
}

/// # Safety
///
/// If `itr` is non-null, it must point to a valid partition iterator.
#[no_mangle]
pub unsafe extern "C" fn pitr_good(itr: *const PartIter) -> c_int {
    if itr.is_null() || unsafe { (*itr).rows < 0 } {
        0
    } else {
        1
    }
}

/// # Safety
///
/// Pointers must refer to valid iterator/vector storage as in upstream
/// `pitr_first`.
#[no_mangle]
pub unsafe extern "C" fn pitr_first(
    itr: *mut PartIter,
    p: *mut IVector,
    mut rows: c_int,
    cols: c_int,
    outer: *mut IVector,
    inner: *mut IVector,
    mut size: c_int,
    opt: c_int,
) -> c_int {
    if itr.is_null() || p.is_null() {
        return -1;
    }
    let use_outer = opt & PITR_USE_OUTER != 0;
    let use_inner = opt & PITR_USE_INNER != 0;
    let use_size = opt & PITR_USE_SIZE != 0;
    unsafe {
        (*itr).part = p;
        (*itr).outer = if use_outer { outer } else { ptr::null_mut() };
        (*itr).inner = if use_inner { inner } else { ptr::null_mut() };
        (*itr).opt = opt;
    }
    if (use_outer && outer.is_null()) || (use_inner && inner.is_null()) {
        unsafe {
            (*itr).rows = -1;
        }
        return 0;
    }

    if cols == 0 {
        rows = 0;
    }
    if use_size && rows > size {
        rows = size;
    }
    if use_outer {
        let outer_len = unsafe { (*outer).length as c_int };
        if rows > outer_len {
            rows = outer_len;
        }
        while rows > 0 && unsafe { part_entry(outer, rows - 1) } == 0 {
            rows -= 1;
        }
    }
    unsafe {
        (*itr).rows = rows;
        (*itr).length = rows;
        iv_set_zero(p);
    }

    if use_inner {
        let inner_len = unsafe { (*inner).length as c_int };
        if inner_len > rows && unsafe { part_entry(inner, rows) } != 0 {
            unsafe {
                (*itr).rows = -1;
            }
            return 0;
        }
        if rows > 0 && cols < unsafe { part_entry(inner, 0) } {
            unsafe {
                (*itr).rows = -1;
            }
            return 0;
        }
    }

    let mut inner_sz = 0;
    if use_size {
        if size > rows.saturating_mul(cols) {
            unsafe {
                (*itr).rows = -1;
            }
            return 0;
        }
        if use_inner {
            inner_sz = unsafe { iv_sum(inner) };
            if size < inner_sz {
                unsafe {
                    (*itr).rows = -1;
                }
                return 0;
            }
        }
    }

    let mut r = 0;
    while r < rows {
        let mut c = cols;
        if use_outer {
            c = c.min(unsafe { part_entry(outer, r) });
        }
        if use_size {
            let mut avail = size;
            if use_inner {
                inner_sz -= unsafe { part_entry(inner, r) };
                avail -= inner_sz;
            }
            if avail == 0 {
                unsafe {
                    (*itr).length = r;
                }
                return 0;
            }
            c = c.min(avail);
            size -= c;
        }
        unsafe {
            *ivector_data_mut(p).add(r as usize) = c;
        }
        r += 1;
    }

    if use_size && size > 0 {
        unsafe {
            (*itr).rows = -1;
        }
        return 0;
    }
    unsafe {
        (*itr).length = r;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn pitr_box_first(
    itr: *mut PartIter,
    p: *mut IVector,
    rows: c_int,
    cols: c_int,
) {
    unsafe { pitr_first(itr, p, rows, cols, ptr::null_mut(), ptr::null_mut(), 0, 0) };
}

#[no_mangle]
pub unsafe extern "C" fn pitr_box_sz_first(
    itr: *mut PartIter,
    p: *mut IVector,
    rows: c_int,
    cols: c_int,
    size: c_int,
) {
    unsafe {
        pitr_first(
            itr,
            p,
            rows,
            cols,
            ptr::null_mut(),
            ptr::null_mut(),
            size,
            PITR_USE_SIZE,
        )
    };
}

#[no_mangle]
pub unsafe extern "C" fn pitr_sub_first(itr: *mut PartIter, p: *mut IVector, outer: *mut IVector) {
    let rows = if outer.is_null() {
        0
    } else {
        unsafe { (*outer).length as c_int }
    };
    let cols = if rows == 0 {
        0
    } else {
        unsafe { part_entry(outer, 0) }
    };
    unsafe {
        pitr_first(
            itr,
            p,
            rows,
            cols,
            outer,
            ptr::null_mut(),
            0,
            PITR_USE_OUTER,
        )
    };
}

#[no_mangle]
pub unsafe extern "C" fn pitr_sub_sz_first(
    itr: *mut PartIter,
    p: *mut IVector,
    outer: *mut IVector,
    size: c_int,
) {
    let rows = if outer.is_null() {
        0
    } else {
        unsafe { (*outer).length as c_int }
    };
    let cols = if rows == 0 {
        0
    } else {
        unsafe { part_entry(outer, 0) }
    };
    unsafe {
        pitr_first(
            itr,
            p,
            rows,
            cols,
            outer,
            ptr::null_mut(),
            size,
            PITR_USE_OUTER | PITR_USE_SIZE,
        )
    };
}

#[no_mangle]
pub unsafe extern "C" fn pitr_between_first(
    itr: *mut PartIter,
    p: *mut IVector,
    outer: *mut IVector,
    inner: *mut IVector,
) {
    let rows = if outer.is_null() {
        0
    } else {
        unsafe { (*outer).length as c_int }
    };
    let cols = if rows == 0 {
        0
    } else {
        unsafe { part_entry(outer, 0) }
    };
    unsafe {
        pitr_first(
            itr,
            p,
            rows,
            cols,
            outer,
            inner,
            0,
            PITR_USE_OUTER | PITR_USE_INNER,
        )
    };
}

#[no_mangle]
pub unsafe extern "C" fn pitr_between_sz_first(
    itr: *mut PartIter,
    p: *mut IVector,
    outer: *mut IVector,
    inner: *mut IVector,
    size: c_int,
) {
    let rows = if outer.is_null() {
        0
    } else {
        unsafe { (*outer).length as c_int }
    };
    let cols = if rows == 0 {
        0
    } else {
        unsafe { part_entry(outer, 0) }
    };
    unsafe {
        pitr_first(
            itr,
            p,
            rows,
            cols,
            outer,
            inner,
            size,
            PITR_USE_OUTER | PITR_USE_INNER | PITR_USE_SIZE,
        )
    };
}

/// # Safety
///
/// `itr` must point to an initialized partition iterator.
#[no_mangle]
pub unsafe extern "C" fn pitr_next(itr: *mut PartIter) {
    if itr.is_null() || unsafe { (*itr).rows < 0 } {
        return;
    }
    let p = unsafe { (*itr).part };
    let outer = unsafe { (*itr).outer };
    let inner = unsafe { (*itr).inner };
    let rows = unsafe { (*itr).rows };
    let opt = unsafe { (*itr).opt };
    let use_outer = opt & PITR_USE_OUTER != 0;
    let use_inner = opt & PITR_USE_INNER != 0;
    let use_size = opt & PITR_USE_SIZE != 0;

    let mut outer_row = rows;
    let mut size = 0;
    let mut inner_sz = 0;
    let mut outer_sz = 0;

    let mut r = unsafe { (*itr).length } - 1;
    while r >= 0 {
        if use_size {
            size += unsafe { part_entry(p, r) };
        }
        if use_size && use_inner {
            inner_sz += unsafe { part_entry(inner, r) };
        }

        let mut c = unsafe { part_entry(p, r) } - 1;
        if use_inner && c < unsafe { part_entry(inner, r) } {
            r -= 1;
            continue;
        }

        if use_size && use_outer {
            while outer_row > 0 && unsafe { part_entry(outer, outer_row - 1) } < c {
                outer_row -= 1;
                outer_sz += unsafe { part_entry(outer, outer_row) };
            }
        }

        if use_size && size > c.saturating_mul(outer_row - r).saturating_add(outer_sz) {
            r -= 1;
            continue;
        }

        if c == 0 {
            unsafe {
                *ivector_data_mut(p).add(r as usize) = 0;
                (*itr).length = r;
            }
            return;
        }

        unsafe {
            (*itr).length = rows;
        }
        let mut rr = r;
        while rr < outer_row {
            if !use_size && use_outer && c > unsafe { part_entry(outer, rr) } {
                break;
            }
            if use_size {
                let mut avail = size;
                if use_inner {
                    inner_sz -= unsafe { part_entry(inner, rr) };
                    avail -= inner_sz;
                }
                if avail == 0 {
                    break;
                }
                c = c.min(avail);
                size -= c;
            }
            unsafe {
                *ivector_data_mut(p).add(rr as usize) = c;
            }
            rr += 1;
        }
        if use_outer {
            while rr < rows {
                c = unsafe { part_entry(outer, rr) };
                if use_size {
                    let mut avail = size;
                    if use_inner {
                        inner_sz -= unsafe { part_entry(inner, rr) };
                        avail -= inner_sz;
                    }
                    if avail == 0 {
                        break;
                    }
                    c = c.min(avail);
                    size -= c;
                }
                unsafe {
                    *ivector_data_mut(p).add(rr as usize) = c;
                }
                rr += 1;
            }
        }
        let old_length = unsafe { (*itr).length };
        let mut j = rr;
        while j < old_length {
            unsafe {
                *ivector_data_mut(p).add(j as usize) = 0;
            }
            j += 1;
        }
        unsafe {
            (*itr).length = rr;
        }
        return;
    }
    unsafe {
        (*itr).rows = -1;
    }
}

/// # Safety
///
/// Non-null vector pointers must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn lrit_new(
    outer: *const IVector,
    inner: *const IVector,
    content: *const IVector,
    maxrows: c_int,
    maxcols: c_int,
    partsz: c_int,
) -> *mut LrTabIter {
    if outer.is_null() {
        return ptr::null_mut();
    }
    let outer_values = unsafe { ivector_values(outer) };
    let inner_values = if inner.is_null() {
        &[][..]
    } else {
        unsafe { ivector_values(inner) }
    };
    let content_values = if content.is_null() {
        &[][..]
    } else {
        unsafe { ivector_values(content) }
    };
    if !abi_valid_partition(outer_values)
        || (!inner.is_null() && !abi_valid_partition(inner_values))
        || (!content.is_null() && !abi_valid_partition(content_values))
    {
        return ptr::null_mut();
    }
    if !inner.is_null() && !abi_partition_leq(inner_values, outer_values) {
        return unsafe { lrit_empty() };
    }

    let len = abi_part_length(outer_values);
    let ilen = if inner.is_null() {
        0
    } else {
        unsafe { (*inner).length as usize }.min(len)
    };
    let clen = if content.is_null() {
        0
    } else {
        abi_part_length(content_values)
    };
    let out0 = if len == 0 { 0 } else { outer_values[0] };

    let mut size = 0i32;
    let mut maxdepth = match i32::try_from(clen) {
        Ok(value) => value,
        Err(_) => return ptr::null_mut(),
    };
    for row in 0..len {
        let inn_r = if row < ilen { inner_values[row] } else { 0 };
        let rowsz = outer_values[row] - inn_r;
        if rowsz < 0 {
            return unsafe { lrit_empty() };
        }
        size = match size.checked_add(rowsz) {
            Some(value) => value,
            None => return ptr::null_mut(),
        };
        if rowsz > 0 {
            maxdepth += 1;
        }
    }
    let mut maxrows = if maxrows < 0 || maxrows > maxdepth {
        maxdepth
    } else {
        maxrows
    };

    let mut array_len = match usize::try_from(size) {
        Ok(value) => value.saturating_add(2),
        Err(_) => return ptr::null_mut(),
    };
    if maxcols >= 0 {
        let clim = maxcols - out0;
        let mut c1 = 0;
        for row in (0..clen).rev() {
            let c0 = content_values[row];
            if c1 < c0 && c1 < maxcols && c0 > clim {
                array_len = match array_len.checked_add(1) {
                    Some(value) => value,
                    None => return ptr::null_mut(),
                };
            }
            c1 = c0;
        }
        if c1 >= maxcols {
            array_len = array_len.saturating_sub(1);
        }
    }

    let lrit = unsafe { lrit_alloc(array_len) };
    if lrit.is_null() {
        return ptr::null_mut();
    }
    if partsz < maxrows {
        maxrows = maxrows.max(0);
    }
    let cont_len_i32 = partsz.max(maxrows);
    let Ok(cont_len) = u32::try_from(cont_len_i32) else {
        unsafe { libc::free(lrit.cast::<c_void>()) };
        return ptr::null_mut();
    };
    let cont = iv_new(cont_len);
    if cont.is_null() {
        unsafe { libc::free(lrit.cast::<c_void>()) };
        return ptr::null_mut();
    }
    unsafe {
        (*lrit).cont = cont;
        (*lrit).size = -1;
    }

    if maxrows < clen as i32 {
        return lrit;
    }
    {
        let cont_values_mut = unsafe { ivector_values_mut(cont) };
        for row in 0..clen {
            cont_values_mut[row] = content_values[row];
        }
        for value in cont_values_mut.iter_mut().skip(clen) {
            *value = 0;
        }
    }
    if maxcols >= 0 && clen > 0 && content_values[0] > maxcols {
        return lrit;
    }
    if maxcols >= 0 && out0 > maxcols {
        return lrit;
    }

    let size_usize = match usize::try_from(size) {
        Ok(value) => value,
        Err(_) => {
            unsafe { lrit_free(lrit) };
            return ptr::null_mut();
        }
    };
    let array = unsafe { lrit_array_mut(lrit) };
    let mut s = 0usize;
    let mut out1 = 0;
    let mut inn0 = if len == 0 {
        out0
    } else if len <= ilen {
        inner_values[len - 1]
    } else {
        0
    };
    for row in (0..len).rev() {
        let out2 = out1;
        let inn1 = inn0;
        out1 = outer_values[row];
        inn0 = if row == 0 {
            out0
        } else if row <= ilen {
            inner_values[row - 1]
        } else {
            0
        };
        if inn1 < out1 {
            maxdepth -= 1;
        }
        for col in inn1..out1 {
            let right = if col + 1 < out1 { s + 1 } else { array_len - 1 };
            let above = if col >= inn0 {
                let delta = match usize::try_from(out1 - inn0) {
                    Ok(value) => value,
                    Err(_) => {
                        unsafe { lrit_free(lrit) };
                        return ptr::null_mut();
                    }
                };
                s + delta
            } else {
                size_usize
            };
            let max = if col < out2 {
                let source = match usize::try_from(s as i64 - i64::from(out2) + i64::from(inn1)) {
                    Ok(value) => value,
                    Err(_) => {
                        unsafe { lrit_free(lrit) };
                        return ptr::null_mut();
                    }
                };
                array[source].max - 1
            } else {
                maxrows - 1
            }
            .min(maxdepth);
            let Ok(above) = c_int::try_from(above) else {
                unsafe { lrit_free(lrit) };
                return ptr::null_mut();
            };
            let Ok(right) = c_int::try_from(right) else {
                unsafe { lrit_free(lrit) };
                return ptr::null_mut();
            };
            array[s] = LritBox {
                value: 0,
                max,
                above,
                right,
            };
            s += 1;
        }
    }

    array[array_len - 1].value = maxrows - 1;
    array[size_usize].value = -1;
    if maxcols >= 0 {
        let clim = maxcols - out0;
        let mut c1 = 0;
        let mut extra = array_len.saturating_sub(2);
        let mut col = out0;
        for row in (0..clen).rev() {
            let c0 = content_values[row];
            if c1 < c0 && c1 < maxcols && c0 > clim {
                array[extra].value = row as i32;
                while col > maxcols - c0 && col > 0 {
                    col -= 1;
                    let index_i64 = i64::from(size) - i64::from(out0) + i64::from(col);
                    let Ok(index) = usize::try_from(index_i64) else {
                        unsafe { lrit_free(lrit) };
                        return ptr::null_mut();
                    };
                    let Ok(extra_i32) = c_int::try_from(extra) else {
                        unsafe { lrit_free(lrit) };
                        return ptr::null_mut();
                    };
                    array[index].above = extra_i32;
                }
                extra = extra.saturating_sub(1);
            }
            c1 = c0;
        }
    }

    for index in (0..size_usize).rev() {
        let above = match usize::try_from(array[index].above) {
            Ok(value) => value,
            Err(_) => {
                unsafe { lrit_free(lrit) };
                return ptr::null_mut();
            }
        };
        let x = array[above].value + 1;
        if x > array[index].max {
            return lrit;
        }
        array[index].value = x;
        let Ok(x_index) = usize::try_from(x) else {
            unsafe { lrit_free(lrit) };
            return ptr::null_mut();
        };
        unsafe {
            let cont_values_mut = ivector_values_mut(cont);
            cont_values_mut[x_index] += 1;
        }
    }

    unsafe {
        (*lrit).size = size;
    }
    lrit
}

/// # Safety
///
/// `lrit` must be null or a pointer returned by `lrit_new`.
#[no_mangle]
pub unsafe extern "C" fn lrit_free(lrit: *mut LrTabIter) {
    if lrit.is_null() {
        return;
    }
    unsafe {
        iv_free((*lrit).cont);
        libc::free(lrit.cast::<c_void>());
    }
}

/// # Safety
///
/// `lrit` must point to a valid LR tableau iterator.
#[no_mangle]
pub unsafe extern "C" fn lrit_good(lrit: *const LrTabIter) -> c_int {
    if lrit.is_null() {
        0
    } else if unsafe { (*lrit).size >= 0 } {
        1
    } else {
        0
    }
}

/// # Safety
///
/// `lrit` must point to a valid LR tableau iterator.
#[no_mangle]
pub unsafe extern "C" fn lrit_next(lrit: *mut LrTabIter) {
    if lrit.is_null() || unsafe { (*lrit).size < 0 } {
        return;
    }
    let cont = unsafe { (*lrit).cont };
    if cont.is_null() {
        unsafe {
            (*lrit).size = -1;
        }
        return;
    }

    let size = unsafe { (*lrit).size as usize };
    let cont_values = unsafe { ivector_values_mut(cont) };
    let array = unsafe { lrit_array_mut(lrit) };
    for index in 0..size {
        let right = match usize::try_from(array[index].right) {
            Ok(value) => value,
            Err(_) => {
                unsafe { (*lrit).size = -1 };
                return;
            }
        };
        let mut max = array[right].value;
        if max > array[index].max {
            max = array[index].max;
        }
        let mut x = array[index].value;
        let Ok(x_index) = usize::try_from(x) else {
            unsafe { (*lrit).size = -1 };
            return;
        };
        cont_values[x_index] -= 1;
        x += 1;
        while x <= max {
            let Ok(current) = usize::try_from(x) else {
                unsafe { (*lrit).size = -1 };
                return;
            };
            if current == 0 || cont_values[current] != cont_values[current - 1] {
                break;
            }
            x += 1;
        }
        if x > max {
            continue;
        }

        array[index].value = x;
        let Ok(x_index) = usize::try_from(x) else {
            unsafe { (*lrit).size = -1 };
            return;
        };
        cont_values[x_index] += 1;
        let mut fill = index;
        while fill != 0 {
            fill -= 1;
            let above = match usize::try_from(array[fill].above) {
                Ok(value) => value,
                Err(_) => {
                    unsafe { (*lrit).size = -1 };
                    return;
                }
            };
            let x = array[above].value + 1;
            array[fill].value = x;
            let Ok(x_index) = usize::try_from(x) else {
                unsafe { (*lrit).size = -1 };
                return;
            };
            cont_values[x_index] += 1;
        }
        return;
    }

    unsafe {
        (*lrit).size = -1;
    }
}

/// # Safety
///
/// `lrit` must point to a valid LR tableau iterator.
#[no_mangle]
pub unsafe extern "C" fn lrit_count(lrit: *mut LrTabIter) -> *mut IvLinComb {
    if lrit.is_null() {
        return ptr::null_mut();
    }
    let lc = ivlc_new(IVLC_HASHTABLE_SZ, IVLC_ARRAY_SZ);
    if lc.is_null() {
        return ptr::null_mut();
    }
    while unsafe { lrit_good(lrit) } != 0 {
        let cont = unsafe { (*lrit).cont };
        if cont.is_null()
            || unsafe { ivlc_add_element(lc, 1, cont, iv_hash(cont) as u32, LC_COPY_KEY) } != 0
        {
            unsafe { ivlc_free_all(lc) };
            return ptr::null_mut();
        }
        unsafe { lrit_next(lrit) };
    }
    lc
}

/// # Safety
///
/// Vector pointers must be valid as for `lrit_new`.
#[no_mangle]
pub unsafe extern "C" fn lrit_expand(
    outer: *const IVector,
    inner: *const IVector,
    content: *const IVector,
    maxrows: c_int,
    maxcols: c_int,
    partsz: c_int,
) -> *mut IvLinComb {
    let lrit = unsafe { lrit_new(outer, inner, content, maxrows, maxcols, partsz) };
    if lrit.is_null() {
        return ptr::null_mut();
    }
    let lc = unsafe { lrit_count(lrit) };
    unsafe { lrit_free(lrit) };
    lc
}

/// # Safety
///
/// `lrit` must point to a valid LR tableau iterator.
#[no_mangle]
pub unsafe extern "C" fn lrit_print_skewtab(
    lrit: *mut LrTabIter,
    outer: *const IVector,
    inner: *const IVector,
) {
    if lrit.is_null() || unsafe { (*lrit).size < 0 } {
        return;
    }
    let size = unsafe { (*lrit).size as usize };
    let array = unsafe { lrit_array_mut(lrit) };
    let outer_values = if outer.is_null() {
        &[][..]
    } else {
        unsafe { ivector_values(outer) }
    };
    let inner_values = if inner.is_null() {
        &[][..]
    } else {
        unsafe { ivector_values(inner) }
    };
    let ilen = inner_values.len();
    let mut len = abi_part_length(outer_values);
    if len <= ilen {
        while len > 0 && inner_values[len - 1] == outer_values[len - 1] {
            len -= 1;
        }
    }
    if len == 0 {
        return;
    }
    let col_first = if ilen < len { 0 } else { inner_values[len - 1] };
    let mut row = 0;
    while row < ilen && inner_values[row] == outer_values[row] {
        row += 1;
    }
    let mut pos = size;
    while row < len {
        let inner_part = inner_values.get(row).copied().unwrap_or(0).max(0) as usize;
        let outer_part = outer_values[row].max(0) as usize;
        let row_size = outer_part.saturating_sub(inner_part);
        pos = pos.saturating_sub(row_size);
        for _ in col_first.max(0) as usize..inner_part {
            print!("  ");
        }
        for col in 0..row_size {
            print!("{:2}", array[pos + col].value);
        }
        println!();
        row += 1;
    }
}

/// # Safety
///
/// `lrit` must point to a valid LR tableau iterator.
#[no_mangle]
pub unsafe extern "C" fn lrit_dump(lrit: *mut LrTabIter) {
    if lrit.is_null() {
        return;
    }
    println!("size={} array_len={}", unsafe { (*lrit).size }, unsafe {
        (*lrit).array_len
    });
    let array_len = unsafe { (*lrit).array_len.max(0) as usize };
    let array = unsafe { lrit_array_mut(lrit) };
    for (index, entry) in array.iter().take(array_len).enumerate() {
        println!(
            "{index:02}: value={} max={} right={} above={}",
            entry.value, entry.max, entry.right, entry.above
        );
    }
}

/// # Safety
///
/// `lrit` must point to a valid LR tableau iterator.
#[no_mangle]
pub unsafe extern "C" fn lrit_dump_skew(
    lrit: *mut LrTabIter,
    _outer: *const IVector,
    _inner: *const IVector,
) {
    unsafe { lrit_dump(lrit) };
}

/// # Safety
///
/// `ht` must point to writable storage for an `IvLinComb`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_init(ht: *mut IvLinComb, tabsz: u32, eltsz: u32) -> c_int {
    if ht.is_null() {
        return -1;
    }
    let table_sz = tabsz.max(1);
    let elts_sz = eltsz.max(1);
    let Some(table_bytes) = ivlc_table_alloc_size(table_sz) else {
        return -1;
    };
    let Some(elts_bytes) = ivlc_elts_alloc_size(elts_sz) else {
        return -1;
    };
    let table = unsafe { libc::calloc(table_sz as usize, mem::size_of::<u32>()) }.cast::<u32>();
    if table.is_null() {
        return -1;
    }
    let elts = unsafe { libc::malloc(elts_bytes) }.cast::<IvlcKeyVal>();
    if elts.is_null() {
        unsafe { libc::free(table.cast::<c_void>()) };
        return -1;
    }
    debug_assert_eq!(table_bytes, table_sz as usize * mem::size_of::<u32>());
    unsafe {
        *ht = IvLinComb {
            table,
            elts,
            card: 0,
            free_elts: 0,
            elts_len: 1,
            elts_sz,
            table_sz,
        };
    }
    0
}

#[no_mangle]
pub extern "C" fn ivlc_new(tabsz: u32, eltsz: u32) -> *mut IvLinComb {
    let ht = unsafe { libc::malloc(mem::size_of::<IvLinComb>()) }.cast::<IvLinComb>();
    if ht.is_null() {
        return ptr::null_mut();
    }
    if unsafe { ivlc_init(ht, tabsz, eltsz) } != 0 {
        unsafe { libc::free(ht.cast::<c_void>()) };
        return ptr::null_mut();
    }
    ht
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid `IvLinComb`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_card(ht: *const IvLinComb) -> u32 {
    if ht.is_null() {
        0
    } else {
        unsafe { (*ht).card }
    }
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid `IvLinComb`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_dealloc(ht: *mut IvLinComb) {
    if ht.is_null() {
        return;
    }
    unsafe {
        libc::free((*ht).table.cast::<c_void>());
        libc::free((*ht).elts.cast::<c_void>());
        (*ht).table = ptr::null_mut();
        (*ht).elts = ptr::null_mut();
        (*ht).card = 0;
        (*ht).free_elts = 0;
        (*ht).elts_len = 0;
        (*ht).elts_sz = 0;
        (*ht).table_sz = 0;
    }
}

/// # Safety
///
/// `ht` must be null or a pointer returned by `ivlc_new`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_free(ht: *mut IvLinComb) {
    if ht.is_null() {
        return;
    }
    unsafe {
        ivlc_dealloc(ht);
        libc::free(ht.cast::<c_void>());
    }
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid `IvLinComb`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_reset(ht: *mut IvLinComb) {
    if ht.is_null() {
        return;
    }
    unsafe {
        ptr::write_bytes((*ht).table, 0, (*ht).table_sz as usize);
        (*ht).card = 0;
        (*ht).free_elts = 0;
        (*ht).elts_len = 1;
    }
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid `IvLinComb`.
#[no_mangle]
pub unsafe extern "C" fn ivlc__grow_table(ht: *mut IvLinComb, sz: u32) -> c_int {
    if ht.is_null() {
        return -1;
    }
    let Some(new_sz) = ivlc_new_table_size(sz) else {
        return -1;
    };
    let new_table = unsafe { libc::calloc(new_sz as usize, mem::size_of::<u32>()) }.cast::<u32>();
    if new_table.is_null() {
        return -1;
    }
    unsafe {
        for old_index in 0..(*ht).table_sz {
            let mut i = *(*ht).table.add(old_index as usize);
            while i != 0 {
                let next = (*(*ht).elts.add(i as usize)).next;
                let new_index = (*(*ht).elts.add(i as usize)).hash % new_sz;
                (*(*ht).elts.add(i as usize)).next = *new_table.add(new_index as usize);
                *new_table.add(new_index as usize) = i;
                i = next;
            }
        }
        libc::free((*ht).table.cast::<c_void>());
        (*ht).table = new_table;
        (*ht).table_sz = new_sz;
    }
    0
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid `IvLinComb`.
#[no_mangle]
pub unsafe extern "C" fn ivlc__grow_elts(ht: *mut IvLinComb, sz: u32) -> c_int {
    if ht.is_null() {
        return -1;
    }
    let Some(new_sz) = sz.checked_mul(2) else {
        return -1;
    };
    let Some(new_bytes) = ivlc_elts_alloc_size(new_sz) else {
        return -1;
    };
    let elts =
        unsafe { libc::realloc((*ht).elts.cast::<c_void>(), new_bytes) }.cast::<IvlcKeyVal>();
    if elts.is_null() {
        return -1;
    }
    unsafe {
        (*ht).elts = elts;
        (*ht).elts_sz = new_sz;
    }
    0
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid `IvLinComb`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_makeroom(ht: *mut IvLinComb, sz: u32) -> c_int {
    unsafe { ivlc_makeroom_inner(ht, sz) }
}

/// # Safety
///
/// Pointers must refer to a valid linear combination and vector key.
#[no_mangle]
pub unsafe extern "C" fn ivlc_lookup(
    ht: *mut IvLinComb,
    key: *const IVector,
    hash: u32,
) -> *mut IvlcKeyVal {
    let index = unsafe { ivlc_lookup_index(ht, key, hash) };
    if index == 0 {
        ptr::null_mut()
    } else {
        unsafe { (*ht).elts.add(index as usize) }
    }
}

/// # Safety
///
/// `ht` must be a valid linear combination and `key` must be a valid vector.
#[no_mangle]
pub unsafe extern "C" fn ivlc_insert(
    ht: *mut IvLinComb,
    key: *mut IVector,
    hash: u32,
    value: i32,
) -> *mut IvlcKeyVal {
    if ht.is_null() || key.is_null() {
        return ptr::null_mut();
    }
    let Some(new_card) = unsafe { (*ht).card }.checked_add(1) else {
        return ptr::null_mut();
    };
    if unsafe { ivlc_makeroom_inner(ht, new_card) } != 0 {
        return ptr::null_mut();
    }
    let i = unsafe {
        if (*ht).free_elts != 0 {
            let i = (*ht).free_elts;
            (*ht).free_elts = (*(*ht).elts.add(i as usize)).next;
            i
        } else {
            let i = (*ht).elts_len;
            (*ht).elts_len += 1;
            i
        }
    };
    unsafe {
        (*ht).card = new_card;
        let table_index = ivlc_table_index(hash, (*ht).table_sz);
        let kv = (*ht).elts.add(i as usize);
        (*kv).key = key;
        (*kv).value = value;
        (*kv).hash = hash;
        (*kv).next = *(*ht).table.add(table_index as usize);
        *(*ht).table.add(table_index as usize) = i;
        kv
    }
}

/// # Safety
///
/// Pointers must refer to a valid linear combination and vector key.
#[no_mangle]
pub unsafe extern "C" fn ivlc_remove(
    ht: *mut IvLinComb,
    key: *const IVector,
    hash: u32,
) -> *mut IvlcKeyVal {
    if ht.is_null() || key.is_null() || unsafe { (*ht).table_sz } == 0 {
        return ptr::null_mut();
    }
    unsafe {
        let table_index = ivlc_table_index(hash, (*ht).table_sz);
        let mut link = (*ht).table.add(table_index as usize);
        let mut i = *link;
        while i != 0 {
            let kv = (*ht).elts.add(i as usize);
            if iv_cmp(key, (*kv).key) == 0 {
                *link = (*kv).next;
                (*kv).next = (*ht).free_elts;
                (*ht).free_elts = i;
                (*ht).card -= 1;
                return kv;
            }
            link = ptr::addr_of_mut!((*kv).next);
            i = *link;
        }
    }
    ptr::null_mut()
}

/// # Safety
///
/// `itr` must point to an iterator initialized by `ivlc_first`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_good(itr: *const IvlcIter) -> c_int {
    if itr.is_null() || unsafe { (*itr).i } == 0 {
        0
    } else {
        1
    }
}

/// # Safety
///
/// `ht` must be a valid linear combination and `itr` writable iterator storage.
#[no_mangle]
pub unsafe extern "C" fn ivlc_first(ht: *mut IvLinComb, itr: *mut IvlcIter) {
    if itr.is_null() {
        return;
    }
    unsafe {
        (*itr).ht = ht;
        (*itr).index = 0;
        (*itr).i = 0;
        if ht.is_null() {
            return;
        }
        let mut index = 0;
        while index < (*ht).table_sz && *(*ht).table.add(index as usize) == 0 {
            index += 1;
        }
        if index == (*ht).table_sz {
            return;
        }
        (*itr).index = index;
        (*itr).i = *(*ht).table.add(index as usize);
    }
}

/// # Safety
///
/// `itr` must point to an iterator initialized by `ivlc_first`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_next(itr: *mut IvlcIter) {
    if itr.is_null() || unsafe { (*itr).ht }.is_null() || unsafe { (*itr).i } == 0 {
        return;
    }
    unsafe {
        let ht = (*itr).ht;
        let current = (*ht).elts.add((*itr).i as usize);
        if (*current).next != 0 {
            (*itr).i = (*current).next;
            return;
        }
        let mut index = (*itr).index + 1;
        while index < (*ht).table_sz && *(*ht).table.add(index as usize) == 0 {
            index += 1;
        }
        if index == (*ht).table_sz {
            (*itr).i = 0;
            return;
        }
        (*itr).index = index;
        (*itr).i = *(*ht).table.add(index as usize);
    }
}

/// # Safety
///
/// `itr` must point to a good iterator.
#[no_mangle]
pub unsafe extern "C" fn ivlc_key(itr: *const IvlcIter) -> *mut IVector {
    if unsafe { ivlc_good(itr) } == 0 {
        return ptr::null_mut();
    }
    unsafe { (*(*(*itr).ht).elts.add((*itr).i as usize)).key }
}

/// # Safety
///
/// `itr` must point to a good iterator.
#[no_mangle]
pub unsafe extern "C" fn ivlc_value(itr: *const IvlcIter) -> i32 {
    if unsafe { ivlc_good(itr) } == 0 {
        return 0;
    }
    unsafe { (*(*(*itr).ht).elts.add((*itr).i as usize)).value }
}

/// # Safety
///
/// `itr` must point to a good iterator.
#[no_mangle]
pub unsafe extern "C" fn ivlc_keyval(itr: *const IvlcIter) -> *mut IvlcKeyVal {
    if unsafe { ivlc_good(itr) } == 0 {
        return ptr::null_mut();
    }
    unsafe { (*(*itr).ht).elts.add((*itr).i as usize) }
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid `IvLinComb` whose keys are
/// owned by the table.
#[no_mangle]
pub unsafe extern "C" fn ivlc_dealloc_refs(ht: *mut IvLinComb) {
    if ht.is_null() {
        return;
    }
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(ht, &mut itr);
        while ivlc_good(&itr) != 0 {
            let key = ivlc_key(&itr);
            iv_free(key);
            ivlc_next(&mut itr);
        }
    }
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid `IvLinComb`.
#[no_mangle]
pub unsafe extern "C" fn ivlc_dealloc_all(ht: *mut IvLinComb) {
    if ht.is_null() {
        return;
    }
    unsafe {
        ivlc_dealloc_refs(ht);
        ivlc_dealloc(ht);
    }
}

/// # Safety
///
/// `ht` must be null or a pointer returned by `ivlc_new`; keys must be owned by
/// the table.
#[no_mangle]
pub unsafe extern "C" fn ivlc_free_all(ht: *mut IvLinComb) {
    if ht.is_null() {
        return;
    }
    unsafe {
        ivlc_dealloc_all(ht);
        libc::free(ht.cast::<c_void>());
    }
}

/// # Safety
///
/// `ht` must be a valid linear combination. `key` must be valid; ownership is
/// copied when `LC_COPY_KEY` is set and otherwise transferred.
#[no_mangle]
pub unsafe extern "C" fn ivlc_add_element(
    ht: *mut IvLinComb,
    c: i32,
    key: *mut IVector,
    hash: u32,
    opt: c_int,
) -> c_int {
    if ht.is_null() || key.is_null() {
        return -1;
    }
    if c == 0 {
        if opt & LC_COPY_KEY == 0 {
            unsafe { iv_free(key) };
        }
        return 0;
    }

    let existing = unsafe { ivlc_lookup(ht, key, hash) };
    if !existing.is_null() {
        if opt & LC_COPY_KEY == 0 {
            unsafe { iv_free(key) };
        }
        unsafe {
            (*existing).value = (*existing).value.wrapping_add(c);
            if (*existing).value == 0 && opt & LC_FREE_ZERO != 0 {
                let existing_key = (*existing).key;
                let existing_hash = (*existing).hash;
                ivlc_remove(ht, existing_key, existing_hash);
                iv_free(existing_key);
            }
        }
        return 0;
    }

    let Some(new_card) = unsafe { (*ht).card }.checked_add(1) else {
        if opt & LC_COPY_KEY == 0 {
            unsafe { iv_free(key) };
        }
        return -1;
    };
    if unsafe { ivlc_makeroom_inner(ht, new_card) } != 0 {
        if opt & LC_COPY_KEY == 0 {
            unsafe { iv_free(key) };
        }
        return -1;
    }
    let stored_key = if opt & LC_COPY_KEY != 0 {
        let copied = unsafe { iv_new_copy(key) };
        if copied.is_null() {
            return -1;
        }
        copied
    } else {
        key
    };
    if unsafe { ivlc_insert(ht, stored_key, hash, c) }.is_null() {
        unsafe { iv_free(stored_key) };
        return -1;
    }
    0
}

/// # Safety
///
/// `dst` and `src` must be valid linear combinations.
#[no_mangle]
pub unsafe extern "C" fn ivlc_add_multiple(
    dst: *mut IvLinComb,
    c: i32,
    src: *mut IvLinComb,
    opt: c_int,
) -> c_int {
    if dst.is_null() || src.is_null() {
        return -1;
    }
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(src, &mut itr);
        while ivlc_good(&itr) != 0 {
            let kv = ivlc_keyval(&itr);
            let value = c.wrapping_mul((*kv).value);
            if ivlc_add_element(dst, value, (*kv).key, (*kv).hash, opt) != 0 {
                return -1;
            }
            ivlc_next(&mut itr);
        }
    }
    0
}

/// # Safety
///
/// `ht1` and `ht2` must point to valid linear combinations.
#[no_mangle]
pub unsafe extern "C" fn ivlc_equals(
    ht1: *mut IvLinComb,
    ht2: *mut IvLinComb,
    opt_zero: c_int,
) -> c_int {
    if ht1.is_null() || ht2.is_null() {
        return 0;
    }
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(ht1, &mut itr);
        while ivlc_good(&itr) != 0 {
            let kv = ivlc_keyval(&itr);
            if (*kv).value != 0 || opt_zero != 0 {
                let other = ivlc_lookup(ht2, (*kv).key, (*kv).hash);
                if other.is_null() || (*other).value != (*kv).value {
                    return 0;
                }
            }
            ivlc_next(&mut itr);
        }
        ivlc_first(ht2, &mut itr);
        while ivlc_good(&itr) != 0 {
            let kv = ivlc_keyval(&itr);
            if (*kv).value != 0 || opt_zero != 0 {
                let other = ivlc_lookup(ht1, (*kv).key, (*kv).hash);
                if other.is_null() || (*other).value != (*kv).value {
                    return 0;
                }
            }
            ivlc_next(&mut itr);
        }
    }
    1
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid linear combination.
#[no_mangle]
pub unsafe extern "C" fn ivlc_print(ht: *mut IvLinComb, opt_zero: c_int) {
    if ht.is_null() {
        return;
    }
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(ht, &mut itr);
        while ivlc_good(&itr) != 0 {
            let value = ivlc_value(&itr);
            if value != 0 || opt_zero != 0 {
                print!("{value}  ");
                iv_print(ivlc_key(&itr));
                println!();
            }
            ivlc_next(&mut itr);
        }
    }
}

/// # Safety
///
/// If `ht` is non-null, it must point to a valid linear combination.
#[no_mangle]
pub unsafe extern "C" fn ivlc_print_stat(ht: *mut IvLinComb) {
    if ht.is_null() {
        return;
    }
    let range = 20usize;
    let mut stat = vec![0usize; range];
    let mut used = 0usize;
    let mut compares = 0usize;
    unsafe {
        for index in 0..(*ht).table_sz as usize {
            let mut i = *(*ht).table.add(index);
            if i == 0 {
                continue;
            }
            used += 1;
            let mut count = 0usize;
            while i != 0 {
                count += 1;
                i = (*(*ht).elts.add(i as usize)).next;
            }
            compares += (count + 1) * count / 2;
            let bucket = count.min(range).saturating_sub(1);
            stat[bucket] += count;
        }
        println!("Hash table size: {}", (*ht).table_sz);
        println!("Hash table used: {used}");
        println!("Total elements: {}", (*ht).card);
        if (*ht).card != 0 {
            println!(
                "Average compares: {}",
                compares as f64 / f64::from((*ht).card)
            );
        }
    }
    print!("Table distribution:");
    for value in stat {
        print!(" {value}");
    }
    println!();
}

/// # Safety
///
/// If `lc` is non-null, it must point to a valid linear combination.
#[no_mangle]
pub unsafe extern "C" fn part_print_lincomb(lc: *mut IvLinComb) {
    if lc.is_null() {
        return;
    }
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(lc, &mut itr);
        while ivlc_good(&itr) != 0 {
            let value = ivlc_value(&itr);
            if value != 0 {
                print!("{value}  ");
                part_printnl(ivlc_key(&itr));
            }
            ivlc_next(&mut itr);
        }
    }
}

/// # Safety
///
/// If `lc` is non-null, it must point to a valid linear combination.
#[no_mangle]
pub unsafe extern "C" fn part_qprint_lincomb(lc: *mut IvLinComb, level: c_int) {
    if lc.is_null() {
        return;
    }
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(lc, &mut itr);
        while ivlc_good(&itr) != 0 {
            let value = ivlc_value(&itr);
            if value != 0 {
                print!("{value}  ");
                part_qprintnl(ivlc_key(&itr), level);
            }
            ivlc_next(&mut itr);
        }
    }
}

unsafe fn maple_term(c: c_int, v: *const IVector, letter: *const c_char, nz: bool) {
    let sign = if c < 0 { '-' } else { '+' };
    let coeff = c.wrapping_abs();
    let letter = if letter.is_null() {
        "s".into()
    } else {
        unsafe { CStr::from_ptr(letter) }.to_string_lossy()
    };
    print!("{sign}{coeff}*{letter}[");
    if !v.is_null() {
        for (index, &value) in unsafe { ivector_values(v) }.iter().enumerate() {
            if nz && value == 0 {
                break;
            }
            if index != 0 {
                print!(",");
            }
            print!("{value}");
        }
    }
    print!("]");
}

/// # Safety
///
/// `ht` must point to a valid linear combination; `letter` must be null or a
/// valid NUL-terminated C string.
#[no_mangle]
pub unsafe extern "C" fn maple_print_lincomb(ht: *mut IvLinComb, letter: *const c_char, nz: c_int) {
    print!("0");
    if ht.is_null() {
        println!();
        return;
    }
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(ht, &mut itr);
        while ivlc_good(&itr) != 0 {
            let value = ivlc_value(&itr);
            if value != 0 {
                maple_term(value, ivlc_key(&itr), letter, nz != 0);
            }
            ivlc_next(&mut itr);
        }
    }
    println!();
}

/// # Safety
///
/// `lc` must point to a valid linear combination; `letter` must be null or a
/// valid NUL-terminated C string.
#[no_mangle]
pub unsafe extern "C" fn maple_qprint_lincomb(
    lc: *mut IvLinComb,
    level: c_int,
    letter: *const c_char,
) {
    print!("0");
    if lc.is_null() {
        println!();
        return;
    }
    let letter = if letter.is_null() {
        "s".into()
    } else {
        unsafe { CStr::from_ptr(letter) }.to_string_lossy()
    };
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(lc, &mut itr);
        while ivlc_good(&itr) != 0 {
            let value = ivlc_value(&itr);
            if value != 0 {
                let key = ivlc_key(&itr);
                let degree = part_qdegree(key, level);
                let coeff = value.wrapping_abs();
                let sign = if value < 0 { '-' } else { '+' };
                print!("{sign}{coeff}*q^{degree}*{letter}[");
                let values = ivector_values(key);
                for index in 0..values.len() {
                    let x = part_qentry(key, index as c_int, degree, level);
                    if x == 0 {
                        break;
                    }
                    if index != 0 {
                        print!(",");
                    }
                    print!("{x}");
                }
                print!("]");
            }
            ivlc_next(&mut itr);
        }
    }
    println!();
}

/// # Safety
///
/// `w` must point to a valid `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn trans(w: *const IVector, vars: c_int) -> *mut IvLinComb {
    if w.is_null() {
        return ptr::null_mut();
    }
    let values = unsafe { ivector_values(w) };
    let terms = match trans_polynomial(values, vars) {
        Ok(terms) => terms,
        Err(_) => return ptr::null_mut(),
    };
    unsafe { ivlc_from_i32_terms(&terms) }
}

/// # Safety
///
/// `slc` must point to a valid linear combination.
#[no_mangle]
pub unsafe extern "C" fn monk(i: c_int, slc: *mut IvLinComb, rank: c_int) -> *mut IvLinComb {
    if slc.is_null() {
        return ptr::null_mut();
    }
    let source = unsafe { ivlc_collect_i32_terms(slc) };
    let terms = match monk_product(i, &source, rank) {
        Ok(terms) => terms,
        Err(_) => return ptr::null_mut(),
    };
    unsafe { ivlc_from_i32_terms(&terms) }
}

/// # Safety
///
/// `poly` must point to a valid linear combination whose keys are owned by the
/// table. `perm` must point to a valid `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn mult_poly_schubert(
    poly: *mut IvLinComb,
    perm: *const IVector,
    rank: c_int,
) -> *mut IvLinComb {
    if poly.is_null() || perm.is_null() {
        return ptr::null_mut();
    }
    let old_keys = unsafe { ivlc_collect_owned_keys(poly) };
    let source = unsafe { ivlc_collect_i32_terms(poly) };
    unsafe {
        ivlc_reset(poly);
        for key in old_keys {
            iv_free(key);
        }
    }
    let perm_values = unsafe { ivector_values(perm) };
    let terms = match multiply_poly_schubert(&source, perm_values, rank) {
        Ok(terms) => terms,
        Err(_) => {
            unsafe { ivlc_free_all(poly) };
            return ptr::null_mut();
        }
    };
    if unsafe { ivlc_fill_i32_terms(poly, &terms) } != 0 {
        unsafe { ivlc_free_all(poly) };
        return ptr::null_mut();
    }
    poly
}

/// # Safety
///
/// `w1` and `w2` must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn mult_schubert(
    w1: *const IVector,
    w2: *const IVector,
    rank: c_int,
) -> *mut IvLinComb {
    if w1.is_null() || w2.is_null() {
        return ptr::null_mut();
    }
    let left = unsafe { ivector_values(w1) };
    let right = unsafe { ivector_values(w2) };
    let terms = match multiply_schubert(left, right, rank) {
        Ok(terms) => terms,
        Err(_) => return ptr::null_mut(),
    };
    unsafe { ivlc_from_i32_terms(&terms) }
}

/// # Safety
///
/// `str1` and `str2` must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn mult_schubert_str(
    str1: *const IVector,
    str2: *const IVector,
) -> *mut IvLinComb {
    if str1.is_null() || str2.is_null() {
        return ptr::null_mut();
    }
    let left = unsafe { ivector_values(str1) };
    let right = unsafe { ivector_values(str2) };
    let terms = match multiply_schubert_strings(left, right) {
        Ok(terms) => terms,
        Err(_) => return ptr::null_mut(),
    };
    unsafe { ivlc_from_i32_terms(&terms) }
}

/// # Safety
///
/// Non-null pointers must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn lrcoef_count(
    outer: *const IVector,
    inner: *const IVector,
    content: *const IVector,
) -> c_longlong {
    if outer.is_null() || inner.is_null() || content.is_null() {
        return -1;
    }
    let outer = unsafe { ivector_values(outer) };
    let inner = unsafe { ivector_values(inner) };
    let content = unsafe { ivector_values(content) };
    lrcoef_i64(outer, inner, content).unwrap_or(-1)
}

/// # Safety
///
/// Non-null pointers must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn schur_lrcoef(
    outer: *const IVector,
    inner1: *const IVector,
    inner2: *const IVector,
) -> c_longlong {
    if outer.is_null() || inner1.is_null() || inner2.is_null() {
        return -1;
    }
    let outer = unsafe { ivector_values(outer) };
    let inner1 = unsafe { ivector_values(inner1) };
    let inner2 = unsafe { ivector_values(inner2) };
    lrcoef_i64(outer, inner1, inner2).unwrap_or(-1)
}

/// # Safety
///
/// Non-null pointers must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn schur_mult(
    sh1: *const IVector,
    sh2: *const IVector,
    rows: c_int,
    cols: c_int,
    partsz: c_int,
) -> *mut IvLinComb {
    if sh1.is_null() || sh2.is_null() {
        return ptr::null_mut();
    }
    let sh1_values = unsafe { ivector_values(sh1) };
    let sh2_values = unsafe { ivector_values(sh2) };
    let terms = match schur_product_expansion(sh1_values, sh2_values, rows, cols) {
        Ok(terms) => terms,
        Err(_) => return ptr::null_mut(),
    };
    let key_len = default_product_key_len(sh1_values, sh2_values, rows, partsz);
    unsafe { ivlc_from_terms(&terms, key_len) }
}

/// # Safety
///
/// `la` and `tmp` must point to valid `IVector` allocations of the same length.
#[no_mangle]
pub unsafe extern "C" fn fusion_reduce(la: *mut IVector, level: c_int, tmp: *mut IVector) -> c_int {
    if la.is_null() || tmp.is_null() || unsafe { (*la).length != (*tmp).length } {
        return 0;
    }
    let values = unsafe { ivector_values(la) }.to_vec();
    let (reduced, sign) = match fusion_reduce_values(&values, level) {
        Ok(Some(reduced)) => reduced,
        Ok(None) | Err(_) => return 0,
    };
    unsafe {
        ivector_values_mut(la).copy_from_slice(&reduced);
    }
    sign
}

/// # Safety
///
/// `lc` must point to a valid linear combination whose keys are owned by the
/// table.
#[no_mangle]
pub unsafe extern "C" fn fusion_reduce_lc(lc: *mut IvLinComb, level: c_int) -> c_int {
    if lc.is_null() {
        return -1;
    }

    let mut terms = Vec::<(*mut IVector, i32)>::new();
    let mut itr = IvlcIter {
        ht: ptr::null_mut(),
        index: 0,
        i: 0,
    };
    unsafe {
        ivlc_first(lc, &mut itr);
        while ivlc_good(&itr) != 0 {
            terms.push((ivlc_key(&itr), ivlc_value(&itr)));
            ivlc_next(&mut itr);
        }
        ivlc_reset(lc);
    }

    while let Some((key, value)) = terms.pop() {
        let values = unsafe { ivector_values(key) }.to_vec();
        let reduced = match fusion_reduce_values(&values, level) {
            Ok(reduced) => reduced,
            Err(_) => {
                unsafe { iv_free(key) };
                for (remaining_key, _) in terms {
                    unsafe { iv_free(remaining_key) };
                }
                return -1;
            }
        };
        let coefficient = match reduced {
            Some((reduced, sign)) => {
                unsafe {
                    ivector_values_mut(key).copy_from_slice(&reduced);
                }
                match value.checked_mul(sign) {
                    Some(coefficient) => coefficient,
                    None => {
                        unsafe { iv_free(key) };
                        for (remaining_key, _) in terms {
                            unsafe { iv_free(remaining_key) };
                        }
                        return -1;
                    }
                }
            }
            None => 0,
        };
        let hash = unsafe { iv_hash(key) } as u32;
        if unsafe { ivlc_add_element(lc, coefficient, key, hash, LC_FREE_ZERO) } != 0 {
            for (remaining_key, _) in terms {
                unsafe { iv_free(remaining_key) };
            }
            return -1;
        }
    }

    0
}

/// # Safety
///
/// Non-null pointers must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn schur_mult_fusion(
    sh1: *const IVector,
    sh2: *const IVector,
    rows: c_int,
    level: c_int,
) -> *mut IvLinComb {
    if sh1.is_null() || sh2.is_null() {
        return ptr::null_mut();
    }
    let sh1_values = unsafe { ivector_values(sh1) };
    let sh2_values = unsafe { ivector_values(sh2) };
    let terms = match schur_product_fusion_expansion(sh1_values, sh2_values, rows, level) {
        Ok(terms) => terms,
        Err(_) => return ptr::null_mut(),
    };
    let key_len = usize::try_from(rows).unwrap_or(0);
    unsafe { ivlc_from_signed_terms(&terms, key_len) }
}

/// # Safety
///
/// Non-null pointers must point to valid `IVector` allocations.
#[no_mangle]
pub unsafe extern "C" fn schur_skew(
    outer: *const IVector,
    inner: *const IVector,
    rows: c_int,
    partsz: c_int,
) -> *mut IvLinComb {
    if outer.is_null() || inner.is_null() {
        return ptr::null_mut();
    }
    let outer_values = unsafe { ivector_values(outer) };
    let inner_values = unsafe { ivector_values(inner) };
    let key_len = default_skew_key_len(outer_values, inner_values, rows, partsz);
    unsafe { ivlc_from_skew_expansion_direct(outer_values, inner_values, rows, key_len) }
}

/// # Safety
///
/// `sh` must point to a valid `IVector` allocation.
#[no_mangle]
pub unsafe extern "C" fn schur_coprod(
    sh: *const IVector,
    rows: c_int,
    cols: c_int,
    partsz: c_int,
    all: c_int,
) -> *mut IvLinComb {
    if sh.is_null() {
        return ptr::null_mut();
    }
    let sh_values = unsafe { ivector_values(sh) };
    let terms = match schur_coproduct_expansion(sh_values, rows, cols, all != 0) {
        Ok(terms) => terms,
        Err(_) => return ptr::null_mut(),
    };
    let key_len = default_coproduct_key_len(sh_values, rows, partsz);
    unsafe { ivlc_from_terms(&terms, key_len) }
}

unsafe fn zero_skew_shape(ss: *mut SkewShapeAbi) {
    if !ss.is_null() {
        unsafe {
            (*ss).outer = ptr::null_mut();
            (*ss).inner = ptr::null_mut();
            (*ss).cont = ptr::null_mut();
            (*ss).sign = 0;
        }
    }
}

unsafe fn fill_skew_shape(
    ss: *mut SkewShapeAbi,
    outer: &[i32],
    inner: Option<&[i32]>,
    cont: Option<&[i32]>,
    sign: c_int,
) -> c_int {
    if ss.is_null() {
        return -1;
    }
    unsafe { zero_skew_shape(ss) };
    let outer_ptr = unsafe { ivector_from_partition(outer, outer.len()) };
    if outer_ptr.is_null() {
        return -1;
    }
    let inner_ptr = if let Some(inner) = inner {
        let ptr = unsafe { ivector_from_partition(inner, inner.len()) };
        if ptr.is_null() {
            unsafe { iv_free(outer_ptr) };
            return -1;
        }
        ptr
    } else {
        ptr::null_mut()
    };
    let cont_ptr = if let Some(cont) = cont {
        let ptr = unsafe { ivector_from_partition(cont, cont.len()) };
        if ptr.is_null() {
            unsafe {
                iv_free(outer_ptr);
                iv_free(inner_ptr);
            }
            return -1;
        }
        ptr
    } else {
        ptr::null_mut()
    };
    unsafe {
        (*ss).outer = outer_ptr;
        (*ss).inner = inner_ptr;
        (*ss).cont = cont_ptr;
        (*ss).sign = sign;
    }
    0
}

/// # Safety
///
/// If `ss` is non-null, it must point to a `skew_shape` allocated by the
/// caller.
#[no_mangle]
pub unsafe extern "C" fn sksh_dealloc(ss: *mut SkewShapeAbi) {
    if ss.is_null() {
        return;
    }
    unsafe {
        iv_free((*ss).outer);
        iv_free((*ss).inner);
        iv_free((*ss).cont);
        zero_skew_shape(ss);
    }
}

/// # Safety
///
/// Non-null pointers must point to valid vectors.
#[no_mangle]
pub unsafe extern "C" fn sksh_print(
    outer: *const IVector,
    inner: *const IVector,
    cont: *const IVector,
) {
    if outer.is_null() {
        return;
    }
    let outer_values = unsafe { ivector_values(outer) };
    let inner_values = if inner.is_null() {
        &[][..]
    } else {
        unsafe { ivector_values(inner) }
    };
    let cont_values = if cont.is_null() {
        &[][..]
    } else {
        unsafe { ivector_values(cont) }
    };
    let mut len = abi_part_length(outer_values);
    let mut ilen = inner_values.len();
    if len <= ilen {
        while len > 0 && inner_values[len - 1] == outer_values[len - 1] {
            len -= 1;
        }
        ilen = len;
    }
    let mut row0 = 0;
    while row0 < ilen && inner_values[row0] == outer_values[row0] {
        row0 += 1;
    }
    let left = if len == 0 || ilen < len {
        0
    } else {
        inner_values[len - 1]
    };
    let right = if len == 0 { 0 } else { outer_values[0] };
    for &part in cont_values.iter().take(abi_part_length(cont_values)) {
        for _ in left..right {
            print!(" ");
        }
        for _ in 0..part {
            print!("c");
        }
        println!();
    }
    for row in row0..len {
        let inn = inner_values.get(row).copied().unwrap_or(0);
        let out = outer_values[row];
        for _ in 0..inn {
            print!(" ");
        }
        for _ in inn..out {
            print!("s");
        }
        println!();
    }
}

/// # Safety
///
/// `ss`, `sh1`, and optional `sh2` must be valid pointers.
#[no_mangle]
pub unsafe extern "C" fn optim_mult(
    ss: *mut SkewShapeAbi,
    sh1: *const IVector,
    sh2: *const IVector,
    maxrows: c_int,
    maxcols: c_int,
) -> c_int {
    if ss.is_null() || sh1.is_null() {
        return -1;
    }
    unsafe { zero_skew_shape(ss) };
    if unsafe { part_valid(sh1) } == 0 || (!sh2.is_null() && unsafe { part_valid(sh2) } == 0) {
        return -1;
    }
    let mut left = unsafe { ivector_values(sh1) };
    let mut right = if sh2.is_null() {
        &[][..]
    } else {
        unsafe { ivector_values(sh2) }
    };
    let mut len1 = abi_part_length(left);
    let mut len2 = abi_part_length(right);
    let mut width1 = if len1 == 0 { 0 } else { left[0] };
    let mut width2 = if len2 == 0 { 0 } else { right[0] };
    if maxrows >= 0 && (len1 > maxrows as usize || len2 > maxrows as usize) {
        return 0;
    }
    if maxcols >= 0 && (width1 > maxcols || width2 > maxcols) {
        return 0;
    }
    if maxrows >= 0 && maxcols >= 0 {
        let start = if len1 + len2 < maxrows as usize {
            len2
        } else {
            (maxrows as usize).saturating_sub(len1)
        };
        for r in start..len2 {
            let left_index = maxrows as usize - r - 1;
            if left[left_index] + right[r] > maxcols {
                return 0;
            }
        }
    }

    let mut fc1 = if maxrows >= 0 && len1 == maxrows as usize && len1 > 0 {
        left[len1 - 1]
    } else {
        0
    };
    let mut fr1 = 0usize;
    while maxcols >= 0 && fr1 < len1 && left[fr1] == maxcols {
        fr1 += 1;
    }
    let mut fc2 = if maxrows >= 0 && len2 == maxrows as usize && len2 > 0 {
        right[len2 - 1]
    } else {
        0
    };
    let mut fr2 = 0usize;
    while maxcols >= 0 && fr2 < len2 && right[fr2] == maxcols {
        fr2 += 1;
    }

    let size1 = left[fr1..len1]
        .iter()
        .copied()
        .sum::<i32>()
        .saturating_sub((len1 - fr1) as i32 * fc1);
    let size2 = right[fr2..len2]
        .iter()
        .copied()
        .sum::<i32>()
        .saturating_sub((len2 - fr2) as i32 * fc2);
    if size1 > size2 {
        std::mem::swap(&mut left, &mut right);
        std::mem::swap(&mut len1, &mut len2);
        std::mem::swap(&mut width1, &mut width2);
        std::mem::swap(&mut fc1, &mut fc2);
        std::mem::swap(&mut fr1, &mut fr2);
    }
    let outer = left[fr1..len1]
        .iter()
        .map(|part| part - fc1)
        .collect::<Vec<_>>();
    let clen = if fc1 + fc2 > 0 {
        maxrows.max(0) as usize
    } else {
        len2 + fr1
    };
    let mut cont = vec![fc1; clen];
    for entry in cont.iter_mut().take(fr1) {
        *entry = maxcols;
    }
    for r in 0..len2 {
        if fr1 + r < cont.len() {
            cont[fr1 + r] = right[r] + fc1;
        }
    }
    unsafe { fill_skew_shape(ss, &outer, None, Some(&cont), 1) }
}

/// # Safety
///
/// Pointers must be valid as in upstream `optim_fusion`.
#[no_mangle]
pub unsafe extern "C" fn optim_fusion(
    ss: *mut SkewShapeAbi,
    sh1: *const IVector,
    sh2: *const IVector,
    rows: c_int,
    level: c_int,
) -> c_int {
    if ss.is_null() || sh1.is_null() || sh2.is_null() || rows < 0 {
        return -1;
    }
    unsafe { zero_skew_shape(ss) };
    if unsafe { part_length(sh1) } > rows || unsafe { part_length(sh2) } > rows {
        return 0;
    }
    let mut left = unsafe { ivector_values(sh1) };
    let mut right = unsafe { ivector_values(sh2) };
    let mut d1 = 0;
    let mut d2 = 0;
    let mut s1 = rows.saturating_mul(level);
    let mut s2 = s1;
    for d in 1..=rows {
        let s = (rows - d).saturating_mul(level)
            - rows.saturating_mul(unsafe { part_entry(sh1, d - 1) });
        if s < s1 {
            d1 = d;
            s1 = s;
        }
        let s = (rows - d).saturating_mul(level)
            - rows.saturating_mul(unsafe { part_entry(sh2, d - 1) });
        if s < s2 {
            d2 = d;
            s2 = s;
        }
    }
    if s1 > s2 {
        std::mem::swap(&mut left, &mut right);
        d1 = d2;
    }
    let d = d1;
    let sh1d = left.get((d - 1) as usize).copied().unwrap_or(0);
    let rows_usize = rows as usize;
    let d_usize = d as usize;
    let mut nsh1 = vec![0; rows_usize];
    let mut nsh2 = vec![0; rows_usize];
    for i in 0..rows_usize.saturating_sub(d_usize) {
        nsh1[i] = left.get(d_usize + i).copied().unwrap_or(0) - sh1d + level;
    }
    for i in 0..d_usize {
        nsh1[rows_usize - d_usize + i] = left.get(i).copied().unwrap_or(0) - sh1d;
    }
    for i in 0..d_usize {
        nsh2[i] = right
            .get(rows_usize.saturating_sub(d_usize) + i)
            .copied()
            .unwrap_or(0)
            + sh1d;
    }
    for i in 0..rows_usize.saturating_sub(d_usize) {
        nsh2[d_usize + i] = right.get(i).copied().unwrap_or(0) + sh1d - level;
    }
    unsafe { fill_skew_shape(ss, &nsh1, None, Some(&nsh2), 1) }
}

/// # Safety
///
/// Pointers must be valid as in upstream `optim_skew`.
#[no_mangle]
pub unsafe extern "C" fn optim_skew(
    ss: *mut SkewShapeAbi,
    outer: *const IVector,
    inner: *const IVector,
    content: *const IVector,
    maxrows: c_int,
) -> c_int {
    if ss.is_null() || outer.is_null() {
        return -1;
    }
    unsafe { zero_skew_shape(ss) };
    if inner.is_null() {
        return unsafe { optim_mult(ss, outer, content, maxrows, -1) };
    }
    if unsafe { part_valid(outer) } == 0
        || unsafe { part_valid(inner) } == 0
        || (!content.is_null() && unsafe { part_valid(content) } == 0)
    {
        return -1;
    }
    if unsafe { part_leq(inner, outer) } == 0 {
        return 0;
    }
    if maxrows >= 0 && !content.is_null() && unsafe { part_length(content) } > maxrows {
        return 0;
    }
    let outer_values = unsafe { ivector_values(outer) };
    let inner_values = unsafe { ivector_values(inner) };
    let content_values = if content.is_null() {
        &[][..]
    } else {
        unsafe { ivector_values(content) }
    };
    unsafe {
        fill_skew_shape(
            ss,
            &outer_values[..abi_part_length(outer_values)],
            Some(&inner_values[..abi_part_length(inner_values)]),
            Some(&content_values[..abi_part_length(content_values)]),
            1,
        )
    }
}

/// # Safety
///
/// Pointers must be valid as in upstream `optim_coef`.
#[no_mangle]
pub unsafe extern "C" fn optim_coef(
    ss: *mut SkewShapeAbi,
    out: *const IVector,
    sh1: *const IVector,
    sh2: *const IVector,
) -> c_int {
    if ss.is_null() || out.is_null() || sh1.is_null() || sh2.is_null() {
        return -1;
    }
    unsafe { zero_skew_shape(ss) };
    let out_values = unsafe { ivector_values(out) };
    let sh1_values = unsafe { ivector_values(sh1) };
    let sh2_values = unsafe { ivector_values(sh2) };
    match native_optim_coef(out_values, sh1_values, sh2_values) {
        Ok(OptimizedCoef::Zero) => 0,
        Ok(OptimizedCoef::One) => {
            unsafe {
                (*ss).sign = 1;
            }
            0
        }
        Ok(OptimizedCoef::Count(shape)) => unsafe {
            fill_skew_shape(
                ss,
                &shape.outer,
                Some(&shape.inner),
                Some(&shape.content),
                2,
            )
        },
        Err(_) => -1,
    }
}

/// # Safety
///
/// `av` must be a valid `argv` array with `ac` elements. This follows
/// upstream's `optind`-based parser.
#[no_mangle]
pub unsafe extern "C" fn get_vect_arg(ac: c_int, av: *mut *mut c_char) -> *mut IVector {
    if av.is_null() || ac < 0 {
        return ptr::null_mut();
    }
    unsafe {
        if optind == ac {
            return ptr::null_mut();
        }
        if optind == 0 {
            optind += 1;
        } else if optind < ac {
            let arg = *av.add(optind as usize);
            if !arg.is_null() {
                let first = *arg;
                let second = *arg.add(1);
                if (first == b'-' as c_char || first == b'/' as c_char) && second == 0 {
                    optind += 1;
                }
            }
        }
    }

    let mut values = Vec::new();
    unsafe {
        while optind < ac {
            let arg = *av.add(optind as usize);
            if arg.is_null() {
                break;
            }
            let text = CStr::from_ptr(arg).to_string_lossy();
            let Ok(value) = text.parse::<i32>() else {
                break;
            };
            values.push(value);
            optind += 1;
        }
    }
    if values.is_empty() {
        return ptr::null_mut();
    }
    unsafe { ivector_from_partition(&values, values.len()) }
}

#[no_mangle]
pub extern "C" fn lrcalc_new_abi_version() -> u32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    unsafe fn vector_from_values(values: &[i32]) -> *mut IVector {
        let v = iv_new_zero(values.len() as u32);
        assert!(!v.is_null());
        unsafe { ivector_values_mut(v) }.copy_from_slice(values);
        v
    }

    unsafe fn collect_lc_trimmed(lc: *mut IvLinComb) -> Vec<(Vec<i32>, i32)> {
        let mut terms = Vec::new();
        let mut itr = IvlcIter {
            ht: ptr::null_mut(),
            index: 0,
            i: 0,
        };
        unsafe {
            ivlc_first(lc, &mut itr);
            while ivlc_good(&itr) != 0 {
                let key = ivlc_key(&itr);
                let mut partition = ivector_values(key).to_vec();
                while partition.last() == Some(&0) {
                    partition.pop();
                }
                terms.push((partition, ivlc_value(&itr)));
                ivlc_next(&mut itr);
            }
        }
        terms.sort();
        terms
    }

    unsafe fn collect_lc_exact(lc: *mut IvLinComb) -> Vec<(Vec<i32>, i32)> {
        let mut terms = Vec::new();
        let mut itr = IvlcIter {
            ht: ptr::null_mut(),
            index: 0,
            i: 0,
        };
        unsafe {
            ivlc_first(lc, &mut itr);
            while ivlc_good(&itr) != 0 {
                let key = ivlc_key(&itr);
                terms.push((ivector_values(key).to_vec(), ivlc_value(&itr)));
                ivlc_next(&mut itr);
            }
        }
        terms.sort();
        terms
    }

    unsafe fn collect_lc_lengths(lc: *mut IvLinComb) -> Vec<(Vec<i32>, u32, i32)> {
        let mut terms = Vec::new();
        let mut itr = IvlcIter {
            ht: ptr::null_mut(),
            index: 0,
            i: 0,
        };
        unsafe {
            ivlc_first(lc, &mut itr);
            while ivlc_good(&itr) != 0 {
                let key = ivlc_key(&itr);
                let mut partition = ivector_values(key).to_vec();
                while partition.last() == Some(&0) {
                    partition.pop();
                }
                terms.push((partition, (*key).length, ivlc_value(&itr)));
                ivlc_next(&mut itr);
            }
        }
        terms.sort();
        terms
    }

    unsafe fn collect_lrit_contents(lrit: *mut LrTabIter) -> Vec<Vec<i32>> {
        let mut contents = Vec::new();
        unsafe {
            while lrit_good(lrit) != 0 {
                let mut content = ivector_values((*lrit).cont).to_vec();
                while content.last() == Some(&0) {
                    content.pop();
                }
                contents.push(content);
                lrit_next(lrit);
            }
        }
        contents.sort();
        contents
    }

    #[test]
    fn ivector_allocation_copy_hash_and_sum() {
        unsafe {
            let v = iv_new_zero(3);
            assert!(!v.is_null());
            (*v).array[0] = 2;
            *ptr::addr_of_mut!((*v).array).cast::<i32>().add(1) = 1;
            *ptr::addr_of_mut!((*v).array).cast::<i32>().add(2) = 0;
            assert_eq!(iv_sum(v), 3);

            let copy = iv_new_copy(v);
            assert_eq!(iv_cmp(v, copy), 0);
            assert_eq!(iv_hash(v), iv_hash(copy));

            iv_free(copy);
            iv_free(v);
        }
    }

    #[test]
    fn vector_partition_and_list_helpers_match_c_surface() {
        unsafe {
            let v = vector_from_values(&[4, 2, 0]);
            let w = vector_from_values(&[2, 1, 0]);
            let dst = iv_new_zero(3);
            assert_eq!(part_valid(v), 1);
            assert_eq!(part_decr(v), 1);
            assert_eq!(part_length(v), 2);
            assert_eq!(part_entry(v, 7), 0);
            assert_eq!(part_leq(w, v), 1);
            assert_eq!(iv_lesseq(w, v), 1);

            iv_mult(dst, 2, w);
            assert_eq!(ivector_values(dst), &[4, 2, 0]);
            iv_reverse(dst, v);
            assert_eq!(ivector_values(dst), &[0, 2, 4]);
            assert_eq!(iv_max(v), 4);
            assert_eq!(iv_min(v), 0);
            let gcd_input = vector_from_values(&[6, 9, 15]);
            assert_eq!(iv_gcd(gcd_input), 3);
            iv_free(gcd_input);
            let initialized = iv_new_init(3, 5, 4, 3, 0, 0, 0, 0, 0);
            assert_eq!(ivector_values(initialized), &[5, 4, 3]);

            let conj = part_conj(v);
            assert!(!conj.is_null());
            assert_eq!(ivector_values(conj), &[2, 2, 1, 1]);

            let list = il_new(1);
            assert!(!list.is_null());
            assert_eq!(il_append(list, 7), 0);
            assert_eq!(il_insert(list, 0, 3), 0);
            assert_eq!(il_poplast(list), 7);
            assert_eq!(il_poplast(list), 3);
            let initialized_list = il_new_init(2, 2, 11, 13, 0, 0, 0, 0, 0, 0);
            assert!(!initialized_list.is_null());
            assert_eq!((*initialized_list).length, 2);
            assert_eq!(*(*initialized_list).array.add(0), 11);
            assert_eq!(*(*initialized_list).array.add(1), 13);

            iv_free(conj);
            il_free(initialized_list);
            il_free(list);
            iv_free(initialized);
            iv_free(dst);
            iv_free(w);
            iv_free(v);
        }
    }

    #[test]
    fn permutation_string_and_partition_iterator_helpers_work() {
        unsafe {
            let perm = vector_from_values(&[2, 1, 3]);
            assert_eq!(perm_valid(perm), 1);
            assert_eq!(perm_length(perm), 1);
            assert_eq!(perm_group(perm), 2);

            let string = vector_from_values(&[0, 1, 0]);
            let dimvec = str2dimvec(string);
            assert!(!dimvec.is_null());
            assert_eq!(ivector_values(dimvec), &[2, 3]);

            let as_perm = string2perm(string);
            assert!(!as_perm.is_null());
            assert_eq!(ivector_values(as_perm), &[1, 3, 2]);
            let back = perm2string(as_perm, dimvec);
            assert!(!back.is_null());
            assert_eq!(ivector_values(back), &[0, 1, 0]);

            let strings = all_strings(dimvec);
            assert!(!strings.is_null());
            assert_eq!((*strings).length, 3);

            let p = iv_new_zero(2);
            let mut itr = PartIter {
                part: ptr::null_mut(),
                outer: ptr::null_mut(),
                inner: ptr::null_mut(),
                length: 0,
                rows: 0,
                opt: 0,
            };
            pitr_box_first(&mut itr, p, 2, 2);
            let mut parts = Vec::new();
            while pitr_good(&itr) != 0 {
                parts.push(ivector_values(p).to_vec());
                pitr_next(&mut itr);
            }
            assert_eq!(
                parts,
                vec![
                    vec![2, 2],
                    vec![2, 1],
                    vec![2, 0],
                    vec![1, 1],
                    vec![1, 0],
                    vec![0, 0]
                ]
            );

            iv_free(p);
            ivl_free_all(strings);
            iv_free(back);
            iv_free(as_perm);
            iv_free(dimvec);
            iv_free(string);
            iv_free(perm);
        }
    }

    #[test]
    fn part_quantum_helpers_match_upstream_formula() {
        unsafe {
            let p = vector_from_values(&[3, 2, 1]);
            let degree = part_qdegree(p, 2);
            assert_eq!(degree, 1);
            assert_eq!(part_qentry(p, 0, degree, 2), 1);
            assert_eq!(part_qentry(p, 1, degree, 2), 0);

            let rectangle = vector_from_values(&[2, 2, 2]);
            assert_eq!(part_qdegree(rectangle, 2), 0);
            assert_eq!(part_qentry(rectangle, 0, 0, 2), 2);
            assert_eq!(part_qentry(rectangle, 2, 0, 2), 2);

            iv_free(rectangle);
            iv_free(p);
        }
    }

    #[test]
    fn lrit_iterator_lists_small_skew_contents() {
        unsafe {
            let outer = vector_from_values(&[2, 1]);
            let inner = vector_from_values(&[1]);
            let lrit = lrit_new(outer, inner, ptr::null(), -1, -1, -1);
            assert!(!lrit.is_null());
            assert_eq!(collect_lrit_contents(lrit), vec![vec![1, 1], vec![2]]);

            lrit_free(lrit);
            iv_free(inner);
            iv_free(outer);
        }
    }

    #[test]
    fn lrit_count_expand_and_lrcoef_count_return_abi_results() {
        unsafe {
            let outer = vector_from_values(&[2, 1]);
            let inner = vector_from_values(&[1]);
            let lrit = lrit_new(outer, inner, ptr::null(), -1, -1, -1);
            assert!(!lrit.is_null());
            let counted = lrit_count(lrit);
            assert!(!counted.is_null());
            assert_eq!(
                collect_lc_trimmed(counted),
                vec![(vec![1, 1], 1), (vec![2], 1)]
            );
            lrit_free(lrit);

            let expanded = lrit_expand(outer, inner, ptr::null(), -1, -1, -1);
            assert!(!expanded.is_null());
            assert_eq!(
                collect_lc_trimmed(expanded),
                vec![(vec![1, 1], 1), (vec![2], 1)]
            );

            let content = vector_from_values(&[1, 1]);
            assert_eq!(lrcoef_count(outer, inner, content), 1);

            iv_free(content);
            ivlc_free_all(expanded);
            ivlc_free_all(counted);
            iv_free(inner);
            iv_free(outer);
        }
    }

    #[test]
    fn lrit_iterator_respects_row_bound_and_empty_shapes() {
        unsafe {
            let outer = vector_from_values(&[2, 1]);
            let inner = vector_from_values(&[1]);
            let row_lrit = lrit_new(outer, inner, ptr::null(), 1, -1, -1);
            assert!(!row_lrit.is_null());
            assert_eq!(collect_lrit_contents(row_lrit), vec![vec![2]]);
            lrit_free(row_lrit);

            let bad_inner = vector_from_values(&[2, 2]);
            let empty_lrit = lrit_new(outer, bad_inner, ptr::null(), -1, -1, -1);
            assert!(!empty_lrit.is_null());
            assert!(collect_lrit_contents(empty_lrit).is_empty());
            lrit_free(empty_lrit);

            iv_free(bad_inner);
            iv_free(inner);
            iv_free(outer);
        }
    }

    #[test]
    fn ivlincomb_add_iterates_combines_and_removes_zero() {
        unsafe {
            let lc = ivlc_new(5, 2);
            assert!(!lc.is_null());
            let key = vector_from_values(&[2, 1]);
            let hash = iv_hash(key) as u32;

            assert_eq!(ivlc_add_element(lc, 2, key, hash, LC_COPY_KEY), 0);
            assert_eq!(ivlc_add_element(lc, 3, key, hash, LC_COPY_KEY), 0);
            assert_eq!(ivlc_card(lc), 1);
            assert_eq!(collect_lc_trimmed(lc), vec![(vec![2, 1], 5)]);

            assert_eq!(
                ivlc_add_element(lc, -5, key, hash, LC_COPY_KEY | LC_FREE_ZERO),
                0
            );
            assert_eq!(ivlc_card(lc), 0);
            assert!(collect_lc_trimmed(lc).is_empty());

            iv_free(key);
            ivlc_free_all(lc);
        }
    }

    #[test]
    fn ivlincomb_add_multiple_scales_terms() {
        unsafe {
            let src = ivlc_new(5, 2);
            let dst = ivlc_new(5, 2);
            assert!(!src.is_null());
            assert!(!dst.is_null());
            let key = vector_from_values(&[2]);
            let hash = iv_hash(key) as u32;

            assert_eq!(ivlc_add_element(src, 4, key, hash, LC_COPY_KEY), 0);
            assert_eq!(ivlc_add_multiple(dst, 3, src, LC_COPY_KEY), 0);
            assert_eq!(collect_lc_trimmed(dst), vec![(vec![2], 12)]);

            iv_free(key);
            ivlc_free_all(dst);
            ivlc_free_all(src);
        }
    }

    #[test]
    fn ivlincomb_transfer_key_ownership_path() {
        unsafe {
            let lc = ivlc_new(5, 2);
            assert!(!lc.is_null());
            let key = vector_from_values(&[3, 1]);
            let duplicate = vector_from_values(&[3, 1]);
            let hash = iv_hash(key) as u32;

            assert_eq!(ivlc_add_element(lc, 7, key, hash, LC_FREE_ZERO), 0);
            assert_eq!(ivlc_add_element(lc, -2, duplicate, hash, LC_FREE_ZERO), 0);
            assert_eq!(collect_lc_trimmed(lc), vec![(vec![3, 1], 5)]);

            ivlc_free_all(lc);
        }
    }

    #[test]
    fn schubert_trans_and_monk_work_through_abi() {
        unsafe {
            let permutation = vector_from_values(&[2, 1]);
            let poly = trans(permutation, 0);
            assert!(!poly.is_null());
            assert_eq!(collect_lc_exact(poly), vec![(vec![1], 1)]);

            let slc = ivlc_new(5, 2);
            assert!(!slc.is_null());
            let identity = vector_from_values(&[]);
            let hash = iv_hash(identity) as u32;
            assert_eq!(ivlc_add_element(slc, 1, identity, hash, LC_FREE_ZERO), 0);

            let product = monk(1, slc, 0);
            assert!(!product.is_null());
            assert_eq!(collect_lc_trimmed(product), vec![(vec![2, 1], 1)]);

            ivlc_free_all(product);
            ivlc_free_all(slc);
            ivlc_free_all(poly);
            iv_free(permutation);
        }
    }

    #[test]
    fn schubert_multiply_polynomial_mutates_input_lc() {
        unsafe {
            let poly = ivlc_new(5, 2);
            assert!(!poly.is_null());
            let monomial = vector_from_values(&[1]);
            let hash = iv_hash(monomial) as u32;
            assert_eq!(ivlc_add_element(poly, 1, monomial, hash, LC_FREE_ZERO), 0);

            let permutation = vector_from_values(&[2, 1]);
            let product = mult_poly_schubert(poly, permutation, 0);
            assert_eq!(product, poly);
            assert_eq!(collect_lc_trimmed(product), vec![(vec![3, 1, 2], 1)]);

            ivlc_free_all(product);
            iv_free(permutation);
        }
    }

    #[test]
    fn schubert_multiply_permutations_through_abi() {
        unsafe {
            let left = vector_from_values(&[2, 1]);
            let right = vector_from_values(&[2, 1]);
            let product = mult_schubert(left, right, 0);
            assert!(!product.is_null());
            assert_eq!(collect_lc_trimmed(product), vec![(vec![3, 1, 2], 1)]);

            ivlc_free_all(product);
            iv_free(right);
            iv_free(left);
        }
    }

    #[test]
    fn schubert_multiply_strings_through_abi() {
        unsafe {
            let left = vector_from_values(&[0, 1]);
            let right = vector_from_values(&[1, 0]);
            let product = mult_schubert_str(left, right);
            assert!(!product.is_null());
            assert_eq!(collect_lc_exact(product), vec![(vec![1, 0], 1)]);

            ivlc_free_all(product);
            iv_free(right);
            iv_free(left);
        }
    }

    #[test]
    fn schur_lrcoef_computes_upstream_example() {
        unsafe {
            let outer = iv_new_zero(3);
            let inner1 = iv_new_zero(2);
            let inner2 = iv_new_zero(2);
            assert!(!outer.is_null());
            assert!(!inner1.is_null());
            assert!(!inner2.is_null());

            (*outer).array[0] = 3;
            *ptr::addr_of_mut!((*outer).array).cast::<i32>().add(1) = 2;
            *ptr::addr_of_mut!((*outer).array).cast::<i32>().add(2) = 1;
            (*inner1).array[0] = 2;
            *ptr::addr_of_mut!((*inner1).array).cast::<i32>().add(1) = 1;
            (*inner2).array[0] = 2;
            *ptr::addr_of_mut!((*inner2).array).cast::<i32>().add(1) = 1;

            assert_eq!(schur_lrcoef(outer, inner1, inner2), 2);

            iv_free(inner2);
            iv_free(inner1);
            iv_free(outer);
        }
    }

    #[test]
    fn optimization_wrappers_fill_skew_shape_records() {
        unsafe {
            let sh1 = vector_from_values(&[1]);
            let sh2 = vector_from_values(&[1]);
            let mut ss = SkewShapeAbi {
                outer: ptr::null_mut(),
                inner: ptr::null_mut(),
                cont: ptr::null_mut(),
                sign: 0,
            };
            assert_eq!(optim_mult(&mut ss, sh1, sh2, -1, -1), 0);
            assert_eq!(ss.sign, 1);
            assert_eq!(ivector_values(ss.outer), &[1]);
            assert_eq!(ivector_values(ss.cont), &[1]);
            sksh_dealloc(&mut ss);

            let out = vector_from_values(&[3, 2, 1]);
            let coef_inner = vector_from_values(&[2, 1]);
            assert_eq!(optim_coef(&mut ss, out, coef_inner, coef_inner), 0);
            assert_eq!(ss.sign, 2);
            assert!(!ss.outer.is_null());
            assert!(!ss.inner.is_null());
            assert!(!ss.cont.is_null());
            sksh_dealloc(&mut ss);

            iv_free(coef_inner);
            iv_free(out);
            iv_free(sh2);
            iv_free(sh1);
        }
    }

    #[test]
    fn schur_mult_returns_linear_combination() {
        unsafe {
            let sh1 = vector_from_values(&[1]);
            let sh2 = vector_from_values(&[1]);
            let lc = schur_mult(sh1, sh2, -1, -1, -1);
            assert!(!lc.is_null());
            assert_eq!(collect_lc_trimmed(lc), vec![(vec![1, 1], 1), (vec![2], 1)]);

            ivlc_free_all(lc);
            iv_free(sh2);
            iv_free(sh1);
        }
    }

    #[test]
    fn schur_mult_respects_bounds_through_abi() {
        unsafe {
            let sh1 = vector_from_values(&[2, 1]);
            let sh2 = vector_from_values(&[2, 1]);

            let row_lc = schur_mult(sh1, sh2, 3, -1, -1);
            assert!(!row_lc.is_null());
            assert_eq!(
                collect_lc_trimmed(row_lc),
                vec![
                    (vec![2, 2, 2], 1),
                    (vec![3, 2, 1], 2),
                    (vec![3, 3], 1),
                    (vec![4, 1, 1], 1),
                    (vec![4, 2], 1),
                ]
            );

            let col_lc = schur_mult(sh1, sh2, -1, 2, -1);
            assert!(!col_lc.is_null());
            assert_eq!(
                collect_lc_trimmed(col_lc),
                vec![(vec![2, 2, 1, 1], 1), (vec![2, 2, 2], 1)]
            );

            ivlc_free_all(col_lc);
            ivlc_free_all(row_lc);
            iv_free(sh2);
            iv_free(sh1);
        }
    }

    #[test]
    fn schur_mult_honors_part_size_padding() {
        unsafe {
            let sh1 = vector_from_values(&[1]);
            let sh2 = vector_from_values(&[1]);
            let lc = schur_mult(sh1, sh2, -1, -1, 4);
            assert!(!lc.is_null());
            assert_eq!(
                collect_lc_lengths(lc),
                vec![(vec![1, 1], 4, 1), (vec![2], 4, 1)]
            );

            ivlc_free_all(lc);
            iv_free(sh2);
            iv_free(sh1);
        }
    }

    #[test]
    fn fusion_reduce_abi_mutates_vector() {
        unsafe {
            let la = vector_from_values(&[2, 0]);
            let tmp = iv_new_zero(2);
            assert!(!tmp.is_null());

            assert_eq!(fusion_reduce(la, 0, tmp), -1);
            assert_eq!(ivector_values(la), &[1, 1]);

            iv_free(tmp);
            iv_free(la);
        }
    }

    #[test]
    fn fusion_reduce_lc_cancels_zero_terms() {
        unsafe {
            let lc = ivlc_new(5, 2);
            assert!(!lc.is_null());
            let first = vector_from_values(&[2, 0]);
            let first_hash = iv_hash(first) as u32;
            let second = vector_from_values(&[1, 1]);
            let second_hash = iv_hash(second) as u32;

            assert_eq!(ivlc_add_element(lc, 1, first, first_hash, LC_FREE_ZERO), 0);
            assert_eq!(
                ivlc_add_element(lc, 1, second, second_hash, LC_FREE_ZERO),
                0
            );
            assert_eq!(fusion_reduce_lc(lc, 0), 0);
            assert!(collect_lc_trimmed(lc).is_empty());

            ivlc_free_all(lc);
        }
    }

    #[test]
    fn schur_mult_fusion_returns_reduced_product() {
        unsafe {
            let sh1 = vector_from_values(&[2, 1]);
            let sh2 = vector_from_values(&[2, 1]);
            let lc = schur_mult_fusion(sh1, sh2, 3, 2);
            assert!(!lc.is_null());
            assert_eq!(
                collect_lc_trimmed(lc),
                vec![(vec![2, 2, 2], 1), (vec![3, 2, 1], 1)]
            );

            ivlc_free_all(lc);
            iv_free(sh2);
            iv_free(sh1);
        }
    }

    #[test]
    fn schur_mult_fusion_uses_row_length_keys() {
        unsafe {
            let sh1 = vector_from_values(&[1]);
            let sh2 = vector_from_values(&[1]);
            let lc = schur_mult_fusion(sh1, sh2, 2, 1);
            assert!(!lc.is_null());
            assert_eq!(collect_lc_lengths(lc), vec![(vec![1, 1], 2, 1)]);

            ivlc_free_all(lc);
            iv_free(sh2);
            iv_free(sh1);
        }
    }

    #[test]
    fn schur_skew_returns_linear_combination() {
        unsafe {
            let outer = vector_from_values(&[2, 1]);
            let inner = vector_from_values(&[1]);
            let lc = schur_skew(outer, inner, -1, -1);
            assert!(!lc.is_null());
            assert_eq!(collect_lc_trimmed(lc), vec![(vec![1, 1], 1), (vec![2], 1)]);

            ivlc_free_all(lc);
            iv_free(inner);
            iv_free(outer);
        }
    }

    #[test]
    fn schur_skew_respects_row_bound_through_abi() {
        unsafe {
            let outer = vector_from_values(&[3, 2, 1]);
            let inner = vector_from_values(&[2, 1]);
            let lc = schur_skew(outer, inner, 2, -1);
            assert!(!lc.is_null());
            assert_eq!(collect_lc_trimmed(lc), vec![(vec![2, 1], 2), (vec![3], 1)]);

            ivlc_free_all(lc);
            iv_free(inner);
            iv_free(outer);
        }
    }

    #[test]
    fn schur_skew_honors_part_size_padding() {
        unsafe {
            let outer = vector_from_values(&[2, 1]);
            let inner = vector_from_values(&[1]);
            let lc = schur_skew(outer, inner, -1, 3);
            assert!(!lc.is_null());
            assert_eq!(
                collect_lc_lengths(lc),
                vec![(vec![1, 1], 3, 1), (vec![2], 3, 1)]
            );

            ivlc_free_all(lc);
            iv_free(inner);
            iv_free(outer);
        }
    }

    #[test]
    fn schur_coprod_returns_linear_combination() {
        unsafe {
            let sh = vector_from_values(&[1]);
            let lc = schur_coprod(sh, 1, 1, -1, 0);
            assert!(!lc.is_null());
            assert_eq!(collect_lc_trimmed(lc), vec![(vec![2], 1)]);

            ivlc_free_all(lc);
            iv_free(sh);
        }
    }

    #[test]
    fn schur_coprod_all_keeps_redundant_terms() {
        unsafe {
            let sh = vector_from_values(&[1]);
            let lc = schur_coprod(sh, 1, 1, -1, 1);
            assert!(!lc.is_null());
            assert_eq!(collect_lc_trimmed(lc), vec![(vec![1, 1], 1), (vec![2], 1)]);

            ivlc_free_all(lc);
            iv_free(sh);
        }
    }

    #[test]
    fn schur_coprod_honors_part_size_padding() {
        unsafe {
            let sh = vector_from_values(&[1]);
            let lc = schur_coprod(sh, 1, 1, 4, 0);
            assert!(!lc.is_null());
            assert_eq!(collect_lc_lengths(lc), vec![(vec![2], 4, 1)]);

            ivlc_free_all(lc);
            iv_free(sh);
        }
    }
}
