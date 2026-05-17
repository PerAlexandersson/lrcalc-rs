//! C ABI types and functions matching the original `lrcalc` headers.

use crate::lrcoef::lrcoef_i64;
use crate::schubert::{
    monk_product, multiply_poly_schubert, multiply_schubert, multiply_schubert_strings,
    trans_polynomial, LinearCombination,
};
use crate::schur::{
    fusion_reduce_values, schur_coproduct_expansion, schur_product_expansion,
    schur_product_fusion_expansion, schur_skew_expansion, SchurTerm, SignedSchurTerm,
};
use libc::{c_int, c_longlong, c_void};
use std::mem;
use std::ptr;
use std::slice;

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

unsafe fn ivector_values<'a>(v: *const IVector) -> &'a [i32] {
    let length = unsafe { (*v).length as usize };
    let data = unsafe { ptr::addr_of!((*v).array).cast::<i32>() };
    unsafe { slice::from_raw_parts(data, length) }
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
    let terms = match schur_skew_expansion(outer_values, inner_values, rows) {
        Ok(terms) => terms,
        Err(_) => return ptr::null_mut(),
    };
    let key_len = default_skew_key_len(outer_values, inner_values, rows, partsz);
    unsafe { ivlc_from_terms(&terms, key_len) }
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
