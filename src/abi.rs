//! C ABI types and functions matching the original `lrcalc` headers.

use crate::lrcoef::lrcoef_i64;
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
