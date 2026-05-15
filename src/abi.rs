//! C ABI types and functions matching the original `lrcalc` headers.

use crate::lrcoef::lrcoef_i64;
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
pub struct IvLinComb {
    _private: [u8; 0],
}

#[repr(C)]
pub struct IvlcIter {
    pub ht: *mut IvLinComb,
    pub index: u32,
    pub i: u32,
}

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

#[no_mangle]
pub extern "C" fn lrcalc_new_abi_version() -> u32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
