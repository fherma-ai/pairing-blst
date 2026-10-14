//! The optimal ate pairing on BLS12-381, through blst's assembly.
//!
//! The library is a C API over hand-written assembly, so this is the one
//! answer in the family that is mostly `unsafe`: every call takes raw
//! pointers and every field element is held in Montgomery form, which is why
//! the limbs are converted in and out rather than copied.
//!
//! The flattening reads blst's own tower — `fp12.fp6[i].fp2[j].fp[k]` — in
//! the order the specification names, which is the order arkworks' tower
//! gives as well. Whether the two libraries agree element for element is not
//! assumed here: the specification is judged against a reference answer, and
//! a disagreement would be the measurement's first finding rather than a
//! silent difference.
use blst::{
    blst_final_exp, blst_fp, blst_fp12, blst_fp2, blst_fp_from_uint64, blst_miller_loop,
    blst_p1_affine, blst_p2_affine, blst_uint64_from_fp,
};

use crate::fherma::{Inputs, Outputs, Point, Tensor};

/// The `i`-th row of a limb tensor as a base-field element, Montgomery form.
fn fp(data: &[u64], i: usize, limbs: usize) -> blst_fp {
    let mut words = [0u64; 6];
    words[..limbs].copy_from_slice(&data[i * limbs..(i + 1) * limbs]);
    let mut out = blst_fp::default();
    unsafe { blst_fp_from_uint64(&mut out, words.as_ptr()) };
    out
}

/// Rows `i` and `i+1` as an element of Fp2, c0 first.
fn fp2(data: &[u64], i: usize, limbs: usize) -> blst_fp2 {
    blst_fp2 { fp: [fp(data, i, limbs), fp(data, i + 1, limbs)] }
}

/// Out of Montgomery form and into the specification's twelve elements.
fn flatten(e: blst_fp12, limbs: usize) -> Tensor<u64> {
    let mut data = Vec::with_capacity(12 * limbs);
    for six in e.fp6 {
        for two in six.fp2 {
            for one in two.fp {
                let mut words = [0u64; 6];
                unsafe { blst_uint64_from_fp(words.as_mut_ptr(), &one) };
                data.extend_from_slice(&words[..limbs]);
            }
        }
    }
    Tensor { shape: vec![12, limbs as i64], data }
}

#[derive(Default)]
pub struct State {
    limbs: usize,
}

pub fn init(p: &Point) -> State {
    State { limbs: p.L as usize }
}

pub fn run(state: &mut State, inp: &Inputs) -> Outputs {
    let limbs = state.limbs;

    let p = blst_p1_affine {
        x: fp(&inp.px.data, 0, limbs),
        y: fp(&inp.py.data, 0, limbs),
    };
    let q = blst_p2_affine {
        x: fp2(&inp.qx.data, 0, limbs),
        y: fp2(&inp.qy.data, 0, limbs),
    };

    // The two halves the specification names, as this library exposes them:
    // the Miller loop into Fp12, then the final exponentiation.
    let mut loops = blst_fp12::default();
    let mut answer = blst_fp12::default();
    unsafe {
        blst_miller_loop(&mut loops, &q, &p);
        blst_final_exp(&mut answer, &loops);
    }

    Outputs { e: flatten(answer, limbs) }
}

pub fn free(state: State) {
    drop(state);
}
