# The pairing on BLS12-381, in assembly

An answer to `pairing/bls12381@1.0.0` on FHERMA with
[blst](https://github.com/supranational/blst) 0.3.17 — the BLS12-381 library
Ethereum's clients sign and verify with, hand-written assembly behind a C API.

It exists to be read against the arkworks answer to the same specification:
one library whose field arithmetic is portable Rust, one whose arithmetic is
assembly, on one curve, one machine and the same cases.

`src/solve.rs` is the whole answer. Two calls, which are the two halves the
specification names:

```rust
blst_miller_loop(&mut loops, &q, &p);
blst_final_exp(&mut answer, &loops);
```

Everything around them is form rather than arithmetic: blst holds a field
element in Montgomery form, so the case's limbs are converted in with
`blst_fp_from_uint64` and the answer's out with `blst_uint64_from_fp`, and the
twelve elements are read off blst's own tower — `fp12.fp6[i].fp2[j].fp[k]` —
in the order the specification fixes.

That the two libraries agree was checked before either was measured: both were
run over the same two cases on one machine, and their twelve field elements
were identical byte for byte. The specification is judged against the arkworks
answer exactly, so a disagreement would have been the first thing a
measurement said rather than a silent difference between curves' conventions.
