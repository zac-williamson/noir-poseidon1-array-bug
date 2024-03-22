// Instantiations of Poseidon constants, permutations and sponge for prime field of the same order as BN254
mod perm;
mod consts;
mod optconsts;

use crate::poseidon::PoseidonOptConfig;
use crate::poseidon::apply_opt_matrix;
use crate::poseidon::apply_cheap_matrix;

#[field(bn254)]
fn Sigma<O>(x: [Field; O]) -> [Field; O] {
    let mut y = x;
    for i in 0..O {
        let t = y[i];
        let tt = t * t;
        let tttt = tt * tt;
        y[i] *= tttt;
    }
    y
}

#[field(bn254)]
fn Ark<O, N>(C: [Field; N], r: u64, x: [Field; O]) -> [Field; O] {
    let mut y = x;
    for i in 0..O {
        y[i] += C[i + r];
    }
    y
}

#[field(bn254)]
fn Mix<O, M>(Mat: [[Field; M]; M], x: [Field; O]) -> [Field; O] {
    let mut y = x;
    y[0] *= Mat[0][0];
    y[0] += x[1] * Mat[0][1] + x[2] * Mat[0][2];
    y[1] *= Mat[1][1];
    y[1] += x[0] * Mat[1][0];
    y[1] += x[2] * Mat[1][2];
    y[2] *= Mat[2][2];
    y[2] += x[0] * Mat[2][0];
    y[2] += x[1] * Mat[2][1];
    // for i in 0..O {
    //     let mut lc = 0;
    //     for j in 0..O {
    //         lc += x[j] * Mat[i][j];
    //     }
    //     y[i] = lc;
    // }
    y
}

#[field(bn254)]
fn MixLast<O, M>(Mat: [[Field; M]; M], x: [Field; O], s: u64) -> Field {
    let mut lc = 0;
    for j in 0..O {
        lc += Mat[j][s] * x[j];
    }
    lc
}

// A*x where A is an n x n matrix in row-major order and x an n-vector
fn MixSOld<N, M>(s: [Field; M], x: [Field; N], r: u64) -> [Field; N] {
    let mut y = x;
    for i in 0..x.len() {
        y[0] += (x[i]* s[(x.len()*2-1)*r+i]);
    }
    for i in 1..x.len() {
        y[i] = x[i] + (y[0] * s[(x.len()*2-1)*r + x.len() + i -1]); //  * x[i];
    }

    y
}

fn MixS<O, M>(s: [Field; M], x: [Field; O], r: u64) -> [Field; O] {
    let mut y = x;
    y[0] += x[0] * s[(x.len()*2-1)*r+0];
    y[0] += x[1] * s[(x.len()*2-1)*r+1];
    y[0] += x[2] * s[(x.len()*2-1)*r+2];

    y[1] = x[1] + (y[0] * s[(x.len()*2-1)*r + x.len() + 0 -1]); //  * x[i];
    y[2] = x[2] + (y[0] * s[(x.len()*2-1)*r + x.len() + 1 -1]); //  * x[i];

    // for i in 0..x.len() {
    //     y[0] += (x[i] * s[(x.len()*2-1)*r+i]);
    // }
    // for i in 1..x.len() {
    //     y[i] = x[i] + (y[0] * s[(x.len()*2-1)*r + x.len() + i -1]); //  * x[i];
    // }

    y
}

// Optimised permutation for this particular field; uses hardcoded rf and rp values,
// which should agree with those in pos_conf.
#[field(bn254)]
pub fn permute<M, N, X, O>(pos_conf: PoseidonOptConfig<M, N, X>, mut state: [Field; O]) -> [Field; O] {
    let PoseidonOptConfig {t, rf: config_rf, rp: config_rp, alpha, C, Mat, S } = pos_conf;
    let rf: u8 = 8;
    let rp: u8 = [56, 57, 56, 60, 60, 63, 64, 63, 60, 66, 60, 65, 70, 60, 64, 68][state.len() - 2];
    let mut s0 = state[0];
    let mut s1 = state[1];
    let mut s2 = state[2];
    for _r in 0..31 {
        s0 += C[((rf as u64 / 2)+1)*state.len() + _r as u64];
        let t = s0;
        let tt = t * t;
        let tttt = tt * tt;
        s0 *= tttt;
        s0 += s0 * S[(5)*_r ];
        s0 += s1 * S[(5)*_r +1];
        s0 += s2 * S[(5)*_r +2];
        s1 += s0 * S[5 * _r +2];
        s2 += s0 * S[5*_r  + 3];
    }
    state[0] = s0;
    state[1] = s1;
    state[2] = s2;
    state
}

#[field(bn254)]
pub fn permutebad<M, N, X, O>(pos_conf: PoseidonOptConfig<M, N, X>, mut state: [Field; O]) -> [Field; O] {
    let PoseidonOptConfig {t, rf: config_rf, rp: config_rp, alpha, C, Mat, S } = pos_conf;
    let rf: u8 = 8;
    let rp: u8 = [56, 57, 56, 60, 60, 63, 64, 63, 60, 66, 60, 65, 70, 60, 64, 68][state.len() - 2];

    // assert(t == state.len() as Field);
    // assert(rf == config_rf);
    // assert(rp == config_rp);

    //   state = Ark(C, 0, state);

    // for _r in 0..rf / 2 {
    //     state = Sigma(state);
    //     //     state = Ark(C, (_r as u64 + 1) * state.len(), state);
    //     //   state = Mix(Mat, state);
    // }

    let mut s0 = state[0];
    let mut s1 = state[1];
    let mut s2 = state[2];
    for _r in 0..31 {
        s0 += C[((rf as u64 / 2)+1)*state.len() + _r as u64];
        let t = s0;
        let tt = t * t;
        let tttt = tt * tt;
        s0 *= tttt;
        // if (_r > 0) {
        // let mut y = state;
        s0 += s0 * S[(5)*_r ];
        s0 += s1 * S[(5)*_r +1];
        s0 += s2 * S[(5)*_r +2];
        s1 +=  s0 * S[5 * _r  + 2];
        s2 += s0 * S[5*_r  + 3]; //  * x[i];
        //}
    }
    state[0] = s0;
    state[1] = s1;
    state[2] = s2;
    // let count : u64 = (rf as u64 / 2 + 1) * state.len() + rp as u64;
    for _r in 0..rf / 2 - 1 {
        // let b : u64 = rp as u64;
        // let c : u64 = _r as u64 * state.len();
        //    state = Ark(C, (rf as u64 / 2 + 1) * state.len() + rp as u64 + (_r as u64 * state.len()), state);
        state = Sigma(state);
        //   state = Mix(Mat, state);
    }
    state = Sigma(state);
    //   state = Mix(Mat, state);

    // let mut count = 0;
    // // First half of full rounds
    // for _r in 0..rf / 2 {
    //     for i in 0..state.len() {
    //         state[i] = state[i] + C[count + i];
    //     } // Shift by round constants
    //     for i in 0..state.len() {
    //         state[i] = state[i].pow_32(alpha);
    //     }

    //     state = apply_opt_matrix(Mat, state); // Apply MDS matrix
    //     count = count + t as u64;
    // }
    // // Partial rounds
    // for _r in 0..rp {
    //     state[0] = state[0] + C[count];
    //     let x = state[0];
    //     state[0] *= state[0];
    //     state[0] *= state[0];
    //     state[0] *= x;
    //     //  state[0] = state[0].pow_32(alpha);
    //     if (_r > 0) {
    //         state = apply_cheap_matrix(S, state, _r as u64); // Apply MDS matrix
    //     }
    //     count = count + t as u64;
    // }
    // // Second half of full rounds
    // for _r in 0..rf / 2 {
    //     for i in 0..state.len() {
    //         state[i] = state[i] + C[count + i];
    //     } // Shift by round constants
    //     for i in 0..state.len() {
    //         state[i] = state[i].pow_32(alpha);
    //     }

    //     state = apply_opt_matrix(Mat, state); // Apply MDS matrix
    //     count = count + t as u64;
    // }

    state
}

// Corresponding absorption.
#[field(bn254)]
fn absorb<M, N, X, O, P>(
    pos_conf: PoseidonOptConfig<M, N, X>,
    mut state: [Field; O], // Initial state; usually [0; O]
    rate: Field, // Rate
    capacity: Field, // Capacity; usually 1
    msg: [Field; P] // Arbitrary length message
) -> [Field; O] {
    assert(pos_conf.t == rate + capacity);

    let mut i = 0;

    for k in 0..msg.len() {
        // Add current block to state
        state[capacity + i] += msg[k];
        i = i+1;
        // Enough to absorb
        if i == rate {
            state = permute(pos_conf, state);
            i = 0;
        }
    }
    // If we have one more block to permute
    if i != 0 {
        state = permute(pos_conf, state);
    }

    state
}
// Variable-length Poseidon-128 sponge as suggested in second bullet point of §3 of https://eprint.iacr.org/2019/458.pdf
#[field(bn254)]
pub fn sponge<N>(msg: [Field; N]) -> Field {
    absorb(optconsts::x5_5_opt_config(), [0; 5], 4, 1, msg)[1]
}
// Various instances of the Poseidon hash function
// Consistent with Circom's implementation
pub fn hash_1(input: [Field; 1]) -> Field {
    let mut state = [0; 2];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_2(state)[0]
}

pub fn hash_2(input: [Field; 2]) -> Field {
    let mut state = [0; 3];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_3(state)[0]
}

pub fn hash_3(input: [Field; 3]) -> Field {
    let mut state = [0; 4];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_4(state)[0]
}

pub fn hash_4(input: [Field; 4]) -> Field {
    let mut state = [0; 5];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_5(state)[0]
}

pub fn hash_5(input: [Field; 5]) -> Field {
    let mut state = [0; 6];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_6(state)[0]
}

pub fn hash_6(input: [Field; 6]) -> Field {
    let mut state = [0; 7];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_7(state)[0]
}

pub fn hash_7(input: [Field; 7]) -> Field {
    let mut state = [0; 8];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_8(state)[0]
}

pub fn hash_8(input: [Field; 8]) -> Field {
    let mut state = [0; 9];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_9(state)[0]
}

pub fn hash_9(input: [Field; 9]) -> Field {
    let mut state = [0; 10];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_10(state)[0]
}

pub fn hash_10(input: [Field; 10]) -> Field {
    let mut state = [0; 11];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_11(state)[0]
}

pub fn hash_11(input: [Field; 11]) -> Field {
    let mut state = [0; 12];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_12(state)[0]
}

pub fn hash_12(input: [Field; 12]) -> Field {
    let mut state = [0; 13];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_13(state)[0]
}

pub fn hash_13(input: [Field; 13]) -> Field {
    let mut state = [0; 14];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_14(state)[0]
}

pub fn hash_14(input: [Field; 14]) -> Field {
    let mut state = [0; 15];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_15(state)[0]
}

pub fn hash_15(input: [Field; 15]) -> Field {
    let mut state = [0; 16];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_16(state)[0]
}

pub fn hash_16(input: [Field; 16]) -> Field {
    let mut state = [0; 17];
    for i in 0..input.len() {
        state[i+1] = input[i];
    }

    perm::x5_17(state)[0]
}
