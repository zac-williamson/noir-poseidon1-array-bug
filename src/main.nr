mod poseidon;

fn main(x: Field, y: pub Field) {
    let state = poseidon::bn254::permute(poseidon::bn254::optconsts::x5_3_opt_config(), [x, y, x + y]);
    // let z = poseidon::bn254::hash_3([x, y, x + y]);
    assert(state[0] != y);
}

#[test]
fn test_main() {
    main(1, 2);
    // Uncomment to make test fail
    // main(1, 1);
}
