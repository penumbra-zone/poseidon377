use ark_ed_on_bls12_377::{Fq, FqParameters};
use ark_ff::FpParameters;
use poseidon_paramgen::poseidon_build;
use std::{
    env, fs,
    io::{BufWriter, Write},
    path::PathBuf,
};

fn main() {
    // We use the default `OUT_DIR` set by Cargo when a build script exists.
    let output_location: PathBuf =
        PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR environmental variable must be set"))
            .join("params.rs");

    let security_level = 128;
    // Recall t = rate + capacity, so t=2 is rate=1, capacity=1 (i.e. 1:1 hash)
    let t_values = vec![2, 3, 4, 5, 6, 7, 8];
    let params_codegen =
        poseidon_build::compile::<Fq>(security_level, t_values, FqParameters::MODULUS, true);

    let fh = fs::File::create(output_location).expect("can create source file");
    let mut f = BufWriter::new(fh);
    f.write_all(params_codegen.as_bytes())
        .expect("can write parameters to file");
}
