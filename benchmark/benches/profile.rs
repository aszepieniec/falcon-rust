// use falcon_rust::math::ntru_gen;
use falcon_rust::falcon1024::keygen;
use falcon_rust::falcon1024::sign;
use falcon_rust::falcon1024::verify;
use falcon_rust::profiling;
use rand::{rng, Rng};

// Run with:
// `cargo bench --bench profile --features profiling`

fn main() {
    let mut rng = rng();

    #[cfg(feature = "profiling")]
    {
        println!("\n=== Falcon Key Gen ===");
        let (sk, pk) = keygen(rng.random());
        profiling::print_summary(
            6,
            &["ffldl", "gram_schmidt_norm_squared", "babai_reduce_i32"],
        );
        profiling::reset();

        println!("\n=== Falcon Sign ===");
        let m = rng.random::<[u8; 15]>();
        let sig = sign(&m, &sk);
        profiling::print_summary(5, &["ffsampling"]);
        profiling::reset();

        println!("\n=== Falcon Verify ===");
        verify(&m, &sig, &pk);
        profiling::print_summary(5, &[]);
        profiling::reset();
    }

    #[cfg(not(feature = "profiling"))]
    {
        println!("Profiling feature not enabled. Run with: cargo bench --bench profile --features profiling");
    }
}
