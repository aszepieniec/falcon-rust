use rand::{rng, RngExt};

// Run with:
// `cargo bench --bench profile --features profiling`

fn main() {
    let mut rng = rng();

    #[cfg(feature = "profiling")]
    {
        println!("\n=== Falcon Key Gen ===");
        let (sk, pk) = falcon_rust::falcon1024::keygen(rng.random());
        falcon_rust::profiling::print_summary(
            6,
            &["ffldl", "gram_schmidt_norm_squared", "babai_reduce_i32"],
        );
        falcon_rust::profiling::reset();

        println!("\n=== Falcon Sign ===");
        let m = rng.random::<[u8; 15]>();
        let sig = falcon_rust::falcon1024::sign(&m, &sk);
        falcon_rust::profiling::print_summary(5, &["ffsampling"]);
        falcon_rust::profiling::reset();

        println!("\n=== Falcon Verify ===");
        falcon_rust::falcon1024::verify(&m, &sig, &pk);
        falcon_rust::profiling::print_summary(5, &[]);
        falcon_rust::profiling::reset();
    }

    #[cfg(not(feature = "profiling"))]
    {
        println!("Profiling feature not enabled. Run with: cargo bench --bench profile --features profiling");
        let _ = rng.random_bool(0.5); // avoid unused variable warning
    }
}
