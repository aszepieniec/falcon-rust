//! Deterministic keygen profiler: mean span breakdown over many distinct seeds.
//!
//! Run: cargo run --release --features profiling --example profile_keygen [N]
//!
//! Flattens each per-keygen span tree by function name (aggregating recursive
//! nodes across depths) and reports mean ms/keygen + % of total, sorted. This
//! identifies the dominant *function*; depths are conflated on purpose.

use std::collections::HashMap;

use falcon_rust::falcon1024;
use falcon_rust::profiling;

fn flatten(node: &profiling::Span, acc: &mut HashMap<String, (f64, usize)>) {
    let ms = node.duration.as_secs_f64() * 1000.0;
    let e = acc.entry(node.name.clone()).or_insert((0.0, 0));
    e.0 += ms;
    e.1 += node.call_count;
    for c in &node.children {
        flatten(c, acc);
    }
}

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(60);

    let mut agg: HashMap<String, (f64, usize)> = HashMap::new();
    let mut keygen_total_ms = 0.0;

    for i in 0u64..n as u64 {
        let mut seed = [0u8; 32];
        seed[0..8].copy_from_slice(&i.to_le_bytes());

        profiling::reset();
        let _ = falcon1024::keygen(seed);

        profiling::BUILDER.with(|b| {
            let b = b.borrow();
            for root in &b.roots {
                if root.name == "keygen" {
                    keygen_total_ms += root.duration.as_secs_f64() * 1000.0;
                }
                flatten(root, &mut agg);
            }
        });
    }

    let mean_keygen = keygen_total_ms / n as f64;
    println!("\nkeygen(1024) mean over {n} distinct seeds: {mean_keygen:.2} ms\n");

    let mut rows: Vec<(String, f64, usize)> = agg
        .into_iter()
        .map(|(name, (ms, calls))| (name, ms / n as f64, calls))
        .collect();
    rows.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    println!(
        "{:<34} {:>10} {:>8} {:>10}",
        "function (flattened by name)", "ms/keygen", "%total", "calls/kg"
    );
    println!("{}", "-".repeat(66));
    for (name, ms, calls) in rows.iter().take(30) {
        let pct = ms / mean_keygen * 100.0;
        let calls_per = *calls as f64 / n as f64;
        println!("{name:<34} {ms:>10.3} {pct:>7.1}% {calls_per:>10.1}");
    }
}
