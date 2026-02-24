use std::cell::RefCell;
use std::time::Duration;
use std::time::Instant;

// The "Node" in the profile tree
#[derive(Debug, Clone)]
pub struct Span {
    pub name: String,
    pub duration: Duration,
    pub children: Vec<Span>,
    pub call_count: usize,
}

// The "State" that tracks the current nesting
pub struct ProfileBuilder {
    stack: Vec<(String, Instant, Vec<Span>, usize)>,
    pub roots: Vec<Span>,
}

thread_local! {
    // Each thread gets its own independent tree builder
    pub static BUILDER: RefCell<ProfileBuilder> = const { RefCell::new(ProfileBuilder {
        stack: Vec::new(),
        roots: Vec::new(),
    }) };
}

pub fn push(name: &str) {
    BUILDER.with(|b| {
        let mut b = b.borrow_mut();
        // Check if this function is already anywhere in the current stack
        let is_recursive = b.stack.iter().any(|(n, _, _, _)| n == name);

        if is_recursive {
            // Increment the recursion depth of the existing entry
            // We just push a "marker" to keep the stack balanced with the pops
            b.stack
                .push((name.to_string(), Instant::now(), Vec::new(), 1));
        } else {
            // First time seeing this function in this branch
            b.stack
                .push((name.to_string(), Instant::now(), Vec::new(), 0));
        }
    });
}

pub fn pop() {
    BUILDER.with(|b| {
        let mut b = b.borrow_mut();
        if let Some((name, start, children, depth)) = b.stack.pop() {
            let elapsed = start.elapsed();

            if depth > 0 {
                // This was a recursive call.
                // We find the "original" caller of this name and add our time/count to it.
                if let Some(_original) = b
                    .stack
                    .iter_mut()
                    .rev()
                    .find(|(n, _, _, d)| n == &name && *d == 0)
                {
                    // Note: We don't add elapsed time here because the "original"
                    // timer is still running and will encompass this time.
                    // We can, however, track that a recursive call occurred.
                }
                // We also need to move any children collected here to the parent
                // so we don't lose spans triggered inside the recursion
                if let Some(parent) = b.stack.last_mut() {
                    parent.2.extend(children);
                }
            } else {
                // This was the base caller.
                let span = Span {
                    name,
                    duration: elapsed,
                    children,
                    call_count: 1,
                };

                if let Some(parent) = b.stack.last_mut() {
                    parent.2.push(span);
                } else {
                    b.roots.push(span);
                }
            }
        }
    });
}

pub fn reset() {
    BUILDER.with(|b| {
        let mut b = b.borrow_mut();
        b.stack.clear();
        b.roots.clear();
    });
}

// --- The Core Logic ---

/// Calculates how many visual columns a string occupies.
/// This ignores the fact that a '├' is 3 bytes and counts it as 1 column.
fn visual_len(s: &str) -> usize {
    s.chars().count()
}

/// Helper to format duration with aligned decimals and ASCII-safe units.
/// We use 'us' instead of 'µs' to ensure 1 char = 1 column across all terminals.
fn format_duration(duration: Duration) -> (String, String, &'static str) {
    let secs = duration.as_secs_f64();
    let (val, unit) = if secs >= 1.0 {
        (secs, "s")
    } else if secs >= 0.001 {
        (secs * 1000.0, "ms")
    } else if secs >= 0.000_001 {
        (secs * 1_000_000.0, "us")
    } else {
        (secs * 1_000_000_000.0, "ns")
    };

    let s = format!("{:.2}", val);
    let parts: Vec<&str> = s.split('.').collect();
    (parts[0].to_string(), parts[1].to_string(), unit)
}

pub struct PrintConfig<'a> {
    pub max_depth: usize,
    pub collapse_names: &'a [&'a str],
    pub total_ms: f64,
    pub max_name_width: usize,
}

fn render_node(
    node: &Span,
    prefix: &str,
    is_last: bool,
    config: &PrintConfig,
    current_depth: usize,
) {
    let marker = if is_last { "└── " } else { "├── " };
    let name_line = format!("{}{}{}", prefix, marker, node.name);
    let name_line = if node.call_count > 1 {
        format!("{} (x{})", name_line, node.call_count)
    } else {
        name_line
    };

    // 1. Calculate manual padding
    // We count characters, not bytes.
    let v_len = visual_len(&name_line);
    let padding_count = config.max_name_width.saturating_sub(v_len);
    let padding = " ".repeat(padding_count);

    // 2. Format Time Parts
    let (whole, decimal, unit) = format_duration(node.duration);
    let node_ms = node.duration.as_secs_f64() * 1000.0;
    let percentage = if config.total_ms > 0.0 {
        (node_ms / config.total_ms) * 100.0
    } else {
        0.0
    };

    // 3. Construct the line.
    // We use NO width modifiers on the strings containing tree markers.
    println!(
        "{}{}  {:>6}.{:<2} {:<3} {:>8.1}%",
        name_line, padding, whole, decimal, unit, percentage
    );

    // --- Recursion ---
    let is_collapsed = config.collapse_names.contains(&node.name.as_str());
    if !node.children.is_empty() {
        let next_prefix = format!("{}{}", prefix, if is_last { "    " } else { "│   " });
        if current_depth < config.max_depth && !is_collapsed {
            let len = node.children.len();
            for (i, child) in node.children.iter().enumerate() {
                render_node(child, &next_prefix, i == len - 1, config, current_depth + 1);
            }
        } else {
            let reason = if is_collapsed {
                "collapsed"
            } else {
                "depth limit"
            };
            let hint = format!(
                "{}└── ... ({} nodes hidden via {})",
                next_prefix,
                node.children.len(),
                reason
            );
            println!("{}", hint);
        }
    }
}

// --- Measurement Pass ---

fn calculate_max_width(
    node: &Span,
    prefix_len: usize,
    max_w: &mut usize,
    current_depth: usize,
    max_depth: usize,
    collapse_names: &[&str],
) {
    // 4 is the length of "├── " or "└── "
    let mut current_line_width = prefix_len + 4 + node.name.chars().count();
    if node.call_count > 1 {
        current_line_width += format!(" (x{})", node.call_count).chars().count();
    }

    *max_w = std::cmp::max(*max_w, current_line_width);

    let is_collapsed = collapse_names.contains(&node.name.as_str());
    if current_depth < max_depth && !is_collapsed {
        for child in &node.children {
            calculate_max_width(
                child,
                prefix_len + 4,
                max_w,
                current_depth + 1,
                max_depth,
                collapse_names,
            );
        }
    }
}

pub fn print_summary(max_depth: usize, collapse_names: &[&str]) {
    BUILDER.with(|b| {
        let b = b.borrow();
        if b.roots.is_empty() {
            return;
        }

        let mut max_name_width = 0;
        for root in &b.roots {
            calculate_max_width(root, 0, &mut max_name_width, 0, max_depth, collapse_names);
        }

        let total_time: Duration = b.roots.iter().map(|r| r.duration).sum();
        let config = PrintConfig {
            max_depth,
            collapse_names,
            total_ms: total_time.as_secs_f64() * 1000.0,
            max_name_width,
        };

        for (i, root) in b.roots.iter().enumerate() {
            render_node(root, "", i == b.roots.len() - 1, &config, 0);
        }
    });
}
