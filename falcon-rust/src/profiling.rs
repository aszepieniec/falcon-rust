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

pub fn print_summary() {
    BUILDER.with(|b| {
        let b = b.borrow();
        if b.roots.is_empty() {
            return;
        }

        // Pass 1: Calculate the maximum width of the name/prefix column
        let mut max_w = 0;
        for root in &b.roots {
            calculate_width(root, 0, &mut max_w);
        }
        // Add a small buffer for the markers (├── )
        max_w += 4;

        // Pass 2: Calculate total time for percentages
        let total_time: Duration = b.roots.iter().map(|r| r.duration).sum();
        let total_ms = total_time.as_secs_f64() * 1000.0;

        // Pass 3: Render
        for (i, root) in b.roots.iter().enumerate() {
            let is_last = i == b.roots.len() - 1;
            render_node(root, "", is_last, total_ms, max_w);
        }
    });
}

fn calculate_width(node: &Span, current_indent: usize, max_w: &mut usize) {
    let mut name_len = node.name.len();
    if node.call_count > 1 {
        name_len += format!(" (x{})", node.call_count).len();
    }

    let total_row_width = current_indent + name_len;
    if total_row_width > *max_w {
        *max_w = total_row_width;
    }

    for child in &node.children {
        // Each level adds 4 characters of indentation ("│   ")
        calculate_width(child, current_indent + 4, max_w);
    }
}

fn format_duration(duration: Duration) -> (String, &'static str) {
    let secs = duration.as_secs_f64();
    if secs >= 1.0 {
        (format!("{:.2}", secs), "s")
    } else if secs >= 0.001 {
        (format!("{:.2}", secs * 1000.0), "ms")
    } else if secs >= 0.000_001 {
        (format!("{:.2}", secs * 1_000_000.0), "µs")
    } else {
        (format!("{:.2}", secs * 1_000_000_000.0), "ns")
    }
}

fn render_node(node: &Span, prefix: &str, is_last: bool, total_ms: f64, max_w: usize) {
    let marker = if is_last { "└── " } else { "├── " };

    let mut name_display = format!("{}{}{}", prefix, marker, node.name);
    if node.call_count > 1 {
        name_display.push_str(&format!(" (x{})", node.call_count));
    }

    // 1. Format duration into (value_string, unit)
    let (val_str, unit) = format_duration(node.duration);

    // 2. Split the value into whole and decimal for alignment
    // We expect "123.45". We want to pad the whole number part.
    let parts: Vec<&str> = val_str.split('.').collect();
    let whole = parts[0];
    let decimal = parts[1];

    // 3. Calculate percentages
    let node_ms = node.duration.as_secs_f64() * 1000.0;
    let percentage = if total_ms > 0.0 {
        (node_ms / total_ms) * 100.0
    } else {
        0.0
    };

    // 4. Print with specific alignment:
    // {:<width$} -> The Tree/Name column
    // {:>6}      -> The whole number part (padded to 6 digits)
    // .{:>2}     -> The decimal part
    // {:<3}      -> The unit with a space
    // {:>8.1}%   -> The percentage
    println!(
        "{:<width$} {:>6}.{:>2} {:<3} {:>8.1}%",
        name_display,
        whole,
        decimal,
        unit,
        percentage,
        width = max_w
    );

    let new_prefix = format!("{}{}", prefix, if is_last { "    " } else { "│   " });
    let len = node.children.len();
    for (i, child) in node.children.iter().enumerate() {
        let child_is_last = i == len - 1;
        render_node(child, &new_prefix, child_is_last, total_ms, max_w);
    }
}

pub fn reset() {
    BUILDER.with(|b| {
        let mut b = b.borrow_mut();
        b.stack.clear();
        b.roots.clear();
    });
}
