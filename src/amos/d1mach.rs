pub(crate) fn d1mach(i: i32) -> f64 {
    match i {
        1 => f64::MIN_POSITIVE,              // Smallest positive normalized number (~2.225e-308)
        2 => f64::MAX,                       // Largest finite number (~1.797e308)
        3 => f64::EPSILON / 2.0,             // Smallest relative spacing (~1.110e-16)
        4 => f64::EPSILON,                   // Machine epsilon / largest relative spacing (~2.220e-16)
        5 => std::f64::consts::LOG10_2,      // Log10 of the radix (approx 0.30103)
        _ => panic!("D1MACH -- I out of bounds: {}", i),
    }
}