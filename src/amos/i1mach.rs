pub(crate) fn i1mach(i: i32) -> i32 {
    match i {
        1 => 5,             // Standard input unit
        2 => 6,             // Standard output unit
        3 => 7,             // Standard punch unit (often 0 or 7)
        4 => 6,             // Standard error message unit
        5 => 32,            // Number of bits per integer storage unit
        6 => 4,             // Number of characters per integer storage unit
        7 => 2,             // Base of integer arithmetic (A)
        8 => 31,            // Number of base-A digits (S)
        9 => i32::MAX,      // Largest integer magnitude (2^31 - 1)
        10 => 2,            // Base of floating-point arithmetic (B)
        11 => 24,           // Number of base-B digits for single precision (T)
        12 => -125,         // Smallest exponent for single precision (EMIN)
        13 => 128,          // Largest exponent for single precision (EMAX)
        14 => 53,           // Number of base-B digits for double precision (T)
        15 => -1021,        // Smallest exponent for double precision (EMIN)
        16 => 1024,         // Largest exponent for double precision (EMAX)
        _ => panic!("I1MACH -- I out of bounds: {}", i),
    }
}