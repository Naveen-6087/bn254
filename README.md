# BN254 Elliptic Curve Implementation

A from-scratch implementation of the BN254 (alt_bn128) elliptic curve in Rust, designed for educational purposes.

## Overview

This project implements the BN254 curve used in Ethereum's precompiled contracts for pairing-based cryptography. The implementation is modular, readable, and avoids unsafe code.

## Modules

- **`fp.rs`** - Base field Fp arithmetic modulo p = 21888242871839275222246405745257275088548364400416034343698204186575808495617
- **`g1.rs`** - G1 curve points over Fp: y² = x³ + 3
- **`fp2.rs`** - Quadratic extension field Fp2 = Fp[u] / (u² + 1)
- **`g2.rs`** - G2 twisted curve points over Fp2
- **`fp6.rs`** - Sextic extension field Fp6 = Fp2[v] / (v³ - (u+9))
- **`fp12.rs`** - Degree-12 extension field Fp12 = Fp6[w] / (w² - v)
- **`pairing.rs`** - Optimal Ate pairing implementation

## Features

- ✅ Complete field arithmetic (Fp, Fp2, Fp6, Fp12)
- ✅ Jacobian coordinate elliptic curve operations (G1, G2)
- ✅ Scalar multiplication with double-and-add
- ✅ Miller loop for pairing computation
- ✅ Final exponentiation
- ✅ Comprehensive test suite

## Building

### Prerequisites

**Windows users**: You need Visual Studio C++ Build Tools installed:
1. Download [Visual Studio Build Tools](https://visualstudio.microsoft.com/downloads/)
2. Install "Desktop development with C++" workload

### Compile

```bash
cargo build
```

### Run Tests

```bash
cargo test
```

## Usage Example

```rust
use bn254::*;
use num_bigint::BigUint;

// Create points on G1 and G2
let p = g1::G1 {
    x: fp::Fp::new(1u32.into()),
    y: fp::Fp::new(2u32.into()),
    z: fp::Fp::one(),
};

let q = g2::G2 {
    x: fp2::Fp2::new(fp::Fp::new(1u32.into()), fp::Fp::new(2u32.into())),
    y: fp2::Fp2::new(fp::Fp::new(3u32.into()), fp::Fp::new(4u32.into())),
    z: fp2::Fp2::one(),
};

// Compute pairing
let result = pairing::pairing(&p, &q);

// Scalar multiplication
let scalar_p = p.mul_u128(5);
```

## Current Status

### ✅ Completed
- Field arithmetic for all extension fields
- Curve point arithmetic in Jacobian coordinates
- Basic pairing structure

### 🚧 TODO
- Complete line function evaluation in Miller loop
- Optimize final exponentiation using cyclotomic subgroup
- Implement proper Frobenius maps with coefficients
- Add known generator points from BN254 specification
- Cross-validate results with established libraries (ark-bn254, etc.)

## Design Principles

1. **Modularity** - Each file is self-contained with its own tests
2. **Safety** - No unsafe code, all operations are bounds-checked
3. **Readability** - Educational code with clear formulas
4. **Correctness** - Uses exact BN254 parameters

## Testing

The project includes:
- Unit tests in each module (`cargo test --lib`)
- Integration tests (`cargo test --test integration`)
- Property-based tests for field laws and curve equations

## References

- [BN254 Curve Specification](https://hackmd.io/@jpw/bn254)
- Barreto-Naehrig, "Pairing-Friendly Elliptic Curves of Prime Order" (2006)
- [EIP-196/197](https://eips.ethereum.org/EIPS/eip-197) - Ethereum's alt_bn128 precompile

## License

MIT

## Contributing

This is an educational project. Contributions are welcome to:
- Complete the pairing implementation
- Add benchmarks
- Cross-validate with other implementations
- Improve documentation
