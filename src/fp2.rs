use crate::fp::Fp;
use std::ops::{Add, Mul, Neg, Sub};

/// Fp2 represents the quadratic extension field Fp2 = Fp[u] / (u² + 1)
/// where u² = -1
/// An element is represented as c0 + c1*u
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fp2 {
    pub c0: Fp,
    pub c1: Fp,
}

impl Fp2 {
    pub fn new(c0: Fp, c1: Fp) -> Self {
        Fp2 { c0, c1 }
    }

    pub fn zero() -> Self {
        Fp2 {
            c0: Fp::zero(),
            c1: Fp::zero(),
        }
    }

    pub fn one() -> Self {
        Fp2 {
            c0: Fp::one(),
            c1: Fp::zero(),
        }
    }

    /// Conjugate: (a + bu)* = a - bu
    pub fn conjugate(&self) -> Self {
        Fp2 {
            c0: self.c0.clone(),
            c1: -self.c1.clone(),
        }
    }

    /// Inverse: (a + bu)^(-1) = (a - bu) / (a² + b²)
    /// since u² = -1, norm = a² - b²u² = a² + b²
    pub fn inv(&self) -> Self {
        let norm = self.c0.clone() * self.c0.clone() + self.c1.clone() * self.c1.clone();
        let norm_inv = norm.inv();
        Fp2 {
            c0: self.c0.clone() * norm_inv.clone(),
            c1: -self.c1.clone() * norm_inv,
        }
    }
}

/// Addition: (a + bu) + (c + du) = (a + c) + (b + d)u
impl Add for Fp2 {
    type Output = Fp2;
    fn add(self, rhs: Fp2) -> Fp2 {
        Fp2 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}

impl<'a, 'b> Add<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn add(self, rhs: &'b Fp2) -> Fp2 {
        Fp2 {
            c0: self.c0.clone() + rhs.c0.clone(),
            c1: self.c1.clone() + rhs.c1.clone(),
        }
    }
}

/// Subtraction: (a + bu) - (c + du) = (a - c) + (b - d)u
impl Sub for Fp2 {
    type Output = Fp2;
    fn sub(self, rhs: Fp2) -> Fp2 {
        Fp2 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}

impl<'a, 'b> Sub<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn sub(self, rhs: &'b Fp2) -> Fp2 {
        Fp2 {
            c0: self.c0.clone() - rhs.c0.clone(),
            c1: self.c1.clone() - rhs.c1.clone(),
        }
    }
}

/// Multiplication: (a + bu)(c + du) = (ac - bd) + (ad + bc)u
/// using u² = -1
impl Mul for Fp2 {
    type Output = Fp2;
    fn mul(self, rhs: Fp2) -> Fp2 {
        let ac = self.c0.clone() * rhs.c0.clone();
        let bd = self.c1.clone() * rhs.c1.clone();
        let ad_plus_bc = (self.c0.clone() + self.c1.clone()) * (rhs.c0.clone() + rhs.c1.clone())
            - ac.clone()
            - bd.clone();
        Fp2 {
            c0: ac - bd,
            c1: ad_plus_bc,
        }
    }
}

impl<'a, 'b> Mul<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;
    fn mul(self, rhs: &'b Fp2) -> Fp2 {
        let ac = self.c0.clone() * rhs.c0.clone();
        let bd = self.c1.clone() * rhs.c1.clone();
        let ad_plus_bc = (self.c0.clone() + self.c1.clone()) * (rhs.c0.clone() + rhs.c1.clone())
            - ac.clone()
            - bd.clone();
        Fp2 {
            c0: ac - bd,
            c1: ad_plus_bc,
        }
    }
}

/// Negation: -(a + bu) = -a - bu
impl Neg for Fp2 {
    type Output = Fp2;
    fn neg(self) -> Fp2 {
        Fp2 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigUint;
    use rand::Rng;

    #[test]
    fn test_basic_ops() {
        let a = Fp2::new(Fp::new(3u32.into()), Fp::new(5u32.into()));
        let b = Fp2::new(Fp::new(7u32.into()), Fp::new(11u32.into()));

        // Test addition
        let sum = &a + &b;
        assert_eq!(sum.c0, Fp::new(10u32.into()));
        assert_eq!(sum.c1, Fp::new(16u32.into()));

        // Test subtraction
        let diff = &b - &a;
        assert_eq!(diff.c0, Fp::new(4u32.into()));
        assert_eq!(diff.c1, Fp::new(6u32.into()));
    }

    #[test]
    fn test_multiplication() {
        let a = Fp2::new(Fp::new(3u32.into()), Fp::new(5u32.into()));
        let b = Fp2::new(Fp::new(7u32.into()), Fp::new(11u32.into()));

        // (3 + 5u)(7 + 11u) = 21 + 33u + 35u + 55u²
        // = 21 + 68u - 55 (since u² = -1)
        // = -34 + 68u
        let prod = &a * &b;
        let expected_c0 = Fp::new(3u32.into()) * Fp::new(7u32.into())
            - Fp::new(5u32.into()) * Fp::new(11u32.into());
        let expected_c1 = Fp::new(3u32.into()) * Fp::new(11u32.into())
            + Fp::new(5u32.into()) * Fp::new(7u32.into());
        assert_eq!(prod.c0, expected_c0);
        assert_eq!(prod.c1, expected_c1);
    }

    #[test]
    fn test_inverse() {
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let a_val: u64 = rng.gen_range(1..1000);
            let b_val: u64 = rng.gen_range(1..1000);
            let a = Fp2::new(
                Fp::new(a_val.to_biguint().unwrap()),
                Fp::new(b_val.to_biguint().unwrap()),
            );
            let inv = a.inv();
            let prod = &a * &inv;
            assert_eq!(prod.c0, Fp::one());
            assert_eq!(prod.c1, Fp::zero());
        }
    }

    #[test]
    fn test_conjugate() {
        let a = Fp2::new(Fp::new(3u32.into()), Fp::new(5u32.into()));
        let conj = a.conjugate();
        assert_eq!(conj.c0, Fp::new(3u32.into()));
        assert_eq!(conj.c1, -Fp::new(5u32.into()));
    }

    #[test]
    fn test_field_laws() {
        let a = Fp2::new(Fp::new(3u32.into()), Fp::new(5u32.into()));
        let b = Fp2::new(Fp::new(7u32.into()), Fp::new(11u32.into()));
        let c = Fp2::new(Fp::new(13u32.into()), Fp::new(17u32.into()));

        // Associativity of addition
        let assoc1 = &(&a + &b) + &c;
        let assoc2 = &a + &(&b + &c);
        assert_eq!(assoc1, assoc2);

        // Commutativity of addition
        let comm1 = &a + &b;
        let comm2 = &b + &a;
        assert_eq!(comm1, comm2);

        // Associativity of multiplication
        let assoc_mul1 = &(&a * &b) * &c;
        let assoc_mul2 = &a * &(&b * &c);
        assert_eq!(assoc_mul1, assoc_mul2);

        // Commutativity of multiplication
        let comm_mul1 = &a * &b;
        let comm_mul2 = &b * &a;
        assert_eq!(comm_mul1, comm_mul2);

        // Distributivity
        let dist1 = &a * &(&b + &c);
        let dist2 = &(&a * &b) + &(&a * &c);
        assert_eq!(dist1, dist2);
    }
}
