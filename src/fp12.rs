use crate::fp2::Fp2;
use crate::fp6::Fp6;
use num_bigint::BigUint;
use num_traits::Zero;
use std::ops::{Add, Mul, Neg, Sub};

/// Fp12 represents the degree-12 extension Fp12 = Fp6[w] / (w² - v)
/// An element is represented as c0 + c1*w
/// where w² = v (a non-residue in Fp6)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
}

impl Fp12 {
    pub fn new(c0: Fp6, c1: Fp6) -> Self {
        Fp12 { c0, c1 }
    }

    pub fn zero() -> Self {
        Fp12 {
            c0: Fp6::zero(),
            c1: Fp6::zero(),
        }
    }

    pub fn one() -> Self {
        Fp12 {
            c0: Fp6::one(),
            c1: Fp6::zero(),
        }
    }

    /// Non-residue for Fp12: v in Fp6 = (0, 1, 0)
    #[allow(dead_code)]
    fn non_residue() -> Fp6 {
        Fp6::new(Fp2::zero(), Fp2::one(), Fp2::zero())
    }

    /// Multiply by non-residue (multiply by v)
    fn mul_by_non_residue(a: &Fp6) -> Fp6 {
        // (c0, c1, c2) * v = (c2*ξ, c0, c1) where ξ = u+9
        Fp6::new(
            Fp6::mul_by_non_residue(&a.c2),
            a.c0.clone(),
            a.c1.clone(),
        )
    }

    pub fn inv(&self) -> Self {
        // (c0 + c1*w)^(-1) = (c0 - c1*w) / (c0² - c1²*v)
        let c0_sq = &self.c0 * &self.c0;
        let c1_sq = &self.c1 * &self.c1;
        let t = &c0_sq - &Self::mul_by_non_residue(&c1_sq);
        let t_inv = t.inv();

        Fp12 {
            c0: &self.c0 * &t_inv,
            c1: &(-self.c1.clone()) * &t_inv,
        }
    }

    /// Frobenius endomorphism
    pub fn frobenius_map(&self, _power: usize) -> Self {
        // TODO: Implement proper Frobenius map with frobenius coefficients
        // For now, this is a placeholder
        self.clone()
    }

    /// Exponentiation
    pub fn pow(&self, exp: &BigUint) -> Self {
        let mut res = Self::one();
        let mut base = self.clone();
        let mut e = exp.clone();

        while !e.is_zero() {
            if e.bit(0) {
                res = &res * &base;
            }
            base = &base * &base;
            e >>= 1;
        }

        res
    }
}

impl Add for Fp12 {
    type Output = Fp12;
    fn add(self, rhs: Fp12) -> Fp12 {
        Fp12 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}

impl<'a, 'b> Add<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;
    fn add(self, rhs: &'b Fp12) -> Fp12 {
        Fp12 {
            c0: &self.c0 + &rhs.c0,
            c1: &self.c1 + &rhs.c1,
        }
    }
}

impl Sub for Fp12 {
    type Output = Fp12;
    fn sub(self, rhs: Fp12) -> Fp12 {
        Fp12 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}

impl<'a, 'b> Sub<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;
    fn sub(self, rhs: &'b Fp12) -> Fp12 {
        Fp12 {
            c0: &self.c0 - &rhs.c0,
            c1: &self.c1 - &rhs.c1,
        }
    }
}

impl Mul for Fp12 {
    type Output = Fp12;
    fn mul(self, rhs: Fp12) -> Fp12 {
        // (a0 + a1*w)(b0 + b1*w) = (a0*b0 + a1*b1*v) + (a0*b1 + a1*b0)*w
        let a0_b0 = &self.c0 * &rhs.c0;
        let a1_b1 = &self.c1 * &rhs.c1;
        let a0_b1 = &self.c0 * &rhs.c1;
        let a1_b0 = &self.c1 * &rhs.c0;

        Fp12 {
            c0: &a0_b0 + &Fp12::mul_by_non_residue(&a1_b1),
            c1: &a0_b1 + &a1_b0,
        }
    }
}

impl<'a, 'b> Mul<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;
    fn mul(self, rhs: &'b Fp12) -> Fp12 {
        // (a0 + a1*w)(b0 + b1*w) = (a0*b0 + a1*b1*v) + (a0*b1 + a1*b0)*w
        let a0_b0 = &self.c0 * &rhs.c0;
        let a1_b1 = &self.c1 * &rhs.c1;
        let a0_b1 = &self.c0 * &rhs.c1;
        let a1_b0 = &self.c1 * &rhs.c0;

        Fp12 {
            c0: &a0_b0 + &Fp12::mul_by_non_residue(&a1_b1),
            c1: &a0_b1 + &a1_b0,
        }
    }
}

impl Neg for Fp12 {
    type Output = Fp12;
    fn neg(self) -> Fp12 {
        Fp12 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp::Fp;

    #[test]
    fn test_basic_ops() {
        let a = Fp12::new(
            Fp6::new(
                Fp2::new(Fp::new(1u32.into()), Fp::new(2u32.into())),
                Fp2::new(Fp::new(3u32.into()), Fp::new(4u32.into())),
                Fp2::new(Fp::new(5u32.into()), Fp::new(6u32.into())),
            ),
            Fp6::new(
                Fp2::new(Fp::new(7u32.into()), Fp::new(8u32.into())),
                Fp2::new(Fp::new(9u32.into()), Fp::new(10u32.into())),
                Fp2::new(Fp::new(11u32.into()), Fp::new(12u32.into())),
            ),
        );

        let b = Fp12::one();
        let _sum = &a + &b;
        let prod = &a * &b;
        assert_eq!(prod, a);
    }

    #[test]
    fn test_inverse() {
        let a = Fp12::new(
            Fp6::new(
                Fp2::new(Fp::new(1u32.into()), Fp::new(2u32.into())),
                Fp2::new(Fp::new(3u32.into()), Fp::new(4u32.into())),
                Fp2::new(Fp::new(5u32.into()), Fp::new(6u32.into())),
            ),
            Fp6::new(
                Fp2::new(Fp::new(7u32.into()), Fp::new(8u32.into())),
                Fp2::new(Fp::new(9u32.into()), Fp::new(10u32.into())),
                Fp2::new(Fp::new(11u32.into()), Fp::new(12u32.into())),
            ),
        );
        let inv = a.inv();
        let prod = &a * &inv;
        assert_eq!(prod, Fp12::one());
    }
}
