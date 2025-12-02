use crate::fp2::Fp2;
use std::ops::{Add, Mul, Neg, Sub};

/// Fp6 represents the cubic extension Fp6 = Fp2[v] / (v³ - (u+9))
/// An element is represented as c0 + c1*v + c2*v²
/// where v³ = u+9 (the non-residue in Fp2)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fp6 {
    pub c0: Fp2,
    pub c1: Fp2,
    pub c2: Fp2,
}

impl Fp6 {
    pub fn new(c0: Fp2, c1: Fp2, c2: Fp2) -> Self {
        Fp6 { c0, c1, c2 }
    }

    pub fn zero() -> Self {
        Fp6 {
            c0: Fp2::zero(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    pub fn one() -> Self {
        Fp6 {
            c0: Fp2::one(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    /// Non-residue: u+9 in Fp2
    fn non_residue() -> Fp2 {
        use crate::fp::Fp;
        Fp2::new(Fp::new(9u32.into()), Fp::new(1u32.into()))
    }

    /// Multiply by non-residue
    pub fn mul_by_non_residue(a: &Fp2) -> Fp2 {
        a * &Self::non_residue()
    }

    pub fn inv(&self) -> Self {
        // Using the formula from "Implementing Cryptographic Pairings"
        let _nr = Self::non_residue();

        let c0 = &(&self.c0 * &self.c0) - &Self::mul_by_non_residue(&(&self.c1 * &self.c2));
        let c1 = &Self::mul_by_non_residue(&(&self.c2 * &self.c2)) - &(&self.c0 * &self.c1);
        let c2 = &(&self.c1 * &self.c1) - &(&self.c0 * &self.c2);

        let t = &(&self.c2 * &Self::mul_by_non_residue(&c1))
            + &(&self.c1 * &Self::mul_by_non_residue(&c2))
            + (&self.c0 * &c0);
        let t_inv = t.inv();

        Fp6 {
            c0: &c0 * &t_inv,
            c1: &c1 * &t_inv,
            c2: &c2 * &t_inv,
        }
    }
}

impl Add for Fp6 {
    type Output = Fp6;
    fn add(self, rhs: Fp6) -> Fp6 {
        Fp6 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
            c2: self.c2 + rhs.c2,
        }
    }
}

impl<'a, 'b> Add<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;
    fn add(self, rhs: &'b Fp6) -> Fp6 {
        Fp6 {
            c0: &self.c0 + &rhs.c0,
            c1: &self.c1 + &rhs.c1,
            c2: &self.c2 + &rhs.c2,
        }
    }
}

impl Sub for Fp6 {
    type Output = Fp6;
    fn sub(self, rhs: Fp6) -> Fp6 {
        Fp6 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
            c2: self.c2 - rhs.c2,
        }
    }
}

impl<'a, 'b> Sub<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;
    fn sub(self, rhs: &'b Fp6) -> Fp6 {
        Fp6 {
            c0: &self.c0 - &rhs.c0,
            c1: &self.c1 - &rhs.c1,
            c2: &self.c2 - &rhs.c2,
        }
    }
}

impl Mul for Fp6 {
    type Output = Fp6;
    fn mul(self, rhs: Fp6) -> Fp6 {
        // Karatsuba multiplication
        let a_a = &self.c0 * &rhs.c0;
        let b_b = &self.c1 * &rhs.c1;
        let c_c = &self.c2 * &rhs.c2;

        let t0 = &(&self.c1 + &self.c2) * &(&rhs.c1 + &rhs.c2);
        let t1 = &(&self.c0 + &self.c1) * &(&rhs.c0 + &rhs.c1);
        let t2 = &(&self.c0 + &self.c2) * &(&rhs.c0 + &rhs.c2);

        let c0 = &a_a + &Fp6::mul_by_non_residue(&(&t0 - &b_b - c_c.clone()));
        let c1 = &t1 - &a_a - b_b.clone() + Fp6::mul_by_non_residue(&c_c);
        let c2 = &t2 - &a_a - c_c + b_b;

        Fp6 { c0, c1, c2 }
    }
}

impl<'a, 'b> Mul<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;
    fn mul(self, rhs: &'b Fp6) -> Fp6 {
        // Karatsuba multiplication
        let a_a = &self.c0 * &rhs.c0;
        let b_b = &self.c1 * &rhs.c1;
        let c_c = &self.c2 * &rhs.c2;

        let t0 = &(&self.c1 + &self.c2) * &(&rhs.c1 + &rhs.c2);
        let t1 = &(&self.c0 + &self.c1) * &(&rhs.c0 + &rhs.c1);
        let t2 = &(&self.c0 + &self.c2) * &(&rhs.c0 + &rhs.c2);

        let c0 = &a_a + &Fp6::mul_by_non_residue(&(&t0 - &b_b - c_c.clone()));
        let c1 = &t1 - &a_a - b_b.clone() + Fp6::mul_by_non_residue(&c_c);
        let c2 = &t2 - &a_a - c_c + b_b;

        Fp6 { c0, c1, c2 }
    }
}

impl Neg for Fp6 {
    type Output = Fp6;
    fn neg(self) -> Fp6 {
        Fp6 {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp::Fp;

    #[test]
    fn test_basic_ops() {
        let a = Fp6::new(
            Fp2::new(Fp::new(1u32.into()), Fp::new(2u32.into())),
            Fp2::new(Fp::new(3u32.into()), Fp::new(4u32.into())),
            Fp2::new(Fp::new(5u32.into()), Fp::new(6u32.into())),
        );
        let b = Fp6::new(
            Fp2::new(Fp::new(7u32.into()), Fp::new(8u32.into())),
            Fp2::new(Fp::new(9u32.into()), Fp::new(10u32.into())),
            Fp2::new(Fp::new(11u32.into()), Fp::new(12u32.into())),
        );

        let sum = &a + &b;
        assert_eq!(sum.c0, Fp2::new(Fp::new(8u32.into()), Fp::new(10u32.into())));
    }

    #[test]
    fn test_inverse() {
        let a = Fp6::new(
            Fp2::new(Fp::new(1u32.into()), Fp::new(2u32.into())),
            Fp2::new(Fp::new(3u32.into()), Fp::new(4u32.into())),
            Fp2::new(Fp::new(5u32.into()), Fp::new(6u32.into())),
        );
        let inv = a.inv();
        let prod = &a * &inv;
        assert_eq!(prod.c0, Fp2::one());
        assert_eq!(prod.c1, Fp2::zero());
        assert_eq!(prod.c2, Fp2::zero());
    }
}
