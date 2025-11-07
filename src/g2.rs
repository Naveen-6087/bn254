use crate::fp::Fp;
use crate::fp2::Fp2;
use num_bigint::BigUint;
use num_traits::Zero;

/// G2 is the twisted curve over Fp2
/// Twist curve equation: y² = x³ + 3/(u + 9)
/// We use the isomorphic curve: y² = x³ + 3*(u+9)
/// in Jacobian coordinates (X:Y:Z) where x = X/Z², y = Y/Z³
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct G2 {
    pub x: Fp2,
    pub y: Fp2,
    pub z: Fp2,
}

impl G2 {
    /// Returns the point at infinity
    pub fn infinity() -> Self {
        Self {
            x: Fp2::zero(),
            y: Fp2::one(),
            z: Fp2::zero(),
        }
    }

    /// Check if this point is the point at infinity
    pub fn is_infinity(&self) -> bool {
        self.z.c0.n.is_zero() && self.z.c1.n.is_zero()
    }

    /// Convert from Jacobian to affine coordinates
    pub fn to_affine(&self) -> (Fp2, Fp2) {
        if self.is_infinity() {
            return (Fp2::zero(), Fp2::zero());
        }
        let z_inv = self.z.inv();
        let z2_inv = &z_inv * &z_inv;
        let z3_inv = &z2_inv * &z_inv;
        let x_aff = &self.x * &z2_inv;
        let y_aff = &self.y * &z3_inv;
        (x_aff, y_aff)
    }

    /// Get the curve coefficient b' = 3/(u+9)
    /// For the twist, we use b' = 3/(9+u)
    fn get_b() -> Fp2 {
        // b' = 3/(9+u)
        // Compute inverse of (9+u)
        let nine_plus_u = Fp2::new(Fp::new(9u32.into()), Fp::new(1u32.into()));
        let inv = nine_plus_u.inv();
        let three = Fp2::new(Fp::new(3u32.into()), Fp::zero());
        &three * &inv
    }

    /// Check if the point is on the curve
    pub fn is_on_curve(&self) -> bool {
        if self.is_infinity() {
            return true;
        }
        let (x, y) = self.to_affine();
        let b = Self::get_b();
        // y² = x³ + b
        let y2 = &y * &y;
        let x3 = &(&x * &x) * &x;
        let rhs = &x3 + &b;
        y2 == rhs
    }

    /// Point doubling in Jacobian coordinates
    /// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
    pub fn double(&self) -> Self {
        if self.is_infinity() {
            return self.clone();
        }

        // A = X1²
        let a = &self.x * &self.x;
        // B = Y1²
        let b = &self.y * &self.y;
        // C = B²
        let c = &b * &b;
        // D = 2*((X1+B)² - A - C)
        let x1_plus_b = &self.x + &b;
        let d_temp = &(&x1_plus_b * &x1_plus_b) - &a;
        let d_temp = &d_temp - &c;
        let d = &d_temp + &d_temp;
        // E = 3*A
        let e = &(&a + &a) + &a;
        // F = E²
        let f = &e * &e;
        // X3 = F - 2*D
        let two_d = &d + &d;
        let x3 = &f - &two_d;
        // Y3 = E*(D - X3) - 8*C
        let eight_c = &(&(&c + &c) + &(&c + &c)) + &(&(&c + &c) + &(&c + &c));
        let y3 = &(&e * &(&d - &x3)) - &eight_c;
        // Z3 = 2*Y1*Z1
        let z3 = &(&self.y * &self.z) + &(&self.y * &self.z);

        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Point addition in Jacobian coordinates
    /// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
    pub fn add(&self, other: &Self) -> Self {
        if self.is_infinity() {
            return other.clone();
        }
        if other.is_infinity() {
            return self.clone();
        }

        // Z1Z1 = Z1²
        let z1z1 = &self.z * &self.z;
        // Z2Z2 = Z2²
        let z2z2 = &other.z * &other.z;
        // U1 = X1*Z2Z2
        let u1 = &self.x * &z2z2;
        // U2 = X2*Z1Z1
        let u2 = &other.x * &z1z1;
        // S1 = Y1*Z2*Z2Z2
        let s1 = &(&self.y * &other.z) * &z2z2;
        // S2 = Y2*Z1*Z1Z1
        let s2 = &(&other.y * &self.z) * &z1z1;

        if u1 == u2 {
            if s1 == s2 {
                return self.double();
            } else {
                return Self::infinity();
            }
        }

        // H = U2 - U1
        let h = &u2 - &u1;
        // I = (2*H)²
        let two_h = &h + &h;
        let i = &two_h * &two_h;
        // J = H*I
        let j = &h * &i;
        // r = 2*(S2 - S1)
        let r = &(&s2 - &s1) + &(&s2 - &s1);
        // V = U1*I
        let v = &u1 * &i;
        // X3 = r² - J - 2*V
        let two_v = &v + &v;
        let x3 = &(&(&r * &r) - &j) - &two_v;
        // Y3 = r*(V - X3) - 2*S1*J
        let two_s1_j = &(&s1 * &j) + &(&s1 * &j);
        let y3 = &(&r * &(&v - &x3)) - &two_s1_j;
        // Z3 = ((Z1+Z2)² - Z1Z1 - Z2Z2)*H
        let z1_plus_z2 = &self.z + &other.z;
        let z3_temp = &(&z1_plus_z2 * &z1_plus_z2) - &z1z1;
        let z3_temp = &z3_temp - &z2z2;
        let z3 = &z3_temp * &h;

        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Scalar multiplication using double-and-add
    pub fn mul_scalar(&self, scalar: &BigUint) -> Self {
        let mut res = Self::infinity();
        let mut base = self.clone();
        let mut s = scalar.clone();

        while !s.is_zero() {
            if s.bit(0) {
                res = res.add(&base);
            }
            base = base.double();
            s >>= 1;
        }

        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigUint;
    use num_traits::{One, Zero};

    #[test]
    fn test_infinity() {
        let inf = G2::infinity();
        assert!(inf.is_infinity());
    }

    #[test]
    fn test_double_vs_add() {
        // Create a point (we'll use a simple construction)
        let x = Fp2::new(Fp::new(1u32.into()), Fp::new(2u32.into()));
        let y = Fp2::new(Fp::new(3u32.into()), Fp::new(4u32.into()));
        let p = G2 {
            x,
            y,
            z: Fp2::one(),
        };

        let doubled = p.double();
        let added = p.add(&p);
        assert_eq!(doubled.to_affine(), added.to_affine());
    }

    #[test]
    fn test_addition_commutative() {
        let p1 = G2 {
            x: Fp2::new(Fp::new(1u32.into()), Fp::new(2u32.into())),
            y: Fp2::new(Fp::new(3u32.into()), Fp::new(4u32.into())),
            z: Fp2::one(),
        };
        let p2 = G2 {
            x: Fp2::new(Fp::new(5u32.into()), Fp::new(6u32.into())),
            y: Fp2::new(Fp::new(7u32.into()), Fp::new(8u32.into())),
            z: Fp2::one(),
        };

        let sum1 = p1.add(&p2);
        let sum2 = p2.add(&p1);
        assert_eq!(sum1.to_affine(), sum2.to_affine());
    }

    #[test]
    fn test_scalar_mul() {
        let p = G2 {
            x: Fp2::new(Fp::new(1u32.into()), Fp::new(2u32.into())),
            y: Fp2::new(Fp::new(3u32.into()), Fp::new(4u32.into())),
            z: Fp2::one(),
        };

        let res0 = p.mul_scalar(&BigUint::zero());
        assert!(res0.is_infinity());

        let res1 = p.mul_scalar(&BigUint::one());
        assert_eq!(res1.to_affine(), p.to_affine());

        let res2 = p.mul_scalar(&2u32.to_biguint().unwrap());
        assert_eq!(res2.to_affine(), p.add(&p).to_affine());

        let res3 = p.mul_scalar(&3u32.to_biguint().unwrap());
        assert_eq!(res3.to_affine(), p.add(&p).add(&p).to_affine());
    }
}
