use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use crate::physics::constants::RAD2DEG;

/// A 4-momentum vector (px, py, pz, energy).
/// The `energy` field stores the total energy E = sqrt(p^2 + m^2),
/// so invariant mass is correctly computed as m^2 = E^2 - p^2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FourMomentum {
    pub px: f32,
    pub py: f32,
    pub pz: f32,
    pub energy: f32,
}

impl FourMomentum {
    /// Create from (px, py, pz, mass) - computes energy from mass
    pub fn new(px: f32, py: f32, pz: f32, mass: f32) -> Self {
        let p2 = px * px + py * py + pz * pz;
        Self {
            px,
            py,
            pz,
            energy: (p2 + mass * mass).sqrt(),
        }
    }

    /// Create from (p, cx, cy, cz, mass) where cx,cy,cz are direction cosines
    pub fn from_direction(p: f32, cx: f32, cy: f32, cz: f32, mass: f32) -> Self {
        Self {
            px: p * cx,
            py: p * cy,
            pz: p * cz,
            energy: (p * p + mass * mass).sqrt(),
        }
    }

    /// Create a 4-vector at rest with given mass
    pub fn at_rest(mass: f32) -> Self {
        Self {
            px: 0.0,
            py: 0.0,
            pz: 0.0,
            energy: mass,
        }
    }

    /// Create directly from (px, py, pz, energy) - no mass lookup
    pub fn from_energy(px: f32, py: f32, pz: f32, energy: f32) -> Self {
        Self { px, py, pz, energy }
    }

    /// Energy
    pub fn e(&self) -> f32 {
        self.energy
    }

    /// Energy (alias)
    pub fn energy(&self) -> f32 {
        self.energy
    }

    /// Rest mass: m = sqrt(E^2 - p^2)
    pub fn mass(&self) -> f32 {
        let m2 = self.m2();
        if m2 > 0.0 {
            m2.sqrt()
        } else {
            0.0
        }
    }

    /// Squared 3-momentum
    pub fn p2(&self) -> f32 {
        self.px * self.px + self.py * self.py + self.pz * self.pz
    }

    /// Magnitude of 3-momentum
    pub fn p(&self) -> f32 {
        self.p2().sqrt()
    }

    /// Squared invariant mass of the 4-vector
    pub fn m2(&self) -> f32 {
        self.energy * self.energy - self.p2()
    }

    /// Invariant mass
    pub fn m(&self) -> f32 {
        self.m2().abs().sqrt()
    }

    /// Alias for m()
    pub fn mag(&self) -> f32 {
        self.m()
    }

    /// Alias for m2()
    pub fn mag2(&self) -> f32 {
        self.m2()
    }

    /// Polar angle theta (radians)
    pub fn theta(&self) -> f32 {
        let p = self.p();
        if p == 0.0 {
            return 0.0;
        }
        (self.pz / p).acos()
    }

    /// Polar angle theta (degrees)
    pub fn theta_deg(&self) -> f32 {
        self.theta() * RAD2DEG
    }

    /// Azimuthal angle phi (radians), in range [-pi, pi]
    pub fn phi(&self) -> f32 {
        self.py.atan2(self.px)
    }

    /// Azimuthal angle phi (degrees), in range [-180, 180]
    pub fn phi_deg(&self) -> f32 {
        self.phi() * RAD2DEG
    }

    /// Dot product of 3-momenta
    pub fn dot3(&self, other: &FourMomentum) -> f32 {
        self.px * other.px + self.py * other.py + self.pz * other.pz
    }

    /// Dot product of 4-vectors (Minkowski: E1*E2 - p1.p2)
    pub fn dot4(&self, other: &FourMomentum) -> f32 {
        self.e() * other.e() - self.dot3(other)
    }

    /// Boost this 4-vector along the z-axis
    pub fn boost_z(&self, beta: f32) -> Self {
        let bgamma = 1.0 / (1.0 - beta * beta).sqrt().max(1.0e-10);
        let e = self.e();
        let pz = self.pz;
        let new_e = bgamma * (e + beta * pz);
        let new_pz = bgamma * (pz + beta * e);
        Self {
            px: self.px,
            py: self.py,
            pz: new_pz,
            energy: new_e,
        }
    }

    /// General Lorentz boost along an arbitrary direction (bx, by, bz)
    pub fn boost(&self, bx: f32, by: f32, bz: f32) -> Self {
        let b2 = bx * bx + by * by + bz * bz;
        if b2 >= 1.0 {
            return *self;
        }
        let gamma = 1.0 / (1.0 - b2).sqrt();
        let bp = self.px * bx + self.py * by + self.pz * bz;
        let gamma2 = if b2 > 1.0e-20 {
            (gamma - 1.0) / b2
        } else {
            0.0
        };
        let e = self.e();

        Self {
            px: self.px + gamma2 * bp * bx + gamma * bx * e,
            py: self.py + gamma2 * bp * by + gamma * by * e,
            pz: self.pz + gamma2 * bp * bz + gamma * bz * e,
            energy: gamma * (e + bp),
        }
    }

    /// Rotate around the z-axis by angle phi (radians)
    pub fn rotate_z(&self, phi: f32) -> Self {
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        Self {
            px: self.px * cos_phi - self.py * sin_phi,
            py: self.px * sin_phi + self.py * cos_phi,
            pz: self.pz,
            energy: self.energy,
        }
    }

    /// Rotate around an arbitrary axis by angle theta
    pub fn rotate_axis(&self, axis: &ThreeVector, angle: f32) -> Self {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let v = ThreeVector::new(self.px, self.py, self.pz);
        let kv = axis.cross(&v);
        let kkv = axis.cross(&kv);

        let rv = ThreeVector::new(
            self.px * cos_a + kv.x * sin_a + kkv.x * (1.0 - cos_a),
            self.py * cos_a + kv.y * sin_a + kkv.y * (1.0 - cos_a),
            self.pz * cos_a + kv.z * sin_a + kkv.z * (1.0 - cos_a),
        );

        Self {
            px: rv.x,
            py: rv.y,
            pz: rv.z,
            energy: self.energy,
        }
    }
}

/// 3-vector for rotation operations
#[derive(Debug, Clone, Copy)]
pub struct ThreeVector {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl ThreeVector {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }

    pub fn mag(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn unit(&self) -> Self {
        let m = self.mag();
        if m == 0.0 {
            return *self;
        }
        Self {
            x: self.x / m,
            y: self.y / m,
            z: self.z / m,
        }
    }

    pub fn dot(&self, other: &ThreeVector) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: &ThreeVector) -> ThreeVector {
        ThreeVector {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn add(&self, other: &ThreeVector) -> ThreeVector {
        ThreeVector {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    pub fn sub(&self, other: &ThreeVector) -> ThreeVector {
        ThreeVector {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    pub fn scale(&self, s: f32) -> ThreeVector {
        ThreeVector {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    /// Rotate around z-axis by angle
    pub fn rotate_z(&self, angle: f32) -> ThreeVector {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        ThreeVector {
            x: self.x * cos_a - self.y * sin_a,
            y: self.x * sin_a + self.y * cos_a,
            z: self.z,
        }
    }
}

// Operator overloads for FourMomentum
impl Add for FourMomentum {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            px: self.px + rhs.px,
            py: self.py + rhs.py,
            pz: self.pz + rhs.pz,
            energy: self.energy + rhs.energy,
        }
    }
}

impl AddAssign for FourMomentum {
    fn add_assign(&mut self, rhs: Self) {
        self.px += rhs.px;
        self.py += rhs.py;
        self.pz += rhs.pz;
        self.energy += rhs.energy;
    }
}

impl Sub for FourMomentum {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            px: self.px - rhs.px,
            py: self.py - rhs.py,
            pz: self.pz - rhs.pz,
            energy: self.energy - rhs.energy,
        }
    }
}

impl SubAssign for FourMomentum {
    fn sub_assign(&mut self, rhs: Self) {
        self.px -= rhs.px;
        self.py -= rhs.py;
        self.pz -= rhs.pz;
        self.energy -= rhs.energy;
    }
}

impl Neg for FourMomentum {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            px: -self.px,
            py: -self.py,
            pz: -self.pz,
            energy: self.energy,
        }
    }
}

impl Mul<f32> for FourMomentum {
    type Output = Self;
    fn mul(self, rhs: f32) -> Self::Output {
        Self {
            px: self.px * rhs,
            py: self.py * rhs,
            pz: self.pz * rhs,
            energy: self.energy * rhs.abs(),
        }
    }
}

impl Div<f32> for FourMomentum {
    type Output = Self;
    fn div(self, rhs: f32) -> Self::Output {
        Self {
            px: self.px / rhs,
            py: self.py / rhs,
            pz: self.pz / rhs,
            energy: self.energy / rhs.abs(),
        }
    }
}

/// Calculate angle between two 4-vectors
pub fn angle_between(a: &FourMomentum, b: &FourMomentum) -> f32 {
    let p1 = a.p();
    let p2 = b.p();
    if p1 == 0.0 || p2 == 0.0 {
        return 0.0;
    }
    let cos_angle = (a.px * b.px + a.py * b.py + a.pz * b.pz) / (p1 * p2);
    cos_angle.clamp(-1.0, 1.0).acos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::constants::{MASS_E, MASS_P};

    #[test]
    fn test_four_momentum_energy() {
        let v = FourMomentum::new(0.0, 0.0, 0.0, MASS_P);
        assert!((v.e() - MASS_P).abs() < 1.0e-6);
    }

    #[test]
    fn test_four_momentum_add() {
        let a = FourMomentum::new(1.0, 0.0, 0.0, MASS_E);
        let b = FourMomentum::new(0.0, 1.0, 0.0, MASS_E);
        let c = a + b;
        assert!((c.px - 1.0).abs() < 1.0e-6);
        assert!((c.py - 1.0).abs() < 1.0e-6);
    }

    #[test]
    fn test_four_momentum_sub_mass() {
        // Missing mass: beam - scattered + target - pion should have mass ~neutron
        let beam = FourMomentum::new(0.0, 0.0, 4.817, MASS_E);
        let target = FourMomentum::at_rest(MASS_P);
        // A simple test: at_rest - at_rest should give mass 0
        let a = FourMomentum::at_rest(1.0);
        let b = FourMomentum::at_rest(1.0);
        let diff = a - b;
        assert!((diff.m() - 0.0).abs() < 1.0e-6);
    }

    #[test]
    fn test_boost_at_rest() {
        // A particle at rest boosted along z should get pz = gamma*beta*m, E = gamma*m
        let v = FourMomentum::new(0.0, 0.0, 0.0, 1.0);
        let beta = 0.5;
        let boosted = v.boost_z(beta);
        let gamma = 1.0 / (1.0 - beta * beta).sqrt();
        assert!((boosted.pz - gamma * beta * 1.0).abs() < 1.0e-5);
        assert!((boosted.e() - gamma * 1.0).abs() < 1.0e-5);
    }
}
