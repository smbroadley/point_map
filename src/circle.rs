use crate::point::Point;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Circle {
    pub center: Point,
    pub radius: f32,
}

impl Circle {
    pub fn intersects(&self, other: &Circle) -> bool {
        // NOTE: we are using distance-squared comparisons,
        //       which benchmarks have show are much quicker.
        //
        let dist = self.center.distance_sq_to(other.center);
        let rsum = self.radius + other.radius;

        dist < (rsum * rsum)
    }

    pub fn new(center: Point, radius: f32) -> Self {
        Self { center, radius }
    }
}
