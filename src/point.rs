#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Point {
    pub x: f32,
    pub y: f32,
}

#[allow(dead_code)]
impl Point {
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    pub fn distance_sq_to(&self, p: Point) -> f32 {
        let dx = (self.x - p.x).abs();
        let dy = (self.y - p.y).abs();

        (dx * dx) + (dy * dy)
    }

    pub fn distance_to(&self, p: Point) -> f32 {
        self.distance_sq_to(p).sqrt()
    }
}
