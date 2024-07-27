use rand::{thread_rng, Rng};

#[derive(Debug, Copy, Clone, PartialEq)]
struct Point {
    x: f32,
    y: f32,
}

impl Point {
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct Bounds<T> {
    left: T,
    right: T,
    top: T,
    bottom: T,
}

// NOTE: templated in case I want to use this for PointMap internals,
// but it is looking like I will not, so it may get removed.
//
impl<T> Bounds<T> {
    pub fn new(left: T, top: T, right: T, bottom: T) -> Bounds<T> {
        Self {
            left,
            top,
            right,
            bottom,
        }
    }
}

struct PointMap {
    width: u32,
    height: u32,
    points: Vec<Point>,
    bounds: Bounds<f32>, // bounds of the points
    mapping: fn(f32, f32) -> (u32, u32),
}

impl PointMap {
    pub fn new(points: Vec<Point>) -> Self {
        let bounds = points_bounds(&points);

        // to make storage as effecient as possible, we re-order the
        // points in the vector to match the order they would appear
        // in the integer-indexed internal map.
        //
        // this allows us to store slices into the vector; we use our
        // own 'offset + size' slice format for further optimization
        // instead of the Rust slice which uses more storage.
        //

        Self {
            width: 0,
            height: 0,
            points,
            bounds,
            mapping: |x, y| (0, 0),
        }
    }

    pub fn nearest(&self, count: i32, sample: Point) -> Vec<Point> {
        todo!()
    }
}

fn points_bounds(points: &Vec<Point>) -> Bounds<f32> {
    let mut bounds = Bounds::<f32>::new(f32::MAX, f32::MAX, f32::MIN, f32::MIN);

    for p in points {
        if p.x < bounds.left {
            bounds.left = p.x
        }
        if p.y < bounds.top {
            bounds.top = p.y
        }
        if p.x > bounds.right {
            bounds.right = p.x
        }
        if p.y > bounds.bottom {
            bounds.bottom = p.y
        }
    }

    bounds
}

fn gen_random_points(count: u32, bounds: Bounds<f32>) -> Vec<Point> {
    let mut rng = thread_rng();
    let mut res = Vec::new();

    let w = bounds.right - bounds.left;
    let h = bounds.bottom - bounds.top;

    for _ in 0..count {
        let x = bounds.left + (rng.gen::<f32>() * w);
        let y = bounds.top + (rng.gen::<f32>() * h);
        res.push(Point::new(x, y))
    }

    res
}

fn test() {
    let bounds = Bounds::new(0.0, 0.0, 100.0, 200.0);
    let points = gen_random_points(1000, bounds);

    let map = PointMap::new(points);

    // NOTE: we could provide a map.nearest(&mut out[0..4])->u32
    //       API if we wanted to remove the allocation for the
    //       returned vector.
    //
    let near = map.nearest(4, Point::new(50.0, 100.0));
}

fn main() {
    println!("{}", thread_rng().gen::<f32>());
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_points_bounds_single() {
        let points = vec![Point::new(10.0, 20.0)];
        assert_eq!(points_bounds(&points), Bounds::new(10.0, 20.0, 10.0, 20.0));
    }

    #[test]
    fn test_points_bounds_two() {
        let points = vec![Point::new(10.0, 20.0), Point::new(30.0, 40.0)];
        assert_eq!(points_bounds(&points), Bounds::new(10.0, 20.0, 30.0, 40.0));
    }
}
