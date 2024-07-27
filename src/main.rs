use rand::{thread_rng, Rng};

struct Point {
    x: f32,
    y: f32,
}

impl Point {
    fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }
}

struct Bounds<T> {
    left: T,
    right: T,
    top: T,
    bottom: T,
}

impl<T> Bounds<T> {
    fn new(left: T, top: T, right: T, bottom: T) -> Bounds<T> {
        Self {
            left, top, right, bottom
        }
    }
}

fn gen_random_points(count: u32, bounds: Bounds<f32>) -> Vec<Point> {
    let mut rng = thread_rng();
    let mut res = Vec::new();

    let bounds_w = bounds.right - bounds.left;
    let bounds_h = bounds.bottom - bounds.top;

    for _ in 0..count {
        let x = bounds.left + (rng.gen::<f32>() * bounds_w);
        let y = bounds.top + (rng.gen::<f32>() * bounds_h);
        res.push(Point::new(x, y))
    }
    res
}

fn main() {
    let points = gen_random_points(1000, Bounds::new(0.0, 0.0, 100.0, 200.0) )
    println!("Hello, world!");
}
