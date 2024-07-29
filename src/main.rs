use rand::Rng;
use std::time::Instant;

mod point;
use crate::point::Point;

mod circle;
use crate::circle::Circle;

mod bounds;
use crate::bounds::BoundsT;

mod point_map;
use point_map::PointMap;

type Bounds = BoundsT<f32>;

fn main() {
    let bounds = Bounds::new(0.0, 0.0, 100.0, 200.0);
    let points = gen_random_points(1_000_000, bounds);

    let map = PointMap::new(points, 100);

    println!("{:?}", map);

    // NOTE: we could provide a map.nearest(&mut out[0..4])->u32
    //       API if we wanted to remove the allocation for the
    //       returned vector.
    //
    let c = Circle::new(Point::new(50.0, 100.0), 40.0);

    let timer = Instant::now();
    let near = map.nearest(4, &c);

    print!(
        "Found nearest {} points within {} units of {:?} in ",
        near.len(),
        c.radius,
        c.center
    );

    let elapsed = timer.elapsed();

    if elapsed.as_millis() > 0 {
        println!("{} ms", elapsed.as_millis());
    } else {
        println!("{} μs", elapsed.as_micros());
    }

    if near.len() > 0 {
        println!("");
        println!("  distance   point");
        println!("╭────────────────────────────────────────────────╮");
        for p in &near {
            println!("│ {:#8.2}   {:>8.2?}  │", p.0, p.1);
        }
        println!("╰────────────────────────────────────────────────╯");
    }
}

fn gen_random_points(count: u32, bounds: Bounds) -> Vec<Point> {
    let mut rng = rand::thread_rng();
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
