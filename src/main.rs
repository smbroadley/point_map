use std::{mem::size_of_val, ops::Range};

use rand::Rng;

#[derive(Debug, Copy, Clone, PartialEq)]
struct Point {
    x: f32,
    y: f32,
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

#[derive(Debug, Copy, Clone, PartialEq)]
struct Circle {
    center: Point,
    radius: f32,
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

    fn new(center: Point, radius: f32) -> Self {
        Self { center, radius }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct Bounds {
    left: f32,
    right: f32,
    top: f32,
    bottom: f32,
}

impl Bounds {
    pub fn new(left: f32, top: f32, right: f32, bottom: f32) -> Bounds {
        Self {
            left,
            top,
            right,
            bottom,
        }
    }

    fn width(&self) -> f32 {
        self.right - self.left
    }

    fn height(&self) -> f32 {
        self.bottom - self.top
    }
}

#[derive(Debug)]
struct PointMapSpan {
    offset: u32,
    size: u32,
}

impl PointMapSpan {
    #[allow(dead_code)]
    pub fn new(offset: u32, count: u32) -> Self {
        Self {
            offset,
            size: count,
        }
    }

    pub fn empty_at(offset: u32) -> Self {
        Self { offset, size: 0 }
    }

    pub fn range(&self) -> Range<usize> {
        self.offset as usize..(self.offset + self.size) as usize
    }
}

struct PointMap {
    size: u32,
    points: Vec<Point>,
    bounds: Bounds, // bounds of the points
    factor: f32,
    index: Vec<PointMapSpan>,
}

impl PointMap {
    pub fn new(mut points: Vec<Point>, size: u32) -> Self {
        let bounds = points_bounds(&points);

        // NOTE: we must ensure that the mapped space is orthonormal(?)
        //       to the input space, so that we can perform simple
        //       circle-circle intersections in our search algorithm.
        //
        let extent = bounds.width().max(bounds.height());
        let factor = size as f32 / extent;

        // to make storage as effecient as possible, we re-order the
        // points in the vector to match the order they would appear
        // in the integer-indexed internal map.
        //
        // this allows us to store slices into the vector; we use our
        // own 'offset + size' slice format for further optimization
        // instead of the Rust slice which uses more storage.
        //
        let cell_index_of = move |p: &Point| {
            let max_idx = size - 1;
            let cx = max_idx.min(((p.x - bounds.left) * factor) as u32);
            let cy = max_idx.min(((p.y - bounds.top) * factor) as u32);

            cx + cy * size
        };

        points.sort_unstable_by_key(cell_index_of);

        // perform indexing on the sorted list of points
        // NOTE: we could just walk the list and build this,
        //       but for simplicity we use a map to gather results.
        //
        let index_size = (size * size) as usize;
        let mut index = Vec::<PointMapSpan>::with_capacity(index_size);

        index.resize_with(index_size, || PointMapSpan::empty_at(0));

        for (i, p) in points.iter().enumerate() {
            let cell = cell_index_of(p);
            let span = &mut index[cell as usize];

            // if the 'size' is 0, we know we are just now initializing
            // this span, and we now set the 'offset' value now, and
            // never again.
            //
            if span.size == 0 {
                span.offset = i as u32;
            }

            span.size += 1;
        }

        Self {
            size,
            points,
            bounds,
            factor,
            index,
        }
    }

    fn point_to_cell(&self, p: &Point) -> (u32, u32) {
        let bounds = self.bounds;
        let factor = self.factor;

        let mut cell_x = ((p.x - bounds.left) * factor) as u32;
        let mut cell_y = ((p.y - bounds.top) * factor) as u32;

        // make sure we don't go out of bounds!
        //
        let max_idx = self.size - 1;

        cell_x = cell_x.clamp(0, max_idx);
        cell_y = cell_y.clamp(0, max_idx);

        (cell_x, cell_y)
    }

    fn cell_to_point(&self, x: u32, y: u32) -> Point {
        let cell_width = 1.0 / self.factor;
        let cell_halfwidth = cell_width / 2.0;

        let px = self.bounds.left + cell_halfwidth + cell_width * x as f32;
        let py = self.bounds.top + cell_halfwidth + cell_width * y as f32;

        Point::new(px, py)
    }

    fn cell_bounds(&self, c: &Circle) -> ((u32, u32), (u32, u32)) {
        let tl = Point::new(c.center.x - c.radius, c.center.y - c.radius);
        let br = Point::new(c.center.x + c.radius, c.center.y + c.radius);

        (self.point_to_cell(&tl), self.point_to_cell(&br))
    }

    /// # Examples
    ///
    /// ```
    /// asserteq!(true, false);
    /// ```
    fn get_points(&self, x: u32, y: u32) -> &[Point] {
        let span = &self.index[(y * self.size + x) as usize];

        // println!("span: {:?}", span);

        &self.points[span.range()]
    }

    pub fn nearest(&self, count: usize, c: &Circle) -> Vec<(f32, Point)> {
        // calculate the top-left, and bottom-right cell
        // indecies as our initial set of cells to consider
        //
        let (tl, br) = self.cell_bounds(c);
        let cap_dist_sq = c.radius * c.radius;

        let mut results = Vec::<_>::with_capacity(count);

        for y in tl.1..=br.1 {
            for x in tl.0..=br.0 {
                // check to see if a circle around the sample point
                // intersects the circle surrounding the cell
                //
                let p = self.cell_to_point(x, y);
                let cr = 1.0 / self.factor;

                if Circle::new(p, cr).intersects(c) {
                    for p in self.get_points(x, y) {
                        let dist_sq = p.distance_sq_to(c.center);
                        if dist_sq <= cap_dist_sq {
                            results.push((dist_sq, *p));
                        }
                    }
                }
            }
        }

        // of course, there is no sorting by f32 or f64 in Rust,
        // even though there is an absolute sorting defined in
        // IEE754... fun.
        //
        // results.sort_by_key(|p| p.distance_sq_to(c.center));
        //
        // we do it with a partial_cmp call...
        //
        results.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        if results.len() > count {
            results.resize(count, (0.0, Point::new(0.0, 0.0)));
        }

        // NOTE: we just compute disatnce-squared to improce
        //       performance, so now we have to .sqrt() them.
        //
        for (dist, _) in results.iter_mut() {
            *dist = dist.sqrt();
        }

        results
    }
}

impl std::fmt::Debug for PointMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "PointMap:")?;
        writeln!(f, "    size: {} x {}", self.size, self.size)?;
        writeln!(f, "    index size: {}", self.index.len())?;
        writeln!(f, "    index bytes: {}", size_of_val(&self.index[..]))?;
        writeln!(f, "    points: {}", self.points.len())?;
        writeln!(f, "    bounds: {:?}", self.bounds)?;

        Ok(())
    }
}

fn points_bounds(points: &Vec<Point>) -> Bounds {
    let mut bounds = Bounds::new(f32::MAX, f32::MAX, f32::MIN, f32::MIN);

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

fn test() {
    let bounds = Bounds::new(0.0, 0.0, 100.0, 200.0);
    let points = gen_random_points(1_000, bounds);

    let map = PointMap::new(points, 100);

    println!("{:?}", map);

    // NOTE: we could provide a map.nearest(&mut out[0..4])->u32
    //       API if we wanted to remove the allocation for the
    //       returned vector.
    //
    let c = Circle::new(Point::new(50.0, 100.0), 4.0);
    let near = map.nearest(4, &c);

    println!(
        "Found {} points within {} units of {:?}",
        near.len(),
        c.radius,
        c.center
    );

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

fn main() {
    test();
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

    #[test]
    fn test_point_map_index_consistency() {
        let points = vec![Point::new(10.0, 20.0), Point::new(30.0, 40.0)];
        let map = PointMap::new(points, 1);

        // validate the points ended up in the same span in the index
        //
        let span = &map.index[0];

        assert_eq!(span.offset, 0);
        assert_eq!(span.size, 2);
    }

    #[test]
    fn test_accuracy() {
        let points = {
            let n = |x: u32, y: u32| Point::new(x as f32, y as f32);
            vec![
                // set extents
                n(0, 0),
                n(4, 4),
                // inner points that sit in centre of cells
                n(1, 1),
                n(2, 1),
                n(3, 1),
                n(1, 2),
                n(2, 2),
                n(3, 2),
                n(1, 3),
                n(2, 3),
                n(3, 3),
            ]
        };

        let map = PointMap::new(points, 2);

        assert_eq!(
            map.bounds,
            Bounds::new(0.0, 0.0, 4.0, 4.0),
            "Bounds should be (0,0)->(5,5)"
        );

        assert_eq!(map.factor, 0.5, "Factor should be 0.5");

        // we have this map of cells:
        // ╭───────────────────────╮
        // │       │       │       │
        // │   A   │   B   │   C   │
        // │       │       │       │
        // │───────┼───────┼───────│
        // │       │       │       │
        // │   D   │   E   │   F   │
        // │       │       │       │
        // │───────┼───────┼───────│
        // │       │       │       │
        // │   G   │   H   │   I   │
        // │       │       │       │
        // ╰───────────────────────╯
        //
        // where the center (E) is physically at point (2,2)
        // and cells are 2 units wide.
        //
        assert_eq!(1, 1);
    }
}
