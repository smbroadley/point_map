use rand::Rng;
use std::time::Instant;
use std::{mem::size_of_val, ops::Range};

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
struct BoundsT<T> {
    left: T,
    right: T,
    top: T,
    bottom: T,
}

use std::ops::{Mul, RangeInclusive, Sub};

impl<T> BoundsT<T>
where
    T: std::ops::Sub + Copy,
    <T as std::ops::Sub>::Output: std::ops::Mul,
{
    pub fn new(left: T, top: T, right: T, bottom: T) -> Self {
        Self {
            left,
            top,
            right,
            bottom,
        }
    }

    fn width(&self) -> <T as Sub>::Output {
        self.right - self.left
    }

    fn height(&self) -> <T as Sub>::Output {
        self.bottom - self.top
    }

    fn area(&self) -> <<T as Sub>::Output as Mul>::Output {
        (self.right - self.left) * (self.bottom - self.top)
    }

    fn x_range(&self) -> RangeInclusive<T> {
        self.left..=self.right
    }

    fn y_range(&self) -> RangeInclusive<T> {
        self.top..=self.bottom
    }
}

type Bounds = BoundsT<f32>;
type CellBounds = BoundsT<u32>;

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
    /// Creates a new [`PointMap`].
    pub fn new(mut points: Vec<Point>, size: u32) -> Self {
        let bounds = points_bounds(&points);

        // we must ensure that cell-space is a simple scale-
        // mappint to and from the point-space, so calulate
        // the extents, and take the largest axis as the
        // one to base the scale factor on
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

    fn cell_bounds(&self, c: &Circle) -> CellBounds {
        let tl = Point::new(c.center.x - c.radius, c.center.y - c.radius);
        let br = Point::new(c.center.x + c.radius, c.center.y + c.radius);

        let ctl = self.point_to_cell(&tl);
        let cbr = self.point_to_cell(&br);

        CellBounds::new(ctl.0, ctl.1, cbr.0, cbr.1)
    }

    /// # Examples
    ///
    /// ```
    /// asserteq!(true, false);
    /// ```
    fn get_points(&self, index: u32) -> &[Point] {
        let span = &self.index[index as usize];

        &self.points[span.range()]
    }

    pub fn nearest(&self, count: usize, c: &Circle) -> Vec<(f32, Point)> {
        // calculate the top-left, and bottom-right cell
        // indecies as our initial set of cells to consider
        //
        let cell_bounds = self.cell_bounds(c);

        // we do not need to return the *real* distance to any
        // particular particle, so we can stay in "squared-
        // distance-space" for better performance when testing
        // for circle intersections
        //
        let search_dist_sq = c.radius * c.radius;

        // create a vec to contain the cells we need to check
        // for particles; it's maximum capacity is calculatable
        // from the cells bounds' area
        //
        let cells_area = cell_bounds.area();
        let mut cells = Vec::<_>::with_capacity(cells_area as usize);

        // iterate over the cells' coordinates in the bounds,
        // testing to see if a circle placed at the center of
        // the cell would intersect our search circle argument.
        //
        for y in cell_bounds.y_range() {
            for x in cell_bounds.x_range() {
                // calculate cell-point (cp) and cell-radius (cr)
                //
                let cp = self.cell_to_point(x, y);
                let cr = (1.0 / self.factor) / 2.0 * std::f32::consts::SQRT_2;

                // does it intersect with our search bounds (c)
                //
                if Circle::new(cp, cr).intersects(c) {
                    let dist = cp.distance_sq_to(c.center); // TODO: get this from previous call
                    let index = x + y * self.size;

                    // add this to the list of potential cells
                    //
                    cells.push((dist, index));
                }
            }
        }

        cells.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // iterate through the cells indexes and add points
        // to the result set...
        //
        let mut results = Vec::<_>::with_capacity(count);

        // create an ever increasing radius of a circle that
        // would entirely cover a cell; that is a circle of
        // radius 'cell-center to cell-corner'
        //
        let cell_inner_radius = (1.0 / self.factor) / 2.0;
        let cell_outer_radius = cell_inner_radius * std::f32::consts::SQRT_2;

        // we expand the radius check each iteration whilds
        // adding points to the potential result set, starting
        // with the smallest radius check possible
        //
        let mut radius_check = cell_outer_radius;

        for &(dist, idx) in &cells {
            // if we find a distance (in our distance-sorted list)
            // that is bigger than our current checkin radius, we
            // can early out if we have ennough points in the result
            // set; otherwise, increase the search radius, and gather
            // more points
            //
            if dist > radius_check {
                if results.len() >= count {
                    break; // bingo!
                }

                radius_check += cell_outer_radius;
            } else {
                for p in self.get_points(idx) {
                    let dist_sq = p.distance_sq_to(c.center);
                    if dist_sq <= search_dist_sq {
                        results.push((dist_sq, *p));
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
        let num_pts = self.points.len();
        let avg_occ = self.points.len() as f32 / (self.size * self.size) as f32;

        writeln!(f, "PointMap:")?;
        writeln!(f, "    points: {}", num_pts)?;
        writeln!(f, "     cells: {} x {}", self.size, self.size)?;
        writeln!(f, "       mem: {} KB", size_of_val(&self.index[..]) / 1024)?;
        writeln!(f, " occupancy: {} points per cell (est)", avg_occ)?;

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
