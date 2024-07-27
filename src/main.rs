use rand::Rng;

#[derive(Debug, Copy, Clone, PartialEq)]
struct Point {
    x: f32,
    y: f32,
}

impl Point {
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    pub fn distance_to(&self, p: Point) -> f32 {
        let dx = (self.x - p.x).abs();
        let dy = (self.y - p.y).abs();

        ((dx * dx) + (dy * dy)).sqrt()
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct Circle {
    center: Point,
    radius: f32,
}

impl Circle {
    pub fn intersects(&self, other: &Circle) -> bool {
        let dist = self.center.distance_to(other.center);

        dist < (self.radius + other.radius)
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

    fn empty_at(offset: u32) -> Self {
        Self { offset, size: 0 }
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
        let mut span_groups = std::collections::HashMap::<u32, PointMapSpan>::new();

        for (i, p) in points.iter().enumerate() {
            let cell = cell_index_of(p);

            // add this cell-index as a single-element span to the
            // PointMap index, or update the span count if it already
            // exists (as we know they are all contiguous from the sort)
            //
            span_groups
                .entry(cell)
                .or_insert_with(|| PointMapSpan::empty_at(i as u32))
                .size += 1;
        }

        let mut index: Vec<PointMapSpan> = span_groups.into_values().collect();

        // NOTE I think we should already be sorted, but I need to check
        //      this... sorting a sorted vec shoudl be super quick, in any
        //      case.
        //
        index.sort_unstable_by_key(|k| k.offset);

        Self {
            size,
            points,
            bounds,
            factor,
            index,
        }
    }

    fn point_to_cell(&self, p: &Point) -> (u32, u32) {
        let max_idx = self.size - 1;
        let bounds = self.bounds;
        let factor = self.factor;

        let cell_x = max_idx.min(((p.x - bounds.left) / factor) as u32);
        let cell_y = max_idx.min(((p.y - bounds.top) / factor) as u32);

        (cell_x, cell_y)
    }

    fn cell_to_point(&self, x: u32, y: u32) -> Point {
        let cell_width = 1.0 / self.factor;
        let cell_halfwidth = cell_width / 2.0;

        let px = self.bounds.left + cell_halfwidth + cell_width * x as f32;
        let py = self.bounds.top + cell_halfwidth + cell_width * y as f32;

        Point::new(px, py)
    }

    fn get_cell_bounds(&self, c: &Circle) -> ((u32, u32), (u32, u32)) {
        let tl = self.point_to_cell(&Point::new(c.center.x - c.radius, c.center.y - c.radius));
        let br = self.point_to_cell(&Point::new(c.center.x + c.radius, c.center.y + c.radius));

        (tl, br)
    }

    pub fn nearest(&self, count: i32, c: &Circle) -> Vec<Point> {
        // calculate the top-left, and bottom-right cell
        // indecies as our initial set of cells to consider
        //
        let (tl, br) = self.get_cell_bounds(c);

        println!("tl: {:?}   br: {:?}", tl, br);

        let mut results = Vec::<Point>::with_capacity(count as usize);

        for y in tl.1..=br.1 {
            for x in tl.0..=br.0 {
                // check to see if a circle around the sample point
                // intersects the circle surrounding the cell
                //
                let p = self.cell_to_point(x, y);
                let cr = 1.0 / self.factor;

                println!(
                    "searching: cell({}, {}) point({:?}) radius({})",
                    x, y, p, cr
                );

                if Circle::new(p, cr).intersects(c) {
                    // we include this cell's points span in
                    // the results
                    //
                    let span = &self.index[(y * self.size + x) as usize];

                    // move all the points into the result set
                    //
                    for point_idx in span.offset..span.size + span.offset {
                        results.push(self.points[point_idx as usize]);
                    }
                }
            }
        }

        results
    }
}

impl std::fmt::Debug for PointMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "PointMap:")?;
        writeln!(f, "    size: {}", self.size)?;
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
    let points = gen_random_points(10, bounds);

    let map = PointMap::new(points, 2);

    println!("{:?}", map);

    // NOTE: we could provide a map.nearest(&mut out[0..4])->u32
    //       API if we wanted to remove the allocation for the
    //       returned vector.
    //
    let c = Circle::new(Point::new(50.0, 100.0), 10.0);
    let near = map.nearest(4, &c);

    println!("Found {} points", near.len());

    for p in &near {
        println!("point: {:?}", p);
    }
}

fn main() {
    test();
    // let points = vec![Point::new(10.0, 20.0), Point::new(30.0, 40.0)];
    // let map = PointMap::new(points, 1);

    // for span in &map.index {
    //     println!("o: {} s: {}", span.offset, span.size);
    // }
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
}
