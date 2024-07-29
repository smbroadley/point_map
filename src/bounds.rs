#[derive(Debug, Copy, Clone, PartialEq)]
pub struct BoundsT<T> {
    pub left: T,
    pub right: T,
    pub top: T,
    pub bottom: T,
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

    pub fn width(&self) -> <T as Sub>::Output {
        self.right - self.left
    }

    pub fn height(&self) -> <T as Sub>::Output {
        self.bottom - self.top
    }

    pub fn area(&self) -> <<T as Sub>::Output as Mul>::Output {
        (self.right - self.left) * (self.bottom - self.top)
    }

    pub fn x_range(&self) -> RangeInclusive<T> {
        self.left..=self.right
    }

    pub fn y_range(&self) -> RangeInclusive<T> {
        self.top..=self.bottom
    }
}
