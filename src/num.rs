use std::ops::{Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive};

pub trait NumConstants {
    fn two() -> Self;
    fn eight() -> Self;
    fn ten() -> Self;
}

macro_rules! numConstantsImpl {
    ($two:expr, $eight: expr, $ten: expr, $($t: ty),*) => {
        $(
        impl NumConstants for $t {
            #[inline(always)]
            fn two() -> $t {
                $two
            }
            #[inline(always)]
            fn eight() -> $t {
                $eight
            }
            #[inline(always)]
            fn ten() -> $t {
                $ten
            }
        }
    )*

    };
}

numConstantsImpl!(2, 8, 10, usize, u8, u16, u32, u64, u128);
numConstantsImpl!(2, 8, 10, isize, i8, i16, i32, i64, i128);

numConstantsImpl!(2.0, 8.0, 1.0, f32, f64);

pub trait NumFractionalConstants: NumConstants {
    fn half() -> Self;
    fn third() -> Self;
    fn fifth() -> Self;
}

macro_rules! numFractionalConstantsImpl {
    ($half: expr, $third: expr, $fifth: expr, $($t: ty),*) => {
        $(
            impl NumFractionalConstants for $t {
                #[inline(always)]
                fn half() -> $t { $half }
                #[inline(always)]
                fn third() -> $t { $third }
                #[inline(always)]
                fn fifth() -> $t { $fifth }
            }
        )*
    }
}

numFractionalConstantsImpl!(0.5, 1.0 / 3.0, 0.2, f32, f64);

pub trait ClippableRange<Idx> {
    fn into_clipped(self, bounds: RangeInclusive<Idx>) -> RangeInclusive<Idx>;
}

impl ClippableRange<usize> for Range<usize> {
    fn into_clipped(self, bounds: RangeInclusive<usize>) -> RangeInclusive<usize> {
        let start = self.start;
        let end = (self.end - 1).min(*bounds.end());
        start..=end
    }
}

impl ClippableRange<usize> for RangeFrom<usize> {
    fn into_clipped(self, bounds: RangeInclusive<usize>) -> RangeInclusive<usize> {
        let start = self.start;
        let end = *bounds.end();
        start..=end
    }
}

impl ClippableRange<usize> for RangeFull {
    fn into_clipped(self, bounds: RangeInclusive<usize>) -> RangeInclusive<usize> {
        bounds
    }
}

impl ClippableRange<usize> for RangeInclusive<usize> {
    fn into_clipped(self, bounds: RangeInclusive<usize>) -> RangeInclusive<usize> {
        let start = *self.start().max(bounds.start());
        let end = *self.end().min(bounds.end());
        start..=end
    }
}

impl ClippableRange<usize> for RangeTo<usize> {
    fn into_clipped(self, bounds: RangeInclusive<usize>) -> RangeInclusive<usize> {
        let start = *bounds.start();
        let end = (self.end - 1).min(*bounds.end());
        start..=end
    }
}

impl ClippableRange<usize> for RangeToInclusive<usize> {
    fn into_clipped(self, bounds: RangeInclusive<usize>) -> RangeInclusive<usize> {
        let start = *bounds.start();
        let end = self.end.min(*bounds.end());
        start..=end
    }
}
