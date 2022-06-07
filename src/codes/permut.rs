//! Lazy permutations.

use std::fmt::Debug;
use std::ops::Index;

///
/// The result of one step on a lazy permutation of a given slice.
///
pub struct Permutation<'a, T> {
    slice: &'a [T],
    ordering: Vec<usize>,
}

impl<'a, T> Index<usize> for Permutation<'a, T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.slice[self.ordering[index]]
    }
}

impl<'a, T: Debug> Debug for Permutation<'a, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list()
            .entries((0..self.slice.len()).map(|i| &self[i]))
            .finish()
    }
}

///
/// A lazy iterator that permutates of a given slice.
///
pub struct Permut<'a, T> {
    elements: &'a [T],

    call_stack: Vec<(usize, usize)>,

    rem: Vec<usize>,
    rem_index: usize,

    so_far: Vec<usize>,
}

impl<'a, T> Permut<'a, T> {
    ///
    /// Creates a new Permuation iterator from a slice.
    ///
    pub fn new(elements: &'a [T]) -> Self {
        Self {
            call_stack: Vec::with_capacity(elements.len()),

            so_far: Vec::with_capacity(elements.len()),

            rem: (0..elements.len()).collect(),
            rem_index: 0,

            elements,
        }
    }
}

impl<'a, T> Iterator for Permut<'a, T> {
    type Item = Permutation<'a, T>;

    fn next(&mut self) -> Option<Self::Item> {
        // Try to go downwards first
        loop {
            if self.rem.is_empty() {
                let current_permut = self.so_far.clone();

                // Bump callstack + move one rem element forward
                let (rem_element, super_pos) = self.call_stack.pop().unwrap();
                self.rem.insert(super_pos, rem_element);
                self.rem_index = super_pos + 1;
                self.so_far.pop();

                // yield element
                return Some(Permutation {
                    slice: self.elements,
                    ordering: current_permut,
                });
            } else {
                // There is some rem element to step
                if self.rem_index >= self.rem.len() {
                    // Rem call ended
                    // bump call stack
                    let (rem_element, super_pos) = match self.call_stack.pop() {
                        Some(v) => v,
                        None => break,
                    };
                    self.rem.insert(super_pos, rem_element);
                    self.rem_index = super_pos + 1;
                    self.so_far.pop();

                    continue;
                }

                // There is a relevant rem element
                // Create new call frame
                let element = self.rem[self.rem_index];
                self.call_stack.push((element, self.rem_index));

                self.so_far.push(element);
                self.rem.remove(self.rem_index);
                self.rem_index = 0;
            }
        }

        None
    }
}

///
/// An extension trait to add permut as an associated funtion to slices.
///
pub trait Permutable {
    type Item;
    fn permut(&self) -> Permut<'_, Self::Item>;
}

impl<T> Permutable for [T] {
    type Item = T;
    fn permut(&self) -> Permut<'_, Self::Item> {
        Permut::new(self)
    }
}

///
/// An alternative callback based permutation.
///

pub fn permut<T>(mut elements: Vec<T>, mut callback: impl FnMut(&[T]))
where
    T: Clone + PartialEq + Debug,
{
    fn inner<T>(so_far: &mut Vec<T>, rem: &mut Vec<T>, callback: &mut impl FnMut(&[T]))
    where
        T: Clone + PartialEq + Debug,
    {
        if rem.is_empty() {
            callback(so_far);
        } else {
            for i in 0..rem.len() {
                let element = rem.remove(i);
                let cp = element.clone(); // DEBUG
                so_far.push(element);

                inner(so_far, rem, callback);

                let element = so_far.pop().unwrap();
                assert_eq!(element, cp); // DEBUG
                rem.insert(i, element);
            }
        }
    }

    inner(
        &mut Vec::with_capacity(elements.len()),
        &mut elements,
        &mut callback,
    );
}
