#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Vec<T, const N: usize> {
    data: [Option<T>; N],
    length: usize,
}

impl<T: Copy, const N: usize> Vec<T, N> {
    // Create a new Vec with a fixed capacity
    pub fn new() -> Self {
        Self {
            data: [None; N],
            length: 0,
        }
    }

    // Add an element to the end of the collection
    pub fn push(&mut self, element: T) -> Result<(), &'static str> {
        if self.length >= N {
            Err("Vec is full")
        } else {
            self.data[self.length] = Some(element);
            self.length += 1;
            Ok(())
        }
    }

    // Extend the collection with the elements from a slice
    pub fn extend_from_slice(&mut self, slice: &[T]) -> Result<(), &'static str> {
        if self.length + slice.len() > N {
            return Err("Vec does not have enough capacity");
        }

        for &item in slice {
            if self.length < N {
                self.data[self.length] = Some(item);
                self.length += 1;
            }
        }

        Ok(())
    }

    // Helper method to retrieve the current length of the collection
    pub fn len(&self) -> usize {
        self.length
    }

    // Check if the collection is empty
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }
}
