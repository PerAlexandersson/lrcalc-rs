#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Partition {
    parts: Vec<i32>,
}

impl Partition {
    pub fn new(parts: Vec<i32>) -> Option<Self> {
        let partition = Self { parts };
        partition.is_valid().then_some(partition)
    }

    pub fn from_slice(parts: &[i32]) -> Option<Self> {
        Self::new(parts.to_vec())
    }

    pub fn parts(&self) -> &[i32] {
        &self.parts
    }

    pub fn size(&self) -> i32 {
        self.parts.iter().sum()
    }

    pub fn length_without_trailing_zeroes(&self) -> usize {
        self.parts
            .iter()
            .rposition(|&part| part != 0)
            .map_or(0, |index| index + 1)
    }

    pub fn is_valid(&self) -> bool {
        self.parts.windows(2).all(|window| window[0] >= window[1])
            && self.parts.iter().all(|&part| part >= 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validates_partitions() {
        assert!(Partition::from_slice(&[3, 2, 1]).is_some());
        assert!(Partition::from_slice(&[3, 0, 0]).is_some());
        assert!(Partition::from_slice(&[2, 3]).is_none());
        assert!(Partition::from_slice(&[2, -1]).is_none());
    }
}
