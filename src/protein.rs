/// The indices here should match the enum values in AminoAcid
const NCBISTDAA_TO_IUPACAA: [u8; 32] = [
    0xff, b'A', b'B', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q',
    b'R', b'S', b'T', b'V', b'W', b'X', b'Y', b'Z', b'U', 0xff, b'O', b'J', 0xff, 0xff, 0xff, 0xff,
];

const fn generate_iupacaa_to_ncbistdaa() -> [u8; 256] {
    let mut r = [0xffu8; 256];
    let mut i = 0;
    while i < NCBISTDAA_TO_IUPACAA.len() {
        if NCBISTDAA_TO_IUPACAA[i] != 0xff {
            r[NCBISTDAA_TO_IUPACAA[i] as usize] = i as u8;
            r[NCBISTDAA_TO_IUPACAA[i].to_ascii_lowercase() as usize] = i as u8;
        }
        i += 1;
    }
    r
}
const IUPACAA_TO_NCBISTDAA: [u8; 256] = generate_iupacaa_to_ncbistdaa();

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinConversionError {
    InvalidBase,
}

impl std::fmt::Display for ProteinConversionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for ProteinConversionError {}

/// Single byte base representation that uses 5 bits per base.
#[repr(transparent)]
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct NcbiStdaaBase(u8);

impl TryFrom<u8> for NcbiStdaaBase {
    type Error = ProteinConversionError;

    #[inline]
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        let base = IUPACAA_TO_NCBISTDAA[value as usize];
        if base != u8::MAX {
            Ok(NcbiStdaaBase(base))
        } else {
            Err(ProteinConversionError::InvalidBase)
        }
    }
}

impl From<NcbiStdaaBase> for u8 {
    #[inline]
    fn from(value: NcbiStdaaBase) -> Self {
        NCBISTDAA_TO_IUPACAA[value.0 as usize]
    }
}

impl From<NcbiStdaaBase> for char {
    #[inline]
    fn from(value: NcbiStdaaBase) -> Self {
        u8::from(value).into()
    }
}

/// BLAST database representation of a protein sequence.
#[derive(Clone, PartialEq, PartialOrd, Eq, Ord, Debug)]
pub struct ProteinSequence {
    seq: Vec<u8>,
}

impl ProteinSequence {
    /// Returns an iterator over the IUPAC amino acid representation fo the sequence.
    pub fn iter(&self) -> impl Iterator<Item = NcbiStdaaBase> + '_ {
        self.seq.iter().copied().map(NcbiStdaaBase)
    }

    /// Returns the number of bases in the sequence
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the raw sequence as bytes.
    pub fn sequence_bytes(&self) -> &[u8] {
        &self.seq[..]
    }
}

impl FromIterator<NcbiStdaaBase> for ProteinSequence {
    fn from_iter<T: IntoIterator<Item = NcbiStdaaBase>>(iter: T) -> Self {
        ProteinSequence {
            seq: iter.into_iter().map(|b| b.0).collect(),
        }
    }
}

/// Create a `ProteinSequence` from an IUPAC amino acid sequence.
/// Accepts both upper and lowercase versions of the bases.
impl std::str::FromStr for ProteinSequence {
    type Err = ProteinConversionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.as_bytes()
            .iter()
            .map(|b| NcbiStdaaBase::try_from(*b))
            .collect::<Result<_, _>>()
    }
}

/// Render the sequence as an IUPAC amino acid sequence.
impl std::fmt::Display for ProteinSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.iter().map(char::from).collect::<String>())
    }
}

#[cfg(test)]
mod test {
    use crate::protein::ProteinConversionError;

    use super::ProteinSequence;
    use std::str::FromStr;

    #[test]
    fn from_iupac() {
        let seq = ProteinSequence::from_str("AARDVARK").unwrap();
        assert_eq!(seq.len(), 8);
        assert_eq!(seq.to_string(), "AARDVARK");
        assert_eq!(seq.sequence_bytes(), [1, 1, 16, 4, 19, 1, 16, 10]);
    }

    #[test]
    fn from_iupac_invalid() {
        assert_eq!(
            ProteinSequence::from_str("AARDVARK1"),
            Err(ProteinConversionError::InvalidBase)
        );
    }
}
