use core::ops::Range;
use std::convert::{From, TryFrom};
use std::iter::{FromIterator, IntoIterator};

use rand::rngs::SmallRng;
use rand::{RngCore, SeedableRng};

/// Base conversion errors for nucleotide sequence conversions.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum NucleotideConversionError {
    /// Input is an invalid base. The attached u8 is the IUPACNA coding of the base.
    InvalidBase(u8),
    /// Ambiguous base that cannot be represented in the target encoding.
    AmbiguousBase,
}

impl std::fmt::Display for NucleotideConversionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for NucleotideConversionError {}

const NCBI8NA_TO_IUPACNA: [u8; 16] = [
    u8::MAX,
    b'A',
    b'C',
    b'M',
    b'G',
    b'R',
    b'S',
    b'V',
    b'T',
    b'W',
    b'Y',
    b'H',
    b'K',
    b'D',
    b'B',
    b'N',
];

const fn generate_iupacna_to_ncbi8na() -> [u8; 256] {
    let mut table = [u8::MAX; 256];
    let mut i = 1;
    while i < NCBI8NA_TO_IUPACNA.len() {
        table[NCBI8NA_TO_IUPACNA[i] as usize] = i as u8;
        table[NCBI8NA_TO_IUPACNA[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    table
}
const IUPACNA_TO_NCBI8NA: [u8; 256] = generate_iupacna_to_ncbi8na();

const fn generate_ncbi8na_to_ncbi2na() -> [u8; 16] {
    let mut table = [u8::MAX; 16];
    let mut i = 0usize;
    while i < 4 {
        table[1 << i] = i as u8;
        i += 1;
    }
    table
}
/// Table for resolve ncbi8na->ncbi2na mapping. Ambiguous mappings return u8::MAX, unambiguous
/// mappings may return any other value. This is only intended to be used with Ncbi8naBase values
/// which guarantees that the values would be within table bounds.
const NCBI8NA_TO_NCBI2NA: [u8; 16] = generate_ncbi8na_to_ncbi2na();

/// Single byte base representation where each of 4 bits represents an unambiguous base.
///
/// `count_zeros()` of each unambiguous base (ACGT) is 1; ambiguous bases are represented by setting
/// bits corresponding to each of the unambiguous bases that it may be (N = 0xf). Note that in
/// this representation 0 is invalid.
#[repr(transparent)]
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct Ncbi8naBase(u8);

impl std::fmt::Display for Ncbi8naBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Convert IUPACNA bases to Ncbi8na format.
impl TryFrom<u8> for Ncbi8naBase {
    type Error = NucleotideConversionError;

    #[inline]
    fn try_from(value: u8) -> Result<Ncbi8naBase, Self::Error> {
        let b = IUPACNA_TO_NCBI8NA[value as usize];
        if b != u8::MAX {
            Ok(Ncbi8naBase(b))
        } else {
            Err(NucleotideConversionError::InvalidBase(value))
        }
    }
}

/// Convert Ncbi8na bases to IUPACNA format.
impl From<Ncbi8naBase> for u8 {
    #[inline]
    fn from(value: Ncbi8naBase) -> Self {
        NCBI8NA_TO_IUPACNA[value.0 as usize]
    }
}

impl From<Ncbi8naBase> for char {
    #[inline]
    fn from(value: Ncbi8naBase) -> Self {
        char::from(u8::from(value))
    }
}

/// Single byte base representation that only consumes 2 bits per base.
///
/// This may only represent unambiguous bases (A, C, G, T/U).
#[repr(transparent)]
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct Ncbi2naBase(u8);

impl TryFrom<Ncbi8naBase> for Ncbi2naBase {
    type Error = NucleotideConversionError;

    #[inline]
    fn try_from(value: Ncbi8naBase) -> Result<Self, Self::Error> {
        let b2 = NCBI8NA_TO_NCBI2NA[value.0 as usize];
        if b2 != u8::MAX {
            Ok(Ncbi2naBase(b2))
        } else {
            Err(NucleotideConversionError::AmbiguousBase)
        }
    }
}

impl From<Ncbi2naBase> for Ncbi8naBase {
    #[inline]
    fn from(value: Ncbi2naBase) -> Self {
        Ncbi8naBase(1 << value.0)
    }
}

impl From<Ncbi2naBase> for u8 {
    #[inline]
    fn from(value: Ncbi2naBase) -> Self {
        u8::from(Ncbi8naBase::from(value))
    }
}

impl From<Ncbi2naBase> for char {
    #[inline]
    fn from(value: Ncbi2naBase) -> Self {
        char::from(Ncbi8naBase::from(value))
    }
}

/// Single byte base representation that packs the unambiguous bases in values 0-3.
///
/// For unambiguous bases this is the same as `Ncbi2naBase`.
#[repr(transparent)]
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct Blast8naBase(u8);

const NCBI8NA_TO_BLAST8NA: [u8; 16] = [
    15, // invalid
    0,  // A
    1,  // C
    6,  // M
    2,  // G
    4,  // R
    9,  // S
    13, // V
    3,  // T
    8,  // W
    5,  // Y
    12, // H
    7,  // K
    11, // D
    10, // B
    14, // N
];
const fn generate_blastna8_to_ncbi8na() -> [u8; 16] {
    let mut table = [0u8; 16];
    let mut i = 0usize;
    while i < table.len() {
        table[NCBI8NA_TO_BLAST8NA[i] as usize] = i as u8;
        i += 1;
    }
    table
}
const BLASTNA8_TO_NCBI8NA: [u8; 16] = generate_blastna8_to_ncbi8na();

impl From<Ncbi8naBase> for Blast8naBase {
    fn from(value: Ncbi8naBase) -> Self {
        Blast8naBase(NCBI8NA_TO_BLAST8NA[value.0 as usize])
    }
}

impl From<Blast8naBase> for Ncbi8naBase {
    fn from(value: Blast8naBase) -> Self {
        Ncbi8naBase(BLASTNA8_TO_NCBI8NA[value.0 as usize])
    }
}

const OLD_MAX_REGION_LEN: usize = 0xf;
const OLD_MAX_OFFSET: usize = 0xffffff;
const NEW_MAX_REGION_LEN: usize = 0xfff;

/// An ambiguous nucleotide region.
#[derive(Clone, Debug)]
struct AmbiguousRegion {
    range: Range<usize>,
    base: Ncbi8naBase,
}

impl AmbiguousRegion {
    fn to_old_format(&self) -> u32 {
        debug_assert!(self.range.start <= OLD_MAX_OFFSET);
        debug_assert!(self.range.len() <= OLD_MAX_REGION_LEN);
        self.range.start as u32 | ((self.range.len() as u32) << 24) | ((self.base.0 as u32) << 28)
    }

    fn from_old_format(w: u32) -> Self {
        let start = w as usize & OLD_MAX_OFFSET;
        let end = start + ((w as usize >> 24) & OLD_MAX_REGION_LEN);
        AmbiguousRegion {
            range: start..end,
            base: Ncbi8naBase((w >> 28) as u8),
        }
    }

    fn to_new_format(&self) -> u64 {
        debug_assert!(self.range.start <= u32::MAX as usize);
        debug_assert!(self.range.len() <= NEW_MAX_REGION_LEN);
        self.range.start as u64 | ((self.range.len() as u64) << 48) | ((self.base.0 as u64) << 60)
    }

    fn from_new_format(w0: u32, w1: u32) -> Self {
        let start = w1 as usize;
        let end = start + ((w0 as usize >> 16) & NEW_MAX_REGION_LEN);
        AmbiguousRegion {
            range: start..end,
            base: Ncbi8naBase((w0 >> 28) as u8),
        }
    }

    fn end_region(len: usize) -> Self {
        AmbiguousRegion {
            range: len..len,
            base: Ncbi8naBase(0xf),
        }
    }
}

/// Encodes ambiguous regions of a nucleotide base sequence.
///
/// Bases are accepted in `Ncbi8naBase` format and each call to add returns an `Ncbi2naBase` that
/// can be encoded in a compressed sequence. If the input base is unambiguous this returns the
/// expected `Ncbi2naBase`; if the input base is ambiguous a random value will be chosen from the
/// possible bases.
///
/// When all bases have been added to the encoder use Vec<u8>::from(AmbiguityEncoder) to
/// obtain the encoded sequence data.
struct AmbiguityEncoder {
    regions: Vec<AmbiguousRegion>,
    offset: usize,
    longest_region: usize,
    rng: SmallRng,
}

impl AmbiguityEncoder {
    pub fn new(size_hint: usize) -> Self {
        AmbiguityEncoder {
            regions: vec![],
            offset: 0,
            longest_region: 0,
            rng: SmallRng::seed_from_u64(size_hint as u64),
        }
    }

    pub fn add(&mut self, base: Ncbi8naBase) -> Ncbi2naBase {
        let mut base2 = NCBI8NA_TO_NCBI2NA[base.0 as usize];
        if base2 == u8::MAX {
            base2 = self.add_ambiguous(base);
        }
        self.offset += 1;
        Ncbi2naBase(base2)
    }

    fn add_ambiguous(&mut self, base: Ncbi8naBase) -> u8 {
        if let Some(last) = self
            .regions
            .last_mut()
            .filter(|last| last.range.end == self.offset && last.base == base)
        {
            last.range.end += 1;
            self.longest_region = std::cmp::max(self.longest_region, last.range.len());
        } else {
            self.regions.push(AmbiguousRegion {
                range: self.offset..(self.offset + 1),
                base,
            })
        }

        // If the base value is 0 then this isn't a valid Ncbi8naBase value.
        debug_assert_ne!(base.0, 0);
        if base.0 == 0xf {
            (self.rng.next_u32() & 0x3) as u8
        } else {
            let mut b = base.0;
            let mut select = self.rng.next_u32() % base.0.count_ones();
            while select > 0 {
                b &= (1 << b.trailing_zeros()) - 1;
                select -= 1;
            }
            b.trailing_zeros() as u8
        }
    }

    fn use_new_format(&self) -> bool {
        self.offset > OLD_MAX_OFFSET || self.longest_region > OLD_MAX_REGION_LEN
    }

    fn into_new_format(self) -> Vec<u8> {
        let mut encoded =
            Vec::with_capacity((self.regions.len() * 2 + 1) * std::mem::size_of::<u32>());
        let hdr = ((self.regions.len() as u32) * 2) | (1u32 << 31);
        encoded.extend_from_slice(&hdr.to_be_bytes());
        for mut region in self.regions {
            while region.range.len() > NEW_MAX_REGION_LEN {
                let subr = AmbiguousRegion {
                    range: region.range.start..(region.range.start + NEW_MAX_REGION_LEN),
                    base: region.base,
                };
                encoded.extend_from_slice(&subr.to_new_format().to_be_bytes());
                region.range.start += NEW_MAX_REGION_LEN;
            }
            encoded.extend_from_slice(&region.to_new_format().to_be_bytes())
        }
        encoded
    }

    fn into_old_format(self) -> Vec<u8> {
        let mut encoded = Vec::with_capacity((self.regions.len() + 1) * std::mem::size_of::<u32>());
        encoded.extend_from_slice(&(self.regions.len() as u32).to_be_bytes());
        for region in self.regions {
            encoded.extend_from_slice(&region.to_old_format().to_be_bytes())
        }
        encoded
    }
}

impl From<AmbiguityEncoder> for Vec<u8> {
    fn from(value: AmbiguityEncoder) -> Vec<u8> {
        if value.regions.is_empty() {
            vec![]
        } else if value.use_new_format() {
            value.into_new_format()
        } else {
            value.into_old_format()
        }
    }
}

struct AmbiguityIterator<'a> {
    amb: &'a [u8],
    new_format: bool,
}

impl<'a> AmbiguityIterator<'a> {
    fn new(mut amb: &'a [u8]) -> Self {
        // round down len of amb to the nearest value divisible by size_of::<u32>()
        amb = &amb[..(amb.len() & !0x3)];
        let mut it = Self {
            amb,
            new_format: false,
        };
        let new_format = it.read_u32().filter(|h| h & (1 << 31) > 0).is_some();
        it.new_format = new_format;
        it
    }

    #[inline]
    fn read_u32(&mut self) -> Option<u32> {
        if !self.amb.is_empty() {
            let (u32_bytes, rest) = self.amb.split_at(std::mem::size_of::<u32>());
            self.amb = rest;
            Some(u32::from_be_bytes(u32_bytes.try_into().unwrap()))
        } else {
            None
        }
    }
}

impl<'a> Iterator for AmbiguityIterator<'a> {
    type Item = AmbiguousRegion;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let w0 = self.read_u32()?;
        if self.new_format {
            let w1 = self.read_u32()?;
            Some(AmbiguousRegion::from_new_format(w0, w1))
        } else {
            Some(AmbiguousRegion::from_old_format(w0))
        }
    }
}

/// BLAST database representation of a nucleotide sequence.
///
/// This representation is divided into two different parts:
/// * Sequence encoding packed down to 2 bits per base. Ambiguous bases are resolved at random.
/// * Description of the ranges containing ambiguous bases.
///
/// This can easily converted to and from a FASTA string:
/// ```
/// use blastdb_sequence_util::{Ncbi8naBase, NucleotideSequence};
/// let fasta = "acgtACGT";
/// let seq: NucleotideSequence = fasta.as_bytes().iter().map(|b| Ncbi8naBase::try_from(*b)).collect::<Result<_, _>>().unwrap();
/// let normalized_fasta: String = seq.iter().map(|b| char::from(b)).collect();  // ACGTACGT
/// ```
pub struct NucleotideSequence {
    seq: Vec<u8>,
    amb: Vec<u8>,
}

impl NucleotideSequence {
    /// Create a `NucleotideSequence` from sequence and ambiguity streams.
    pub fn new(seq: Vec<u8>, amb: Vec<u8>) -> NucleotideSequence {
        Self { seq, amb }
    }

    /// Return the number of bases in the sequence.
    pub fn len(&self) -> usize {
        if let Some(b) = self.seq.last() {
            (self.seq.len() - 1) * 4 + (*b as usize & 0x3)
        } else {
            0
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate over Ncbi8na representation of sequence data.
    pub fn iter(&self) -> SequenceIter<'_, '_> {
        SequenceIter::new(&self.seq, &self.amb)
    }

    /// Return the raw sequence bytes for storage.
    pub fn sequence_bytes(&self) -> &[u8] {
        &self.seq
    }

    /// Return the raw ambiguity bytes for storage.
    pub fn ambiguity_bytes(&self) -> &[u8] {
        &self.amb
    }
}

const NCBI2NA_SHIFT: [u8; 4] = [6, 4, 2, 0];

impl FromIterator<Ncbi8naBase> for NucleotideSequence {
    fn from_iter<T>(into_iter: T) -> Self
    where
        T: IntoIterator<Item = Ncbi8naBase>,
    {
        let iter = into_iter.into_iter();
        let (lower, upper) = iter.size_hint();
        let hint = upper.unwrap_or(lower);
        let mut seq = Vec::with_capacity((hint / 4) + 1);
        let mut amb = AmbiguityEncoder::new(hint);
        let (mut buf, mut nbuf) = (0u8, 0usize);
        for base in iter {
            let base2 = amb.add(base);
            buf |= base2.0 << NCBI2NA_SHIFT[nbuf & 0x3];
            nbuf += 1;
            if nbuf == 4 {
                seq.push(buf);
                buf = 0;
                nbuf = 0;
            }
        }
        buf |= (nbuf & 0x3) as u8;
        seq.push(buf);
        NucleotideSequence {
            seq,
            amb: amb.into(),
        }
    }
}

impl std::str::FromStr for NucleotideSequence {
    type Err = NucleotideConversionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.as_bytes()
            .iter()
            .map(|b| Ncbi8naBase::try_from(*b))
            .collect::<Result<_, _>>()
    }
}

struct CompressedSequenceIter<'a> {
    seq: &'a [u8],
    len: usize,
    i: usize,
}

impl<'a> CompressedSequenceIter<'a> {
    fn new(seq: &'a [u8]) -> Self {
        let rem = seq.last().copied().unwrap_or(0) as usize & 0x3;
        let len = (seq.len() - 1) * 4 + rem;
        Self { seq, len, i: 0 }
    }
}

impl<'a> Iterator for CompressedSequenceIter<'a> {
    type Item = (usize, Ncbi2naBase);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.len {
            let i = self.i;
            let raw_base = (self.seq[i / 4] >> NCBI2NA_SHIFT[i & 3]) & 3;
            self.i += 1;
            Some((i, Ncbi2naBase(raw_base)))
        } else {
            None
        }
    }

    #[inline]
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.i += n;
        self.next()
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.len, Some(self.len))
    }
}

impl<'a> ExactSizeIterator for CompressedSequenceIter<'a> {}

pub struct SequenceIter<'a, 'b> {
    seq_it: CompressedSequenceIter<'a>,
    amb_it: AmbiguityIterator<'b>,
    amb_region: AmbiguousRegion,
}

impl<'a, 'b> SequenceIter<'a, 'b> {
    fn new(seq: &'a [u8], amb: &'b [u8]) -> Self {
        let seq_it = CompressedSequenceIter::new(seq);
        let mut amb_it = AmbiguityIterator::new(amb);
        let amb_region = amb_it
            .next()
            .unwrap_or_else(|| AmbiguousRegion::end_region(seq_it.len()));
        Self {
            seq_it,
            amb_it,
            amb_region,
        }
    }

    fn next_region(&mut self) -> AmbiguousRegion {
        self.amb_it
            .next()
            .unwrap_or_else(|| AmbiguousRegion::end_region(self.seq_it.len()))
    }
}

impl<'a, 'b> Iterator for SequenceIter<'a, 'b> {
    type Item = Ncbi8naBase;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        #[allow(clippy::iter_nth_zero)]
        self.nth(0)
    }

    #[inline]
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        let (i, base2) = self.seq_it.nth(n)?;
        while i >= self.amb_region.range.end {
            self.amb_region = self.next_region();
        }
        if i < self.amb_region.range.start {
            Some(base2.into())
        } else {
            Some(self.amb_region.base)
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.seq_it.size_hint()
    }
}

impl<'a, 'b> ExactSizeIterator for SequenceIter<'a, 'b> {}

/// Render the sequence as an IUPAC nucleic acide sequence.
impl std::fmt::Display for NucleotideSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.iter().map(char::from).collect::<String>())
    }
}

#[cfg(test)]
mod test {
    use super::{Ncbi2naBase, Ncbi8naBase, NucleotideConversionError, NucleotideSequence};

    #[test]
    fn convert_iupacna_ncbi8na() {
        assert_eq!(u8::from(Ncbi8naBase::try_from(b'A').unwrap()), b'A');
        assert_eq!(u8::from(Ncbi8naBase::try_from(b'A').unwrap()), b'A');
        assert_eq!(u8::from(Ncbi8naBase::try_from(b'G').unwrap()), b'G');
        assert_eq!(u8::from(Ncbi8naBase::try_from(b'N').unwrap()), b'N');
        assert_eq!(
            Ncbi8naBase::try_from(b'E'),
            Err(NucleotideConversionError::InvalidBase(b'E'))
        );
    }

    #[test]
    fn convert_iupacna_ncbi2na() {
        assert_eq!(
            u8::from(Ncbi2naBase::try_from(Ncbi8naBase::try_from(b'A').unwrap()).unwrap()),
            b'A'
        );
        assert_eq!(
            u8::from(Ncbi2naBase::try_from(Ncbi8naBase::try_from(b'C').unwrap()).unwrap()),
            b'C'
        );
        assert_eq!(
            u8::from(Ncbi2naBase::try_from(Ncbi8naBase::try_from(b'G').unwrap()).unwrap()),
            b'G'
        );
        assert_eq!(
            u8::from(Ncbi2naBase::try_from(Ncbi8naBase::try_from(b'T').unwrap()).unwrap()),
            b'T'
        );
        assert_eq!(
            Ncbi2naBase::try_from(Ncbi8naBase::try_from(b'D').unwrap()),
            Err(NucleotideConversionError::AmbiguousBase)
        );
    }

    fn test_sequence_conversion(
        iupac_seq: &str,
        expected_seq_bytes: usize,
        expected_amb_words: usize,
    ) {
        let seq = iupac_seq
            .as_bytes()
            .into_iter()
            .map(|b| Ncbi8naBase::try_from(*b))
            .collect::<Result<NucleotideSequence, _>>()
            .unwrap();
        assert_eq!(seq.len(), iupac_seq.len());
        assert_eq!(seq.sequence_bytes().len(), expected_seq_bytes);
        assert_eq!(
            seq.ambiguity_bytes().len() / std::mem::size_of::<u32>(),
            expected_amb_words
        );
        assert_eq!(seq.to_string(), iupac_seq);
    }

    #[test]
    fn simple_sequence() {
        test_sequence_conversion("ACGTACGT", 3, 0);
    }

    #[test]
    fn odd_base_sequence() {
        test_sequence_conversion("ACGTA", 2, 0);
    }

    #[test]
    fn completely_ambiguous_old_format() {
        test_sequence_conversion("NNNNNNNN", 3, 2);
    }

    #[test]
    fn completely_ambiguous_new_format() {
        test_sequence_conversion("NNNNNNNNNNNNNNNN", 5, 3);
    }

    #[test]
    fn some_ambiguous() {
        test_sequence_conversion("AAGNCAGTNNGGCNTA", 5, 4);
    }

    #[test]
    fn mixed_ambiguous() {
        test_sequence_conversion("AAGYCAGTNMGGCNTA", 5, 5);
    }
}
