//! Huffman encoding and decoding, as well as source sampeling and bitstream manipulation.

use std::io::{Read, Write};

mod bitstream;
mod source;

pub use bitstream::*;
pub use source::*;

///
/// Computes the blueprint for a huffman tree in index form
///
pub fn compute_blueprint(p: &mut [usize], n: usize) -> (Vec<usize>, Vec<bool>) {
    let mut pred = vec![0; 2 * n - 2];
    let mut mark = vec![false; 2 * n - 2];

    // [0 .. n-1] leafs [ n .. 2n-1 ] inners

    pred[0] = n;
    pred[1] = n;
    mark[0] = false /*0*/ ;
    mark[1] = true;

    p[n] = p[0] + p[1];

    let mut k = 2;
    let mut h = n;

    for b in (n + 1)..(2 * n - 1) {
        let i;
        let j;

        if (k < n) && p[k] <= p[h] {
            i = k;
            k += 1;
        } else {
            i = h;
            h += 1;
        }

        if (k < n) && (b == h || p[k] <= p[h]) {
            j = k;
            k += 1;
        } else {
            j = h;
            h += 1;
        }

        //

        assert_ne!(i, j);
        pred[i] = b;
        pred[j] = b;
        mark[i] = false;
        mark[j] = true;
        p[b] = p[i] + p[j];
    }

    (pred, mark)
}

///
/// Computes a huffman endcoding tree using the given source distribution.
///
pub fn compute_tree(distr: &FullSourceDistr) -> Tree {
    // (char, count)
    let mut alph = distr.alph_with_p().collect::<Vec<(usize, usize)>>();
    alph.sort_by(|lhs, rhs| lhs.1.cmp(&rhs.1));

    let n = alph.len();

    let mut p = alph.iter().map(|e| e.1).collect::<Vec<usize>>();
    p.resize(2 * n - 1, 0);

    // index --> value
    let label = alph.iter().map(|e| e.0).collect::<Vec<usize>>();
    // value --> index
    let label_index = vec![usize::MAX; 256]
        .into_iter()
        .enumerate()
        .map(|(val, _v)| {
            if let Some((i, _)) = label.iter().enumerate().find(|(_, e)| **e == val) {
                i
            } else {
                _v
            }
        })
        .collect::<Vec<usize>>();

    let (pred, mark) = compute_blueprint(&mut p, n);

    let leafs = vec![Box::new(TreeNode::default()); n];
    let mut inners = vec![Some(Box::new(TreeNode::default())); 2 * n - 1];

    for (i, mut leaf) in leafs.into_iter().enumerate() {
        leaf.value = label[i] as u8;
        if !mark[i] {
            inners[pred[i]].as_mut().unwrap().left = Some(leaf);
        } else {
            inners[pred[i]].as_mut().unwrap().right = Some(leaf);
        }
    }

    for i in n..(2 * n - 2) {
        if !mark[i] {
            inners[pred[i]].as_mut().unwrap().left = Some(inners[i].take().unwrap());
        } else {
            inners[pred[i]].as_mut().unwrap().right = Some(inners[i].take().unwrap());
        }
    }

    let root = inners[2 * n - 2].take().unwrap();

    Tree {
        root,
        mark,
        pred,
        label_index,
        n,
    }
}

///
/// Encodes a text using a given huffman tree as encoder.
///
/// This object additionally provides a cache for efficency of bigger datasets.
///
pub struct Encoder<'a> {
    pub(crate) tree: &'a Tree,
    pub(crate) cache: [Vec<bool>; 256],
}

impl<'a> Encoder<'a> {
    pub fn new(tree: &'a Tree) -> Self {
        Self {
            tree,
            cache: [(); 256].map(|_| Vec::new()),
        }
    }

    #[allow(unused)]
    pub fn fill_cache(&mut self) {
        for i in 0..256 {
            self.cache[i].clear();
            let byte = i as u8;

            let index = self.tree.label_index[byte as usize];
            if index == usize::MAX {
                // byte not encoded
                continue;
            }

            fn iternal_encode(this: &mut Encoder<'_>, index: usize, i: usize) {
                let next = this.tree.pred[index];
                if next != 2 * this.tree.n - 2 {
                    iternal_encode(this, next, i)
                }
                this.cache[i].push(this.tree.mark[index])
            }

            iternal_encode(self, index, i)
        }
    }

    pub fn encode_bytes_cached<W>(
        &self,
        text: impl Iterator<Item = u8>,
        stream: &mut BitWriter<W>,
    ) -> Result<(), std::io::Error>
    where
        W: Write,
    {
        for byte in text {
            self.encode_byte_cached(byte, stream)?
        }

        Ok(())
    }

    pub fn encode_byte_cached<W>(
        &self,
        byte: u8,
        stream: &mut BitWriter<W>,
    ) -> Result<(), std::io::Error>
    where
        W: Write,
    {
        let bits = &self.cache[byte as usize];
        debug_assert!(!bits.is_empty());
        for bit in bits {
            stream.write_bit(*bit)?
        }

        Ok(())
    }
}

///
/// Decodes a given text using a decoding tree, and a reader,
/// putting the results onto an output stream.
///
pub fn decode_text<R>(
    tree: &TreeNode,
    stream: &mut BitReader<R>,
    out: &mut impl Write,
) -> std::io::Result<()>
where
    R: Read,
{
    while !stream.is_empty() {
        out.write_all(&[decode_byte(&tree, stream)])?
    }

    Ok(())
}

///
/// Decodes a single byte using a decoding tree, from a given bit stream.
///
pub fn decode_byte<R>(tree: &TreeNode, stream: &mut BitReader<R>) -> u8
where
    R: Read,
{
    if tree.is_leaf() {
        tree.value
    } else {
        let bit = stream.read_bit().unwrap();
        if bit {
            decode_byte(tree.right.as_ref().unwrap(), stream)
        } else {
            decode_byte(tree.left.as_ref().unwrap(), stream)
        }
    }
}

///
/// A huffman tree defining the word-cword mapping and some metadata.
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Tree {
    root: Box<TreeNode>,
    mark: Vec<bool>,
    pred: Vec<usize>,
    label_index: Vec<usize>,
    n: usize,
}

impl Tree {
    ///
    /// Returns the encoding/decoding tree.
    ///
    pub fn tree(&self) -> &TreeNode {
        &self.root
    }

    ///
    /// Computes the average code word length under the distribution
    /// of source characters provided.
    ///
    pub fn avg_code_word_size(&self, distr: &FullSourceDistr) -> f64 {
        let mut sum = 0.0;
        let mut encoder = Encoder::new(self);
        encoder.fill_cache();

        for (c, p) in distr.alph_with_p() {
            let c = c as u8;

            let mut stream = BitWriter::new(Vec::new());
            encoder.encode_byte_cached(c, &mut stream).unwrap();
            let bytes = stream.finish();

            assert!(bytes.len() >= 2);

            let len = (bytes.len() - 2) * 8 + (*bytes.last().unwrap() as usize);
            sum += len as f64 * (p as f64 / distr.len() as f64);
        }
        sum
    }
}

///
/// A part of an encoding/decoding tree.
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TreeNode {
    left: Option<Box<TreeNode>>,
    right: Option<Box<TreeNode>>,
    value: u8,
}

impl TreeNode {
    ///
    /// The number of leafs under this node.
    ///
    pub fn len(&self) -> usize {
        if self.is_leaf() {
            1
        } else {
            self.left.as_ref().unwrap().len() + self.right.as_ref().unwrap().len()
        }
    }

    ///
    /// Indicates whether this element is a leaf.
    ///
    pub fn is_leaf(&self) -> bool {
        self.left.is_none()
    }

    ///
    /// Encodes the tree itself onto a bitstream.
    ///
    pub fn encode<W>(&self, stream: &mut BitWriter<W>) -> std::io::Result<()>
    where
        W: Write,
    {
        if self.is_leaf() {
            stream.write_bit(true)?;
            stream.write_byte(self.value)
        } else {
            stream.write_bit(false)?;
            self.left.as_ref().unwrap().encode(stream)?;
            self.right.as_ref().unwrap().encode(stream)
        }
    }

    ///
    /// Decodes the tree itself from a bitstream.
    ///
    pub fn decode<R>(stream: &mut BitReader<R>) -> std::io::Result<Box<Self>>
    where
        R: Read,
    {
        let is_leaf = stream.read_bit()?;
        if is_leaf {
            Ok(Box::new(TreeNode {
                left: None,
                right: None,
                value: stream.read_byte()?,
            }))
        } else {
            let left = Self::decode(stream)?;
            let right = Self::decode(stream)?;

            Ok(Box::new(TreeNode {
                left: Some(left),
                right: Some(right),
                value: 0,
            }))
        }
    }
}

impl Default for TreeNode {
    fn default() -> Self {
        Self {
            left: None,
            right: None,
            value: 0,
        }
    }
}
