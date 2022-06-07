// https://crates.io/crates/bitstream-rs

use std::io::{Read, Write};

///
/// An object for writing bit to a bytewise target..
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BitWriter<W>
where
    W: Write,
{
    inner: W,
    byte: u8,
    offset: u8,
}

#[allow(unused)]
impl<W> BitWriter<W>
where
    W: Write,
{
    ///
    /// Creates a new bit-writer uverlaying another writer (e.g. std::fs::File).
    ///
    pub fn new(writer: W) -> Self {
        Self {
            inner: writer,
            byte: 0,
            offset: 0,
        }
    }

    ///
    /// Returns a reference to the contained writer.
    ///
    pub fn inner(&self) -> &W {
        &self.inner
    }

    ///
    /// Writes a byte.
    ///
    pub fn write_byte(&mut self, byte: u8) -> std::io::Result<()> {
        let mut mask = 128u8;
        for _ in 0..8 {
            let bit = (byte & mask) != 0;
            self.write_bit(bit)?;

            mask = mask >> 1;
        }
        Ok(())
    }

    ///
    /// Writes a bit. The bit may not nessecaryly be flushed.
    ///
    pub fn write_bit(&mut self, bit: bool) -> std::io::Result<()> {
        // print!("{}", if bit { "1" } else { "0" });
        if bit {
            let mask = 128_u8 >> self.offset;
            self.byte |= mask;
        }

        self.offset += 1;
        if self.offset == 8 {
            self.inner.write_all(&[self.byte])?;

            self.byte = 0;
            self.offset = 0;
        }

        Ok(())
    }

    ///
    /// Flushes and closes the writer, returning the contained handle.
    ///
    pub fn finish(mut self) -> W {
        if self.offset != 0 {
            self.inner.write_all(&[self.byte, self.offset]).unwrap();
        } else {
            self.inner.write_all(&[8u8]).unwrap();
        }

        self.inner
    }
}

// TODO: Drop Impl

///
/// A bit reader over a bytewise reader.
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BitReader<R>
where
    R: Read,
{
    inner: R,

    fin: bool,
    // num valid bytes in buf
    fill: u8,
    // bit idx in current byte
    current_bit: u8,
    // next three bytes including current
    buf: [u8; 3],
    // ???
    byte_fill: u8,
}
impl<R> BitReader<R>
where
    R: Read,
{
    ///
    /// Indicates whether the reader has any more bits to read.
    ///
    pub fn is_empty(&self) -> bool {
        self.fin && self.fill == 1 && (self.current_bit == self.byte_fill || self.byte_fill == 0)
    }

    ///
    /// Creates a new bit reader using a handle to a readable object.
    ///
    pub fn new(reader: R) -> Self {
        Self {
            inner: reader,
            fin: false,
            fill: 0,
            buf: [0; 3],
            current_bit: 0,
            byte_fill: 8,
        }
    }

    fn fill_buffer(&mut self) -> std::io::Result<()> {
        while !self.fin && self.fill != 3 {
            // Read up to 3 bytes
            match self.inner.read(&mut self.buf[self.fill as usize..]) {
                // EOF
                Ok(0) => {
                    // println!("Finshed filling buffers: {:?} fill {}", self.buf, self.fill);
                    self.fin = true;
                    self.fill -= 1;
                    // TODO: CHECK
                    self.byte_fill = self.buf[self.fill as usize];
                }
                // Normal
                Ok(n) => self.fill += n as u8, /* n = 3 */
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }

    ///
    /// Reads a byte.
    ///
    pub fn read_byte(&mut self) -> std::io::Result<u8> {
        let mut result = 0;

        for _ in 0..8 {
            result = result << 1;
            let bit = self.read_bit()?;
            if bit {
                result |= 1;
            }
        }
        Ok(result)
    }

    ///
    /// Reads a bit.
    ///
    pub fn read_bit(&mut self) -> std::io::Result<bool> {
        // if valid bytes && bit idx == byte_fill == 8
        if self.fill > 0 && self.current_bit == self.byte_fill {
            self.buf = [self.buf[1], self.buf[2], 0];
            self.current_bit = 0;
            self.fill -= 1;
        }

        self.fill_buffer()?;
        if self.fill > 0 {
            // dbg!(self.fill, self.byte_fill, self.current_bit);
            if 128u8.checked_shr(self.current_bit as u32).is_none() {
                dbg!(
                    self.fin,
                    self.fill,
                    self.buf,
                    self.current_bit,
                    self.byte_fill
                );
            }

            let res = (self.buf[0] & (128u8 >> self.current_bit)) != 0;
            // let res = (self.buf[0] & (128u8 >> self.current_bit)) == (128u8 >> self.current_bit);
            self.current_bit += 1;
            // print!("{}", if res { "1" } else { "0" });
            Ok(res)
        } else {
            unreachable!()
        }
    }
}
