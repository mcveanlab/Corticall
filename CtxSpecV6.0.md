# Cortex Graph (.ctx) File Specification

This is for the version 6.0 graph format.  Types are C99.  Types are little-endian with the exception of the binary kmer, which is big-endian.  This spec follows the [example found here](https://hackworthy.blogspot.se/2013/01/how-to-write-binary-file-format.html).

## Header

#### Table 1. Cortex header
| Offset        | Size (Bytes) | Type/Contents | Description                                           |
|---------------|--------------|---------------|-------------------------------------------------------|
| 0             | 6            | "CORTEX"      | Signature                                             |
| 6             | 4            | `0x6`         | The graph format version                              |
| 10            | 4            | uint32_t      | The kmer size, `k`                                    |
| 14            | 4            | uint32_t      | The kmer container size in uint64_t's, `s`            |
| 18            | 4            | uint32_t      | Number of colors, `c`, in the graph                   |
| 22            | 4*`c`        | uint32_t[`c`] | The mean read length for each color                   |
| 22+(4*`c`)    | 8*`c`        | uint64_t[`c`] | Total sequence for each color (definition is unclear) |
| 22+(12*`c`)   | ?            | Block[]       | `c` color name blocks (see Table 2)                   |
| 22+(12*`c`)+? | 16           | char[16]      | Error rate (currently unused)                         |
| 38+(12*`c`)+? | ?            | Block[]       | `c` color information blocks (see Table 3)            |
| 38+(12*`c`)+? | 6            | "CORTEX"      | Trailer                                               |

#### Table 2. Color name block: color name lengths and color names for each color
| Offset | Size (Bytes) | Type/Contents | Description               |
|--------|--------------|---------------|---------------------------|
| 0      | 4            | uint32_t      | Length of color name, `L` |
| 4      | `L`          | char[`L`]     | Color name                |

#### Table 3. Color information block
| Offset | Size (Bytes) | Type/Contents | Description                                                           |
|--------|--------------|---------------|-----------------------------------------------------------------------|
| 0      | 1            | _Bool         | Is tip-clipping applied?                                              |
| 1      | 1            | _Bool         | Have low-coverage unitigs been removed?                               |
| 2      | 1            | _Bool         | Have low-coverage kmers been removed?                                 |
| 3      | 1            | _Bool         | Has this graph been cleaned against another graph?                    |
| 4      | 4            | uint32_t      | Coverage threshold on unitigs                                         |
| 8      | 4            | uint32_t      | Coverage threshold on kmers                                           |
| 12     | 4            | uint32_t      | Length `G` of name of graph against which this graph has been cleaned |
| 16     | `G`          | char[`G`]     | Name of graph against which this graph has been cleaned               |

## Body
Let each kmer record size, `S`, equal `8*s + 5*c`.  Then the number of kmer records in the file is calculated as `N = (file_size - header_size)/S`.

#### Table 4. Records
| Offset | Size (Bytes)      | Type/Contents | Description        |
|--------|-------------------|---------------|--------------------|
| 0      | N*S               | Block[]       | Kmer record blocks |

#### Table 5. Record block
| Offset        | Size (Bytes) | Type/Contents | Description              |
|---------------|--------------|---------------|--------------------------|
| 0             | 8*`s`        | uint64_t[`s`] | Binary kmer (big endian, See table 6) |
| 8*`s`         | 4*`c`        | uint32_t[`c`] | Coverage per color       |
| 8*`s` + 4*`c` | 1*`c`        | uint8_t[`c`]  | Edges per color |

### Binary kmer specification
The binary kmer is a big endian representation of a fixed size string of the letters A, C, G, and T.
Each letter is represented by two bits.  The conversion table is:

#### Table 6
| Bit value | Letter |
| --------- | ------ |
| `0b00` | A |
| `0b01` | C |
| `0b10` | G |
| `0b11` | T |

### Edge specification
An edge is a bit mask that determines the presence of incoming and outgoing edges. 
From highest to lowest order bit (little endian), the edges are `acgtACGT`, 
where lower-case letters represent incoming edges, and upper-case letters represent outgoing edges.

## Pseudocode for decoding a binary kmer
```java
// decode a binary kmer (long[containerSize] into a kmer (char[kmerSize])
char[] decodeBinaryKmer(long[] binaryKmer, int kmerSize, int containerSize) {
    char[] rawKmer = new byte[kmerSize];

    for (int i = 0; i < binaryKmer.length; i++) {
        binaryKmer[i] = reverse(binaryKmer[i]);
    }

    for (int i = kmerSize - 1; i >= 0; i--) {
        rawKmer[i] = binaryNucleotideToChar(binaryKmer[containerSize - 1] & 0x3);

        shiftBinaryKmerByOneBase(binaryKmer, containerSize);
    }

    return rawKmer;
}

// return a long that was originally read in little-endian format as big-endian
long reverse(long x) {    
    ByteBuffer bbuf = ByteBuffer.allocate(8);
    bbuf.order(ByteOrder.BIG_ENDIAN);
    bbuf.putLong(x);
    bbuf.order(ByteOrder.LITTLE_ENDIAN);

    return bbuf.getLong(0);
}

// convert a binary nucleotide (00, 01, 02, 03) as a char (A, C, G, T)
char binaryNucleotideToChar(long nucleotide) {
    switch ((int) nucleotide) {
        case 0: return 'A'; break;
        case 1: return 'C'; break;
        case 2: return 'G'; break;
        case 3: return 'T'; break;
        default:
            throw error
    }
}

// shift the binary kmer array by one nucleotide to prepare to decode the next nucleotide
void shiftBinaryKmerByOneBase(long[] binaryKmer, int containerSize) {
    for(int i = containerSize - 1; i > 0; i--) {
        binaryKmer[i] >>>= 2;  // This is an unsigned right-shift operation.  Rather than
                               // preserving the sign by copying the left-most bit after
                               // the shift, this always pads the left-most bit with 0.
        binaryKmer[i] |= (binaryKmer[i-1] << 62);
    }
    binaryKmer[0] >>>= 2;
}
```

## Pseudocode for decoding a binary edge

```java
// decode binary edge into byte[] array (lowercase for incoming edges, 
// uppercase for outgoing edges, '.' for absent edge)
byte[] decodeBinaryEdge(byte edge) {
    byte[] edges = new byte[8];
    byte[] str = {'a', 'c', 'g', 't', 'A', 'C', 'G', 'T'};

    int left = (edge >> 4);
    int right = (edge & 0xf);

    for (int i = 0; i < 4; i++) {
        int leftEdge = (left & (0x1 << (3 - i)));
        edges[i] = (byte) ((leftEdge != 0) ? str[i] : '.');

        int rightEdge = (right & (0x1 << i));
        edges[i + 4] = (byte) ((rightEdge != 0) ? str[i + 4] : '.');
    }

    return edges;
}
```
