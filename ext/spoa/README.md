# Spoa

[![Latest GitHub release](https://img.shields.io/github/release/rvaser/spoa.svg)](https://github.com/rvaser/spoa/releases/latest)
![Build status for gcc/clang](https://github.com/rvaser/spoa/actions/workflows/spoa.yml/badge.svg)
[![Published in Genome Research](https://img.shields.io/badge/published%20in-Genome%20Research-blue.svg)](https://doi.org/10.1101/gr.214270.116)

Spoa (SIMD POA) is a c++ implementation of the partial order alignment (POA) algorithm (as described in 10.1093/bioinformatics/18.3.452) which is used to generate consensus sequences (as described in 10.1093/bioinformatics/btg109). It supports three alignment modes: local (Smith-Waterman), global (Needleman-Wunsch) and semi-global alignment (overlap), and three gap modes: linear, affine and convex (piecewise affine). It also supports Intel SSE4.1+ and AVX2 vectorization (marginally faster due to high latency shifts), [SIMDe](https://github.com/simd-everywhere/simde) and dispatching.

## Build

### Dependencies

- gcc 7+ | clang 4+
- (spoa_exe)(spoa_test) zlib 1.2.8+

#### Hidden

- \[optional\] USCiLab/cereal 1.3.0
- \[optional\] simd-everywhere/simde 0.7.6
- \[optional\] google/cpu_features 0.6.0
- (spoa_exe)(spoa_test) rvaser/bioparser 3.1.0
- (spoa_exe)(spoa_test) rvaser/biosoup 0.11.0
- (spoa_test) google/googletest 1.10.0

### CMake (3.12+)

```bash
git clone https://github.com/rvaser/spoa && cd spoa
cmake -B build -DCMAKE_BUILD_TYPE=Release
make -C build
```

#### Options

- `spoa_install`: generate library install target
- `spoa_build_exe`: build executable
- `spoa_build_tests`: build unit tests
- `spoa_optimize_for_native`: build with `-march=native`
- `spoa_optimize_for_portability`: build with `-msse4.1`
- `spoa_use_cereal`: use cereal library
- `spoa_use_simde`: build with SIMDe for porting vectorized code
- `spoa_use_simde_nonvec`: use SIMDe library for nonvectorized code
- `spoa_use_simde_openmp`: use SIMDe support for OpenMP SIMD
- `spoa_generate_dispatch`: use SIMDe to generate x86 dispatch

### Meson (0.60.0+)

```bash
git clone https://github.com/rvaser/spoa && cd spoa
meson setup build
ninja -C build
```

#### Options

- `exe`: build executable
- `tests`: build unit tests
- `avx2`: build with `-mavx2`
- `sse41`: build with `-msse4.1`
- `cereal`: build serialization funcitons
- `simde`: build with SIMDe
- `simde_nonvec`: use SIMDe for nonvectorized code
- `simde_openmp`: use SIMDe support for OpenMP SIMD
- `dispatch`: use SIMDe to generate x86 dispatch

## Usage

```bash
usage: spoa [options ...] <sequences>

  # default output is stdout
  <sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    -m <int>
      default: 5
      score for matching bases
    -n <int>
      default: -4
      score for mismatching bases
    -g <int>
      default: -8
      gap opening penalty (must be non-positive)
    -e <int>
      default: -6
      gap extension penalty (must be non-positive)
    -q <int>
      default: -10
      gap opening penalty of the second affine function
      (must be non-positive)
    -c <int>
      default: -4
      gap extension penalty of the second affine function
      (must be non-positive)
    -l, --algorithm <int>
      default: 0
      alignment mode:
        0 - local (Smith-Waterman)
        1 - global (Needleman-Wunsch)
        2 - semi-global
    -r, --result <int> (option can be used multiple times)
      default: 0
      result mode:
        0 - consensus (FASTA)
        1 - multiple sequence alignment (FASTA)
        2 - 0 & 1 (FASTA)
        3 - partial order graph (GFA)
        4 - 0 & 3 (GFA)
    -d, --dot <file>
      output file for the partial order graph in DOT format
    -s, --strand-ambiguous
      for each sequence pick the strand with the better alignment
    --version
      prints the version number
    -h, --help
      prints the usage

  gap mode:
    linear if g >= e
    affine if g <= q or e >= c
    convex otherwise (default)
```

## Examples

```cpp
#include <iostream>

#include "spoa/spoa.hpp"

int main(int argc, char** argv) {

  std::vector<std::string> sequences = {
      "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
      "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
      "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
      "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
      "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
      "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
  };

  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps

  spoa::Graph graph{};

  for (const auto& it : sequences) {
    auto alignment = alignment_engine->Align(it, graph);
    graph.AddAlignment(alignment, it);
  }

  auto consensus = graph.GenerateConsensus();

  std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
            << consensus << std::endl;

  auto msa = graph.GenerateMultipleSequenceAlignment();

  for (const auto& it : msa) {
    std::cerr << it << std::endl;
  }

  return 0;
}
```

## Acknowledgement

This work has been supported in part by Croatian Science Foundation under projects UIP-11-2013-7353 and IP-2018-01-5886.
