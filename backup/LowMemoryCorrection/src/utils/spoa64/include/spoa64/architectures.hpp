// Copyright (c) 2020 Mario Brcic, Robert Vaser

#ifndef SPOA_ARCHITECTURES_HPP_64
#define SPOA_ARCHITECTURES_HPP_64

namespace spoa64 {

enum class Architecture {
  kAVX2,
  kSSE4_1,
  kSSE2,
  kAutomatic
};

}  // namespace spoa

#endif  // SPOA_ARCHITECTURES_HPP_
