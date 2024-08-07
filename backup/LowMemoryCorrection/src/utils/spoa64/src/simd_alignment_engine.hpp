// Copyright (c) 2020 Robert Vaser

#ifndef SIMD_ALIGNMENT_ENGINE_HPP_
#define SIMD_ALIGNMENT_ENGINE_HPP_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "spoa64/alignment_engine.hpp"
#include "spoa64/architectures.hpp"

namespace spoa64 {

class Graph;

std::unique_ptr<AlignmentEngine> CreateSimdAlignmentEngine(  // for dispatcher
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c);

template<Architecture A>
class SimdAlignmentEngine: public AlignmentEngine {
 public:
  SimdAlignmentEngine(const SimdAlignmentEngine&) = delete;
  SimdAlignmentEngine& operator=(const SimdAlignmentEngine&) = delete;

  SimdAlignmentEngine(SimdAlignmentEngine&&) = default;
  SimdAlignmentEngine& operator=(SimdAlignmentEngine&&) = delete;

  ~SimdAlignmentEngine() = default;

  static std::unique_ptr<AlignmentEngine> Create(
      AlignmentType type,
      AlignmentSubtype subtype,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c);

  //void Prealloc(
  //    std::uint64_t max_sequence_len,
  //    std::uint64_t alphabet_size) override;

  Alignment Align(
      const char* sequence, std::uint64_t sequence_len,
      const Graph& graph,
      std::int64_t* score) override;

  Alignment Align(
      const std::vector<std::uint64_t>& sequence, std::uint64_t sequence_len,
      const Graph& graph,
      std::int64_t* score) override;


  friend std::unique_ptr<AlignmentEngine> CreateSimdAlignmentEngine(
      AlignmentType type,
      AlignmentSubtype subtype,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c);

 private:
  SimdAlignmentEngine(
      AlignmentType type,
      AlignmentSubtype subtype,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c);

  template<typename T>
  Alignment Linear(
      std::uint64_t sequence_len,
      const Graph& graph,
      std::int64_t* score) noexcept;

  template<typename T>
  Alignment Affine(
      std::uint64_t sequence_len,
      const Graph& graph,
      std::int64_t* score) noexcept;

  template<typename T>
  Alignment Convex(
      std::uint64_t sequence_len,
      const Graph& graph,
      std::int64_t* score) noexcept;


  void Realloc(
      std::uint64_t matrix_width,
      std::uint64_t matrix_height,
      std::uint8_t num_codes);

  template<typename T>
  void Initialize(
      const char* sequence,
      const Graph& graph,
      std::uint64_t normal_matrix_width,
      std::uint64_t matrix_width,
      std::uint64_t matrix_height) noexcept;

  template<typename T>
  void Initialize(
      const std::vector<std::uint64_t>& sequence,
      const Graph& graph,
      std::uint64_t normal_matrix_width,
      std::uint64_t matrix_width,
      std::uint64_t matrix_height) noexcept;
      
  struct Implementation;
  std::unique_ptr<Implementation> pimpl_;
};

}  // namespace spoa

#endif  // SIMD_ALIGNMENT_ENGINE_HPP_
