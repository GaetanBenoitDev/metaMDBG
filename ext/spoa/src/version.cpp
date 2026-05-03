// Copyright (c) 2023 Robert Vaser

#include "spoa/version.hpp"

#include "spoa_config.h"

namespace spoa {

std::string Version() {
    return std::string(SPOA_VERSION);
}

}  // namespace spoa
