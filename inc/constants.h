// Common project-wide constants
#pragma once

namespace oca {

// Number of detectors actively used by displays/analysis
// Detectors are indexed 0..N_DETECTORS-1
constexpr int N_DETECTORS = 3;   // D (index 3) not used

// Total number of detectors in data/calibration (Event expects this)
constexpr int N_DETECTORS_TOTAL = 4;

// Number of channels per detector
constexpr int N_CHANNELS  = 384;

}
