#include "pdg_quantum_numbers.hh"

#include <COLA/EventData.hh>

#include <array>
#include <cstddef>
#include <cstdlib>

namespace urqmd_test {

  int HepChgTimes3(int pdg_id) {
    // ichg(1:109) from urqmd-4.0/hepchg.f
    static const std::array<int, 110> hep_chg_table = {
        0, -1, 2, -1, 2, -1, 2, -1, 2, 0, 0, -3, 0, -3, 0, -3, 0, -3, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0,
        0, 0,  0, 0,  0, 0,  3, 0,  0, 3, 0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 6, 3, 6, 0,
        0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0,  0, 0,  0, 0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0};

    const int abs_pdg = std::abs(pdg_id);
    const int million_digit = (abs_pdg / 1000000) % 10;
    const int thousands_digit = (abs_pdg / 1000) % 10;
    const int hundreds_digit = (abs_pdg / 100) % 10;
    const int tens_digit = (abs_pdg / 10) % 10;
    const int ones_digit = abs_pdg % 10;
    const int remnant_below_10000 = abs_pdg % 10000;

    int charge_times3 = 0;

    if (abs_pdg == 0 || abs_pdg >= 10000000) {
      return 0;
    }
    if (abs_pdg <= 100) {
      if (static_cast<std::size_t>(abs_pdg) < hep_chg_table.size()) {
        charge_times3 = hep_chg_table[static_cast<std::size_t>(abs_pdg)];
      }
    } else if (ones_digit == 0) {
      charge_times3 = 0;
    } else if (million_digit > 0 && remnant_below_10000 <= 100) {
      if (static_cast<std::size_t>(remnant_below_10000) < hep_chg_table.size()) {
        charge_times3 = hep_chg_table[static_cast<std::size_t>(remnant_below_10000)];
      }
      if (abs_pdg == 1000017 || abs_pdg == 1000018) {
        charge_times3 = 0;
      }
      if (abs_pdg == 1000034 || abs_pdg == 1000052) {
        charge_times3 = 0;
      }
      if (abs_pdg == 1000053 || abs_pdg == 1000054) {
        charge_times3 = 0;
      }
      if (abs_pdg == 9900061 || abs_pdg == 9900062) {
        charge_times3 = 6;
      }
    } else if (abs_pdg == 9221132) {
      charge_times3 = 3;
    } else if (abs_pdg == 9331122) {
      charge_times3 = -6;
    } else if (thousands_digit == 0) {
      if (static_cast<std::size_t>(hundreds_digit) < hep_chg_table.size() &&
          static_cast<std::size_t>(tens_digit) < hep_chg_table.size()) {
        charge_times3 = hep_chg_table[static_cast<std::size_t>(hundreds_digit)] -
                        hep_chg_table[static_cast<std::size_t>(tens_digit)];
      }
      if (hundreds_digit == 3 || hundreds_digit == 5) {
        if (static_cast<std::size_t>(tens_digit) < hep_chg_table.size() &&
            static_cast<std::size_t>(hundreds_digit) < hep_chg_table.size()) {
          charge_times3 = hep_chg_table[static_cast<std::size_t>(tens_digit)] -
                          hep_chg_table[static_cast<std::size_t>(hundreds_digit)];
        }
      }
    } else if (tens_digit == 0) {
      if (static_cast<std::size_t>(thousands_digit) < hep_chg_table.size() &&
          static_cast<std::size_t>(hundreds_digit) < hep_chg_table.size()) {
        charge_times3 = hep_chg_table[static_cast<std::size_t>(thousands_digit)] +
                        hep_chg_table[static_cast<std::size_t>(hundreds_digit)];
      }
    } else {
      if (static_cast<std::size_t>(thousands_digit) < hep_chg_table.size() &&
          static_cast<std::size_t>(hundreds_digit) < hep_chg_table.size() &&
          static_cast<std::size_t>(tens_digit) < hep_chg_table.size()) {
        charge_times3 = hep_chg_table[static_cast<std::size_t>(thousands_digit)] +
                        hep_chg_table[static_cast<std::size_t>(hundreds_digit)] +
                        hep_chg_table[static_cast<std::size_t>(tens_digit)];
      }
    }

    if (pdg_id < 0 && charge_times3 != 0) {
      charge_times3 = -charge_times3;
    }
    return charge_times3;
  }

  int BaryonNumber(int pdg_id) {
    const int abs_pdg = std::abs(pdg_id);
    if (abs_pdg == 0) {
      return 0;
    }
    const int billion_digit = (abs_pdg / 1000000000) % 10;
    if (abs_pdg >= 10000000 && billion_digit == 1) {
      const auto az = cola::PdgToAZ(pdg_id);
      const int mass_number = static_cast<int>(az.first);
      return pdg_id > 0 ? mass_number : -mass_number;
    }
    const int thousands_digit = (abs_pdg / 1000) % 10;
    const int tens_digit = (abs_pdg / 10) % 10;
    const int ones_digit = abs_pdg % 10;
    if (ones_digit == 0) {
      return 0;
    }
    if (thousands_digit == 0) {
      return 0;
    }
    if (tens_digit == 0) {
      return 0;
    }
    return pdg_id > 0 ? 1 : -1;
  }

}  // namespace urqmd_test
