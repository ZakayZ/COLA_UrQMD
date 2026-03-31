#pragma once

#include <COLA/EventData.hh>

#include <array>

namespace urqmd_test {

int HepChgTimes3(int pdgId) {
    // ichg(1:109) from urqmd-4.0/hepchg.f
    static const std::array<int, 110> hepChgTable = {
        0,
        -1, 2,   -1, 2,   -1, 2,   -1, 2,   0,   0,
        -3, 0,   -3, 0,   -3, 0,   -3, 0,   0,   0,
        0,  0,   0,  3,   0,  0,   0,  0,   0,   0,
        0,  0,   0,  3,   0,  0,   3,  0,   0,   0,
        0,  0,   0,  0,   0,  0,   0,  0,   0,   0,
        0,  6,   3,  6,   0,  0,   0,  0,   0,   0,
        0,  0,   0,  0,   0,  0,   0,  0,   0,   0,
        0,  0,   0,  0,   0,  0,   0,  0,   0,   0,
        0,  0,   0,  0,   0,  0,   0,  0,   0,   0,
        0,  0,   0,  0,   0,  0,   0,  0,   0,   0,
        0,  0,   0,  0,   0,  0,   0,  0,   0};

    int absPdg = std::abs(pdgId);
    int millionDigit = (absPdg / 1000000) % 10;
    int thousandsDigit = (absPdg / 1000) % 10;
    int hundredsDigit = (absPdg / 100) % 10;
    int tensDigit = (absPdg / 10) % 10;
    int onesDigit = absPdg % 10;
    int remnantBelow10000 = absPdg % 10000;

    int chargeTimes3 = 0;

    if (absPdg == 0 || absPdg >= 10000000) {
        return 0;
    }
    if (absPdg <= 100) {
        if (static_cast<std::size_t>(absPdg) < hepChgTable.size()) {
            chargeTimes3 = hepChgTable[static_cast<std::size_t>(absPdg)];
        }
    } else if (onesDigit == 0) {
        chargeTimes3 = 0;
    } else if (millionDigit > 0 && remnantBelow10000 <= 100) {
        if (static_cast<std::size_t>(remnantBelow10000) < hepChgTable.size()) {
            chargeTimes3 = hepChgTable[static_cast<std::size_t>(remnantBelow10000)];
        }
        if (absPdg == 1000017 || absPdg == 1000018) {
            chargeTimes3 = 0;
        }
        if (absPdg == 1000034 || absPdg == 1000052) {
            chargeTimes3 = 0;
        }
        if (absPdg == 1000053 || absPdg == 1000054) {
            chargeTimes3 = 0;
        }
        if (absPdg == 9900061 || absPdg == 9900062) {
            chargeTimes3 = 6;
        }
    } else if (absPdg == 9221132) {
        chargeTimes3 = 3;
    } else if (absPdg == 9331122) {
        chargeTimes3 = -6;
    } else if (thousandsDigit == 0) {
        if (static_cast<std::size_t>(hundredsDigit) < hepChgTable.size() &&
            static_cast<std::size_t>(tensDigit) < hepChgTable.size()) {
            chargeTimes3 = hepChgTable[static_cast<std::size_t>(hundredsDigit)] -
                           hepChgTable[static_cast<std::size_t>(tensDigit)];
        }
        if (hundredsDigit == 3 || hundredsDigit == 5) {
            if (static_cast<std::size_t>(tensDigit) < hepChgTable.size() &&
                static_cast<std::size_t>(hundredsDigit) < hepChgTable.size()) {
                chargeTimes3 = hepChgTable[static_cast<std::size_t>(tensDigit)] -
                               hepChgTable[static_cast<std::size_t>(hundredsDigit)];
            }
        }
    } else if (tensDigit == 0) {
        if (static_cast<std::size_t>(thousandsDigit) < hepChgTable.size() &&
            static_cast<std::size_t>(hundredsDigit) < hepChgTable.size()) {
            chargeTimes3 = hepChgTable[static_cast<std::size_t>(thousandsDigit)] +
                           hepChgTable[static_cast<std::size_t>(hundredsDigit)];
        }
    } else {
        if (static_cast<std::size_t>(thousandsDigit) < hepChgTable.size() &&
            static_cast<std::size_t>(hundredsDigit) < hepChgTable.size() &&
            static_cast<std::size_t>(tensDigit) < hepChgTable.size()) {
            chargeTimes3 = hepChgTable[static_cast<std::size_t>(thousandsDigit)] +
                           hepChgTable[static_cast<std::size_t>(hundredsDigit)] +
                           hepChgTable[static_cast<std::size_t>(tensDigit)];
        }
    }

    if (pdgId < 0 && chargeTimes3 != 0) {
        chargeTimes3 = -chargeTimes3;
    }
    return chargeTimes3;
}

int BaryonNumber(int pdgId) {
    int absPdg = std::abs(pdgId);
    if (absPdg == 0) {
        return 0;
    }
    int billionDigit = (absPdg / 1000000000) % 10;
    if (absPdg >= 10000000 && billionDigit == 1) {
        auto az = cola::PdgToAZ(pdgId);
        int massNumber = static_cast<int>(az.first);
        return pdgId > 0 ? massNumber : -massNumber;
    }
    int thousandsDigit = (absPdg / 1000) % 10;
    int tensDigit = (absPdg / 10) % 10;
    int onesDigit = absPdg % 10;
    if (onesDigit == 0) {
        return 0;
    }
    if (thousandsDigit == 0) {
        return 0;
    }
    if (tensDigit == 0) {
        return 0;
    }
    return pdgId > 0 ? 1 : -1;
}

} // namespace urqmd_test
