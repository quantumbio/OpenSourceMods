// gemmi-convert-proper.cpp

// ================================
// Compilation Instructions
// ================================
// If you have Gemmi cloned locally and built nothing, compile using:
//    $ g++ -std=c++11 gemmi-convert-proper.cpp -o gemmi-convert-proper -lz -I/path/to/gemmi/include
//
// Example:
//    $ g++ -std=c++14 gemmi-convert-proper.cpp     -Igemmi/include     -Lgemmi/build -lgemmi_cpp -lz     -o gemmi-convert-proper
//    $ ./gemmi-convert-proper input.pdb output.cif
//
//
// If you've built libgemmi.a using CMake:
//    $ g++ -std=c++14 gemmi-convert-proper.cpp -o gemmi-convert-proper -lz -lgemmi -Igemmi/include -Lgemmi/build
//

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/to_pdb.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input.pdb output.[pdb|cif]" << std::endl;
        return 1;
    }

    std::string in = argv[1];
    std::string out = argv[2];

    gemmi::Structure st = gemmi::read_pdb_file(in);

    if (out.size() >= 4 && out.compare(out.size() - 4, 4, ".pdb") == 0) {
        std::ofstream fp(out);
        gemmi::PdbWriteOptions opt;
        opt.minimal_file = false
        gemmi::write_pdb(st, fp, opt);
        std::cout << "Wrote PDB → " << out << std::endl;
    }
    else {
        std::cerr << "Unsupported output format. Use .pdb or .cif\n";
        return 1;
    }

    return 0;
}
