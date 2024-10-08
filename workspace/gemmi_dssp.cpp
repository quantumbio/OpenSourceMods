// NOTE: to Compile:
//  $ g++ -std=c++20 gemmi_dssp.cpp -o gemmi_dssp -lcifpp -ldssp -lgemmi_cpp -lz

// On my laptop, with a local build of dssp and gemmi, I use this approach:
//      $ g++ -std=c++20 gemmi_dssp.cpp -o gemmi_dssp -lcifpp -ldssp -lgemmi_cpp -lz -Igemmi/include -Idssp/build/_deps/cifpp-src/include -Idssp/libdssp/include -Ldssp/build/_deps/cifpp-build -Ldssp/build/libdssp -Lgemmi/
//      $ export LIBCIFPP_DATA_DIR=${PWD}/dssp/build/_deps/cifpp-src/rsrc
//      $ ./gemmi_dssp 

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>

#include <gemmi/model.hpp>    // for Structure, impl::find_or_add
#include <gemmi/pdb.hpp>     // to read
#include <gemmi/to_pdb.hpp>     // to write
#include <gemmi/polyheur.hpp>     // to write
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#include <cif++.hpp>
#include <dssp.hpp>

// This function is used to update the SEQRES if no SEQRES is found in the PDB file AND if no FASTA is provided.
void ApplySEQRES (gemmi::Structure& st1)
{
     for (gemmi::Model& model : st1.models)
         for (gemmi::Chain& ch : model.chains)  {
             gemmi::Entity* entity = st1.get_entity_of(ch.get_polymer());
             if (entity && entity->full_sequence.empty() && entity->entity_type == gemmi::EntityType::Polymer)
                 for (std::string& tmpString : ch.get_polymer().extract_sequence())
                    entity->full_sequence.push_back (tmpString);
         }
}

// This function is used to call DSSP on the gemmi::Structure and determine helies and sheets
//      We should call this every time we read in a PDB, mmCIF, or MOL2 file and each time we write a PDB or mmCIF file
//          There is some below noted overhead due to DSSP's PDB->CIF reader which appears to be very slow. 
//          Therefore, we may only want to do this on READ when we think we need the data. If not then skip it. However, we should do it on write for completeness of the PDB/mmCIF file.
void StructureDSSP (gemmi::Structure& st1)
{
    std::stringstream buf;
    std::stringstream obuf;
    
    // Originally, I had explored the use of the DSSP Structure/model data structures, but this will take a lot more work and it isn't clear if it will work as well or at all.
    //      IF we end up finding the DSSP is useful for our predictions, then perhaps it would be useful to explore and/or layer on a proper API.
    
    // I had started out copying CIF-based data strcutures, but it turns out that DSSP's CIF reader doesn't work with Gemmi CIFs.
    //      If this ever works, then we can uncomment these lines and remove those associated with the PDB file.
//    gemmi::cif::Document doc;
//    doc.blocks.resize(1);
//    gemmi::MmcifOutputGroups groups(true);
//    groups.auth_all = true;
//    gemmi::update_mmcif_block(st1, doc.blocks[0], groups);    
//    gemmi::cif::write_cif_to_stream(buf, doc, gemmi::cif::WriteOptions());

    // So I switched to using PDB. I don't like it as much, but it should work.
    gemmi::PdbWriteOptions opt=gemmi::PdbWriteOptions();
    opt.minimal_file = false;
    opt.seqres_records = false;
    opt.ssbond_records = false;
    opt.link_records = false;
    opt.cispep_records = false;
    opt.end_record = false;
    std::vector<std::string> tmp_raw_remarks = st1.raw_remarks;   
    st1.raw_remarks.clear();                                        // DSSP's PDB reader doesn't like Gemmi REMARKS
    write_pdb(st1, buf, opt);
    st1.raw_remarks = tmp_raw_remarks;

    // Now use DSSP's tools to process the PDB file.
        // NOTE: The CIF++ PDB->CIF reader is extremely slow - may want to replace this with the dssp MM structure approach (but it's unclear if that will be better and it's clearly harder to implement)
    cif::file f = cif::pdb::read(buf);   

    dssp dssp(f.front(), 1, 3, true);
    dssp.annotate(f.front(), false, true);
    
    // TODO: f.front() sheets do not have sense set because Gemmi-CIF uses _struct_sheet_order but DSSP doesn't set that. Instead, it uses _dssp_struct_ladder. So to fix this, need to copy _dssp_struct_ladder -> _struct_sheet_order.
    
    // Now read in the CIF-formatted output from DSSP and copy the helix and sheet structures (these should be complete)
   obuf << f.front();
   gemmi::Structure st_tmp = gemmi::make_structure(gemmi::cif::read_string(obuf.str()));
   // TODO: ADD ERROR TRAPS (for disagreement between st1 and st_tmp)
   
   st1.helices.clear();
   st1.sheets.clear();
   st1.helices = st_tmp.helices;
   st1.sheets = st_tmp.sheets;
}


int main(int argc, char *argv[]) {
    std::string infilename;
    int return_status = 1;
    if (argc == 2) { 
        infilename = argv[1];
    
        // Read in the pdb file
        std::cout << "READING old PDBfile: " << infilename << std::endl;
        gemmi::Structure st1 = gemmi::read_pdb_file(infilename);
    //    gemmi::Structure st1 = gemmi::read_pdb_file("4y5u.pdb");
        
        gemmi::setup_entities(st1);
    
        // DEMO: HOW to create some remarks in Gemmi
        std::vector<std::string> new_raw_remarks;
        new_raw_remarks.push_back("REMARK   TEST THIS IS OURS 1");
        new_raw_remarks.push_back("REMARK   TEST THIS IS OURS 2");
        new_raw_remarks.push_back("REMARK   TEST THIS IS OURS 3");
    //    st1.raw_remarks = new_raw_remarks;                            // To stomp on the old REMARKS when we update them.
        
        // Make sure we have standard PDB file information.
        st1.info["_struct_keywords.pdbx_keywords"] = "PROTEIN";     // NOTE: we should update this to be what it is (PROTEIN, DNA, RNA, mixtures)
        if (st1.get_info("_struct.title").empty()) st1.info["_struct.title"] = "QM/MM Refinement using DivCon Suite vXXX"; 
    
        // Use DSSP to obtain secondary structures   
        StructureDSSP (st1);
        
        // Update the secondary structure to be complete
        ApplySEQRES (st1);
        
        std::string outpdbname = std::filesystem::path(infilename).stem().string() + "-dssp.pdb";
        std::cout << "WRITING new PDBfile: " << outpdbname << std::endl;
        std::ofstream of1(outpdbname);
        gemmi::PdbWriteOptions opt=gemmi::PdbWriteOptions();
        write_pdb(st1, of1, opt);
    
        std::string outcifname = std::filesystem::path(infilename).stem().string() + "-dssp.cif";
        std::cout << "WRITING new CIFfile: " << outcifname << std::endl;
        std::ofstream of2(outcifname);
        gemmi::cif::Document doc;
        doc.blocks.resize(1);
        gemmi::MmcifOutputGroups groups(true);
        groups.auth_all = true;
        gemmi::update_mmcif_block(st1, doc.blocks[0], groups);    
        gemmi::cif::write_cif_to_stream(of2, doc, gemmi::cif::WriteOptions());
        return_status = 0;
    }
    else
    {
        std::cout << "ERROR: gemmi_dssp requires a pdbfile and only a pdbfile as input:\n\t$ gemmi_dssp infile.pdb" << std::endl;
    }

    return return_status;
}
