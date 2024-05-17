#include <string>
#include <vector>
#include <iostream>
#include <fstream>

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

void StructureDSSP (gemmi::Structure& st1){
    std::stringstream buf;
    std::stringstream obuf;
    
    // lmw - I had started out copying CIF-based data strcutures, but it turns out that DSSP's CIF reader doesn't work with Gemmi CIF's.
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
    std::vector<std::string> tmp_raw_remarks = st1.raw_remarks;   // DSSP's PDB reader doesn't like Gemmi REMARKS
    st1.raw_remarks.clear();
    write_pdb(st1, buf, opt);
    st1.raw_remarks = tmp_raw_remarks;

    // Now use DSSP's tools to process the PDB file.
    cif::file f = cif::pdb::read(buf);   
    dssp dssp(f.front(), 1, 3, true);
    dssp.annotate(f.front(), false, true);
    
    
    // Now read in the CIF-formatted output from DSSP and copy the helix and sheet structures (these should be complete)
   obuf << f.front();
   gemmi::Structure st_tmp = gemmi::make_structure(gemmi::cif::read_string(obuf.str()));
   st1.helices.clear();
   st1.sheets.clear();
   st1.helices = st_tmp.helices;
   st1.sheets = st_tmp.sheets;
    
    // OLD CODE ALERT: the code below was the piece-wise attempt before we tried the copy above. If the copy doesn't work, we can do the one below.

// std::cout << f.front() << std::endl;
    
//    dssp.write_legacy_output(obuf);
//    std::cout << obuf.str() << std::endl;
//    std::cout << f.front() << std::endl;
    
//    obuf << f.front();
    
//    gemmi::Structure st_tmp = gemmi::make_structure(gemmi::cif::read_string(obuf.str()));

// dssp.write_legacy_output(std::cout);



//     if (dssp.empty())   {
//         std::cout << "No secondary structure information found" << std::endl; 
//     } else {
//     // Go through all of the Helices from DSSP and populate gemmi::Structure Helices
//     	auto st = dssp.begin(), lt = st;
// 		auto lastSS = st->type();
// 
// 		for (auto t = dssp.begin();; lt = t, ++t)
// 		{
// 			bool stop = t == dssp.end();
// 
// 			bool flush = (stop or t->type() != lastSS);
// 
// 			if (flush and lastSS != dssp::structure_type::Loop)
// 			{
// 				auto &rb = *st;
// 				auto &re = *lt;
// 
// 				int id = 0;
// 				switch (lastSS)
// 				{
// 					case dssp::structure_type::Helix_3:
// 						id = gemmi::Helix::HelixClass::R310;    // HELX_RH_3T_P
// 						break;
// 
// 					case dssp::structure_type::Alphahelix:
// 						id = gemmi::Helix::HelixClass::RAlpha;  // HELX_RH_AL_P
// 						break;
// 
// 					case dssp::structure_type::Helix_5:
// 						id = gemmi::Helix::HelixClass::RPi;     // HELX_RH_PI_P
// 						break;
// 
// 					case dssp::structure_type::Helix_PPII:
// 						id = gemmi::Helix::HelixClass::HelixPolyProlineNone;    // HELX_LH_PP_P
// 						break;
// 				}
// 				if (id > 0) {
//                     gemmi::Helix helix;
//                     helix.start.chain_name = rb.auth_asym_id();
//                     helix.start.res_id = gemmi::pdb_impl::read_res_id(std::to_string(rb.auth_seq_id()).c_str(),rb.compound_id().c_str());
//                     helix.end.chain_name = re.auth_asym_id();
//                     helix.end.res_id = gemmi::pdb_impl::read_res_id(std::to_string(re.auth_seq_id()).c_str(),re.compound_id().c_str());
//                     helix.set_helix_class_as_int(id);
//                     st1.helices.emplace_back(helix);
//                 } else {
//                     id = 0;
//                     switch (lastSS)
//                     {
//                         case dssp::structure_type::Strand:
// //                        case dssp::structure_type::Betabridge:
//                             id = 10;
//                             break;
//                     }
//                     if (id > 0) {
//                           gemmi::Sheet& sheet = gemmi::impl::find_or_add(st1.sheets, sheet_id);
//                           sheet.strands.emplace_back();
//                           gemmi::Sheet::Strand& strand = sheet.strands.back();
//                           strand.start.chain_name = rb.auth_asym_id();
//                           strand.start.res_id = gemmi::pdb_impl::read_res_id(std::to_string(rb.auth_seq_id()).c_str(),rb.compound_id().c_str());
//                           strand.end.chain_name = re.auth_asym_id();
//                           strand.end.res_id = gemmi::pdb_impl::read_res_id(std::to_string(re.auth_seq_id()).c_str(),re.compound_id().c_str());
//                           strand.sense = 0;
//                     }                
//                 }
// /*
// 				if (foundTypes.count(id) == 0)
// 				{
// 					structConfType.emplace({ { "id", id },
// 						{ "criteria", "DSSP" } });
// 					foundTypes[id] = 1;
// 				}
// */
// 				st = t;
// 			}
// 
// 			if (stop)
// 				break;
// 
// 			if (lastSS != t->type())
// 			{
// 				st = t;
// 				lastSS = t->type();
// 			}
// 		}
//     
//     // Go through all of the Sheets from DSSP and populate gemmi::Structure Sheets
// 
//     }
// 

    
}


int main() {
    // Read in the pdb file
//    gemmi::Structure st1 = gemmi::read_pdb_file("5db3_out.pdb");
    gemmi::Structure st1 = gemmi::read_pdb_file("4y5u.pdb");
    
    gemmi::setup_entities(st1);

    // Create some remarkas
    std::vector<std::string> new_raw_remarks;
    new_raw_remarks.push_back("REMARK   TEST THIS IS OURS 1");
    new_raw_remarks.push_back("REMARK   TEST THIS IS OURS 2");
    new_raw_remarks.push_back("REMARK   TEST THIS IS OURS 3");
//    st1.raw_remarks = new_raw_remarks;
    
    // Make sure we have standard PDB file information.
//    st1.info["_struct_keywords.pdbx_keywords"] = "PROTEIN";     // NOTE: we should update this to be what it is (PROTEIN, DNA, RNA, mixtures)
//    if (st1.get_info("_struct.title").empty()) st1.info["_struct.title"] = "QM/MM Refinement using DivCon Suite vXXX"; 

//    st1.info["_software.name"] = "DivCon";
//    st1.info["_software.version"] = "DEV.1234";

    // Use DSSP to obtain secondary structures   
//    StructureDSSP (st1);

    // If the SEQRES is not set & no fasta provided, then set it to the sequence we know.
//     for (gemmi::Model& model : st1.models)
//         for (gemmi::Chain& ch : model.chains){
//             gemmi::Entity* entity = st1.get_entity_of(ch.get_polymer());
//             if (entity->full_sequence.empty() && entity->entity_type == gemmi::EntityType::Polymer)
//                 for (std::string& tmpString : ch.get_polymer().extract_sequence())
//                     entity->full_sequence.push_back (tmpString);
//         }
    
    
    // Ok, time to write:
    std::ofstream of1("output.pdb");
    gemmi::PdbWriteOptions opt=gemmi::PdbWriteOptions();
    write_pdb(st1, of1, opt);

    std::ofstream of2("output.cif");
    gemmi::cif::Document doc;
    doc.blocks.resize(1);
    gemmi::MmcifOutputGroups groups(true);
    groups.auth_all = true;
    gemmi::update_mmcif_block(st1, doc.blocks[0], groups);    
    gemmi::cif::write_cif_to_stream(of2, doc, gemmi::cif::WriteOptions());

    

    return 0;
}
