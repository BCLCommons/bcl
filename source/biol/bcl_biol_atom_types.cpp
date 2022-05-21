// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "biol/bcl_biol_atom_types.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief construct all AtomTypes
    AtomTypes::AtomTypes() :
//                                                 ElementType                              in BB  in SC, HELIX                                                                   STRAND
      N(    AddEnum( "N",    AtomTypeData( " N  ", chemistry::GetElementTypes().e_Nitrogen, true,  false, coord::CylinderCoordinates( -0.245, 1.520, -27.600 / 180 * math::g_Pi), coord::CylinderCoordinates( 0.477, 0.391,  0.998)))), //!< Nitrogen from the peptide bond
      CA(   AddEnum( "CA",   AtomTypeData( " CA ", chemistry::GetElementTypes().e_Carbon,   true,  false, coord::CylinderCoordinates(  0.661, 2.260,   0.000 / 180 * math::g_Pi), coord::CylinderCoordinates( 1.746, 0.778, -0.377)))), //!< Carbon alpha backbone
      C(    AddEnum( "C",    AtomTypeData( " C  ", chemistry::GetElementTypes().e_Carbon,   true,  false, coord::CylinderCoordinates(  1.747, 1.678,  26.335 / 180 * math::g_Pi), coord::CylinderCoordinates( 2.949, 0.534,  1.126)))), //!< Carbon from the carboxyl group
      O(    AddEnum( "O",    AtomTypeData( " O  ", chemistry::GetElementTypes().e_Oxygen,   true,  false, coord::CylinderCoordinates(  2.922, 1.950,  20.250 / 180 * math::g_Pi), coord::CylinderCoordinates( 3.001, 1.720,  1.424)))), //!< Oxygen from the carboxyl group
      CB(   AddEnum( "CB",   AtomTypeData( " CB ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates( -0.151, 3.240,  18.680 / 180 * math::g_Pi), coord::CylinderCoordinates( 1.756, 2.275, -0.106)))), //!< Carbon beta first side chain atom
      CG(   AddEnum( "CG",   AtomTypeData( " CG ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon gamma - second side chain atom
      CG1(  AddEnum( "CG1",  AtomTypeData( " CG1", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon gamma - second side chain atom 1 for two second side chain atoms
      CG2(  AddEnum( "CG2",  AtomTypeData( " CG2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon gamma - second side chain atom 2 for two second side chain atoms
      CD(   AddEnum( "CD",   AtomTypeData( " CD ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon delta - third side chain atom
      CD1(  AddEnum( "CD1",  AtomTypeData( " CD1", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon delta - third side chain atom 1 for two third side chain atoms
      CD2(  AddEnum( "CD2",  AtomTypeData( " CD2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon delta - third side chain atom 2 for two third side chain atoms
      CE(   AddEnum( "CE",   AtomTypeData( " CE ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon epsilon - fourth side chain atom
      CE1(  AddEnum( "CE1",  AtomTypeData( " CE1", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon epsilon - fourth side chain atom 1 for two or three fourth side chain atoms
      CE2(  AddEnum( "CE2",  AtomTypeData( " CE2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon epsilon - fourth side chain atom 2 for two or three fourth side chain atoms
      CE3(  AddEnum( "CE3",  AtomTypeData( " CE3", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon epsilon - fourth side chain atom 3 for two or three fourth side chain atoms
      CZ(   AddEnum( "CZ",   AtomTypeData( " CZ ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon zeta - fifth side chain atom
      CZ2(  AddEnum( "CZ2",  AtomTypeData( " CZ2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon zeta - fifth side chain atom 1 for two fifth side chain atoms
      CZ3(  AddEnum( "CZ3",  AtomTypeData( " CZ3", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon zeta - fifth side chain atom 2 for two fifth side chain atoms
      CH2(  AddEnum( "CH2",  AtomTypeData( " CH2", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Carbon eta - sixth side chain aom as in TRP
      ND1(  AddEnum( "ND1",  AtomTypeData( " ND1", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen delta - third position
      ND2(  AddEnum( "ND2",  AtomTypeData( " ND2", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen delta - third position
      NE(   AddEnum( "NE",   AtomTypeData( " NE ", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen epsilon - fourth position nitrogen as in ARG
      NE1(  AddEnum( "NE1",  AtomTypeData( " NE1", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen epsilon - fourth position nitrogen as in TRP
      NE2(  AddEnum( "NE2",  AtomTypeData( " NE2", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen epsilon - fourth position nitrogen as in GLN or HIS
      NZ(   AddEnum( "NZ",   AtomTypeData( " NZ ", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen zeta - fifth side chain atom as in LYS
      NH1(  AddEnum( "NH1",  AtomTypeData( " NH1", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen eta 1 - sixth side chain atom as in ARG
      NH2(  AddEnum( "NH2",  AtomTypeData( " NH2", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Nitrogen eta 2 - sixth side chain atom as in ARG
      OD1(  AddEnum( "OD1",  AtomTypeData( " OD1", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen delta - third side chain atom as in ASP or ASN
      OD2(  AddEnum( "OD2",  AtomTypeData( " OD2", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen delta - third side chain atom as in ASP
      OG(   AddEnum( "OG",   AtomTypeData( " OG ", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen gamma - second side chain atom as in SER
      OG1(  AddEnum( "OG1",  AtomTypeData( " OG1", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen gamma - second side chain atom as in THR
      OE1(  AddEnum( "OE1",  AtomTypeData( " OE1", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen epsilon - fourth side chain atom as in GLN
      OE2(  AddEnum( "OE2",  AtomTypeData( " OE2", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen epsilon - fourth side chain atom as in GLN or GLU
      OH(   AddEnum( "OH",   AtomTypeData( " OH ", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Oxygen eta - sixth side chain atom as in TYR
      SD(   AddEnum( "SD",   AtomTypeData( " SD ", chemistry::GetElementTypes().e_Sulfur,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Suflur on delta carbon - like in methionin
      SE(   AddEnum( "SE",   AtomTypeData( " SE ", chemistry::GetElementTypes().e_Selenium, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Selenium on delta carbon - like in selenomethionine
      SG(   AddEnum( "SG",   AtomTypeData( " SG ", chemistry::GetElementTypes().e_Sulfur,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Sulfur on gamma carbon - like in cystein
      H(    AddEnum( "H",    AtomTypeData( " H  ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on the backbone nitrogen atom
      HA(   AddEnum( "HA",   AtomTypeData( " HA ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CA alpha Carbon
      HA2(  AddEnum( "HA2",  AtomTypeData( " HA2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates( -0.151, 3.240,  18.680 / 180 * math::g_Pi), coord::CylinderCoordinates( 1.756, 2.275, -0.106)))), //!< Hydrogen on CA for GLY in CB position
      HA3(  AddEnum( "HA3",  AtomTypeData( " HA3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CA for GLY in HA position
      HB(   AddEnum( "HB",   AtomTypeData( " HB ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CB for THR, ILE or VAL
      HB1(  AddEnum( "HB1",  AtomTypeData( " HB1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CB for ALA
      HB2(  AddEnum( "HB2",  AtomTypeData( " HB2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CB
      HB3(  AddEnum( "HB3",  AtomTypeData( " HB3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CB
      HG(   AddEnum( "HG",   AtomTypeData( " HG ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen an CG as in LEU, CYS, SER
      HG1(  AddEnum( "HG1",  AtomTypeData( " HG1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 an CG as in THR
      HG2(  AddEnum( "HG2",  AtomTypeData( " HG2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CG
      HG3(  AddEnum( "HG3",  AtomTypeData( " HG3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CG
      HG11( AddEnum( "HG11", AtomTypeData( "HG11", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CG1 as in VAL
      HG12( AddEnum( "HG12", AtomTypeData( "HG12", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CG2 as in ILE, VAL
      HG13( AddEnum( "HG13", AtomTypeData( "HG13", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CG1 as in ILE, VAL
      HG21( AddEnum( "HG21", AtomTypeData( "HG21", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CG2 as in ILE, VAL, THR
      HG22( AddEnum( "HG22", AtomTypeData( "HG22", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CG2 as in ILE, VAL, THR
      HG23( AddEnum( "HG23", AtomTypeData( "HG23", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CG2 as in ILE, VAL, THR
      HD1(  AddEnum( "HD1",  AtomTypeData( " HD1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CD
      HD2(  AddEnum( "HD2",  AtomTypeData( " HD2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CD
      HD3(  AddEnum( "HD3",  AtomTypeData( " HD3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CD
      HD11( AddEnum( "HD11", AtomTypeData( "HD11", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CD1
      HD12( AddEnum( "HD12", AtomTypeData( "HD12", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CD1
      HD13( AddEnum( "HD13", AtomTypeData( "HD13", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CD1
      HD21( AddEnum( "HD21", AtomTypeData( "HD21", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CD2
      HD22( AddEnum( "HD22", AtomTypeData( "HD22", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CD2
      HD23( AddEnum( "HD23", AtomTypeData( "HD23", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CD2
      HE(   AddEnum( "HE",   AtomTypeData( " HE ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CE as in ARG
      HE1(  AddEnum( "HE1",  AtomTypeData( " HE1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CE
      HE2(  AddEnum( "HE2",  AtomTypeData( " HE2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CE
      HE3(  AddEnum( "HE3",  AtomTypeData( " HE3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on CE
      HE21( AddEnum( "HE21", AtomTypeData( "HE21", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on CE2 as in GLN
      HE22( AddEnum( "HE22", AtomTypeData( "HE22", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on CE2 as in GLN
      HZ(   AddEnum( "HZ",   AtomTypeData( " HZ ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CZ as in PHE
      HZ1(  AddEnum( "HZ1",  AtomTypeData( " HZ1", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on Z1 as in LYS (on Nitrogen NZ)
      HZ2(  AddEnum( "HZ2",  AtomTypeData( " HZ2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on Z as in LYS (on NZ), TRP (on CZ2)
      HZ3(  AddEnum( "HZ3",  AtomTypeData( " HZ3", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 3 on Z as in LYS (on NZ), TRP (on CZ3)
      HH(   AddEnum( "HH",   AtomTypeData( " HH ", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CH as in TYR
      HH2(  AddEnum( "HH2",  AtomTypeData( " HH2", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen on CH2 as in TRP
      HH11( AddEnum( "HH11", AtomTypeData( "HH11", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on NH1 as in ARG
      HH12( AddEnum( "HH12", AtomTypeData( "HH12", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on NH1 as in ARG
      HH21( AddEnum( "HH21", AtomTypeData( "HH21", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 1 on NH2 as in ARG
      HH22( AddEnum( "HH22", AtomTypeData( "HH22", chemistry::GetElementTypes().e_Hydrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen 2 on NH2 as in ARG
      // terminal amine
      H1(   AddEnum( "H1",   AtomTypeData( " H1 ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< leaving hydrogen, connected to backbone Nitrogen
      H2(   AddEnum( "H2",   AtomTypeData( " H2 ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< leaving hydrogen, connected to backbone Nitrogen
      H3(   AddEnum( "H3",   AtomTypeData( " H3 ", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< leaving hydrogen, connected to backbone Nitrogen
      // terminal carboxylic acid
      HXT(  AddEnum( "HXT",  AtomTypeData( " HXT", chemistry::GetElementTypes().e_Hydrogen, false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< Hydrogen leaving hydrogen, connected to C if there is no peptide bond
      OXT(  AddEnum( "OXT",  AtomTypeData( " OXT", chemistry::GetElementTypes().e_Oxygen,   false, false, coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< leaving oxygen, connected to C if there is no peptide bond

      C2(   AddEnum( "C2",   AtomTypeData( " C2 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mttsl C next to N-O nearest to linker
      C3(   AddEnum( "C3",   AtomTypeData( " C3 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C in "pentane" ring connected to linker
      C4(   AddEnum( "C4",   AtomTypeData( " C4 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C in pentane ring connected to H and connected to C3 by double bond
      C5(   AddEnum( "C5",   AtomTypeData( " C5 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C next to N-O nearest to linker on side with C=C-H and opposite C2
      C6(   AddEnum( "C6",   AtomTypeData( " C6 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C above C5 when looking at ring from linker to N-O and oriented C=C-H
      C7(   AddEnum( "C7",   AtomTypeData( " C7 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C below C5 when looking at ring from linker to N-O and oriented C=C-H
      C8(   AddEnum( "C8",   AtomTypeData( " C8 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C above C2 when looking at ring from linker to N-O and oriented C=C-H
      C9(   AddEnum( "C9",   AtomTypeData( " C9 ", chemistry::GetElementTypes().e_Carbon,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl C below C2 when looking at ring from linker to N-O and oriented C=C-H
      N1(   AddEnum( "N1",   AtomTypeData( " N1 ", chemistry::GetElementTypes().e_Nitrogen, false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     )))), //!< mtssl N in ring creating N-O nitroxide moiety
      O1(   AddEnum( "O1",   AtomTypeData( " O1 ", chemistry::GetElementTypes().e_Oxygen,   false, true,  coord::CylinderCoordinates(                                          ), coord::CylinderCoordinates(                     ))))  //!< mtssl O in ring creating N-O nitroxide moiety

    {
      // set bond lengths for backbone atom types (taken from rosetta param files)
        N->SetBondLength(   C, 1.329);
        N->SetBondLength(  CA, 1.458);
        N->SetBondLength(   H, 1.03);
       CA->SetBondLength(  HA, 1.09);
       CA->SetBondLength( HA2, 1.09);
       CA->SetBondLength( HA3, 1.09);
       CA->SetBondLength(   N, 1.458);
       CA->SetBondLength(   C, 1.523);
       CA->SetBondLength(  CB, 1.533);
        C->SetBondLength(   N, 1.329);
        C->SetBondLength(   O, 1.231);
        C->SetBondLength(  CA, 1.523);
        O->SetBondLength(   C, 1.231);
       CB->SetBondLength(  CA, 1.533);
        H->SetBondLength(   N, 1.03);
       HA->SetBondLength(  CA, 1.09);
      HA2->SetBondLength(  CA, 1.09);
      HA3->SetBondLength(  CA, 1.09);

      // set obligatory connections for all simple atom types (e.g. connections that will always be made provided that
      // both Atom types are present)
      // The only non-mandatory connection is N <> CD, which only occurs in proline
      // This list was generated (mostly) automatically from the rosetta params files, after verification that all
      // connections other than the proline N <> CD obligatory
      C->SetConnections( storage::Set< AtomType>::Create( CA, O, OXT));
      C->SetDoubleBondConnections( storage::Set< AtomType>::Create( O));
      CA->SetConnections( storage::Set< AtomType>::Create( C, CB, HA, HA2, HA3, N));
      CB->SetConnections( storage::Set< AtomType>::Create( CA, CG, CG1, CG2, HB, HB1, HB2, HB3, OG, OG1, SG));
      CD->SetConnections( storage::Set< AtomType>::Create( CE, CG, HD1, HD2, HD3, NE, NE2, OE1, OE2));
      CD->SetDoubleBondConnections( storage::Set< AtomType>::Create( OE1));
      CD1->SetConnections( storage::Set< AtomType>::Create( CE1, CG, CG1, HD1, HD11, HD12, HD13, NE1));
      CD1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CG));
      CD2->SetConnections( storage::Set< AtomType>::Create( CE2, CE3, CG, HD2, HD21, HD22, HD23, NE2));
      CD2->SetDoubleBondConnections( storage::Set< AtomType>::Create( CE2)); // HIS, TRP, PHE, TYR
      CE->SetConnections( storage::Set< AtomType>::Create( CD, HE1, HE2, HE3, NZ, SD, SE, C3));
      CE1->SetConnections( storage::Set< AtomType>::Create( CD1, CZ, HE1, ND1, NE2));
      CE1->SetDoubleBondConnections( storage::Set< AtomType>::Create( ND1, CD1, CZ));
      CE2->SetConnections( storage::Set< AtomType>::Create( CD2, CZ, CZ2, HE2, NE1));
      CE2->SetDoubleBondConnections( storage::Set< AtomType>::Create( CZ, CD2)); // PHE, TYR
      CE3->SetConnections( storage::Set< AtomType>::Create( CD2, CZ3, HE3));
      CE3->SetDoubleBondConnections( storage::Set< AtomType>::Create( CZ3));
      CG->SetConnections( storage::Set< AtomType>::Create( CB, CD, CD1, CD2, HG, HG1, HG2, HG3, ND1, ND2, OD1, OD2, SD, SE));
      CG->SetDoubleBondConnections( storage::Set< AtomType>::Create( OD1, CD1));
      CG1->SetConnections( storage::Set< AtomType>::Create( CB, CD1, HG11, HG12, HG13));
      CG2->SetConnections( storage::Set< AtomType>::Create( CB, HG21, HG22, HG23));
      CH2->SetConnections( storage::Set< AtomType>::Create( CZ2, CZ3, HH2));
      CH2->SetDoubleBondConnections( storage::Set< AtomType>::Create( CZ2));
      CZ->SetConnections( storage::Set< AtomType>::Create( CE1, CE2, HZ, NE, NH1, NH2, OH));
      CZ->SetDoubleBondConnections( storage::Set< AtomType>::Create( NH1, CE1)); // ARG / PHE, TYR
      CZ2->SetConnections( storage::Set< AtomType>::Create( CE2, CH2, HZ2));
      CZ2->SetDoubleBondConnections( storage::Set< AtomType>::Create( CH2));
      CZ3->SetConnections( storage::Set< AtomType>::Create( CE3, CH2, HZ3));
      CZ3->SetDoubleBondConnections( storage::Set< AtomType>::Create( CE3));
      H->SetConnections( storage::Set< AtomType>::Create( N));
      HA->SetConnections( storage::Set< AtomType>::Create( CA));
      HA2->SetConnections( storage::Set< AtomType>::Create( CA));
      HA3->SetConnections( storage::Set< AtomType>::Create( CA));
      HB->SetConnections( storage::Set< AtomType>::Create( CB));
      HB1->SetConnections( storage::Set< AtomType>::Create( CB));
      HB2->SetConnections( storage::Set< AtomType>::Create( CB));
      HB3->SetConnections( storage::Set< AtomType>::Create( CB));

      // Additional double bonds that are AA-dependent:
      // CG = CD2, not CD1 (Histidine only)

      // Connection ND1 <> HD1 is only present on HIS at pH < 6, but evidently it was decided that pH < 6 is the most
      // important since HIS was given an HD1. The reasons behind this decision should be investigated and documented
      HD1->SetConnections( storage::Set< AtomType>::Create( CD, CD1, ND1));
      HD11->SetConnections( storage::Set< AtomType>::Create( CD1));
      HD12->SetConnections( storage::Set< AtomType>::Create( CD1));
      HD13->SetConnections( storage::Set< AtomType>::Create( CD1));
      HD2->SetConnections( storage::Set< AtomType>::Create( CD, CD2, OD2));
      HD21->SetConnections( storage::Set< AtomType>::Create( CD2, ND2));
      HD22->SetConnections( storage::Set< AtomType>::Create( CD2, ND2));
      HD23->SetConnections( storage::Set< AtomType>::Create( CD2));
      HD3->SetConnections( storage::Set< AtomType>::Create( CD));
      HE->SetConnections( storage::Set< AtomType>::Create( NE));
      HE1->SetConnections( storage::Set< AtomType>::Create( CE, CE1, NE1));
      HE2->SetConnections( storage::Set< AtomType>::Create( CE, CE2, NE2, OE2));
      HE21->SetConnections( storage::Set< AtomType>::Create( NE2));
      HE22->SetConnections( storage::Set< AtomType>::Create( NE2));
      HE3->SetConnections( storage::Set< AtomType>::Create( CE, CE3));
      HG->SetConnections( storage::Set< AtomType>::Create( CG, OG, SG));
      HG1->SetConnections( storage::Set< AtomType>::Create( CG, OG1));
      HG11->SetConnections( storage::Set< AtomType>::Create( CG1));
      HG12->SetConnections( storage::Set< AtomType>::Create( CG1));
      HG13->SetConnections( storage::Set< AtomType>::Create( CG1));
      HG2->SetConnections( storage::Set< AtomType>::Create( CG));
      HG21->SetConnections( storage::Set< AtomType>::Create( CG2));
      HG22->SetConnections( storage::Set< AtomType>::Create( CG2));
      HG23->SetConnections( storage::Set< AtomType>::Create( CG2));
      HG3->SetConnections( storage::Set< AtomType>::Create( CG));
      HH->SetConnections( storage::Set< AtomType>::Create( OH));
      HH11->SetConnections( storage::Set< AtomType>::Create( NH1));
      HH12->SetConnections( storage::Set< AtomType>::Create( NH1));
      HH2->SetConnections( storage::Set< AtomType>::Create( CH2));
      HH21->SetConnections( storage::Set< AtomType>::Create( NH2));
      HH22->SetConnections( storage::Set< AtomType>::Create( NH2));
      HZ->SetConnections( storage::Set< AtomType>::Create( CZ));
      HZ1->SetConnections( storage::Set< AtomType>::Create( NZ));
      HZ2->SetConnections( storage::Set< AtomType>::Create( CZ2, NZ));
      HZ3->SetConnections( storage::Set< AtomType>::Create( CZ3, NZ));
      N->SetConnections( storage::Set< AtomType>::Create( CA, H, H1, H2, H3));
      ND1->SetConnections( storage::Set< AtomType>::Create( CE1, CG, HD1));
      ND1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CE1));
      ND2->SetConnections( storage::Set< AtomType>::Create( CG, HD21, HD22));
      NE->SetConnections( storage::Set< AtomType>::Create( CD, CZ, HE));
      NE1->SetConnections( storage::Set< AtomType>::Create( CD1, CE2, HE1));
      NE2->SetConnections( storage::Set< AtomType>::Create( CD, CD2, CE1, HE2, HE21, HE22));
      NH1->SetConnections( storage::Set< AtomType>::Create( CZ, HH11, HH12));
      NH1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CZ)); // ARG
      NH2->SetConnections( storage::Set< AtomType>::Create( CZ, HH21, HH22));
      NZ->SetConnections( storage::Set< AtomType>::Create( CE, HZ1, HZ2, HZ3));
      O->SetConnections( storage::Set< AtomType>::Create( C));
      O->SetDoubleBondConnections( storage::Set< AtomType>::Create( C)); // Backbone O
      OD1->SetConnections( storage::Set< AtomType>::Create( CG));
      OD1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CG)); // ASP, ASN
      OD2->SetConnections( storage::Set< AtomType>::Create( CG, HD2));
      OE1->SetConnections( storage::Set< AtomType>::Create( CD));
      OE1->SetDoubleBondConnections( storage::Set< AtomType>::Create( CD)); // GLN, GLU
      OE2->SetConnections( storage::Set< AtomType>::Create( CD, HE2));
      OG->SetConnections( storage::Set< AtomType>::Create( CB, HG));
      OG1->SetConnections( storage::Set< AtomType>::Create( CB, HG1));
      OH->SetConnections( storage::Set< AtomType>::Create( CZ, HH));
      SD->SetConnections( storage::Set< AtomType>::Create( CE, CG, SG));
      SE->SetConnections( storage::Set< AtomType>::Create( CG, CE));
      SG->SetConnections( storage::Set< AtomType>::Create( CB, HG, SD));

      // terminal amine
      H1->SetConnections( storage::Set< AtomType>::Create( N));
      H2->SetConnections( storage::Set< AtomType>::Create( N));
      H3->SetConnections( storage::Set< AtomType>::Create( N));

      // terminal carboxylic acid
      HXT->SetConnections( storage::Set< AtomType>::Create( OXT));
      OXT->SetConnections( storage::Set< AtomType>::Create( C, HXT));

      // methanesulfonothioate - specific
      C2->SetConnections( storage::Set< AtomType>::Create( C3, C8, C9, N1));
      C3->SetConnections( storage::Set< AtomType>::Create( C2, C4, CE));
      C3->SetDoubleBondConnections( storage::Set< AtomType>::Create( C4));
      C4->SetConnections( storage::Set< AtomType>::Create( C3, C5));
      C4->SetDoubleBondConnections( storage::Set< AtomType>::Create( C3));
      C5->SetConnections( storage::Set< AtomType>::Create( C4, C6, C7, N1));
      C6->SetConnections( storage::Set< AtomType>::Create( C5));
      C7->SetConnections( storage::Set< AtomType>::Create( C5));
      C8->SetConnections( storage::Set< AtomType>::Create( C2));
      C9->SetConnections( storage::Set< AtomType>::Create( C2));
      N1->SetConnections( storage::Set< AtomType>::Create( C2, C5, O1));
    }

    //! @brief conversion to a string from a Subset
    //! @param SUBSET the subset to get a string for
    //! @return a string representing that subset
    const std::string &AtomTypes::GetSubsetName( const Subset &SUBSET)
    {
      static const std::string s_descriptors[] =
      {
        "All",
        "Backbone",
        "Side-chain"
      };
      return s_descriptors[ size_t( SUBSET)];
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return command line flag for defining the atoms used in the calculation
    //! @return command line flag for defining the the atoms used in the calculation
    util::ShPtr< command::FlagInterface> &AtomTypes::GetFlagAtomTypes()
    {
      // initialize atom types
      storage::Vector< std::string> atom_types
      (
        storage::Vector< std::string>::Create
        (
          GetSubsetName( e_All),
          GetSubsetName( e_Backbone),
          GetSubsetName( e_SideChain)
        )
      );
      atom_types.Append( GetAtomTypes().GetEnumStrings());

      // initialize flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "atoms", "list of atoms to be included in quality calculation, CA by default",
           command::Parameter
          (
            "atom",
            "any atom from the list (backbone and side chain atoms)",
            command::ParameterCheckAllowed( atom_types),
            GetAtomTypes().CA.GetName()
          ),
          0,
          GetAtomTypes().GetEnumCount()
        )
      );

      // end
      return s_flag;
    }

    //! @brief get the set of atoms defined by the flag
    //! @return the set of atoms defined by the flag
    storage::Set< AtomType> AtomTypes::GetCommandLineAtoms()
    {
      // check for "All"
      const storage::Set< std::string> params( GetFlagAtomTypes()->GetObjectSet< std::string>());
      if( params.Contains( GetSubsetName( e_All)))
      {
        // return all atom types
        return storage::Set< AtomType>( GetAtomTypes().Begin(), GetAtomTypes().End());
      }
      else if( params.Contains( GetSubsetName( e_Backbone)))
      {
        // return backbone heavy atom types
        return GetAtomTypes().GetBackBoneAtomTypes();
      }
      else if( params.Contains( GetSubsetName( e_SideChain)))
      {
        // return side chain heavy atom types
        return GetAtomTypes().GetSideChainAtomTypes();
      }

      // get specified atom types
      storage::Set< AtomType> atom_types( GetFlagAtomTypes()->GetObjectSet< AtomType>());

      // add CA if empty
      if( atom_types.IsEmpty())
      {
        atom_types.Insert( GetAtomTypes().CA);
      }

      // end
      return atom_types;
    }

    //! @brief return StorageVector of Backbone Atom Types
    //! @return StorageVector of Backbone Atom Types
    const storage::Set< AtomType> &AtomTypes::GetBackBoneAtomTypes() const
    {
      // initialize static set of AtomTypes from first enum until CB
      static const storage::Set< AtomType> s_backbone_atom_types( Begin(), CB.GetIterator());

      // return
      return s_backbone_atom_types;
    }

    //! @brief return StorageVector of backbone atom names
    //! @return StorageVector of backbone atom names
    const storage::Vector< std::string> &AtomTypes::GetBackBoneAtomNames() const
    {
      // initialize vector of strings for storing atom names
      static const storage::Vector< std::string> s_bb_atom_names( Begin(), CB.GetIterator());

      // return
      return s_bb_atom_names;
    }

    //! @brief return StorageVector of side chain Atom Types
    //! @return StorageVector of side chain Atom Types
    const storage::Set< AtomType> &AtomTypes::GetSideChainAtomTypes() const
    {
      // initialize static set of AtomTypes from CB to SG
      static const storage::Set< AtomType> s_side_chain_atom_types( CB.GetIterator(), H.GetIterator());

      // return
      return s_side_chain_atom_types;
    }

    //! @brief return set of AtomTypes composed of CA and CB
    //! @return set of AtomTypes composed of CA and CB
    const storage::Set< AtomType> &AtomTypes::GetCACB() const
    {
      // initialize a set of AtomTypes composed of CA and CB
      static const storage::Set< AtomType> s_ca_cb( CA, CB);

      // return
      return s_ca_cb;
    }

    //! @brief determines and returns the atom type from the provided pdb atom name
    //! @param PDB_ATOM_NAME AtomName for the pdb atom of interest
    //! @return the atom type from the provided pdb atom name
    AtomType AtomTypes::TypeFromPDBAtomName( const std::string &PDB_ATOM_NAME) const
    {
      static std::map< std::string, AtomType> s_map;
      if( s_map.empty())
      {
        for
        (
          const_iterator type_itr( Begin()), type_itr_end( End());
          type_itr != type_itr_end;
          ++type_itr
        )
        {
          s_map[ type_itr->GetName()] = *type_itr;
        }
      }

      auto ite( s_map.find( util::TrimString( PDB_ATOM_NAME)));
      return ite == s_map.end() ? e_Undefined : ite->second;
    }

    //! @brief access to set of possible first side chain atom
    //! @return set of first sied chain atom types, CB and HA2 (for glycin)
    const storage::Set< AtomType> &AtomTypes::GetFirstSidechainAtomTypes() const
    {
      // set of first side chain atom types
      static const storage::Set< AtomType> s_first_sidechain_atom_types( CB, HA2);

      // end
      return s_first_sidechain_atom_types;
    }

    //! @brief access to map of additional atom types for terminal residues and their PDB atom name defined by PDB
    //! @return Map containing the atoms, only found in the terminal amine and carboxylic acid with their PDB atom NAME (e.g. H1 -> 1H)
    const storage::Map< AtomType, std::string> &AtomTypes::GetTerminalExtraAtomTypes() const
    {
      static const storage::Map< AtomType, std::string> s_terminal_atom_type( TerminalExtraAtomTypes());
      return s_terminal_atom_type;
    }

    //! @brief terminal atomtype from atom name
    //! @param ATOM_NAME name of atom (e.g. H1 or 1H, OXT ...)
    //! @return the terminal atom type for that atom name - undefined if there is non
    AtomType AtomTypes::GetTerminalAtomTypeFromName( const std::string &ATOM_NAME) const
    {
      // find atom type by iterating over mapping
      for
      (
        storage::Map< AtomType, std::string>::const_iterator
          itr( GetTerminalExtraAtomTypes().Begin()), itr_end( GetTerminalExtraAtomTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->second.find( ATOM_NAME) != std::string::npos)
        {
          return itr->first;
        }
      }

      // name can not be found in map
      return TypeFromPDBAtomName( ATOM_NAME);
    }

    //! @brief access to map of additional atom types for terminal residues and their PDB atom name defined by PDB
    //! @return Map containing the atoms, only found in the terminal amine and carboxylic acid with their PDB atom NAME (e.g. H1 -> 1H)
    storage::Map< AtomType, std::string> AtomTypes::TerminalExtraAtomTypes() const
    {
      // set of terminal atom types that only appear at the N and C terminus
      storage::Map< AtomType, std::string> terminal_atom_types;
      terminal_atom_types[  H1] = "1H  ";
      terminal_atom_types[  H2] = "2H  ";
      terminal_atom_types[  H3] = "3H  ";
      terminal_atom_types[ HXT] = " HXT";
      terminal_atom_types[ OXT] = " OXT";

      // end
      return terminal_atom_types;
    }

    //! @brief access to the only instance of AtomTypes
    //! @return reference to only instance of AtomTypes
    const AtomTypes &GetAtomTypes()
    {
      return AtomTypes::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::AtomTypeData, biol::AtomTypes>;

  } // namespace util
} // namespace bcl
