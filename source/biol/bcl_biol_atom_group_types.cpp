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
#include "biol/bcl_biol_atom_group_types.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_atom_group_type_data.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief construct all AtomTypes
    AtomGroupTypes::AtomGroupTypes() :
//                                         Volume   Radius  Bound_Hydrogen,  scattering_length H20, scattering_length D20
//                                         (A^3)    (A)                          (10^-12 cm)             (10^-12 cm)
      H(    AddEnum( "H",   AtomGroupTypeData( GetAtomTypes().H,    5.15, 1.07,  0, -0.3742,  0.6671 ))), //!< Hydrogen
      C(    AddEnum( "C",   AtomGroupTypeData( GetAtomTypes().C,   16.44, 1.58,  0,  0.6651,  0.6651 ))), //!< Carbon
      CH(   AddEnum( "CH",  AtomGroupTypeData( GetAtomTypes().C,   21.59, 1.73,  1,  0.2909,  0.2909 ))), //!< Carbon 1 bound hydrogen
      CH2(  AddEnum( "CH2", AtomGroupTypeData( GetAtomTypes().C,   26.74, 1.85,  2, -0.0833, -0.0833 ))), //!< Carbon 2 bound hydrogens
      CH3(  AddEnum( "CH3", AtomGroupTypeData( GetAtomTypes().C,   31.89, 1.97,  3, -0.4575, -0.4575 ))), //!< Carbon 3 bound hydrogens
      N(    AddEnum( "N",   AtomGroupTypeData( GetAtomTypes().N,    2.49, 0.84,  0,  0.9400,  0.9400 ))), //!< Nitrogen
      NH(   AddEnum( "NH",  AtomGroupTypeData( GetAtomTypes().N,    7.64, 1.22,  1,  0.5658,  1.6071 ))), //!< Nitrogen 1 bound hydrogen
      NH2(  AddEnum( "NH2", AtomGroupTypeData( GetAtomTypes().N,   12.79, 1.45,  2,  0.1916,  2.2742 ))), //!< Nitrogen 2 bound hydrogens
      NH3(  AddEnum( "NH3", AtomGroupTypeData( GetAtomTypes().N,   17.94, 1.62,  3, -0.1826,  2.9413 ))), //!< Nitrogen 3 bound hydrogens
      O(    AddEnum( "O",   AtomGroupTypeData( GetAtomTypes().O,    9.13, 1.30,  0,  0.5804,  0.5804 ))), //!< Oxygen
      OH(   AddEnum( "OH",  AtomGroupTypeData( GetAtomTypes().O,   14.28, 1.50,  1,  0.2062,  1.2475 ))), //!< Oxygen 1 bound hydrogen
      S(    AddEnum( "S",   AtomGroupTypeData( GetAtomTypes().SG,  19.86, 1.68,  0,  0.2847,  0.2847 ))), //!< Sulfur
      SH(   AddEnum( "SH",  AtomGroupTypeData( GetAtomTypes().SG,  25.10, 1.81,  1, -0.0895,  1.2365 ))), //!< Sulfur 1 bound hydrogen
      SE(   AddEnum( "SE",  AtomGroupTypeData( GetAtomTypes().SE,  28.73, 1.90,  0,  0.7971,  0.7971 )))  //!< Selenium
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomGroupTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the only instance of AtomTypes
    //! @return reference to only instance of AtomTypes
    const AtomGroupTypes &GetAtomGroupTypes()
    {
      return AtomGroupTypes::GetEnums();
    }

    //! @brief map to hold Solvent Volume and Radius values for each heavy atom in each amino acid
    //! @param AA The amino acid types
    //! @param ATOMTYPE The atomtype
    //! @return Group Types
    //! The values are taken from Crysol - a program to evaluate x-ray solution scattering of biological
    //! macromolecules from Atomic Coordinates.  D. Svergun et all 1995
    const AtomGroupType &AtomGroupTypes::GetType( const AAType &AA, const AtomType &ATOMTYPE)
    {
      static storage::Map< AAType, storage::Map< AtomType, AtomGroupType> > map;
      if( map.IsEmpty())
      {
        // Compose Map
        map[ GetAATypes().ALA][ GetAtomTypes().N  ] = GetAtomGroupTypes().N;
        map[ GetAATypes().ALA][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ALA][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ALA][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ALA][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().ARG][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ARG][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ARG][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ARG][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ARG][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ARG][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ARG][ GetAtomTypes().CD ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ARG][ GetAtomTypes().NE ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ARG][ GetAtomTypes().CZ ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ARG][ GetAtomTypes().NH1] = GetAtomGroupTypes().NH2;
        map[ GetAATypes().ARG][ GetAtomTypes().NH2] = GetAtomGroupTypes().NH2;
        map[ GetAATypes().ASN][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ASN][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ASN][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ASN][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ASN][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ASN][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ASN][ GetAtomTypes().OD1] = GetAtomGroupTypes().O;
        map[ GetAATypes().ASN][ GetAtomTypes().ND2] = GetAtomGroupTypes().NH2;
        map[ GetAATypes().ASP][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ASP][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ASP][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ASP][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ASP][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ASP][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ASP][ GetAtomTypes().OD1] = GetAtomGroupTypes().O;
        map[ GetAATypes().ASP][ GetAtomTypes().OD2] = GetAtomGroupTypes().OH;
        map[ GetAATypes().CYS][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().CYS][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().CYS][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().CYS][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().CYS][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().CYS][ GetAtomTypes().SG ] = GetAtomGroupTypes().SH;
        map[ GetAATypes().GLN][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().GLN][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().GLN][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLN][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().GLN][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLN][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLN][ GetAtomTypes().CD ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLN][ GetAtomTypes().OE1] = GetAtomGroupTypes().O;
        map[ GetAATypes().GLN][ GetAtomTypes().NE2] = GetAtomGroupTypes().NH2;
        map[ GetAATypes().GLU][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().GLU][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().GLU][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLU][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().GLU][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLU][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLU][ GetAtomTypes().CD ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLU][ GetAtomTypes().OE1] = GetAtomGroupTypes().O;
        map[ GetAATypes().GLU][ GetAtomTypes().OE2] = GetAtomGroupTypes().OH;
        map[ GetAATypes().GLY][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().GLY][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().GLY][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().GLY][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().HIS][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().HIS][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().HIS][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().HIS][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().HIS][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().HIS][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().HIS][ GetAtomTypes().ND1] = GetAtomGroupTypes().NH;
        map[ GetAATypes().HIS][ GetAtomTypes().CD2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().HIS][ GetAtomTypes().CE1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().HIS][ GetAtomTypes().NE2] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ILE][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().ILE][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ILE][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().ILE][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().ILE][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().ILE][ GetAtomTypes().CG1] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().ILE][ GetAtomTypes().CG2] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().ILE][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().LEU][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().LEU][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().LEU][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().LEU][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().LEU][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LEU][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().LEU][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().LEU][ GetAtomTypes().CD2] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().LYS][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().LYS][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().LYS][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().LYS][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().LYS][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LYS][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LYS][ GetAtomTypes().CD ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LYS][ GetAtomTypes().CE ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().LYS][ GetAtomTypes().NZ ] = GetAtomGroupTypes().NH3;
        map[ GetAATypes().MET][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().MET][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().MET][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().MET][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().MET][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().MET][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().MET][ GetAtomTypes().SD ] = GetAtomGroupTypes().S;
        map[ GetAATypes().MET][ GetAtomTypes().CE ] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().MSE][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().MSE][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().MSE][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().MSE][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().MSE][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().MSE][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().MSE][ GetAtomTypes().SE ] = GetAtomGroupTypes().SE;
        map[ GetAATypes().MSE][ GetAtomTypes().CE ] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().PHE][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().PHE][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().PHE][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().PHE][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().PHE][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().PHE][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().CD2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().CE1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().CE2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PHE][ GetAtomTypes().CZ ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PRO][ GetAtomTypes().N  ] = GetAtomGroupTypes().N;
        map[ GetAATypes().PRO][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().PRO][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().PRO][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().PRO][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().PRO][ GetAtomTypes().CG ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().PRO][ GetAtomTypes().CD ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().SER][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().SER][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().SER][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().SER][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().SER][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().SER][ GetAtomTypes().OG ] = GetAtomGroupTypes().OH;
        map[ GetAATypes().THR][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().THR][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().THR][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().THR][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().THR][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().THR][ GetAtomTypes().OG1] = GetAtomGroupTypes().OH;
        map[ GetAATypes().THR][ GetAtomTypes().CG2] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().TRP][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().TRP][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TRP][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().TRP][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().TRP][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TRP][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().CD2] = GetAtomGroupTypes().C;
        map[ GetAATypes().TRP][ GetAtomTypes().NE1] = GetAtomGroupTypes().NH;
        map[ GetAATypes().TRP][ GetAtomTypes().CE2] = GetAtomGroupTypes().C;
        map[ GetAATypes().TRP][ GetAtomTypes().CE3] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().CZ2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().CZ3] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TRP][ GetAtomTypes().CH2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().TYR][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TYR][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().TYR][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH2;
        map[ GetAATypes().TYR][ GetAtomTypes().CG ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TYR][ GetAtomTypes().CD1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().CD2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().CE1] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().CE2] = GetAtomGroupTypes().CH;
        map[ GetAATypes().TYR][ GetAtomTypes().CZ ] = GetAtomGroupTypes().C;
        map[ GetAATypes().TYR][ GetAtomTypes().OH ] = GetAtomGroupTypes().OH;
        map[ GetAATypes().VAL][ GetAtomTypes().N  ] = GetAtomGroupTypes().NH;
        map[ GetAATypes().VAL][ GetAtomTypes().CA ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().VAL][ GetAtomTypes().C  ] = GetAtomGroupTypes().C;
        map[ GetAATypes().VAL][ GetAtomTypes().O  ] = GetAtomGroupTypes().O;
        map[ GetAATypes().VAL][ GetAtomTypes().CB ] = GetAtomGroupTypes().CH;
        map[ GetAATypes().VAL][ GetAtomTypes().CG1] = GetAtomGroupTypes().CH3;
        map[ GetAATypes().VAL][ GetAtomTypes().CG2] = GetAtomGroupTypes().CH3;
      }
      // return group type
      return map.GetValue( AA).GetValue( ATOMTYPE);
    }

    //! @brief map to hold Solvent Volume and Radius values for each heavy atom in each amino acid
    //! @param AA The amino acid types
    //! @param ATOMTYPE The atomtype
    //! @return Group Types
    //! The values are taken from Crysol - a program to evaluate x-ray solution scattering of biological
    //! macromolecules from Atomic Coordinates.  D. Svergun et all 1995
    const std::string &AtomGroupTypes::GetTypeString( const AAType &AA, const AtomType &ATOMTYPE)
    {
      static storage::Map< AAType, storage::Map< AtomType, std::string> > map;
      if( map.IsEmpty())
      {
        // Compose Map
        map[ GetAATypes().ALA][ GetAtomTypes().N  ] = "N";
        map[ GetAATypes().ALA][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ALA][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ALA][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ALA][ GetAtomTypes().CB ] = "CH3";
        map[ GetAATypes().ARG][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().ARG][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ARG][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ARG][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ARG][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().ARG][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().ARG][ GetAtomTypes().CD ] = "CH2";
        map[ GetAATypes().ARG][ GetAtomTypes().NE ] = "NH";
        map[ GetAATypes().ARG][ GetAtomTypes().CZ ] = "C";
        map[ GetAATypes().ARG][ GetAtomTypes().NH1] = "NH2";
        map[ GetAATypes().ARG][ GetAtomTypes().NH2] = "NH2";
        map[ GetAATypes().ASN][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().ASN][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ASN][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ASN][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ASN][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().ASN][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().ASN][ GetAtomTypes().OD1] = "O";
        map[ GetAATypes().ASN][ GetAtomTypes().ND2] = "NH2";
        map[ GetAATypes().ASP][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().ASP][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ASP][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ASP][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ASP][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().ASP][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().ASP][ GetAtomTypes().OD1] = "O";
        map[ GetAATypes().ASP][ GetAtomTypes().OD2] = "OH";
        map[ GetAATypes().CYS][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().CYS][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().CYS][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().CYS][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().CYS][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().CYS][ GetAtomTypes().SG ] = "SH";
        map[ GetAATypes().GLN][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().GLN][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().GLN][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().GLN][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().GLN][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().GLN][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().GLN][ GetAtomTypes().CD ] = "C";
        map[ GetAATypes().GLN][ GetAtomTypes().OE1] = "O";
        map[ GetAATypes().GLN][ GetAtomTypes().NE2] = "NH2";
        map[ GetAATypes().GLU][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().GLU][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().GLU][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().GLU][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().GLU][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().GLU][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().GLU][ GetAtomTypes().CD ] = "C";
        map[ GetAATypes().GLU][ GetAtomTypes().OE1] = "O";
        map[ GetAATypes().GLU][ GetAtomTypes().OE2] = "OH";
        map[ GetAATypes().GLY][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().GLY][ GetAtomTypes().CA ] = "CH2";
        map[ GetAATypes().GLY][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().GLY][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().HIS][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().HIS][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().HIS][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().HIS][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().HIS][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().HIS][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().HIS][ GetAtomTypes().ND1] = "NH";
        map[ GetAATypes().HIS][ GetAtomTypes().CD2] = "CH";
        map[ GetAATypes().HIS][ GetAtomTypes().CE1] = "CH";
        map[ GetAATypes().HIS][ GetAtomTypes().NE2] = "NH";
        map[ GetAATypes().ILE][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().ILE][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().ILE][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().ILE][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().ILE][ GetAtomTypes().CB ] = "CH";
        map[ GetAATypes().ILE][ GetAtomTypes().CG1] = "CH2";
        map[ GetAATypes().ILE][ GetAtomTypes().CG2] = "CH3";
        map[ GetAATypes().ILE][ GetAtomTypes().CD1] = "CH3";
        map[ GetAATypes().LEU][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().LEU][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().LEU][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().LEU][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().LEU][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().LEU][ GetAtomTypes().CG ] = "CH";
        map[ GetAATypes().LEU][ GetAtomTypes().CD1] = "CH3";
        map[ GetAATypes().LEU][ GetAtomTypes().CD2] = "CH3";
        map[ GetAATypes().LYS][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().LYS][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().LYS][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().LYS][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().LYS][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().LYS][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().LYS][ GetAtomTypes().CD ] = "CH2";
        map[ GetAATypes().LYS][ GetAtomTypes().CE ] = "CH2";
        map[ GetAATypes().LYS][ GetAtomTypes().NZ ] = "NH3";
        map[ GetAATypes().MET][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().MET][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().MET][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().MET][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().MET][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().MET][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().MET][ GetAtomTypes().SD ] = "S";
        map[ GetAATypes().MET][ GetAtomTypes().CE ] = "CH3";
        map[ GetAATypes().MSE][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().MSE][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().MSE][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().MSE][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().MSE][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().MSE][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().MSE][ GetAtomTypes().SE ] = "SE";
        map[ GetAATypes().MSE][ GetAtomTypes().CE ] = "CH3";
        map[ GetAATypes().PHE][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().PHE][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().PHE][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().PHE][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().PHE][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().PHE][ GetAtomTypes().CD1] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().CD2] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().CE1] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().CE2] = "CH";
        map[ GetAATypes().PHE][ GetAtomTypes().CZ ] = "CH";
        map[ GetAATypes().PRO][ GetAtomTypes().N  ] = "N";
        map[ GetAATypes().PRO][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().PRO][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().PRO][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().PRO][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().PRO][ GetAtomTypes().CG ] = "CH2";
        map[ GetAATypes().PRO][ GetAtomTypes().CD ] = "CH2";
        map[ GetAATypes().SER][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().SER][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().SER][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().SER][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().SER][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().SER][ GetAtomTypes().OG ] = "OH";
        map[ GetAATypes().THR][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().THR][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().THR][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().THR][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().THR][ GetAtomTypes().CB ] = "CH";
        map[ GetAATypes().THR][ GetAtomTypes().OG1] = "OH";
        map[ GetAATypes().THR][ GetAtomTypes().CG2] = "CH3";
        map[ GetAATypes().TRP][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().TRP][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().TRP][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().TRP][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().TRP][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().TRP][ GetAtomTypes().CD1] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().CD2] = "C";
        map[ GetAATypes().TRP][ GetAtomTypes().NE1] = "NH";
        map[ GetAATypes().TRP][ GetAtomTypes().CE2] = "C";
        map[ GetAATypes().TRP][ GetAtomTypes().CE3] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().CZ2] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().CZ3] = "CH";
        map[ GetAATypes().TRP][ GetAtomTypes().CH2] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().TYR][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().TYR][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().TYR][ GetAtomTypes().CB ] = "CH2";
        map[ GetAATypes().TYR][ GetAtomTypes().CG ] = "C";
        map[ GetAATypes().TYR][ GetAtomTypes().CD1] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().CD2] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().CE1] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().CE2] = "CH";
        map[ GetAATypes().TYR][ GetAtomTypes().CZ ] = "C";
        map[ GetAATypes().TYR][ GetAtomTypes().OH ] = "OH";
        map[ GetAATypes().VAL][ GetAtomTypes().N  ] = "NH";
        map[ GetAATypes().VAL][ GetAtomTypes().CA ] = "CH";
        map[ GetAATypes().VAL][ GetAtomTypes().C  ] = "C";
        map[ GetAATypes().VAL][ GetAtomTypes().O  ] = "O";
        map[ GetAATypes().VAL][ GetAtomTypes().CB ] = "CH";
        map[ GetAATypes().VAL][ GetAtomTypes().CG1] = "CH3";
        map[ GetAATypes().VAL][ GetAtomTypes().CG2] = "CH3";
      }
      // return group type
      return map.GetValue( AA).GetValue( ATOMTYPE);
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::AtomGroupTypeData, biol::AtomGroupTypes>;

  } // namespace util
} // namespace bcl
