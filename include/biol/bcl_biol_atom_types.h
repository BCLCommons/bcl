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

#ifndef BCL_BIOL_ATOM_TYPES_H_
#define BCL_BIOL_ATOM_TYPES_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_atom_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomTypes
    //! @brief Enumerator class to be used for accessing Atom type information in biology related classes.
    //! @details This enumerator class has an individual enum for each Atom type that are relevant to proteins, thus
    //! it is more limited compared to chemistry::AtomTypes.
    //!
    //! @see @link example_biol_atom_types.cpp @endlink
    //! @author karakam, woetzen
    //! @date 11/25/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomTypes :
      public util::Enumerate< AtomTypeData, AtomTypes>
    {
      friend class util::Enumerate< AtomTypeData, AtomTypes>;
    public:

    //////////
    // data //
    //////////

      // declare all atom types
      const AtomType N;    //!< Nitrogen from the peptide bond
      const AtomType CA;   //!< Carbon alpha backbone
      const AtomType C;    //!< Carbon from the carboxyl group
      const AtomType O;    //!< Oxygen from the carboxyl group
      const AtomType CB;   //!< Carbon beta first side chain atom
      const AtomType CG;   //!< Carbon gamma - second side chain atom
      const AtomType CG1;  //!< Carbon gamma - second side chain atom 1 for two second side chain atoms
      const AtomType CG2;  //!< Carbon gamma - second side chain atom 2 for two second side chain atoms
      const AtomType CD;   //!< Carbon delta - third side chain atom
      const AtomType CD1;  //!< Carbon delta - third side chain atom 1 for two third side chain atoms
      const AtomType CD2;  //!< Carbon delta - third side chain atom 2 for two third side chain atoms
      const AtomType CE;   //!< Carbon epsilon - fourth side chain atom
      const AtomType CE1;  //!< Carbon epsilon - fourth side chain atom 1 for two or three fourth side chain atoms
      const AtomType CE2;  //!< Carbon epsilon - fourth side chain atom 2 for two or three fourth side chain atoms
      const AtomType CE3;  //!< Carbon epsilon - fourth side chain atom 3 for two or three fourth side chain atoms
      const AtomType CZ;   //!< Carbon zeta - fifth side chain atom
      const AtomType CZ2;  //!< Carbon zeta - fifth side chain atom 1 for two fifth side chain atoms
      const AtomType CZ3;  //!< Carbon zeta - fifth side chain atom 2 for two fifth side chain atoms
      const AtomType CH2;  //!< Carbon eta - sixth side chain aom as in TRP
      const AtomType ND1;  //!< Nitrogen delta - third position
      const AtomType ND2;  //!< Nitrogen delta - third position
      const AtomType NE;   //!< Nitrogen epsilon - fourth position nitrogen as in ARG
      const AtomType NE1;  //!< Nitrogen epsilon - fourth position nitrogen as in TRP
      const AtomType NE2;  //!< Nitrogen epsilon - fourth position nitrogen as in GLN or HIS
      const AtomType NZ;   //!< Nitrogen zeta - fifth side chain atom as in LYS
      const AtomType NH1;  //!< Nitrogen eta 1 - sixth side chain atom as in ARG
      const AtomType NH2;  //!< Nitrogen eta 2 - sixth side chain atom as in ARG
      const AtomType OD1;  //!< Oxygen delta - third side chain atom as in ASP or ASN
      const AtomType OD2;  //!< Oxygen delta - third side chain atom as in ASP
      const AtomType OG;   //!< Oxygen gamma - second side chain atom as in SER
      const AtomType OG1;  //!< Oxygen gamma - second side chain atom as in THR
      const AtomType OE1;  //!< Oxygen epsilon - fourth side chain atom as in GLN
      const AtomType OE2;  //!< Oxygen epsilon - fourth side chain atom as in GLU
      const AtomType OH;   //!< Oxygen eta - sixth side chain atom as in TYR
      const AtomType SD;   //!< Suflur on delta carbon - like in methionine
      const AtomType SE;   //!< Selenium on delta carbon - like in selenomethionine
      const AtomType SG;   //!< Sulfur on gamma carbon - like in cystein
      const AtomType H;    //!< Hydrogen on the backbone nitrogen atom
      const AtomType HA;   //!< Hydrogen on CA alpha Carbon
      const AtomType HA2;  //!< Hydrogen on CA for GLY in CB position
      const AtomType HA3;  //!< Hydrogen on CA for GLY in HA position
      const AtomType HB;   //!< Hydrogen on CB for THR, ILE or VAL
      const AtomType HB1;  //!< Hydrogen on CB for ALA
      const AtomType HB2;  //!< Hydrogen 2 on CB
      const AtomType HB3;  //!< Hydrogen 3 on CB
      const AtomType HG;   //!< Hydrogen an CG as in LEU, CYS, SER
      const AtomType HG1;  //!< Hydrogen 1 an CG as in THR
      const AtomType HG2;  //!< Hydrogen 2 on CG
      const AtomType HG3;  //!< Hydrogen 3 on CG
      const AtomType HG11; //!< Hydrogen 1 on CG1 as in VAL
      const AtomType HG12; //!< Hydrogen 2 on CG2 as in ILE, VAL
      const AtomType HG13; //!< Hydrogen 3 on CG1 as in ILE, VAL
      const AtomType HG21; //!< Hydrogen 1 on CG2 as in ILE, VAL, THR
      const AtomType HG22; //!< Hydrogen 2 on CG2 as in ILE, VAL, THR
      const AtomType HG23; //!< Hydrogen 3 on CG2 as in ILE, VAL, THR
      const AtomType HD1;  //!< Hydrogen 1 on CD
      const AtomType HD2;  //!< Hydrogen 2 on CD
      const AtomType HD3;  //!< Hydrogen 3 on CD
      const AtomType HD11; //!< Hydrogen 1 on CD1
      const AtomType HD12; //!< Hydrogen 2 on CD1
      const AtomType HD13; //!< Hydrogen 3 on CD1
      const AtomType HD21; //!< Hydrogen 1 on CD2
      const AtomType HD22; //!< Hydrogen 2 on CD2
      const AtomType HD23; //!< Hydrogen 3 on CD2
      const AtomType HE;   //!< Hydrogen on CE as in ARG
      const AtomType HE1;  //!< Hydrogen 1 on CE
      const AtomType HE2;  //!< Hydrogen 2 on CE
      const AtomType HE3;  //!< Hydrogen 3 on CE
      const AtomType HE21; //!< Hydrogen 1 on CE2 as in GLN
      const AtomType HE22; //!< Hydrogen 2 on CE2 as in GLN
      const AtomType HZ;   //!< Hydrogen on CZ as in PHE
      const AtomType HZ1;  //!< Hydrogen on Z1 as in LYS (on Nitrogen NZ)
      const AtomType HZ2;  //!< Hydrogen 2 on Z as in LYS (on NZ), TRP (on CZ2)
      const AtomType HZ3;  //!< Hydrogen 3 on Z as in LYS (on NZ), TRP (on CZ3)
      const AtomType HH;   //!< Hydrogen on CH as in TYR
      const AtomType HH2;  //!< Hydrogen on CH2 as in TRP
      const AtomType HH11; //!< Hydrogen 1 on NH1 as in ARG
      const AtomType HH12; //!< Hydrogen 2 on NH1 as in ARG
      const AtomType HH21; //!< Hydrogen 1 on NH2 as in ARG
      const AtomType HH22; //!< Hydrogen 2 on NH2 as in ARG

      // terminal amine
      const AtomType H1;   //!< Hydrogen 1 on backbone terminal NH3
      const AtomType H2;   //!< Hydrogen 1 on backbone terminal NH3
      const AtomType H3;   //!< Hydrogen 1 on backbone terminal NH3

      // terminal carboxylic acid
      const AtomType HXT;  //!< Hydrogen leaving hydrogen, connected to C if there is no peptide bond
      const AtomType OXT;  //!< leaving oxygen, connected to C if there is no peptide bond

      const AtomType C2;   //!< mttsl C next to N-O nearest to linker
      const AtomType C3;   //!< mtssl C in "pentane" ring connected to linker
      const AtomType C4;   //!< mtssl C in pentane ring connected to H and connected to C3 by double bond
      const AtomType C5;   //!< mtssl C next to N-O nearest to linker on side with C=C-H and opposite C2
      const AtomType C6;   //!< mtssl C above C5 when looking at ring from linker to N-O and oriented C=C-H
      const AtomType C7;   //!< mtssl C below C5 when looking at ring from linker to N-O and oriented C=C-H
      const AtomType C8;   //!< mtssl C above C2 when looking at ring from linker to N-O and oriented C=C-H
      const AtomType C9;   //!< mtssl C below C2 when looking at ring from linker to N-O and oriented C=C-H
      const AtomType N1;   //!< mtssl N in ring creating N-O nitroxide moiety
      const AtomType O1;   //!< mtssl O in ring creating N-O nitroxide moiety

      // enumerator for subsets
      enum Subset
      {
        e_All,
        e_Backbone,
        e_SideChain
      };

      //! @brief conversion to a string from a Subset
      //! @param SUBSET the subset to get a string for
      //! @return a string representing that subset
      static const std::string &GetSubsetName( const Subset &SUBSET);

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all AtomTypes
      AtomTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return command line flag for defining the atoms used in the calculation
      //! @return command line flag for defining the the atoms used in the calculation
      static util::ShPtr< command::FlagInterface> &GetFlagAtomTypes();

      //! @brief get the set of atoms defined by the flag
      //! @return the set of atoms defined by the flag
      static storage::Set< AtomType> GetCommandLineAtoms();

      //! @brief return StorageVector of Backbone Atom Types
      //! @return StorageVector of Backbone Atom Types
      const storage::Set< AtomType> &GetBackBoneAtomTypes() const;

      //! @brief return StorageVector of backbone atom names
      //! @return StorageVector of backbone atom names
      const storage::Vector< std::string> &GetBackBoneAtomNames() const;

      //! @brief return StorageVector of side chain Atom Types
      //! @return StorageVector of side chain Atom Types
      const storage::Set< AtomType> &GetSideChainAtomTypes() const;

      //! @brief return set of AtomTypes composed of CA and CB
      //! @return set of AtomTypes composed of CA and CB
      const storage::Set< AtomType> &GetCACB() const;

      //! @brief determines and returns the atom type from the provided pdb atom name
      //! @param PDB_ATOM_NAME AtomName for the pdb atom of interest
      //! @return the atom type from the provided pdb atom name
      AtomType TypeFromPDBAtomName( const std::string &PDB_ATOM_NAME) const;

      //! @brief access to set of possible first side chain atom
      //! @return set of first sied chain atom types, CB and HA2 (for glycine)
      const storage::Set< AtomType> &GetFirstSidechainAtomTypes() const;

      //! @brief access to map of additional atom types for terminal residues and their PDB atom name defined by PDB
      //! @return Map containing the atoms, only found in the terminal amine and carboxylic acid with their PDB atom NAME (e.g. H1 -> 1H)
      const storage::Map< AtomType, std::string> &GetTerminalExtraAtomTypes() const;

      //! @brief terminal atomtype from atom name
      //! @param ATOM_NAME name of atom (e.g. H1 or 1H, OXT ...)
      //! @return the terminal atom type for that atom name - undefined if there is non
      AtomType GetTerminalAtomTypeFromName( const std::string &ATOM_NAME) const;

    private:

      //! @brief create map of additional atom types for terminal residues and their PDB atom name defined by PDB
      //! @return Map containing the atoms, only found in the terminal amine and carboxylic acid with their PDB atom NAME (e.g. H1 -> 1H)
      storage::Map< AtomType, std::string> TerminalExtraAtomTypes() const;

    }; // class AtomTypes

    //! @brief access to the only instance of AtomTypes
    //! @return reference to only instance of AtomTypes
    BCL_API
    const AtomTypes &GetAtomTypes();

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< biol::AtomTypeData, biol::AtomTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_BIOL_ATOM_TYPES_H_
