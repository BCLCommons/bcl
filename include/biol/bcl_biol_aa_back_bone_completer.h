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

#ifndef BCL_BIOL_AA_BACK_BONE_COMPLETER_H_
#define BCL_BIOL_AA_BACK_BONE_COMPLETER_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_atom_types.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AABackBoneCompleter
    //! @brief Adds missing backbone hydrogens or oxygens to the given model or sequence
    //! @details This class works on a given ProteinModel or AASequence by making a AAComplete copy of the given argument
    //! while inserting any missing backbone hydrogens and oxygens depending on the way the class itself is initialized.
    //!
    //! @see @link example_biol_aa_back_bone_completer.cpp @endlink
    //! @author karakam
    //! @date Mar 24, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AABackBoneCompleter :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! whether or not to add missing amide hydrogens
      bool m_AddAmideHydrogens;

      //! whether or not to add missing carbonyl oxygens
      bool m_AddCarbonylOxygens;

      //! whether or not to add missing HA hydrogens
      bool m_AddHAHydrogens;

    public:

      //! the bond angle between the CA->C->N atoms where N is of the residue following the CA and C atom residue
      static const double s_CACNBondAngle;

      //! the bond angle between the O->C->N atoms where N  is of the residue following the O and C atom residue
      static const double s_OCNAngle;

      //! the bond angle between the C->N->CA atoms where C is of the residue preceding the N and CA atom residue
      static const double s_CNCABondAngle;

      //! the bond angle between the C->N->H atoms where C is of the residue preceding the N and H atom residue
      static const double s_CNHBondAngle;

      //! the ideal omega bond angle CA->C->N->CA
      static const double s_OmegaCACNCA;

      //! the ideal omega bond angle CA->C->N->H
      static const double s_OmegaCACNH;

      //! the ideal omega bond angle O->C->N->CA
      static const double s_OmegaOCNCA;

      //! the ideal omega bond angle O->C->N->H
      static const double s_OmegaOCNH;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from two booleans whether or not to add missing amide hydrogen and carbonyl oxygens
      //! @param ADD_AMIDE_HYDROGENS whether or not to add missing amide hydrogens
      //! @param ADD_CARBONYL_OXYGENS whether or not to add missing carbonyl oxygens
      //! @param ADD_HA_HYDROGENS whether or not to add missing HA hydrogens
      AABackBoneCompleter
      (
        const bool ADD_AMIDE_HYDROGENS = true,
        const bool ADD_CARBONYL_OXYGENS = true,
        const bool ADD_HA_HYDROGENS = true
      );

      //! @brief Clone function
      //! @return pointer to new AABackBoneCompleter
      AABackBoneCompleter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get whether or not to add missing amide hydrogens
      //! @return whether or not to add missing amide hydrogens
      bool GetAddAmideHydrogens() const
      {
        return m_AddAmideHydrogens;
      }

      //! @brief get whether or not to add missing carbonyl oxygens
      //! @return whether or not to add missing carbonyl oxygens
      bool GetAddCarbonylOxygens() const
      {
        return m_AddCarbonylOxygens;
      }

      //! @brief get whether or not to add missing HA hydrogens
      //! @return whether or not to add missing HA hydrogens
      bool GetAddHAHydrogens() const
      {
        return m_AddHAHydrogens;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief function to complete backbones for a given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return ProteinModel with backbone completed and of type AAComplete
      //! @return ShPtr to a new ProteinModel with completed backbone
      util::ShPtr< assemble::ProteinModel> CompleteProteinModel( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief function to complete backbones for a given Chain
      //! @param CHAIN Chain of interest
      //! @return Chain with backbone completed and of type AAComplete
      //! @return ShPtr to a new Chain with completed backbone
      util::ShPtr< assemble::Chain> CompleteChain( const assemble::Chain &CHAIN) const;

      //! @brief function to complete backbones for a given AASequence
      //! @param AA_SEQUENCE AASequence of interest
      //! @param SP_PREV_AA amino acid previous to the sequence
      //! @param SP_NEXT_AA amino acid following the sequence
      //! @return AASequence with backbone completed and of type AAComplete
      //! @return ShPtr to a new AASequence with completed backbone
      util::ShPtr< AASequence> CompleteAASequence
      (
        const AASequence &AA_SEQUENCE,
        const util::SiPtr< const AABase> &SP_PREV_AA,
        const util::SiPtr< const AABase> &SP_NEXT_AA
      ) const;

      //! @brief function to generate a missing hydrogen from amino acid and previous C coordinates
      //! @param AMINO_ACID amino acid of interest
      //! @param PREV_C_COORD Coordinates of C of previous amino acid
      //! @param BOND_LENGTH N-H bond length
      //! @return generated hydrogen atom
      static Atom GenerateHydrogen
      (
        const AABase &AMINO_ACID,
        const linal::Vector3D &PREV_C_COORD,
        const double BOND_LENGTH = GetAtomTypes().N->GetBondLength( GetAtomTypes().H)
      );

      //! @brief function to generate a missing oxygen from a amino acid and next N coordinates
      //! @param AMINO_ACID amino acid of interest
      //! @param NEXT_N_COORD Coordinates of N of next amino acid
      //! @param BOND_LENGTH C=O bond length
      //! @return generated oxygen atom
      static Atom GenerateOxygen
      (
        const AABase &AMINO_ACID,
        const linal::Vector3D &NEXT_N_COORD,
        const double BOND_LENGTH = GetAtomTypes().C->GetBondLength( GetAtomTypes().O)
      );

      //! @brief function to generate a missing HA from a amino acid
      //! @param AMINO_ACID amino acid of interest
      //! @param BOND_LENGTH CA-HA bond length
      //! @return generated oxygen atom
      static Atom GenerateHA
      (
        const AABase &AMINO_ACID,
        const double BOND_LENGTH = GetAtomTypes().CA->GetBondLength( GetAtomTypes().HA)
      );

      //! @brief function to generate the N of the neighboring amino acid in the peptide bond
      //! @param AMINO_ACID amino acid of interest, with psi angle defined (using oxgyen)
      //! @param BOND_LENGTH C-N bond length
      //! @return generated nitrogen (c-terminal direction)
      static Atom GenerateN
      (
        const AABase &AMINO_ACID,
        const double BOND_LENGTH = GetAtomTypes().C->GetBondLength( GetAtomTypes().N)
      );

      //! @brief function to generate the C of the neighboring amino acid in the petide bond
      //! @param AMINO_ACID amino acid of interest
      //! @param PHI the phi angle
      //! @param BOND_LENGTH the C-N bond length
      //! @return generated carbon (n-terminal direction)
      static Atom GenerateC
      (
        const AABase &AMINO_ACID,
        const double PHI,
        const double BOND_LENGTH = GetAtomTypes().C->GetBondLength( GetAtomTypes().N)
      );

      //! @brief function to generate C, O and CA on an N-terminal amino acid
      //! @param AMINO_ACID amino acid of interest
      //! @param PHI the phi angle on the AMINO acid
      //! @return map of atoms with C, O and CA atoms, peptide bond to AMINO_ACID
      static storage::Map< AtomType, Atom> GenerateCOCA
      (
        const AABase &AMINO_ACID,
        const double PHI
      );

      //! @brief function to generate H, N and CA on an C-terminal amino acid
      //! @param AMINO_ACID amino acid of interest
      //! @return map of atoms with H, N and CA atoms, peptide bond to AMINO_ACID
      static storage::Map< AtomType, Atom> GenerateHNCA
      (
        const AABase &AMINO_ACID
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class AABackBoneCompleter

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_AA_BACK_BONE_COMPLETER_H_ 
