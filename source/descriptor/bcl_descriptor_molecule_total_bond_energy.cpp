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
#include "descriptor/bcl_descriptor_molecule_total_bond_energy.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "command/bcl_command_command_state.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeTotalBondEnergy::MoleculeTotalBondEnergy() :
      m_BondEnergiesFilename( "chembl_bdb_constitution.clean.expanded_bond_types.STATS.element_bond_counts.energies.txt.gz")
    {
    }

    //! @brief specify bond energies file constructor
    MoleculeTotalBondEnergy::MoleculeTotalBondEnergy( const std::string &BOND_ENERGIES_FILENAME) :
      m_BondEnergiesFilename( BOND_ENERGIES_FILENAME)
    {
//      // Read in the atom environments file directly
//      io::IFStream( file);
//      io::File::MustOpenIFStream
//      (
//        file,
//        m_BondEnergiesFilename
//      );
//      io::Serialize::Read( m_BondEnergies, file);
//      io::File::CloseClearFStream( file);
//
//      // Make sure the bond propensities file is not empty
//      BCL_Assert( !m_BondEnergies.IsEmpty(), "Statistical bond energies file is empty!");
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeTotalBondEnergy
    MoleculeTotalBondEnergy *MoleculeTotalBondEnergy::Clone() const
    {
      return new MoleculeTotalBondEnergy( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeTotalBondEnergy::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeTotalBondEnergy::GetAlias() const
    {
      static const std::string s_name( "MoleculeTotalBondEnergy");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeTotalBondEnergy::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      // initialize total bond energy term
      math::RunningMinMax< float> minmax_bond_energy;
      math::RunningAverageSD< float> avesd_bond_energy;

      // we deal with fewer bond cases than possible
      chemistry::ConfigurationalBondType canonical_bonds[ 11] =
      {
          chemistry::GetConfigurationalBondTypes().e_Undefined,
          chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedSingleBond,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedTripleBond,
          chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBondInRing,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedSingleBondInRing,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedTripleBondInRing,
          chemistry::GetConfigurationalBondTypes().e_AmideSingleBond,
          chemistry::GetConfigurationalBondTypes().e_AromaticBond
      };
      // go over each atom in molecule
      for
      (
          auto itr_atoms( this->GetCurrentObject()->GetIterator());
          itr_atoms.NotAtEnd();
          ++itr_atoms
      )
      {
        // go over each bond from each atom
        const chemistry::AtomConformationalInterface &atom( *itr_atoms);
        for( auto bond_itr( atom.GetBonds().Begin()), bond_itr_end( atom.GetBonds().End()); bond_itr != bond_itr_end; ++bond_itr)
        {
          // make a key of element_type-element_type-bond_type to match with statistical bond energy
          storage::Triplet< chemistry::ElementType, chemistry::ElementType, chemistry::ConfigurationalBondType> key
          (
            atom.GetElementType(),
            bond_itr->GetTargetAtom().GetElementType(),
            canonical_bonds[ bond_itr->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)]
          );
          // increment total bond energy
          auto value( m_BondEnergies.Find( key));
          if( value == m_BondEnergies.End())
          {
            minmax_bond_energy += m_BondEnergies.ReverseBegin()->second;
            avesd_bond_energy += m_BondEnergies.ReverseBegin()->second;
          }
          else
          {
            minmax_bond_energy += value->second;
            avesd_bond_energy += value->second;
          }
        }
      }
      STORAGE(0) = avesd_bond_energy.GetAverage() * avesd_bond_energy.GetWeight();
      STORAGE(1) = avesd_bond_energy.GetAverage();
      STORAGE(2) = avesd_bond_energy.GetSampleStandardDeviation();
      STORAGE(3) = minmax_bond_energy.GetMax();
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeTotalBondEnergy::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // Read in bond energies file
      io::IFStream( file);
      io::File::MustOpenIFStream
      (
        file,
        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") +
        "atom_environments/statistical_bond_energies/" + ( m_BondEnergiesFilename)
      );
      io::Serialize::Read( m_BondEnergies, file);
      io::File::CloseClearFStream( file);

      // Make sure the bond energies file has shit in it
      BCL_Assert( !m_BondEnergies.IsEmpty(), "Bond energies file is empty! Exiting...");
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeTotalBondEnergy::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the sum of bond energies derived from the statistical bond potential"
      );
//      parameters.AddInitializer
//      (
//        "bond_energy_file",
//        "the file containing the pre-generated statistical bond energies",
//        io::Serialization::GetAgent( &m_BondEnergiesFilename),
//        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "atom_environments/statistical_bond_energies/"
//        + "chembl_bdb_constitution.clean.expanded_bond_types.STATS.element_bond_counts.energies.txt.gz"
//      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
