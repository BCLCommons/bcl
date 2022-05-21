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
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_flag_static.h"

// external includes - sorted alphabetically
//#define BCL_PROFILE_ConformationComparisonInterface
#ifdef BCL_PROFILE_ConformationComparisonInterface
#include "util/bcl_util_stopwatch.h"
#endif

namespace bcl
{
  namespace chemistry
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief prepare the class for comparing conformations in the given ensemble
    //! @param ENSEMBLE the ensemble to prepare to compare
    void ConformationComparisonInterface::PrepareEnsemble( const FragmentEnsemble &ENSEMBLE) const
    {
      for( FragmentEnsemble::const_iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End()); itr != itr_end; ++itr)
      {
        Prepare( *itr);
      }
    }

    //! @brief get a bool that indicates whether to ignore atom / bond types in ConformationsAreComparable, below
    //! @note this is necessary when comparing parts of molecules, rather than complete ligands
    bool &ConformationComparisonInterface::GetIgnoreAtomAndBondTypesWhenDeterminingComparability()
    {
      static bool s_ignore_atom_bond_types( false);
      return s_ignore_atom_bond_types;
    }

    //! @brief Get a flag to turn on/off strict atom type checking
    const util::ShPtr< command::FlagInterface> &ConformationComparisonInterface::GetDisableStrictAtomBondTypeCheckingFlag()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "disable_strict_atom_bond_type_checks",
          "disables strict atom/bond type checking. Use only if you are confident that the atom / bond order is the same but "
          " that the atom/bond types have been messed up by external software (esp conversion to/from pdbs)"
        )
      );
      return s_flag;
    }

    //! @brief determine whether the conformations represent identical, aligned, constitutions
    //! @param MOLECULE_A, MOLECULE_B the conformations to check for being identical and aligned
    //! @return true iff the conformations represent identical, aligned, constitutions
    bool ConformationComparisonInterface::ConformationsAreComparable
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    )
    {
      // Determine if molecules are equal in # atoms
      if( MOLECULE_A.GetNumberAtoms() != MOLECULE_B.GetNumberAtoms())
      {
        // nope, so return an undefined double
        BCL_MessageCrt( "Non-equal # of atoms! " + util::Format()( MOLECULE_A.GetNumberAtoms()) + " vs " + util::Format()( MOLECULE_B.GetNumberAtoms()));
        return false;
      }

      // Determine if molecules are equal in # bonds
      if( MOLECULE_A.GetNumberBonds() != MOLECULE_B.GetNumberBonds())
      {
        // nope, so return an undefined double
        BCL_MessageCrt( "Non-equal # of bonds!");
        return false;
      }

#ifdef BCL_PROFILE_ConformationComparisonInterface
      static util::Stopwatch s_are_comparable( "ConformationsAreComparable", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_are_comparable.Stop();
      s_are_comparable.Start();
#endif
      const bool strict_checking
      (
        !GetIgnoreAtomAndBondTypesWhenDeterminingComparability()
        && !GetDisableStrictAtomBondTypeCheckingFlag()->GetFlag()
      );

      // Determine if molecules have the same atom types, in the same order
      if( strict_checking)
      {
        for( auto itra( MOLECULE_A.GetAtomsIterator()), itrb( MOLECULE_B.GetAtomsIterator()); itra.NotAtEnd(); ++itra, ++itrb)
        {
          if( itra->GetAtomType() != itrb->GetAtomType())
          {
            // nope, so return an undefined double
            BCL_MessageCrt( "Atom types differed\n" + MOLECULE_A.GetAtomTypesString() + "\n" + MOLECULE_B.GetAtomTypesString());
            #ifdef BCL_PROFILE_ConformationComparisonInterface
            s_are_comparable.Stop();
            #endif
            return false;
          }
        }
      }
      else
      {
        for( auto itra( MOLECULE_A.GetAtomsIterator()), itrb( MOLECULE_B.GetAtomsIterator()); itra.NotAtEnd(); ++itra, ++itrb)
        {
          if( itra->GetElementType() != itrb->GetElementType())
          {
            // nope, so return an undefined double
            BCL_MessageCrt
            (
              "Element types differed\n" + MOLECULE_A.GetElementTypesString()
              + "\n" + MOLECULE_B.GetElementTypesString()
            );
            #ifdef BCL_PROFILE_ConformationComparisonInterface
            s_are_comparable.Stop();
            #endif
            return false;
          }
        }
      }

      // Determine if molecules have equivalent connectivity
      const ConfigurationalBondTypeData::DataEnum bond_comparison_level
      (
        !strict_checking
        ? ConfigurationalBondTypeData::e_Identity
        : ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
      );
      for( auto itra( MOLECULE_A.GetAtomsIterator()), itrb( MOLECULE_B.GetAtomsIterator()); itra.NotAtEnd(); ++itra, ++itrb)
      {
        if( itra->GetBonds().GetSize() != itrb->GetBonds().GetSize())
        {
          BCL_MessageCrt( "Connectivity differed");
          #ifdef BCL_PROFILE_ConformationComparisonInterface
          s_are_comparable.Stop();
          #endif
          return false;
        }
        for
        (
          auto itr_bnd_a( itra->GetBonds().Begin()), itr_bnd_a_end( itra->GetBonds().End()), itr_bnd_b( itrb->GetBonds().Begin());
          itr_bnd_a != itr_bnd_a_end;
          ++itr_bnd_a, ++itr_bnd_b
        )
        {
          if
          (
            MOLECULE_A.GetAtomIndex( itr_bnd_a->GetTargetAtom()) != MOLECULE_B.GetAtomIndex( itr_bnd_b->GetTargetAtom())
            || itr_bnd_a->GetBondType()->GetBondData( bond_comparison_level) != itr_bnd_b->GetBondType()->GetBondData( bond_comparison_level)
          )
          {
            BCL_MessageCrt( "Connectivity differed");
            #ifdef BCL_PROFILE_ConformationComparisonInterface
            s_are_comparable.Stop();
            #endif
            return false;
          }
        }
      }

      #ifdef BCL_PROFILE_ConformationComparisonInterface
      s_are_comparable.Stop();
      #endif
      // the molecules are aligned and represent identical constitutions
      return true;
    }

  } // namespace chemistry
} // namespace bcl
