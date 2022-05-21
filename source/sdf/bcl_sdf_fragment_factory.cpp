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
#include "sdf/bcl_sdf_fragment_factory.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"

namespace bcl
{
  namespace sdf
  {

    //! @brief read complete fragment from an mdl file
    //! @param MDL_HANDLER handler that has the information
    //! @return a fragment complete
    chemistry::FragmentComplete
      FragmentFactory::MakeFragment
      (
        const CTabHandler &HANDLER,
        const HydrogenHandlingPref &H_PREF,
        const NeutralizationPref &NEUTRALIZATION_PREF,
        const std::string &MOLECULE_ID
      )
    {
      // create atoms
      chemistry::AtomVector< chemistry::AtomComplete> atoms( HANDLER.GetAtomInfo(), HANDLER.GetBondInfo());

      if( !HANDLER.WereBondTypesRead() || !HANDLER.WereAtomTypesRead())
      {
        // determine atom and bond types
        chemistry::AtomsCompleteStandardizer standardizer( atoms, MOLECULE_ID, false);
      }
      else
      {
        chemistry::AtomsCompleteStandardizer::SetConjugationOfBondTypes( atoms);
      }
      chemistry::AtomsCompleteStandardizer::TryNeutralize( atoms, NEUTRALIZATION_PREF);

      if( !HANDLER.WasDoubleBondIsometryRead())
      {
        // add bond isometry information
        chemistry::BondIsometryHandler::AddIsometryInformation( atoms, true);
      }
      if( !HANDLER.WasChiralityRead())
      {
        // add stereocenter information
        chemistry::StereocentersHandler::AddChiralityFromConformation( atoms);
      }

      // for each atom, put the preference as to whether to add H into a vector
      storage::Vector< size_t> add_hydrogen;
      add_hydrogen.AllocateMemory( HANDLER.GetAtomInfo().GetSize());
      // iterate over all atoms
      for
      (
        storage::Vector< AtomInfo>::const_iterator
          itr( HANDLER.GetAtomInfo().Begin()), itr_end( HANDLER.GetAtomInfo().End());
        itr != itr_end;
        ++itr
      )
      {
        add_hydrogen.PushBack( itr->CanAddH());
      }
      chemistry::HydrogensHandler::HandleHydrogenPref( atoms, H_PREF, add_hydrogen);

      // Give this molecule a blank name for the time being
      return chemistry::FragmentComplete( atoms, MOLECULE_ID);
    }

    //! @brief read complete fragment from an mdl file
    //! @param MDL_HANDLER handler that has the information
    //! @return a fragment complete
    chemistry::FragmentComplete
      FragmentFactory::MakeFragment
      (
        const MolfileHandler &HANDLER,
        const HydrogenHandlingPref &H_PREF,
        const NeutralizationPref &NEUTRALIZATION_PREF,
        const std::string &MOLECULE_ID
      )
    {
      chemistry::FragmentComplete new_fragment( MakeFragment( static_cast< CTabHandler>( HANDLER), H_PREF, NEUTRALIZATION_PREF, MOLECULE_ID));
      new_fragment.SetName( HANDLER.GetDescription());
      return new_fragment;
    }

    //! @brief read complete fragment from an mdl file
    //! @param MDL_HANDLER handler that has the information
    //! @return a fragment complete
    chemistry::FragmentComplete
      FragmentFactory::MakeFragment
      (
        const MdlHandler &HANDLER,
        const HydrogenHandlingPref &H_PREF,
        const NeutralizationPref &NEUTRALIZATION_PREF,
        const std::string &MOLECULE_ID
      )
    {
      chemistry::FragmentComplete new_fragment( MakeFragment( HANDLER.GetMolfile(), H_PREF, NEUTRALIZATION_PREF, MOLECULE_ID));
      const storage::Map< std::string, std::string> &prop_map( HANDLER.GetMiscProperties());
      for
      (
        storage::Map< std::string, std::string>::const_iterator itr_map( prop_map.Begin()), itr_map_end( prop_map.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        new_fragment.StoreProperty( itr_map->first, itr_map->second);
      }
      new_fragment.CacheNumeric( prop_map);
      return new_fragment;
    }

    //! @brief read fragment conformation from MDL file
    //! @param MDL_HANDLER handler that has the conformation information
    //! @return a fragment conformation shared
    chemistry::FragmentConformationShared
      FragmentFactory::MakeConformation( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF)
    {
      return chemistry::FragmentConformationShared( MakeFragment( HANDLER, H_PREF));
    }

    //! @brief read fragment configuration from MDL file
    //! @param MDL_HANDLER handler that has the configuration information
    //! @return a fragment configuration shared
    chemistry::FragmentConfigurationShared
      FragmentFactory::MakeConfiguration( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF)
    {
      return chemistry::FragmentConfigurationShared( MakeFragment( HANDLER, H_PREF));
    }

    //! @brief read fragment constitution from MDL file
    //! @param MDL_HANDLER handler that has the constitution information
    //! @return a fragment constitution shared
    chemistry::FragmentConstitutionShared
      FragmentFactory::MakeConstitution( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF)
    {
      return chemistry::FragmentConstitutionShared( MakeFragment( HANDLER, H_PREF));
    }

  } // namespace sdf
} // namespace bcl

