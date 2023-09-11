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
#include "math/bcl_math_limits.h"
#include "mm/bcl_mm_rdkit_energy_uff.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"

// external includes - sorted alphabetically
#include "ForceField/ForceField.h"
#include "GraphMol/ForceFieldHelpers/UFF/Builder.h"
#include "GraphMol/ForceFieldHelpers/UFF/AtomTyper.h"
#include "GraphMol/RWMol.h"

namespace bcl
{
  namespace mm
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    RDKitEnergyUFF::RDKitEnergyUFF() :
      m_NonbondedThreshold( 100.0),
      m_IgnoreInterFragmentInteractions( true)
    {
    }

    //! full constructor
    RDKitEnergyUFF::RDKitEnergyUFF
    (
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS
    ) :
      m_NonbondedThreshold( NON_BONDED_THRESHOLD),
      m_IgnoreInterFragmentInteractions( IGNORE_INTER_FRAG_INTERACTIONS)
    {
    }

    //! virtual copy constructor
    RDKitEnergyUFF *RDKitEnergyUFF::Clone() const
    {
      return new RDKitEnergyUFF( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RDKitEnergyUFF::GetAlias() const
    {
      static const std::string s_name( "Energy_UFF");
      return s_name;
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RDKitEnergyUFF::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the non-bonded threshold
    double RDKitEnergyUFF::GetNonbondedThreshold() const
    {
      return m_NonbondedThreshold;
    }

    //! @brief returns whether to ignore fragment interactions
    bool RDKitEnergyUFF::GetIgnoreInterFragmentInteractions() const
    {
      return m_IgnoreInterFragmentInteractions;
    }

  ///////////////////
  //   operations  //
  ///////////////////

    //! @brief sets the non-bonded threshold
    void RDKitEnergyUFF::SetNonbondedThreshold( const double THRESHOLD)
    {
      m_NonbondedThreshold = THRESHOLD;
    }

    //! @brief sets whether to ignore fragment interactions
    void RDKitEnergyUFF::SetIgnoreInterFragmentInteractions( const bool IGNORE_INTER_FRAGMENT_INTERACTIONS)
    {
      m_IgnoreInterFragmentInteractions = IGNORE_INTER_FRAGMENT_INTERACTIONS;
    }

    //! @brief splits the molecule according to GetComponentVertices
    //! @param MOLECULE the molecule for which the energy will be computed
    //! @return the energy of MOLECULE
    double RDKitEnergyUFF::CalculateEnergy( const chemistry::FragmentComplete &MOLECULE) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
//      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, GetMMFFVariantAsString());
//      if( !mmff_mol_properties.isValid())
//      {
//        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
//        return util::GetUndefinedDouble();
//      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::UFF::constructForceField( *rdkit_mol, m_NonbondedThreshold, -1, m_IgnoreInterFragmentInteractions);
      ff->initialize();

      // compute energy
      return ff->calcEnergy();
    }

    //! @brief Computes the MMFF94 potential energy of a molecule
    //! @param MOLECULE the molecule for which the energy will be computed
    //! @param MMFF_VARIANT whether to use MMFF94 or MMFF94s
    //! @param NON_BONDED_THRESHOLD the threshold to be used in adding non-bonded terms to the force field.
    //! @param IGNORE_INTER_FRAG_INTERACTIONS If true, nonbonded terms will not be added between fragments
    //! @return the energy of MOLECULE
    double RDKitEnergyUFF::CalculateEnergy
    (
      const chemistry::FragmentComplete &MOLECULE,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS
    )
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
//      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, MMFF_VARIANT);
//      if( !mmff_mol_properties.isValid())
//      {
//        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
//        return util::GetUndefinedDouble();
//      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::UFF::constructForceField( *rdkit_mol, NON_BONDED_THRESHOLD, -1, IGNORE_INTER_FRAG_INTERACTIONS);
      ff->initialize();

      // compute energy
      return ff->calcEnergy();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RDKitEnergyUFF::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes the potential energy of a molecule with the UFF force field.");
      parameters.AddInitializer
      (
        "non_bonded_threshold",
        "The threshold to be used in adding non-bonded terms to the force field. "
        "Any non-bonded contact whose current distance is greater than nonBondedThresh * the minimum value "
        "for that contact will not be included.",
        io::Serialization::GetAgent( &m_NonbondedThreshold),
        "100.0"
      );
      parameters.AddInitializer
      (
        "ignore_inter_fragment_interactions",
        "If true, nonbonded terms will not be added between fragments",
        io::Serialization::GetAgent( &m_IgnoreInterFragmentInteractions),
        "true"
      );
      return parameters;
    }

  } // namespace mm
} // namespace bcl
