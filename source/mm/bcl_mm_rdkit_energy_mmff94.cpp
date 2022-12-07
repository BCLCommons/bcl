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
#include "mm/bcl_mm_rdkit_energy_mmff94.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"

// external includes - sorted alphabetically
#include "ForceField/ForceField.h"
#include "GraphMol/ForceFieldHelpers/MMFF/MMFF.h"
#include "GraphMol/RWMol.h"

namespace bcl
{
  namespace mm
  {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////


    //! default constructor
    RDKitEnergyMMFF94::RDKitEnergyMMFF94() :
      m_MMFFVariant( e_MMFF94s),
      m_NonbondedThreshold( 100.0),
      m_IgnoreInterFragmentInteractions( true)
    {
      if( !m_MMFFVariantString.empty())
      {
        SetMMFFVariantFromString( m_MMFFVariantString);
      }
    }

    //! full constructor
    RDKitEnergyMMFF94::RDKitEnergyMMFF94
    (
      const MMFFVariant &VARIANT,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS
    ) :
      m_MMFFVariant( VARIANT),
      m_NonbondedThreshold( NON_BONDED_THRESHOLD),
      m_IgnoreInterFragmentInteractions( IGNORE_INTER_FRAG_INTERACTIONS)
    {
      if( !m_MMFFVariantString.empty())
      {
        SetMMFFVariantFromString( m_MMFFVariantString);
      }
    }

    //! virtual copy constructor
    RDKitEnergyMMFF94 *RDKitEnergyMMFF94::Clone() const
    {
      return new RDKitEnergyMMFF94( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RDKitEnergyMMFF94::GetAlias() const
    {
      static const std::string s_name( m_MMFFVariant == e_MMFF94 ? "Energy_MMFF94" : "Energy_MMFF94s");
      return s_name;
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RDKitEnergyMMFF94::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the MMFF94 variant
    RDKitEnergyMMFF94::MMFFVariant RDKitEnergyMMFF94::GetMMFFVariant() const
    {
      return m_MMFFVariant;
    }

    //! @brief returns the MMFF94 variant as a string
    std::string RDKitEnergyMMFF94::GetMMFFVariantAsString() const
    {
      return m_MMFFVariant == e_MMFF94 ? "MMFF94" : "MMFF94s";
    }

    //! @brief returns the non-bonded threshold
    double RDKitEnergyMMFF94::GetNonbondedThreshold() const
    {
      return m_NonbondedThreshold;
    }

    //! @brief returns whether to ignore fragment interactions
    bool RDKitEnergyMMFF94::GetIgnoreInterFragmentInteractions() const
    {
      return m_IgnoreInterFragmentInteractions;
    }

    ///////////////////
    //   operations  //
    ///////////////////

    //! @brief sets the MMFF94 variant
    void RDKitEnergyMMFF94::SetMMFFVariant( const MMFFVariant &VARIANT)
    {
      m_MMFFVariant = VARIANT;
    }

    //! @brief sets the MMFF94 variant from a string
    void RDKitEnergyMMFF94::SetMMFFVariantFromString( const std::string &VARIANT)
    {
      if( VARIANT == "MMFF94")
      {
        m_MMFFVariant = e_MMFF94;
      }
      else if( VARIANT == "MMFF94s")
      {
        m_MMFFVariant = e_MMFF94s;
      }
      else
      {
        BCL_Exit("Invalid string; choices are either 'MMFF94' or 'MMFF94s'.", 1);
      }
    }

    //! @brief sets the non-bonded threshold
    void RDKitEnergyMMFF94::SetNonbondedThreshold( const double THRESHOLD)
    {
      m_NonbondedThreshold = THRESHOLD;
    }

    //! @brief sets whether to ignore fragment interactions
    void RDKitEnergyMMFF94::SetIgnoreInterFragmentInteractions( const bool IGNORE_INTER_FRAGMENT_INTERACTIONS)
    {
      m_IgnoreInterFragmentInteractions = IGNORE_INTER_FRAGMENT_INTERACTIONS;
    }

    //! @brief splits the molecule according to GetComponentVertices
    //! @param MOLECULE the molecule for which the energy will be computed
    //! @return the energy of MOLECULE
    double RDKitEnergyMMFF94::CalculateEnergy( const chemistry::FragmentComplete &MOLECULE) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, GetMMFFVariantAsString());
      if( !mmff_mol_properties.isValid())
      {
        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
        return util::GetUndefinedDouble();
      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( *rdkit_mol, m_NonbondedThreshold, -1, m_IgnoreInterFragmentInteractions);
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
    double RDKitEnergyMMFF94::CalculateEnergy
    (
      const chemistry::FragmentComplete &MOLECULE,
      const std::string &MMFF_VARIANT,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS
    )
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, MMFF_VARIANT);
      if( !mmff_mol_properties.isValid())
      {
        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
        return util::GetUndefinedDouble();
      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( *rdkit_mol, NON_BONDED_THRESHOLD, -1, IGNORE_INTER_FRAG_INTERACTIONS);
      ff->initialize();

      // compute energy
      return ff->calcEnergy();
    }

    //////////////////////
    // helper functions //
    //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RDKitEnergyMMFF94::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes the potential energy of a molecule with the MMFF94(s) force field.");
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
      parameters.AddInitializer
      (
        "mmff_variant",
        "MMFF force field variant to use; options are 'MMFF94' or 'MMFF94s'.",
        io::Serialization::GetAgent( &m_MMFFVariantString),
        "MMFF94s"
      );
      return parameters;
    }

  } // namespace mm
} // namespace bcl
