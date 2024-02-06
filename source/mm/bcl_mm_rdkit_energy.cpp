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
#include "mm/bcl_mm_rdkit_energy.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"
#include "math/bcl_math_limits.h"
#include "mm/bcl_mm_rdkit_force_field_utils.h"

// external includes - sorted alphabetically
#include "ForceField/ForceField.h"
#include "GraphMol/RWMol.h"

namespace bcl
{
  namespace mm
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    RDKitEnergy::RDKitEnergy() :
      m_ForceFieldEnum( e_UFF),
      m_NonbondedThreshold( 100.0),
      m_IgnoreInterFragmentInteractions( true)
    {
      m_ForceFieldString = GetRdkitForceFieldsName( e_UFF);
    }

    //! full constructor
    RDKitEnergy::RDKitEnergy
    (
      const RdkitForceFieldsEnum &VARIANT,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS
    ) :
      m_ForceFieldEnum( VARIANT),
      m_NonbondedThreshold( NON_BONDED_THRESHOLD),
      m_IgnoreInterFragmentInteractions( IGNORE_INTER_FRAG_INTERACTIONS)
    {
      m_ForceFieldString = GetRdkitForceFieldsName( VARIANT);
    }

    //! virtual copy constructor
    RDKitEnergy *RDKitEnergy::Clone() const
    {
      return new RDKitEnergy( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RDKitEnergy::GetAlias() const
    {
      static const std::string s_name( "RDKitEnergy");  // TODO: change name based on force field
      return s_name;
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RDKitEnergy::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the force field variant
    RdkitForceFieldsEnum RDKitEnergy::GetForceFieldEnum() const
    {
      return m_ForceFieldEnum;
    }

    //! @brief returns the force field variant as a string
    std::string RDKitEnergy::GetForceFieldString() const
    {
      return m_ForceFieldString;
    }

    //! @brief returns the non-bonded threshold
    double RDKitEnergy::GetNonbondedThreshold() const
    {
      return m_NonbondedThreshold;
    }

    //! @brief returns whether to ignore fragment interactions
    bool RDKitEnergy::GetIgnoreInterFragmentInteractions() const
    {
      return m_IgnoreInterFragmentInteractions;
    }

  ///////////////////
  //   operations  //
  ///////////////////

    //! @brief sets the force field variant
    void RDKitEnergy::SetForceFieldFromEnum( const RdkitForceFieldsEnum &VARIANT)
    {
      m_ForceFieldEnum = VARIANT;
      m_ForceFieldString = GetRdkitForceFieldsName( m_ForceFieldEnum);
    }

    //! @brief sets the force field variant from a string
    void RDKitEnergy::SetForceFieldFromString( const std::string &VARIANT)
    {
      m_ForceFieldEnum = RdkitForceFieldUtils::GetRdkitForceFieldsEnum( VARIANT);
      m_ForceFieldString = GetRdkitForceFieldsName( m_ForceFieldEnum);
    }

    //! @brief sets the force field string
    void RDKitEnergy::SetForceFieldString( const std::string &VARIANT)
    {
      m_ForceFieldString = VARIANT;
      m_ForceFieldEnum = RdkitForceFieldUtils::GetRdkitForceFieldsEnum( m_ForceFieldString);
    }

    //! @brief sets the non-bonded threshold
    void RDKitEnergy::SetNonbondedThreshold( const double THRESHOLD)
    {
      m_NonbondedThreshold = THRESHOLD;
    }

    //! @brief sets whether to ignore fragment interactions
    void RDKitEnergy::SetIgnoreInterFragmentInteractions( const bool IGNORE_INTER_FRAGMENT_INTERACTIONS)
    {
      m_IgnoreInterFragmentInteractions = IGNORE_INTER_FRAGMENT_INTERACTIONS;
    }

    //! @brief splits the molecule according to GetComponentVertices
    //! @param MOLECULE the molecule for which the energy will be computed
    //! @return the energy of MOLECULE
    double RDKitEnergy::CalculateEnergy( const chemistry::FragmentComplete &MOLECULE) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // build force field
      ::ForceFields::ForceField *ff
      (
        RdkitForceFieldUtils::ConstructForceField
        (
          *rdkit_mol,
          GetForceFieldString(),
          m_NonbondedThreshold,
          m_IgnoreInterFragmentInteractions,
          true
        )
      );

      // compute energy
      return ff->calcEnergy();
    }

    //! @brief Computes the force field potential energy of a molecule
    //! @param MOLECULE the molecule for which the energy will be computed
    //! @param VARIANT whether to use UFF, MMFF94, or MMFF94s
    //! @param NON_BONDED_THRESHOLD the threshold to be used in adding non-bonded terms to the force field.
    //! @param IGNORE_INTER_FRAG_INTERACTIONS If true, nonbonded terms will not be added between fragments
    //! @return the energy of MOLECULE
    double RDKitEnergy::CalculateEnergy
    (
      const chemistry::FragmentComplete &MOLECULE,
      const std::string &VARIANT,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS
    )
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // build force field
      ::ForceFields::ForceField *ff
      (
        RdkitForceFieldUtils::ConstructForceField
        (
          *rdkit_mol,
          VARIANT,
          NON_BONDED_THRESHOLD,
          IGNORE_INTER_FRAG_INTERACTIONS,
          true
        )
      );

      // compute energy
      return ff->calcEnergy();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RDKitEnergy::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes the potential energy of a molecule with molecular mechanics force field.");
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
        "force_field",
        "Force field to use; options are 'UFF', 'MMFF94' or 'MMFF94s'.",
        io::Serialization::GetAgent( &m_ForceFieldString),
        "UFF"
      );
      return parameters;
    }

  } // namespace mm
} // namespace bcl
