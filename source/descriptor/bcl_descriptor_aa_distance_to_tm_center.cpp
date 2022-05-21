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
#include "descriptor/bcl_descriptor_aa_distance_to_tm_center.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AADistanceToTMCenter::s_AATMPositionInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AADistanceToTMCenter( false)
      )
    );
    const util::SiPtr< const util::ObjectInterface> AADistanceToTMCenter::s_AATMNormPositionInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AADistanceToTMCenter( true)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param NORM whether to take the membrane norm
    AADistanceToTMCenter::AADistanceToTMCenter( const bool &MEMBRANE_NORM) :
      m_MembraneNorm( MEMBRANE_NORM),
      m_DefaultValue( MEMBRANE_NORM ? 3.0 : 30.0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AADistanceToTMCenter *AADistanceToTMCenter::Clone() const
    {
      return new AADistanceToTMCenter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AADistanceToTMCenter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AADistanceToTMCenter::GetAlias() const
    {
      static const std::string s_name( "AA_DistanceFromMembraneCenter"), s_normname( "AA_NormDistanceFromMembraneCenter");
      return m_MembraneNorm ? s_normname : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AADistanceToTMCenter::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A the element of interest
    //! @param STORAGE storage for the descriptor
    void AADistanceToTMCenter::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if the sequence maps are empty, this indicates that this object is now operating over a new sequence,
      // so it is necessary to reload the files
      if( !m_HaveCheckedForMembrane)
      {
        // Get protein model
        util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

        // Check to make sure there is only a single chain in the given PDB file
        BCL_Assert
        (
          sp_protein_model->GetChains().GetSize() == 1,
          "More than one chain in the given protein model: "
           + util::Format()( sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile))
        );

        m_HaveMembrane
          = sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane).IsDefined();
        if( m_HaveMembrane)
        {
          // check that the membrane itself is defined as well
          util::SiPtr< const biol::Membrane> sp_membrane
          (
            sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
          );
          m_HaveMembrane = sp_membrane.IsDefined() && sp_membrane->IsDefined();
          if( m_HaveMembrane)
          {
            // get the half-width of the membrane
            m_MembraneWidth = sp_membrane->GetThickness( biol::GetEnvironmentTypes().e_MembraneCore);
            m_TransitionWidth
              = sp_membrane->GetThickness( biol::GetEnvironmentTypes().e_GapCoreTransition)
                + sp_membrane->GetThickness( biol::GetEnvironmentTypes().e_Transition);
          }
        }

        m_HaveCheckedForMembrane = true;
      }

      STORAGE( 0) = m_HaveMembrane ? ELEMENT->GetFirstSidechainAtom().GetCoordinates().Z() : m_DefaultValue;
      STORAGE( 0) = math::Absolute( STORAGE( 0));

      // handle normalization
      if( m_MembraneNorm && m_HaveMembrane)
      {
        if( STORAGE( 0) < m_MembraneWidth)
        {
          // in membrane core
          STORAGE( 0) /= m_MembraneWidth;
        }
        else
        {
          // in transition or solution
          STORAGE( 0) = 1.0 + ( STORAGE( 0) - m_MembraneWidth) / m_TransitionWidth;
        }
      }
      STORAGE( 0) = std::min( STORAGE( 0), float( m_DefaultValue));
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AADistanceToTMCenter::SetObjectHook()
    {
      m_HaveCheckedForMembrane = false;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AADistanceToTMCenter::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_MembraneNorm
        ? "Normalized distance to center of membrane; 0-1 -> in membrane core, 1-2 in transition region, "
          "2-3 soluble (near 2 for close to transition region)"
        : "Absolute distance from the center of the membrane, as defined by the PDBTM xml file "
          "(which must be present in same directory as the pdb file for this descriptor to work).  When absent, and for "
          "soluble proteins, uses the defualt value parameter"
      );
      parameters.AddInitializer
      (
        "default",
        "value to use for soluble proteins or when the PDBTM xml file is absent",
        io::Serialization::GetAgent( &m_DefaultValue),
        m_MembraneNorm ? "3.0" : "30.0"
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
