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
#include "descriptor/bcl_descriptor_aa_tm_direction.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"
#include "util/bcl_util_enumerated.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////
  // data //
  //////////

    namespace
    {
      util::ObjectInterface *AddInstances()
      {
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATMDirection( sspred::CIPhiPsi::e_Amphipathic, "Amphipathic"));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATMDirection( sspred::CIPhiPsi::e_EnteringCytosol, "TMSegmentEnteringCytosol"));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATMDirection( sspred::CIPhiPsi::e_Inside, "InsideMembrane"));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATMDirection( sspred::CIPhiPsi::e_LeavingCytosol, "TMSegmentLeavingCytosol"));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATMDirection( sspred::CIPhiPsi::e_NonMembrane, "NonMembrane"));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATMDirection( sspred::CIPhiPsi::e_Outside, "OutsideMembrane"));
        return util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATMDirection( sspred::CIPhiPsi::e_Pore, "Pore"));
      }
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AATMDirection::s_Instance
    (
      AddInstances()
    );
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new AATMDirection
    AATMDirection *AATMDirection::Clone() const
    {
      return new AATMDirection( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AATMDirection::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AATMDirection::GetAlias() const
    {
      return m_Alias;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AATMDirection::GetNormalSizeOfFeatures() const
    {
      // calculate and return
      return 1;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AATMDirection::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // get the prediction map for this AA
      util::SiPtr< const sspred::CIPhiPsi> prediction_method_result
      (
        ELEMENT->GetSSPrediction( sspred::GetMethods().e_CIPhiPsi)
      );

      if( prediction_method_result.IsDefined())
      {
        STORAGE( 0) = prediction_method_result->GetTMDirection() == m_Type ? 1.0 : 0.0;
      }
      else
      {
        STORAGE( 0) = util::GetUndefined< float>();
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AATMDirection::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Uses CiPhiPsi to determine the TM segment direction"
      );
      return parameters;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AATMDirection::SetObjectHook()
    {
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace descriptor
} // namespace bcl
