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
#include "descriptor/bcl_descriptor_aa_pore_orientation.h"

// includes from bcl - sorted alphabetically
#include "sspred/bcl_sspred_ci_phi_psi.h"
#include "sspred/bcl_sspred_mahssmi.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAPoreOrientation::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAPoreOrientation( false))
    );
    const util::SiPtr< const util::ObjectInterface> AAPoreOrientation::s_NInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAPoreOrientation( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new AAPoreOrientation
    AAPoreOrientation *AAPoreOrientation::Clone() const
    {
      return new AAPoreOrientation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAPoreOrientation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAPoreOrientation::GetAlias() const
    {
      static const std::string s_alias( "FacesPore"), s_tm_n_term_cytosol( "TMOriginatesInCytosol");
      return m_TMOrientation ? s_tm_n_term_cytosol : s_alias;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAPoreOrientation::GetNormalSizeOfFeatures() const
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
    void AAPoreOrientation::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // get the prediction map for this AA
      util::SiPtr< const sspred::Mahssmi> prediction_method_result
      (
        ELEMENT->GetSSPrediction( sspred::GetMethods().e_MAHSSMI)
      );

      // handle the case where the prediction is unavailable
      if( !prediction_method_result.IsDefined())
      {
        // use CIPhiPsi as backup
        util::SiPtr< const sspred::CIPhiPsi> prediction_method_result_ciphipsi
        (
          ELEMENT->GetSSPrediction( sspred::GetMethods().e_CIPhiPsi)
        );
        if( !prediction_method_result_ciphipsi.IsDefined())
        {
          STORAGE = util::GetUndefined< float>();
          return;
        }
        else
        {
          STORAGE( 0) = !m_TMOrientation ? prediction_method_result_ciphipsi->DoesBetaBarrelResidueFacePore()
                        : prediction_method_result_ciphipsi->GetTMDirection() == sspred::CIPhiPsi::e_LeavingCytosol;
        }
      }

      STORAGE( 0) = !m_TMOrientation ? prediction_method_result->DoesBetaBarrelResidueFacePore()
                    : prediction_method_result->OriginatesInCytosol();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPoreOrientation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_TMOrientation
        ? "Uses the Membrane Aware Hybrid Secondary Structure & Membrane topology Identification algorithm (unpublished) "
        "to find identify TM regions whose N-terminii are in the cytosol (returns 1). Extracellular N-terminii regions can be identified "
        "by taking Multiply(Not(TMOriginatesInCytosol),MembraneTransitionSolution(method=MAHSSMI/CIPhiPsi))"
        : "Uses the Membrane Aware Hybrid Secondary Structure & Membrane topology Identification algorithm (unpublished) "
        "to find beta barrel pores.  This descriptor returns 1 if the residue's side chain points towards the pore and "
        "the residue is in the membrane"
      );
      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace descriptor
} // namespace bcl
