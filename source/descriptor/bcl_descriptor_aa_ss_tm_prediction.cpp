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
#include "descriptor/bcl_descriptor_aa_ss_tm_prediction.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
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
    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AASSTMPrediction( "MembraneEnvironment", false, false));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AASSTMPrediction( "MembraneTransitionSolution", true, false));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AASSTMPrediction( "HelixStrandCoil", false, true));
        return util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AASSTMPrediction( "SecondaryStructure9D", true, true));
      }
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AASSTMPrediction::s_Instance( AddInstances());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a single Method
    //! @param METHOD Method of interest
    //! @param CONSIDER_MEMBRANE whether to consider the membrane states
    //! @param CONSIDER_SS whether to consider secondary structure states
    AASSTMPrediction::AASSTMPrediction
    (
      const std::string &ALIAS,
      const bool CONSIDER_MEMBRANE,
      const bool CONSIDER_SS
    ) :
      m_ConsiderMembraneStates( CONSIDER_MEMBRANE),
      m_ConsiderSecondaryStructure( CONSIDER_SS),
      m_MinHelixSize( 0),
      m_MinStrandSize( 0),
      m_MinMembraneSpan( 0),
      m_Alias( ALIAS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AASSTMPrediction
    AASSTMPrediction *AASSTMPrediction::Clone() const
    {
      return new AASSTMPrediction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AASSTMPrediction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AASSTMPrediction::GetAlias() const
    {
      return m_Alias;
    }

    //! @brief set the method used for calculation
    //! @param METHOD the new method to use
    void AASSTMPrediction::SetMethod( const sspred::Method &METHOD)
    {
      m_Method = METHOD;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AASSTMPrediction::GetNormalSizeOfFeatures() const
    {
      // calculate and return
      return ( 2 * int( m_ConsiderMembraneStates) + 1) * ( 2 * int( m_ConsiderSecondaryStructure) + 1);
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
    void AASSTMPrediction::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // get the prediction map for this AA
      util::SiPtr< const sspred::MethodInterface> prediction_method_result( ELEMENT->GetSSPrediction( m_Method));

      // handle the case where the prediction is unavailable
      if( !prediction_method_result.IsDefined())
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }

      // if considering membrane
      if( m_ConsiderMembraneStates)
      {
        if( m_ConsiderSecondaryStructure)
        {
          // store the 9 state prediction
          linal::Matrix< double> prediction( prediction_method_result->GetNineStatePrediction());
          // copy the 9 state probabilities into the storage
          STORAGE.CopyValues( linal::Vector< float>( prediction.Begin(), prediction.End()));
        }
        else
        {
          // store the 3 state prediction for membrane prediction
          linal::Vector3D prediction( prediction_method_result->GetThreeStateTMPrediction());
          STORAGE.CopyValues( linal::Vector< float>( prediction.Begin(), prediction.End()));
        }
      }
      else if( m_ConsiderSecondaryStructure)
      {
        // store the 3 state prediction
        linal::Vector3D prediction( prediction_method_result->GetThreeStatePrediction());
        STORAGE.CopyValues( linal::Vector< float>( prediction.Begin(), prediction.End()));
      }
      else
      {
        STORAGE( 0) = prediction_method_result->GetOneStateTMPrediction().GetIndex();
      }

      biol::EnvironmentType env_type( prediction_method_result->GetOneStateTMPrediction());
      // filter by membrane span size
      if( m_MinMembraneSpan && env_type == biol::GetEnvironmentTypes().e_MembraneCore)
      {
        size_t membrane_size( 0);
        util::SiPtr< const sspred::MethodInterface> neighbor_prediction_method_result_type
        (
          prediction_method_result
        );
        // walk to the end of the range
        for
        (
          iterate::Generic< const biol::AABase> itr( ELEMENT);
          itr.NotAtEnd()
          && ( neighbor_prediction_method_result_type = itr->GetSSPrediction( m_Method)).IsDefined()
          && neighbor_prediction_method_result_type->GetOneStateTMPrediction() == env_type
          && membrane_size < m_MinMembraneSpan;
          ++itr, ++membrane_size
        )
        {
        }

        // walk to the beginning of the range; subtract one from membrane size because ELEMENT is counted twice
        iterate::Generic< const biol::AABase> itr( ELEMENT);
        for
        (
          --membrane_size;
          !itr.AtBegin()
          && ( neighbor_prediction_method_result_type = itr->GetSSPrediction( m_Method)).IsDefined()
          && neighbor_prediction_method_result_type->GetOneStateTMPrediction() == env_type
          && membrane_size < m_MinMembraneSpan;
          --itr, ++membrane_size
        )
        {
        }
        if
        (
          membrane_size == m_MinMembraneSpan - 1
          && itr.AtBegin()
          && ( neighbor_prediction_method_result_type = itr->GetSSPrediction( m_Method)).IsDefined()
          && neighbor_prediction_method_result_type->GetOneStateTMPrediction() == env_type
        )
        {
          ++membrane_size;
        }
        if( membrane_size < m_MinMembraneSpan)
        {
          env_type = biol::GetEnvironmentTypes().e_Solution;
          if( m_ConsiderMembraneStates)
          {
            if( m_ConsiderSecondaryStructure)
            {
              // store the 9 state prediction
              linal::Matrix< double> prediction
              (
                prediction_method_result->ConvertThreeStateToNineState
                (
                  prediction_method_result->GetThreeStatePrediction(),
                  env_type
                )
              );
              // copy the 9 state probabilities into the storage
              STORAGE.CopyValues( linal::Vector< float>( prediction.Begin(), prediction.End()));
            }
            else
            {
              // Set solution-only prediction
              STORAGE( 0) = STORAGE( 1) = 0.0;
              STORAGE( 2) = 1.0;
            }
          }
          else
          {
            STORAGE( 0) = 2.0;
          }
        }
      }
      biol::SSType type( prediction_method_result->GetOneStateSSPrediction());
      if
      (
        ( m_MinHelixSize && type == biol::GetSSTypes().HELIX)
        ||
        ( m_MinStrandSize && type == biol::GetSSTypes().STRAND)
      )
      {
        const size_t size_constraint( type == biol::GetSSTypes().HELIX ? m_MinHelixSize : m_MinStrandSize);
        size_t sse_size( 0);
        util::SiPtr< const sspred::MethodInterface> neighbor_prediction_method_result_type
        (
          prediction_method_result
        );
        // walk to the end of the range
        for
        (
          iterate::Generic< const biol::AABase> itr( ELEMENT);
          itr.NotAtEnd()
          && ( neighbor_prediction_method_result_type = itr->GetSSPrediction( m_Method)).IsDefined()
          && neighbor_prediction_method_result_type->GetOneStateSSPrediction() == type
          && sse_size < size_constraint;
          ++itr, ++sse_size
        )
        {
        }

        // walk to the beginning of the range; subtract one from sse size because ELEMENT is counted twice

        iterate::Generic< const biol::AABase> itr( ELEMENT);
        for
        (
          --sse_size;
          !itr.AtBegin()
          && ( neighbor_prediction_method_result_type = itr->GetSSPrediction( m_Method)).IsDefined()
          && neighbor_prediction_method_result_type->GetOneStateSSPrediction() == type
          && sse_size < size_constraint;
          --itr, ++sse_size
        )
        {
        }
        if
        (
          sse_size == size_constraint - 1
          && itr.AtBegin()
          && ( neighbor_prediction_method_result_type = itr->GetSSPrediction( m_Method)).IsDefined()
          && neighbor_prediction_method_result_type->GetOneStateSSPrediction() == type
        )
        {
          ++sse_size;
        }
        if( sse_size < size_constraint)
        {
          // if considering membrane
          if( m_ConsiderMembraneStates)
          {
            // store the 9 state prediction
            sspred::Mahssmi revised_sse_type( biol::GetSSTypes().COIL, env_type, false);
            linal::Matrix< double> prediction( revised_sse_type.GetNineStatePrediction());
            // copy the 9 state probabilities into the storage
            STORAGE.CopyValues( linal::Vector< float>( prediction.Begin(), prediction.End()));
          }
          else if( m_ConsiderSecondaryStructure)
          {
            // Set coil-only prediction
            STORAGE( 0) = STORAGE( 1) = 0.0;
            STORAGE( 2) = 1.0;
          }
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AASSTMPrediction::GetSerializer() const
    {
      io::Serializer parameters;
      std::string description;
      if( m_ConsiderSecondaryStructure)
      {
        if( m_ConsiderMembraneStates)
        {
          description =
            "Prediction of probabilities of each combination of secondary structure\n    "
            "(H=Helix, E=Strand, C=Coil) with environment (M=Membrane, T=Transition, S=Solution)\n    "
            "in the following order: HM EM CM HT ET CT HS ES CS";
        }
        else
        {
          description = "Prediction of secondary structure probabilities in the order: Helix Strand Coil";
        }
      }
      else if( m_ConsiderMembraneStates)
      {
        description = "Prediction of protein environment probabilities in the order: Membrane, Transition, Solution";
      }
      else
      {
        description =
          "Protein environment prediction 1 state\n    "
          "returns 0 if AA is predicted to be membrane, 1 for the transition region, and 2 for solution,\n    ";
      }
      description +=
        "requires a prediction file with the methods normal extension (e.g. JUFO needs a .jufo file, SAM needs a .rdb6 file)\n    "
        "with the same basename as the .fasta or .pdb file in the same directory as the pdb or fasta file";
      parameters.SetClassDescription( description);
      parameters.AddInitializer
      (
        "method",
        "the method to use to perform the prediction",
        io::Serialization::GetAgent( &m_Method)
      );
      if( m_ConsiderMembraneStates)
      {
        parameters.AddOptionalInitializer
        (
          "min membrane size",
          "If fewer than this many consecutive residues are in the membrane, return soluble",
          io::Serialization::GetAgent( &m_MinMembraneSpan)
        );
      }
      if( m_ConsiderSecondaryStructure)
      {
        parameters.AddOptionalInitializer
        (
          "min strand size",
          "If fewer than this many consecutive residues are in a strand, it is reclassified as a coil",
          io::Serialization::GetAgent( &m_MinStrandSize)
        );
        parameters.AddOptionalInitializer
        (
          "min helix size",
          "If fewer than this many consecutive residues are in a helix, it is reclassified as a coil",
          io::Serialization::GetAgent( &m_MinHelixSize)
        );
      }
      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace descriptor
} // namespace bcl
