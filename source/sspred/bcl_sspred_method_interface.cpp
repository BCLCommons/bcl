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
#include "sspred/bcl_sspred_method_interface.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "math/bcl_math_statistics.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically
#include <iostream>
#include <sstream>

namespace bcl
{
  namespace sspred
  {

  //////////
  // data //
  //////////

    //! @brief returns static predictions vector3D set to coil
    //! @return static predictions vector3D set to coil
    const linal::Vector3D &MethodInterface::GetDefaultPredictionVector()
    {
      //! static predictions vector3D set to coil
      static const linal::Vector3D s_default_prediction_vector( 0.0, 0.0, 1.0);

      // edn
      return s_default_prediction_vector;
    }

    //! @brief returns static predictions linal::Matrix set to solution coil
    //! @return static predictions linal::Matrix set to solution  oil
    const linal::Matrix< double> &MethodInterface::GetDefaultPredictionMatrix()
    {
      //! static predictions matrix set to coil
      static const linal::Matrix< double> s_default_prediction_matrix( CreateDefaultPredictionMatrix());

      // edn
      return s_default_prediction_matrix;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment prediction
    //! @return three state, environment prediction
    linal::Vector3D MethodInterface::GetThreeStateTMPrediction() const
    {
      // instantiate empty vector
      linal::Vector3D three_state;

      // get the 9-state prediction tranposed so that SS elements are within rows
      linal::Matrix< double> nine_state_prediction_transposed( GetNineStatePrediction().Transposed());

      // sum up each column
      // the sum of properties are already equal to 1, so no need to normalize
      three_state += nine_state_prediction_transposed.GetRow( biol::GetSSTypes().COIL);
      three_state += nine_state_prediction_transposed.GetRow( biol::GetSSTypes().HELIX);
      three_state += nine_state_prediction_transposed.GetRow( biol::GetSSTypes().STRAND);

      // end
      return three_state;
    }

    //! @brief find the SSType with highest prediction and returns it
    //! @return SSType with highest prediction
    biol::SSType MethodInterface::GetOneStateSSPrediction() const
    {
      const linal::Vector3D three_state( GetThreeStatePrediction());
      if( !three_state.IsDefined())
      {
        return biol::GetSSTypes().e_Undefined;
      }
      // get the index of the maximum prediction
      size_t max_index( std::max_element( three_state.Begin(), three_state.End()) - three_state.Begin());
      return biol::SSType( max_index);
    }

    //! @brief find the biol::EnvironmentType with highest prediction and returns it
    //! @return TMTYpe with highest prediction
    biol::EnvironmentType MethodInterface::GetOneStateTMPrediction() const
    {
      // calculate the sum for each row
      linal::Vector< double> three_tm_state( GetThreeStateTMPrediction());

      if( !three_tm_state.IsDefined())
      {
        return biol::GetEnvironmentTypes().e_Undefined;
      }

      // find the maximum index
      const size_t max_index( math::Statistics::MaximumIndex( three_tm_state.Begin(), three_tm_state.End()));

      // end
      return biol::GetEnvironmentTypes().GetReducedTypes()( max_index);
    }

    //! @brief find the SSType-TMType pair with highest prediction and returns it
    //! @return SSType-TMType pair with highest prediction
    storage::Pair< biol::SSType, biol::EnvironmentType> MethodInterface::GetOneStateSSTMPrediction() const
    {
      // get the 9-state prediction
      linal::Matrix< double> nine_state_prediction( GetNineStatePrediction());

      // store the maximum index
      const size_t max_index
      (
        math::Statistics::MaximumIndex( nine_state_prediction.Begin(), nine_state_prediction.End())
      );

      // construct the pair and return it
      return storage::Pair< biol::SSType, biol::EnvironmentType>
      (
        biol::SSType( max_index % 3),
        biol::GetEnvironmentTypes().GetReducedTypes()( max_index / 3)
      );
    }

    //! @brief converts a three state prediction to a nine state predictions for given TM_TYPE
    //! @param THREE_STATE three state secondary structure prediction for given TM_TYPE
    //! @param TM_TYPE environment type which THREE_STATE predictions are for
    //! @return nine state prediction which has the provided THREE_STATE for given TM_TYPE and 0's for rest
    linal::Matrix< double> MethodInterface::ConvertThreeStateToNineState
    (
      const linal::Vector3D &THREE_STATE,
      const biol::EnvironmentType &TM_TYPE
    )
    {
      // instantiate empty matrix
      linal::Matrix< double> nine_state
      (
        biol::GetEnvironmentTypes().GetNumberReducedTypes(), biol::GetSSTypes().COIL.GetIndex() + 1, 0.0
      );

      // if the TM-type is known
      if( TM_TYPE.IsDefined())
      {
        // replace the row corresponding to given TM_TYPE with given THREE_STATE
        nine_state.ReplaceRow( TM_TYPE->GetReducedIndex(), THREE_STATE);
      }
      else
      {
        // else, assume soluble
        nine_state.ReplaceRow( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex(), THREE_STATE);
      }

      // end
      return nine_state;
    }

    //! @brief converts a nine state prediction to a three state predictions
    //! @param NINE_STATE nine state secondary structure prediction
    //! @return three state prediction constructed from a nine state prediction
    linal::Vector3D MethodInterface::ConvertNineStateToThreeState
    (
      const linal::Matrix< double> &NINE_STATE
    )
    {
      // instantiate empty vector
      linal::Vector3D three_state;

      // sum up each row
      // the sum of properties are already equal to 1, so no need to normalize
      three_state += NINE_STATE.GetRow( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex());
      three_state += NINE_STATE.GetRow( biol::GetEnvironmentTypes().e_Transition->GetReducedIndex());
      three_state += NINE_STATE.GetRow( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex());

      // end
      return three_state;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write secondary structure predictions in three-state form to given OSTREAM
    //! @param OSTREAM output stream
    //! @return std::ostream which was written to
    std::ostream &MethodInterface::WriteThreeStatePredictions( std::ostream &OSTREAM) const
    {
      // get three state predictions
      const linal::Vector3D prediction( GetThreeStatePrediction());

      // write data
      OSTREAM << GetOneStateSSPrediction()->GetOneLetterCode() << ' '
             << util::Format().W( 7).FFP( 3)( prediction( biol::GetSSTypes().COIL))
             << util::Format().W( 7).FFP( 3)( prediction( biol::GetSSTypes().HELIX))
             << util::Format().W( 7).FFP( 3)( prediction( biol::GetSSTypes().STRAND))
             << '\n';

      // return
      return OSTREAM;

    }

    //! @brief write secondary structure predictions in nine-state form to given OSTREAM
    //! @param OSTREAM output stream
    //! @return std::ostream which was written to
    std::ostream &MethodInterface::WriteNineStatePredictions( std::ostream &OSTREAM) const
    {
      // get three state predictions
      linal::Matrix< double> prediction( GetNineStatePrediction());

      // get the one state prediction
      storage::Pair< biol::SSType, biol::EnvironmentType> one_state( GetOneStateSSTMPrediction());

      static const util::Format s_format( util::Format().W( 7).FFP( 3));

      // write data
      OSTREAM << one_state.First()->GetOneLetterCode() << ' '
              << one_state.Second()->GetTwoLetterCode()<< ' '
              << s_format( prediction( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex(), biol::GetSSTypes().COIL))
              << s_format( prediction( biol::GetEnvironmentTypes().e_Transition->GetReducedIndex()  , biol::GetSSTypes().COIL))
              << s_format( prediction( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex()    , biol::GetSSTypes().COIL))
              << s_format( prediction( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex(), biol::GetSSTypes().HELIX))
              << s_format( prediction( biol::GetEnvironmentTypes().e_Transition->GetReducedIndex()  , biol::GetSSTypes().HELIX))
              << s_format( prediction( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex()    , biol::GetSSTypes().HELIX))
              << s_format( prediction( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex(), biol::GetSSTypes().STRAND))
              << s_format( prediction( biol::GetEnvironmentTypes().e_Transition->GetReducedIndex()  , biol::GetSSTypes().STRAND))
              << s_format( prediction( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex()    , biol::GetSSTypes().STRAND));

      // iterate over the matrix to compute confidence
      double highest( -std::numeric_limits< double>::max());
      double second_highest( -std::numeric_limits< double>::max());
      for( double *ptr( prediction.Begin()), *ptr_end( prediction.End()); ptr != ptr_end; ++ptr)
      {
        if( *ptr >= highest)
        {
          second_highest = highest;
          highest = *ptr;
        }
        else if( *ptr >= second_highest)
        {
          second_highest = *ptr;
        }
      }
      OSTREAM << s_format( highest - second_highest) << '\n';

      // return
      return OSTREAM;
    }

    //! @brief write secondary structure predictions to the provided OSTREAM
    //! @param OSTREAM output stream
    //! @return std::ostream which was written to
    std::ostream &MethodInterface::WritePredictions( std::ostream &OSTREAM) const
    {
      return WriteThreeStatePredictions( OSTREAM);
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which secondary structure predictions will be read
    //! @param SS_METHOD Method to be read
    //! @return std::istream which was read from
    std::istream &MethodInterface::ReadStandardPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE,
      const Method &SS_METHOD
    )
    {
      // iterate over all amino acids in the sequence
      for
      (
        biol::AASequence::iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        ISTREAM >> std::skipws;
        if( ISTREAM.good() && !ISTREAM.eof() && ISTREAM.peek() != '#')
        {
          // read the predictions for this amino acid
          MethodHandler::ReadPredictionsForAA( ISTREAM, **aa_itr, SS_METHOD);
        }
      }

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief creates a default prediction matrix with coil,solution set to 1, everything else to 0
    //! @return a default prediction matrix with coil,solution set to 1, everything else to 0
    linal::Matrix< double> MethodInterface::CreateDefaultPredictionMatrix()
    {
      // initialize a default matrix
      linal::Matrix< double> prediction_matrix( 3, 3, double( 0.0));

      // replace solution coil to 1
      prediction_matrix( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex(), biol::GetSSTypes().COIL) = 1.0;

      // end
      return prediction_matrix;
    }

  } // namespace sspred
} // namespace bcl
