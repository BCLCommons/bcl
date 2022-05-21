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

#ifndef BCL_SCORE_LOOP_CLOSURE_H_
#define BCL_SCORE_LOOP_CLOSURE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopClosure
    //! @brief based on statistics, penalizes the distance between the end of two SSEs if they cannot be bridged by the
    //! number of residues that are in between
    //! @details This pseudo-scoring function determines the number of loops in the model that can't be closed
    //!
    //! @see @link example_score_loop_closure.cpp @endlink
    //! @author karakam
    //! @date May 10, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoopClosure :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      //! Number of excluded residues from each end of loop
      size_t m_NrExcludedResidues;

      //! width of the sigmoidal function to be used - use 0 for linear
      double m_SigmoidWidth;

      //! fraction of allowed distance to use from the allowed
      double m_FractionAllowedDistance;

      //! true if coil sses should not be considered in scoring
      bool m_ExcludeCoil;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the name of the object
      //! @return the name of the object
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoopClosure();

      //! @brief constructor from a sigmoid width
      //! @param NR_EXCLUDED_RESIDUES Number of excluded residues from each end of loop
      //! @param SIGMOID_WIDTH width of sigmoid function
      //! @param FRACTION_ALLOWED_DISTANCE the allowed distance is multiplied with the fraction before calculating the penalty
      //! @param EXCLUDE_COIL true if coil sses should not be considered in scoring
      LoopClosure
      (
        const size_t NR_EXCLUDED_RESIDUES,
        const double SIGMOID_WIDTH,
        const double FRACTION_ALLOWED_DISTANCE,
        const bool EXCLUDE_COIL = true
      );

      //! @brief Clone function
      //! @return pointer to new LoopClosure
      LoopClosure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the number of excluded residues
      //! @return the number of excluded residues
      const size_t GetNrExcludedResidues() const
      {
        return m_NrExcludedResidues;
      }

      //! @brief gets the width of the sigmoidal function
      //! @return the width of the sigmoidal function
      const double GetSigmoidWidth() const
      {
        return m_SigmoidWidth;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate Sequence distance between aas and Euclidean distance between last and first amino acids of SSEs
      //! if the NR_EXCLUDED_RESIDUES is >= size of sse, it will be reduced, so that at least 1 amino acid remains
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @param NR_EXCLUDED_RESIDUES number of residues to be excluded from each end
      //! @return a pair which first member is the length of the loop and aminoacids and the second is the Euclidean distance
      static
      storage::Pair< size_t, double>
      SequenceAndEuclideanDistanceWithExclusion
      (
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B,
        const size_t NR_EXCLUDED_RESIDUES
      );

      //! @brief score the loop
      //! @param SEQ_EUC_DISTANCE pair of sequence and Euclidean distance
      //! @return the loop closure score
      double ScoreLoop( const storage::Pair< size_t, double> &SEQ_EUC_DISTANCE) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief gives a penalty for two sses based on the difference between the observed distance and the distance the
      //! missing amino acids could bridge
      //! @param SSE_A the first sse
      //! @param SSE_B the second sse
      //! @return the score for two sses based on the possibility to close teh loop between SSE_A and SSE_B
      double operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

      //! @brief operator that scores the chain
      //! @param CHAIN the chain for which all neighbor scores are calculated
      //! @return score
      double operator()( const assemble::Chain &CHAIN) const;

      //! @brief operator that scores the Protein model
      //! @param PROTEIN_MODEL the protein model for which all neighbor scores are calculated
      //! @return score
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    public:

      //! @brief write the Scheme and the Function value for the ARGUMENT to the STREAM
      //! @param SSE_A the first sse
      //! @param SSE_B the second sse
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::SSE &SSE_A, const assemble::SSE &SSE_B,
        std::ostream &OSTREAM
      ) const;

    }; // class LoopClosure

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_LOOP_CLOSURE_H_
