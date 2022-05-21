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

#ifndef BCL_SCORE_PROTEIN_MODEL_SSE_LINEAR_LOOP_PROXIMITY_H_
#define BCL_SCORE_PROTEIN_MODEL_SSE_LINEAR_LOOP_PROXIMITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelSSELinearLoopProximity
    //! @brief score for proximity of imaginary linear loops between all neighborings sses to other sses in the model
    //!
    //! @see @link example_score_protein_model_sse_linear_loop_proximity.cpp @endlink
    //! @author mendenjl
    //! @date Feb 07, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelSSELinearLoopProximity :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! whether to scale by # of aas
      bool m_Normalize;

      //! bool: weight each close contact by
      //! (x+y+2*footpoint_offset)/(1+2*footpoint_offset), where
      //! x = closest point on sse norm distance from end (0-0.5)
      //! y = closest point on loop norm distance from end (0-0.5)
      bool m_ConsiderDistanceAlongSSE;

      //! foot point offset; see equation for m_ConsiderDistanceAlongSSE
      double m_FootPointOffset;

      //! option: whether to consider virtual loop clashes as well
      bool m_ConsiderVirtualLoopClashes;

      //! Maximum sequence distance between sses to consider; useful during folding and with incomplete models, where
      //! linear loop approximation breaks down, potentially in such a way that it prevents non-consecutive insertion of
      //! SSEs.  Set to undefined or similarly high number to consider all loops
      size_t m_MaximumSequenceSeparation;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct parameters
      //! @param NORMALIZE if true, final score will be normalized by the number of sses/sequences in the protein model
      //! @param CONSIDER_DISTANCE_ALONG_SSE onsider the distance along the sse if true; usually useful, see comment
      //!        for m_ConsiderDistanceAlongSSE
      //! @param FOOTPOINT_OFFSET sets m_FootPointOffset, see comment for m_ConsiderDistanceAlongSSE
      //! @param CONSIDER_VIRTUAL_LOOP_CLASHES whether to examine clashes between the virtual loops
      //! @param MAXIMUM_LINEAR_LOOP_RESIDUES maximum number of residues between adjacent sses for virtual loops; useful
      //!                                     so that incomplete models are not penalized in cases where the linear loop
      //!                                     approximation breaks down
      ProteinModelSSELinearLoopProximity
      (
        const bool NORMALIZE = false,
        const bool CONSIDER_DISTANCE_ALONG_SSE = true,
        const double &FOOTPOINT_OFFSET = 0.125,
        const bool &CONSIDER_VIRTUAL_LOOP_CLASHES = true,
        const size_t &MAXIMUM_LINEAR_LOOP_RESIDUES = 40
      );

      //! virtual copy constructor
      ProteinModelSSELinearLoopProximity *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that scores the chain
      //! @param CHAIN the chain for which all neighbor scores are calculated
      //! @param MISSING_SSES sse pool of missing sses
      //! @return score
      double operator()
      (
        const assemble::Chain &CHAIN,
        const assemble::SSEPool &MISSING_SSES
      ) const;

      //! @brief operator that scores the Protein model
      //! @param PROTEIN_MODEL the protein model for which all neighbor scores are calculated
      //! @return score
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        std::ostream &OSTREAM
      ) const;

    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @param INDENT indentation
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief helper function called by WriteDetailedSchemeAndValues and operator() so that the code remains in-sync
      //! @param CHAIN the chain of interest
      //! @param MISSING_SSES sse pool of missing sses
      //! @param OSTREAM the output stream to write the detailed scheme to for this chain
      //! @param DO_WRITE set to true to actually write to the output stream; otherwise, nothing will be written
      //! @return the final score
      double ScoreLoops
      (
        const assemble::Chain &CHAIN,
        const assemble::SSEPool &MISSING_SSES,
        std::ostream &OSTREAM,
        const bool &DO_WRITE
      ) const;

    }; // class ProteinModelSSELinearLoopProximity

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_PROTEIN_MODEL_SSE_LINEAR_LOOP_PROXIMITY_H_
