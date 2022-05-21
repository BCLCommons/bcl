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

#ifndef BCL_SCORE_PROTEIN_MODEL_DEFINED_LOOPS_H_
#define BCL_SCORE_PROTEIN_MODEL_DEFINED_LOOPS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "biol/bcl_biol_aa_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelDefinedLoops
    //! @brief This class scores the loop completeness of a protein model
    //! @detail The completeness of the protein model is scored by computing the number of loops with at least
    //! partially undefined coordinates.
    //!
    //! @see @link example_score_protein_model_defined_loops.cpp @endlink
    //! @author fischea
    //! @date May 7, 2016
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ProteinModelDefinedLoops :
      public ProteinModel
    {

    //////////
    // data //
    //////////

    private:

      //! ignore terminal loops for the calculation of the completeness
      bool m_IgnoreTermLoops;

      //! scheme of this score
      std::string m_Scheme;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief construct from default values
      //! @param IGNORE_TERM_LOOPS ignore terminal loops for the calculation of the completeness
      //! @param SCHEME scheme of the score
      ProteinModelDefinedLoops( bool IGNORE_TERM_LOOPS = false, const std::string &SCHEME = GetDefaultScheme());

      //! @brief returns a pointer to a new ProteinModelDefinedLoops
      //! @return pointer to a new ProteinModelDefinedLoops
      ProteinModelDefinedLoops *Clone() const;

    /////////////////
    // data access //
    /////////////////
      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the scheme of this score
      //! @return the scheme of this score
      const std::string &GetScheme() const;

      //! @brief returns the default scheme of this score
      //! @return the default scheme of this score
      static const std::string &GetDefaultScheme();

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    ////////////////
    // operations //
    ////////////////

      //! @brief scores the loop completeness of the given protein model
      //! @detail the loop completeness is calculated by computing the number of loops with fully defined
      //! backbone coordinates.
      //! @param PROTEIN_MODEL protein model for which to compute the loop completeness
      //! @return loop completeness of the given protein model with -1.0 being complete and 0.0 being empty
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads members from an input stream
      //! @param ISTREAM input stream to read members from
      //! @return the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into an output stream
      //! @param OSTREAM output stream to write members into
      //! @INDENT number of indentations to use
      //! @return the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief determines if all coordinates in the given SSE are defined
      //! @param SSE SSE for which to determine if all coordinates are defined
      //! @return true, if all coordinates of the given SSE are defined
      bool IsDefined( const assemble::SSE &SSE) const;

    }; // class ProteinModelDefinedLoops

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_DEFINED_LOOPS_H_
