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

#ifndef BCL_SCORE_PROTEIN_MODEL_INVERTED_H_
#define BCL_SCORE_PROTEIN_MODEL_INVERTED_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_protein_model_inverter.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelInverted
    //! @brief wrapper class for scoring an inverted model
    //! @details This wrapper class inverts the given ProteinModel and scores it
    //!
    //! @see @link example_score_protein_model_inverted.cpp @endlink
    //! @author karakam
    //! @date Oct 9, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelInverted :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! ShPtr to scoring function to be used
      util::ShPtr< ProteinModel> m_Score;

      //! ShPtr to the inverter to be used
      util::ShPtr< assemble::ProteinModelInverter> m_Inverter;

      //! scheme
      const std::string m_Scheme;

      //! score type
      ProteinModel::Type m_ScoreType;

      //! readable scheme
      std::string m_ReadableScheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelInverted();

      //! @brief constructor from a score function, a inverter and a scheme
      //! @param SP_SCORE ShPtr to ProteinModel scoring function to be used
      //! @param SP_INVERTER ShPtr to ProteinModelInverter to be used
      //! @param SCHEME Scheme to be used
      //! @param SCORE_TYPE score type
      //! @param READABLE_SCHEME scheme that is more human readable
      ProteinModelInverted
      (
        const util::ShPtr< ProteinModel> &SP_SCORE,
        const util::ShPtr< assemble::ProteinModelInverter> &SP_INVERTER,
        const std::string &SCHEME,
        const ProteinModel::Type &SCORE_TYPE = ProteinModel::e_Undefined,
        const std::string &READABLE_SCHEME = ""
      );

      //! @brief Clone function
      //! @return pointer to new ProteinModelInverted
      ProteinModelInverted *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return scheme
      //! @return scheme
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief get a more readable score scheme
      //! @return a more readable score scheme
      const std::string &GetReadableScheme() const
      {
        return m_ReadableScheme;
      }

      //! @brief get score type
      //! @return score type
      Type GetType() const
      {
        return m_ScoreType;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief invert the given ProteinModel and score it
      //! @param PROTEIN_MODEL ProteinModel to be inverted and scored
      //! @return score of the inverted model
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ProteinModelInverted

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_INVERTED_H_ 
