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

#ifndef BCL_SCORE_SSE_PREDICTIONS_H_
#define BCL_SCORE_SSE_PREDICTIONS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_sse_prediction_interface.h"
#include "biol/bcl_biol_ss_types.h"
#include "math/bcl_math_piecewise_function.h"
#include "sspred/bcl_sspred_methods.h"
#include "storage/bcl_storage_map.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPredictions
    //! @brief score the ss prediction by scoring each individual amino acid prediction
    //! @details Uses predictions from specified sspred method, and score the log likelihodd of each amino acids prediction
    //!
    //! @see @link example_score_sse_predictions.cpp @endlink
    //! @author woetzen
    //! @date Jun 30, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPredictions :
      public SSEPredictionInterface
    {

    private:

    //////////
    // data //
    //////////

      //! sspred method to use in evaluation
      sspred::Method m_Method;

      //! in units of z-score, at which z-score, the confidence is negative
      double m_ConfidenceThreshold;

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! potentials depending on SSType
      storage::Map< biol::SSType, math::PiecewiseFunction> m_Potentials;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief multiple of the confidence interval for the z-score
      static const double s_DefaultConfidenceThreshold;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEPredictions();

      //! @brief constructor from a ssmethod and confidence threshold
      //! @param SS_METHOD ssmethod to use
      //! @param CONFIDENCE_THRESHOLD the threshold in units of z-score, above which the score becomes negative
      //! @param SCHEME scheme to be used
      SSEPredictions
      (
        const sspred::Method &SS_METHOD,
        const double CONFIDENCE_THRESHOLD = s_DefaultConfidenceThreshold,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new SSEPredictions
      SSEPredictions *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief access the potentials
      //! @return map of methods and sstype
      const storage::Map< biol::SSType, math::PiecewiseFunction> &GetEnergyFunctions() const
      {
        return m_Potentials;
      }

      //! @brief set the ss prediction method
      //! @param METHOD the method to use for scoring
      void SetMethod( const sspred::Method &METHOD)
      {
        m_Method = METHOD;
      }

      //! @brief access to the method used in scoring
      //! @return the method used
      const sspred::Method &GetMethod() const
      {
        return m_Method;
      }

      //! @brief set the confidence threshold above which scores become negative
      //! @param THRESHOLD confidence threshold in units of z-score
      void SetConfidenceThreshold( const double THRESHOLD)
      {
        m_ConfidenceThreshold = THRESHOLD;
      }

      //! @brief get the confidence threshold used
      //! @return the confidence threshold for scoring
      double GetConfidenceThreshold() const
      {
        return m_ConfidenceThreshold;
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates the score for a given SSE
      //! @param SSE SSE of interest
      //! @param MEMBRANE membrane object
      //! @return score calculated for the given SSE
      storage::Pair< double, size_t> operator()
      (
        const assemble::SSE &SSE,
        const biol::Membrane &MEMBRANE
      ) const;

      //! @brief operator that calculates the score for a given SSE for a specified SSMethod
      //! @param SSE SSE of interest
      //! @param SS_METHOD SSMethod of interest
      //! @return score calculated for the given SSE
      double operator()
      (
        const assemble::SSE &SSE,
        const sspred::Method &SS_METHOD
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param SSE SSE of interest
      //! @param MEMBRANE membrane object
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::SSE &SSE,
        const biol::Membrane &MEMBRANE,
        std::ostream &OSTREAM
      ) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief score a single amino acid
      //! @param AMIO_ACID amino acid to score
      //! @param SS_TYPE the sstype amino acid prediction will be evaluated for
      //! @param ZSCORE the zscore to use
      //! @return the score
      double ScoreAminoAcid( const biol::AABase &AMINO_ACID, const biol::SSType &SS_TYPE, const math::PiecewiseFunction &ZSCORE) const;

      //! @brief read energy distribution for scoring pairs of sses
      void ReadEnergyVector();

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

    }; // class SSEPredictions

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_SSE_PREDICTIONS_H_
