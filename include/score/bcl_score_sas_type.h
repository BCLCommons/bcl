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

#ifndef BCL_SCORE_SAS_TYPE_H_
#define BCL_SCORE_SAS_TYPE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasType
    //! @brief Calculate similarity score between two SAXS profiles based on provided scoring type
    //! @details Currently set up for chi, cumulative integral, and stovgaard scoring types
    //! @see @link example_score_saxs_type.cpp @endlink
    //! @author putnamdk
    //! @date April 24, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasType :
      public math::FunctionInterfaceSerializable< restraint::SasExperimentalAndCalculatedData, double>
    {

    public:

      enum ScoreFunction
      {
        e_chi,
        e_cumulative,
        e_stovgaard,
        s_NumberScoreTypes
      };

      //! @brief ScoreFunction as string
      //! @param SCORE_FUNCTION the ScoreFunction
      //! @return the string for the ScoreFunction
      static const std::string &GetFunctionDescriptor( const ScoreFunction &SCORE_FUNCTION);

      //! @brief ScoreFunctionEnum is used for I/O of ScoreFunction
      typedef util::WrapperEnum< ScoreFunction, &GetFunctionDescriptor, s_NumberScoreTypes> ScoreFunctionEnum;

    private:

    //////////
    // data //
    //////////

      //! @brief flag to use errors in profile comparison
      bool m_UseErrors;
      ScoreFunctionEnum m_ScoreType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Constructor
      SasType();

      //! @brief Constructor
      SasType( const bool &USE_ERRORS, const ScoreFunctionEnum &SCORE_TYPE);

      //! @brief Clone function
      //! @return pointer to new SasType
      SasType *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns whether errors are considered
      //! @return true if errors are considered
      const bool &GetUseErrors() const
      {
        return m_UseErrors;
      }

      //! @brief returns the score function enum
      //! @return the score function enum
      const ScoreFunctionEnum &GetScoreFunction() const
      {
        return m_ScoreType;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief function to calculate chi score with or without the derivative for different SAXS intensity curves
      //! @param SAXS_DATA experimental and calculated saxs data
      //! @return chi score
      double CalculateChiScore( const restraint::SasExperimentalAndCalculatedData &SAXS_DATA) const;

      //! @brief function to calculate cumulative integral score for different SAXS intensity curves
      //! @param SAXS_DATA experimental and calculated saxs data
      //! @return cumulative integral score
      double CalculateCumulativeIntegralScore( const restraint::SasExperimentalAndCalculatedData &SAXS_DATA) const;

      //! @brief function to calculate Stovgaard score for different SAXS intensity curves
      //! @param SAXS_DATA experimental and calculated saxs data
      //! @return Stovgaard score
      double CalculateStovgaardScore( const restraint::SasExperimentalAndCalculatedData &SAXS_DATA) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief overloaded () operator to calculate RMSD from two SAXS curves
      //! @param SAXS_DATA experimental and calculated saxs data
      //! @return return RMSD for two SAXS curves
      double operator()( const restraint::SasExperimentalAndCalculatedData &SAXS_DATA) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class SasType

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_SAS_TYPE_H_
