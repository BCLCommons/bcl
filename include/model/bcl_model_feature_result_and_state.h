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

#ifndef BCL_MODEL_FEATURE_RESULT_AND_STATE_H_
#define BCL_MODEL_FEATURE_RESULT_AND_STATE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureResultAndState
    //! @brief References a feature, result, and state (for classification)
    //!
    //! @tparam t_DataType can be float, float, int, complex, etc...
    //!
    //! @see @link example_model_feature_result_and_state.cpp @endlink
    //! @author mendenjl
    //! @date May 09, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FeatureResultAndState :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      FeatureReference< float>  m_Feature;     //!< feature reference
      FeatureReference< float>  m_Result;      //!< result reference
      FeatureReference< size_t> m_ResultState; //!< state vector, each value is whether m_Result[x] <= activity_cutoff

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FeatureResultAndState();

      //! @brief constructor from members
      FeatureResultAndState
      (
        const FeatureReference< float>  &FEATURE,
        const FeatureReference< float>  &RESULT,
        const FeatureReference< size_t> &RESULT_STATE
      );

      //! @brief Clone function
      //! @return pointer to new FeatureResultAndState
      FeatureResultAndState *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the features
      //! @return the features
      const FeatureReference< float> &GetFeature() const;

      //! @brief get the result
      //! @return the result
      const FeatureReference< float> &GetResult() const;

      //! @brief get the result state
      //! @return the result state
      const FeatureReference< size_t> &GetResultState() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class FeatureResultAndState

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FEATURE_RESULT_AND_STATE_H_ 
