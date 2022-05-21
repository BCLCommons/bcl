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
#include "model/bcl_model_feature_result_and_state.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FeatureResultAndState::s_Instance
    (
      GetObjectInstances().AddInstance( new FeatureResultAndState())
    );

    //! @brief default constructor
    FeatureResultAndState::FeatureResultAndState() :
      m_Feature(),
      m_Result(),
      m_ResultState()
    {
    }

    //! @brief constructor from members
    FeatureResultAndState::FeatureResultAndState
    (
      const FeatureReference< float>  &FEATURE,
      const FeatureReference< float>  &RESULT,
      const FeatureReference< size_t> &RESULT_STATE
    ) :
      m_Feature( FEATURE),
      m_Result( RESULT),
      m_ResultState( RESULT_STATE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FeatureResultAndState
    FeatureResultAndState *FeatureResultAndState::Clone() const
    {
      return new FeatureResultAndState( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FeatureResultAndState::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the features
    //! @return the features
    const FeatureReference< float> &FeatureResultAndState::GetFeature() const
    {
      return m_Feature;
    }

    //! @brief get the result
    //! @return the result
    const FeatureReference< float> &FeatureResultAndState::GetResult() const
    {
      return m_Result;
    }

    //! @brief get the result state
    //! @return the result state
    const FeatureReference< size_t> &FeatureResultAndState::GetResultState() const
    {
      return m_ResultState;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FeatureResultAndState::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Feature, ISTREAM);
      io::Serialize::Read( m_Result, ISTREAM);
      io::Serialize::Read( m_ResultState, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FeatureResultAndState::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Feature, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Result, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ResultState, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
