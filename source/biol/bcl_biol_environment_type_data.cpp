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
#include "biol/bcl_biol_environment_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_environment_types.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EnvironmentTypeData::s_Instance( GetObjectInstances().AddInstance( new EnvironmentTypeData()));

    //! @brief construct undefined environment type
    EnvironmentTypeData::EnvironmentTypeData() :
      m_TwoLetterCode( ""),
      m_ReducedType( ""),
      m_ReducedIndex( util::GetUndefined< size_t>()),
      m_IsGap( false),
      m_DefaultThickness( util::GetUndefined< double>())
    {
    }

    //! @brief construct EnvironmentTypeData
    //! @param TWO_LETTER_CODE two letter code for this environment type
    //! @param REDUCED_TYPE_NAME two letter code of reduced EnvironmentType
    //! @param REDUCED_INDEX index for the reduced EnvironmentType
    //! @param IS_GAP boolean to indicate whether this type corresponds to a gap region
    //! @param DEFAULT_THICKNESS efault thickness for this environment type
    EnvironmentTypeData::EnvironmentTypeData
    (
      const std::string &TWO_LETTER_CODE,
      const std::string &REDUCED_TYPE_NAME,
      const size_t REDUCED_INDEX,
      const bool IS_GAP,
      const double DEFAULT_THICKNESS
    ) :
      m_TwoLetterCode( TWO_LETTER_CODE),
      m_ReducedType( REDUCED_TYPE_NAME),
      m_ReducedIndex( REDUCED_INDEX),
      m_IsGap( IS_GAP),
      m_DefaultThickness( DEFAULT_THICKNESS)
    {
    }

    //! @brief virtual copy constructor
    EnvironmentTypeData *EnvironmentTypeData::Clone() const
    {
      return new EnvironmentTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EnvironmentTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get two letter code
    //! @return two letter code
    const std::string &EnvironmentTypeData::GetTwoLetterCode() const
    {
      return m_TwoLetterCode;
    }

    //! @brief get index for reduced type
    //! @return index for reduced type
    size_t EnvironmentTypeData::GetReducedIndex() const
    {
      return m_ReducedIndex;
    }

    //! @brief get the reduced environment type
    //! this maps environment types to one of three states: Core, Transition or Solution
    //! @return reduced EnvironmentType
    const EnvironmentType &EnvironmentTypeData::GetReducedType() const
    {
      return GetEnvironmentTypes().EnvironmentTypeFromTwoLetterCode( m_ReducedType);
    }

    //! @brief get the reduced EnvironmentType's name
    //! @return name of reduced EnvironmentType for this type
    const std::string &EnvironmentTypeData::GetReducedTypeString() const
    {
      return m_ReducedType;
    }

    //! @brief get default thickness
    double EnvironmentTypeData::GetDefaultThickness() const
    {
      return m_DefaultThickness;
    }

    //! @brief get whether this is a gap region
    //! @return whether this is a gap region
    bool EnvironmentTypeData::IsGap() const
    {
      return m_IsGap;
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EnvironmentTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_TwoLetterCode, ISTREAM);
      io::Serialize::Read( m_ReducedType, ISTREAM);
      io::Serialize::Read( m_ReducedIndex, ISTREAM);
      io::Serialize::Read( m_IsGap, ISTREAM);
      io::Serialize::Read( m_DefaultThickness, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &EnvironmentTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_TwoLetterCode, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ReducedType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ReducedIndex, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IsGap, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DefaultThickness, OSTREAM, INDENT) << '\n';

      // return
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
