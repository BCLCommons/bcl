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

#ifndef BCL_BIOL_ENVIRONMENT_TYPE_DATA_H_
#define BCL_BIOL_ENVIRONMENT_TYPE_DATA_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnvironmentTypeData
    //! @brief This is a low level helper class to store environment properties
    //! @details This class provides storage for representing an individual EnvironmentType enumerated in the class
    //! EnvironmentTypes
    //!
    //! @see @link example_biol_environment_type_data.cpp @endlink
    //! @author karakam, woetzen
    //! @date 08/27/09
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EnvironmentTypeData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string m_TwoLetterCode;    //!< two  letter code
      std::string m_ReducedType;      //!< two letter code for reduced types (3 states): MembraneCore, Transition or Solution
      size_t      m_ReducedIndex;     //!< index for the reduced type
      bool        m_IsGap;            //!< boolean to indicate whether this type corresponds to a gap region
      double      m_DefaultThickness; //!< default thickness for this environment type

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined environment type
      EnvironmentTypeData();

      //! @brief construct EnvironmentTypeData
      //! @param TWO_LETTER_CODE two letter code for this environment type
      //! @param REDUCED_TYPE_NAME two letter code of reduced EnvironmentType
      //! @param REDUCED_INDEX index for the reduced EnvironmentType
      //! @param IS_GAP boolean to indicate whether this type corresponds to a gap region
      //! @param DEFAULT_THICKNESS efault thickness for this environment type
      EnvironmentTypeData
      (
        const std::string &TWO_LETTER_CODE,
        const std::string &REDUCED_TYPE_NAME,
        const size_t REDUCED_INDEX,
        const bool IS_GAP,
        const double DEFAULT_THICKNESS
      );

      //! @brief virtual copy constructor
      EnvironmentTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get two letter code
      //! @return two letter code
      const std::string &GetTwoLetterCode() const;

      //! @brief get index for reduced type
      //! @return index for reduced type
      size_t GetReducedIndex() const;

      //! @brief get the reduced environment type
      //! this maps environment types to one of three states: Core, Transition or Solution
      //! @return reduced EnvironmentType
      const EnvironmentType &GetReducedType() const;

      //! @brief get the reduced EnvironmentType's name
      //! @return name of reduced EnvironmentType for this type
      const std::string &GetReducedTypeString() const;

      //! @brief get default thickness
      double GetDefaultThickness() const;

      //! @brief get whether this is a gap region
      //! @return whether this is a gap region
      bool IsGap() const;

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; //class EnvironmentTypeData

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_ENVIRONMENT_TYPE_DATA_H_

