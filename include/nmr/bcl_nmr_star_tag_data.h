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

#ifndef BCL_NMR_STAR_TAG_DATA_H_
#define BCL_NMR_STAR_TAG_DATA_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_nmr_star_tag_categories.h"
#include "util/bcl_util_cpp_data_types.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StarTagData
    //! @brief stores data for NMR-Star tags to be used with the Enum class, StarTags
    //!
    //! @see @link example_nmr_star_tag_data.cpp @endlink
    //! @author weinerbe
    //! @date Jun 22, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StarTagData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! category this tag belongs to
      StarTagCategory m_TagCategory;

      //! string used to recognize this tag in an NMR-STAR file
      std::string m_Description;

      //! data type
      util::CPPDataTypes::TypeEnum m_DataType;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from tag category and data type
      //! @param CATEGORY category for this tag
      //! @param DESCRIPTION description string for this tag
      //! @param DATA_TYPE data type for this tag
      StarTagData
      (
        const StarTagCategory &CATEGORY = GetStarTagCategories().e_Undefined,
        const std::string &DESCRIPTION = "",
        const util::CPPDataTypes::Types DATA_TYPE = util::CPPDataTypes::e_Unknown
      );

      //! @brief Clone function
      //! @return pointer to new StarTagData
      StarTagData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the category for this tag
      //! @return the category for this tag
      const StarTagCategory &GetTagCategory() const
      {
        return m_TagCategory;
      }

      //! @brief returns the description for this tag
      //! @return the description for this tag
      const std::string &GetDescription() const
      {
        return m_Description;
      }

      //! @brief returns the data type of a tag
      //! @return the data type of a tag
      const util::CPPDataTypes::TypeEnum &GetDataType() const
      {
        return m_DataType;
      }

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

    }; // class StarTagData

  } // namespace nmr
} // namespace bcl

#endif // BCL_NMR_STAR_TAG_DATA_H_
