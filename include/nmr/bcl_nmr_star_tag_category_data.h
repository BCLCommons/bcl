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

#ifndef BCL_NMR_STAR_TAG_CATEGORY_DATA_H_
#define BCL_NMR_STAR_TAG_CATEGORY_DATA_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StarTagCategoryData
    //! @brief stores NMR-Star tag category data to be used in the StarTagCategoris Enum class
    //!
    //! @see @link example_nmr_star_tag_category_data.cpp @endlink
    //! @author weinerbe
    //! @date Jun 22, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StarTagCategoryData :
      public util::ObjectInterface
    {
    public:
    //////////
    // data //
    //////////

      //! save frames that a tag category can belong to
      enum SaveFrame
      {
        e_DistanceConstraints, //!< distance constraints
        e_RDCConstraints,      //!< RDC constraints
        s_NumberSaveFrames
      };

      //! @brief SaveFrame as string
      //! @param FRAME the save frame
      //! @return the string for FRAME
      static const std::string &GetSaveFrameName( const SaveFrame &FRAME);

      //! @brief SaveFrame enum I/O helper
      typedef util::WrapperEnum< SaveFrame, &GetSaveFrameName, s_NumberSaveFrames> SaveFrameEnum;

    private:

    //////////
    // data //
    //////////

      //! save frame that this tag category belongs to
      SaveFrameEnum m_SaveFrame;

      //! string used to recognize this category in an NMR-STAR file
      std::string m_Description;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from StarSaveFrame
      //! @param SAVE_FRAME save frame to assign to this tag category
      //! @param DESCRIPTION string used to recognize this category in an NMR-STAR file
      StarTagCategoryData
      (
        const SaveFrame &SAVE_FRAME = s_NumberSaveFrames,
        const std::string &DESCRIPTION = ""
      );

      //! @brief Clone function
      //! @return pointer to new StarTagCategoryData
      StarTagCategoryData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the save frame
      //! @return the save frame
      const SaveFrameEnum &GetSaveFrame() const
      {
        return m_SaveFrame;
      }

      //! @brief returns the description
      //! @return the description
      const std::string &GetDescription() const
      {
        return m_Description;
      }

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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class StarTagCategoryData

  } // namespace nmr
} // namespace bcl

#endif // BCL_NMR_STAR_TAG_CATEGORY_DATA_H_
