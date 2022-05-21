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

#ifndef BCL_NMR_STAR_TAG_CATEGORIES_H_
#define BCL_NMR_STAR_TAG_CATEGORIES_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_nmr_star_tag_category_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StarTagCategories
    //! @brief Enum-based class that stores NMR-Star tag categories to be used when reading in STAR
    //!        files using the handler.  See http://www.bmrb.wisc.edu/dictionary/3.1html/SuperGroupPage.html for more
    //!        info.
    //!
    //! @see @link example_nmr_star_tag_categories.cpp @endlink
    //! @author weinerbe
    //! @date Jun 22, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StarTagCategories :
      public util::Enumerate< StarTagCategoryData, StarTagCategories>
    {
      friend class util::Enumerate< StarTagCategoryData, StarTagCategories>;

    public:

    //////////
    // data //
    //////////

      StarTagCategory e_DistConstraint;      //!< distance constraint
      StarTagCategory e_DistConstraintList;  //!< distance constraint list
      StarTagCategory e_DistConstraintTree;  //!< distance constraint tree
      StarTagCategory e_DistConstraintValue; //!< distance constraint value
      StarTagCategory e_RDCConstraint;       //!< RDC constraint
      StarTagCategory e_RDCConstraintList;   //!< RDC constraint list

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      StarTagCategories();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a tag category from a string
      //! @param DESCRIPTION string description for the enum
      //! @return a tag category from a string
      static StarTagCategory GetCategoryFromString( const std::string &DESCRIPTION);

    ////////////////
    // operations //
    ////////////////

    }; // class StarTagCategories

    //! @brief get access to all star tag categories
    BCL_API const StarTagCategories &GetStarTagCategories();

  } // namespace nmr

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< nmr::StarTagCategoryData, nmr::StarTagCategories>;

  } // namespace util
} // namespace bcl

#endif // BCL_NMR_STAR_TAG_CATEGORIES_H_
