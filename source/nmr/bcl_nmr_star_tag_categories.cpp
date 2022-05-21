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
#include "nmr/bcl_nmr_star_tag_categories.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StarTagCategories::StarTagCategories() :
      e_DistConstraint(      AddEnum( "DistConstraint",      StarTagCategoryData( StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint"))),
      e_DistConstraintList(  AddEnum( "DistConstraintList",  StarTagCategoryData( StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint_list"))),
      e_DistConstraintTree(  AddEnum( "DistConstraintTree",  StarTagCategoryData( StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint_tree"))),
      e_DistConstraintValue( AddEnum( "DistConstraintValue", StarTagCategoryData( StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint_value"))),
      e_RDCConstraint(       AddEnum( "RDCConstraint",       StarTagCategoryData( StarTagCategoryData::e_RDCConstraints,      "_RDC_constraint"))),
      e_RDCConstraintList(   AddEnum( "RDCConstraintList",   StarTagCategoryData( StarTagCategoryData::e_RDCConstraints,      "_RDC_constraint_list")))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarTagCategories::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a tag category from a string
    //! @param DESCRIPTION string description for the enum
    //! @return a tag category from a string
    StarTagCategory StarTagCategories::GetCategoryFromString( const std::string &DESCRIPTION)
    {
      // iterate through the tag categories
      for
      (
        StarTagCategories::const_iterator category_itr( GetStarTagCategories().Begin()),
          category_itr_end( GetStarTagCategories().End());
        category_itr != category_itr_end;
        ++category_itr
      )
      {
        // if the category has the same description
        if( ( *category_itr)->GetDescription() == DESCRIPTION)
        {
          // return it
          return *category_itr;
        }
      }

      // no matching description found
      return GetStarTagCategories().e_Undefined;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief get access to all star tag categories
    const StarTagCategories &GetStarTagCategories()
    {
      return StarTagCategories::GetEnums();
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace nmr

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< nmr::StarTagCategoryData, nmr::StarTagCategories>;

  } // namespace util
} // namespace bcl
