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

#ifndef BCL_FOLD_MUTATE_SHEET_DIVIDE_H_
#define BCL_FOLD_MUTATE_SHEET_DIVIDE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_mutate_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSheetDivide
    //! @brief divides a given Sheet at a determined location and moves the newly formed sub-sheet.
    //! @details This mutates finds a breaking pointer in a given beta-sheet using its member values and then divides the given
    //! beta-sheet at the given location and then moves one of the newly formed sub-beta-sheets away from the other
    //! one
    //!
    //! @see @link example_fold_mutate_sheet_divide.cpp @endlink
    //! @author karakam
    //! @date May 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSheetDivide :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

      //! minimum sheet size for to be divided
      size_t m_MinSheetSize;

      //! minimum size of the divided sheet
      size_t m_MinDividedSheetSize;

      //! boolean to whether to form a beta sandwich
      bool m_FormBetaSandwich;

      //! minimum translations
      linal::Vector3D m_MinTranslations;

      //! maximum translations
      linal::Vector3D m_MaxTranslations;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from min sheet and divided sheet sizes, and min and max translations
      //! @param MIN_SHEET_SIZE minimum number of strands in the sheet to be divided
      //! @param MIN_DIVIDED_SHEET_SIZE minimum number of strand in the divided sheet
      //! @param FORM_BETA_SANDWICH whether to form a beta sandwich
      //! @param MIN_TRANSLATIONS minimum translation along x,y and z axis
      //! @param MAX_TRANSLATIONS maximum translation along x,y and z axis
      MutateSheetDivide
      (
        const size_t MIN_SHEET_SIZE,
        const size_t MIN_DIVIDED_SHEET_SIZE,
        const bool FORM_BETA_SANDWICH,
        const linal::Vector3D &MIN_TRANSLATIONS = linal::Vector3D( 0.0, 0.0, 0.0),
        const linal::Vector3D &MAX_TRANSLATIONS = linal::Vector3D( 4.0, 4.0, 2.0)
      );

      //! @brief Clone function
      //! @return pointer to new MutateSheetDivide
      MutateSheetDivide *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return minimum sheet size
      //! @return minimum sheet size
      size_t GetMinSheetSize() const
      {
        return m_MinSheetSize;
      }

      //! @brief return minimum divided sheet size
      //! @return minimum divided sheet size
      size_t GetMinDividedSheetSize() const
      {
        return m_MinDividedSheetSize;
      }

      //! @brief boolean to whether to form a beta sandwich
      //! @return boolean to whether to form a beta sandwich
      bool GetFormBetaSandwich() const
      {
        return m_FormBetaSandwich;
      }

      //! @brief return minimum translations
      //! @return minimum translations
      const linal::Vector3D &GetMinTranslations()
      {
        return m_MinTranslations;
      }

      //! @brief return maximum translations
      //! @return maximum translations
      const linal::Vector3D &GetMaxTranslations()
      {
        return m_MaxTranslations;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that takes a Sheet and return a mutated Sheet
      //! @param SHEET Sheet which will be mutated
      //! @return MutateResult with the mutated Sheet
      math::MutateResult< assemble::Domain> operator()( const assemble::Domain &SHEET) const;

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

    }; // class MutateSheetDivide

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SHEET_DIVIDE_H_ 
