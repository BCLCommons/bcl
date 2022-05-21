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

#ifndef BCL_FOLD_MUTATE_SHEET_REGISTER_FIX_H_
#define BCL_FOLD_MUTATE_SHEET_REGISTER_FIX_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSheetRegisterFix
    //! @brief fixes the register of strands in the given Sheet
    //! @details This Mutate fixes the register of strands in a given Sheet by translating them to ensure proper hydrogen bonding
    //! is achieved.
    //!
    //! @see @link example_fold_mutate_sheet_register_fix.cpp @endlink
    //! @author karakam
    //! @date Mar 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSheetRegisterFix :
      public math::MutateInterface< assemble::Domain>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new MutateSheetRegisterFix
      MutateSheetRegisterFix *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief finds possible hydrogen bonding pair residues from a StrandPacking objects
      //! @param STRAND_PACK SSEGeometryPacking between two strands
      //! @param STRAND_A first strand of interest
      //! @param STRAND_B second strand of interest
      //! @return possible hydrogen bonding pair residues from a StrandPacking objects
      static linal::Vector3D CalculateTranslationForHydrogenBonding
      (
        const assemble::SSEGeometryPacking &STRAND_PACK,
        const assemble::SSE &STRAND_A,
        const assemble::SSE &STRAND_B
      );

    }; // class MutateSheetRegisterFix

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SHEET_REGISTER_FIX_H_ 
