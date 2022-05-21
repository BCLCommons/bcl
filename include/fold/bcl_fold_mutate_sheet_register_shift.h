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

#ifndef BCL_FOLD_MUTATE_SHEET_REGISTER_SHIFT_H_
#define BCL_FOLD_MUTATE_SHEET_REGISTER_SHIFT_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSheetRegisterShift
    //! @brief shifts the register of strands in the given Sheet
    //! @details This Mutate shifts the register of strands in a given Sheet by translating them and ensuring hydrogen bonding
    //! is not broken
    //!
    //! @see @link example_fold_mutate_sheet_register_shift.cpp @endlink
    //! @author karakam
    //! @date Apr 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSheetRegisterShift :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

      //! probabilities of regular shift and flip-shift
      storage::VectorND< 2, double> m_ShiftProbabilities;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief static function to return default shift probabilities
      //! @return default shift probabilities
      static const storage::VectorND< 2, double> &GetDefaultShiftProbabilities();

      //! @brief static function to return CA z-translation for a single residue
      //! @return CA z-translation for a single residue
      static double GetCAShiftAlongZAxis();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a shift probabilities vector
      //! @param SHIFT_PROBABILITIES vector that contains probabilities from regular and flip shift
      MutateSheetRegisterShift
      (
        const storage::VectorND< 2, double> &SHIFT_PROBABILITIES = GetDefaultShiftProbabilities()
      );

      //! @brief Clone function
      //! @return pointer to new MutateSheetRegisterShift
      MutateSheetRegisterShift *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return shift probabilities
      //! @return shift probabilities
      storage::VectorND< 2, double> GetShiftProbabilities() const
      {
        return m_ShiftProbabilities;
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

    }; // class MutateSheetRegisterShift

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SHEET_REGISTER_SHIFT_H_ 
