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

#ifndef BCL_FOLD_MUTATE_SHEET_CYCLE_H_
#define BCL_FOLD_MUTATE_SHEET_CYCLE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSheetCycle
    //! @brief Rotates the order of strands in a sheet similar to a conveyer belt
    //! @details This mutate can take a strand and rotate the relative positions of the strands in the sheet in a certain direction
    //! of a certain number of steps.
    //! The rotations can be done with and without preserving the overall sheet geometry and orientations of individiual
    //! strands. Also the rotation can also occur in a consecutive subset of strands in the sheet
    //!
    //! @see @link example_fold_mutate_sheet_cycle.cpp @endlink
    //! @author karakam
    //! @date Mar 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSheetCycle :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

      //! boolean to whether to keep the sheet geometry intact
      bool m_PreserveSheetGeometry;

      //! boolean to whether to keep the orientations of individiual strands intact
      bool m_PreserveStrandOrientations;

      //! boolean to whether to rotate a subset
      bool m_RotateSubset;

      //! minimum and maximum number of rotations
      math::Range< size_t> m_NumberRotations;

      //! minimum and maximum number of strands in the rotation group, applicable only if subset flag is given
      math::Range< size_t> m_SubsetSize;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from applicable bools and ranges
      //! @param PRESERVE_SHEET_GEOMETRY boolean to whether to keep the sheet geometry intact
      //! @param PRESERVE_STRAND_ORIENTATIONS boolean to whether to keep the orientations of individual strands intact
      //! @param ROTATE_SUBSET boolean to whether to rotate a subset
      //! @param NUMBER_ROTATIONS minimum and maximum number of rotations
      //! @param SUBSET_SIZE minimum and maximum number of strands in the rotation group
      MutateSheetCycle
      (
        const bool PRESERVE_SHEET_GEOMETRY = true,
        const bool PRESERVE_STRAND_ORIENTATIONS = true,
        const bool ROTATE_SUBSET = false,
        const math::Range< size_t> &NUMBER_ROTATIONS = ( math::Range< size_t>( 1, 2)),
        const math::Range< size_t> &SUBSET_SIZE = ( math::Range< size_t>( 2, 3))
      );

      //! @brief Clone function
      //! @return pointer to new MutateSheetCycle
      MutateSheetCycle *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get boolean to whether to keep the sheet geometry intact
      //! @return boolean to whether to keep the sheet geometry intact
      bool GetPreserveSheetGeometry() const
      {
        return m_PreserveSheetGeometry;
      }

      //! @brief get boolean to whether to keep the orientations of individiaul strands intact
      //! @return boolean to whether to keep the orientations of individiaul strands intact
      bool GetPreserveStrandOrientations() const
      {
        return m_PreserveStrandOrientations;
      }

      //! @brief return boolean to whether to rotate a subset
      //! @return boolean to whether to rotate a subset
      bool GetRotateSubset() const
      {
        return m_RotateSubset;
      }

      //! @brief returns number of rotations
      //! @return number of rotations
      const math::Range< size_t> &GetNumberRotations() const
      {
        return m_NumberRotations;
      }

      //! @brief returns size of subset
      //! @return size of subset
      const math::Range< size_t> &GetSubsetSize() const
      {
        return m_SubsetSize;
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

    }; // class MutateSheetCycle

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SHEET_CYCLE_H_ 
