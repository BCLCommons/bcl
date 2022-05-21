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

#ifndef BCL_FOLD_MUTATE_SHEET_TWIST_H_
#define BCL_FOLD_MUTATE_SHEET_TWIST_H_

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
    //! @class MutateSheetTwist
    //! @brief This mutates adjusts the twist angle of the given Sheet
    //! @details This Mutate sets twist angles to random values around the preferred angles
    //!
    //! @see @link example_fold_mutate_sheet_twist.cpp @endlink
    //! @author karakam
    //! @date Mar 16, 2010
    //!
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSheetTwist :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

      //! range for the given twist angles
      math::Range< double> m_AngleRange;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateSheetTwist();

      //! @brief constructor from an angle range
      //! @param ANGLE_RANGE range of changes in the twist angles in radians to be applied
      MutateSheetTwist( const math::Range< double> ANGLE_RANGE);

      //! @brief Clone function
      //! @return pointer to new MutateSheetTwist
      MutateSheetTwist *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns angle range
      //! @return angle range
      math::Range< double> GetAngleRange() const
      {
        return m_AngleRange;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

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

    }; // class MutateSheetTwist

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SHEET_TWIST_H_ 
