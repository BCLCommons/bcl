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

#ifndef BCL_SCORE_ACCESSIBILITY_HYDROPHOBIC_MOMENT_MAGNITUDE_H_
#define BCL_SCORE_ACCESSIBILITY_HYDROPHOBIC_MOMENT_MAGNITUDE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_accessibility_hydrophobic_moment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AccessibilityHydrophobicMomentMagnitude
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_score_accessibility_hydrophobic_moment_magnitude.cpp @endlink
    //! @author alexanns
    //! @date Feb 7, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AccessibilityHydrophobicMomentMagnitude :
      public AccessibilityHydrophobicMoment
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AccessibilityHydrophobicMomentMagnitude();

      //! @brief constructor taking member variable parameters
      //! @param ENVIRONMENT the type of environment the accessibility was measured in that should be scored
      //! @param WINDOW_SIZES for sstype, number of restraints included in each window the moment will be calculated for
      AccessibilityHydrophobicMomentMagnitude
      (
        const restraint::AccessibilityAA::EnvironmentType &ENVIRONMENT,
        const storage::Map< biol::SSType, size_t> WINDOW_SIZES
      );

      //! @brief Clone function
      //! @return pointer to new AccessibilityHydrophobicMomentMagnitude
      AccessibilityHydrophobicMomentMagnitude *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an accessibility profile that has been assigned and scores it
      //! @param ASSIGNMENT the profile assignment that is going to be scored
      //! @return double that is the score of the accessibility profile assignment
      double operator()( const restraint::AccessibilityProfileAssignment &ASSIGNMENT) const;

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

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
       std::ostream &WriteDetailedSchemeAndValues
      (
        const restraint::AccessibilityProfileAssignment &ASSIGNMENT,
        std::ostream &OSTREAM
      ) const;

    }; // class AccessibilityHydrophobicMomentMagnitude

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_ACCESSIBILITY_HYDROPHOBIC_MOMENT_MAGNITUDE_H_ 
