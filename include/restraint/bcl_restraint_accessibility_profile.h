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

#ifndef BCL_RESTRAINT_ACCESSIBILITY_PROFILE_H_
#define BCL_RESTRAINT_ACCESSIBILITY_PROFILE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_accessibility_aa.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AccessibilityProfile
    //! @brief restraint for multiple accessibilities for residues
    //! @details Restraint involving a series of accessibility measurements. Typically accessibility measurements are
    //!          conducted over many residues in an amino acid sequence
    //!
    //! @see @link example_restraint_accessibility_profile.cpp @endlink
    //! @author alexanns
    //! @date Apr 7, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AccessibilityProfile :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! list of accessibilities indicating a series of accessibility measurements
      storage::List< AccessibilityAA> m_Accessibilities;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AccessibilityProfile();

      //! @brief constructor taking member variable
      //! @param ACCESSIBILITIES the list of accessibilities indicating a series of accessibilit measurements
      AccessibilityProfile( const storage::List< AccessibilityAA> &ACCESSIBILITIES);

      //! @brief Clone function
      //! @return pointer to new AccessibilityProfile
      AccessibilityProfile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the list of accessibilities in this profile
      //! @return the list of accessibilities in this profile
      const storage::List< AccessibilityAA> &GetAccessibilities() const
      {
        return m_Accessibilities;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief function for creating an assignment of the profile to residues in a protein model
      //! @param PROTEIN_MODEL model from which the assignment will be created
      //! @return AccessibilityProfileAssignment which assigns the profile to the residues in PROTEIN_MODEL
      AccessibilityProfileAssignment GenerateAssignment( const assemble::ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

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

    private:

    }; // class AccessibilityProfile

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ACCESSIBILITY_PROFILE_H_ 
