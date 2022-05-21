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

#ifndef BCL_RESTRAINT_ACCESSIBILITY_AA_ASSIGNMENT_H_
#define BCL_RESTRAINT_ACCESSIBILITY_AA_ASSIGNMENT_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_accessibility_aa.h"
#include "assemble/bcl_assemble_aa_exposure_interface.h"
#include "assemble/bcl_assemble_aa_neighbor_list.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AccessibilityAAAssignment
    //! @brief assignment for an accessibility restraint
    //! @details contains the AANeighborList for residue of interest as well as experimental accessibilities for
    //!          comparison to the exposure of the residue
    //!
    //! @see @link example_restraint_accessibility_aa_assignment.cpp @endlink
    //! @author alexanns
    //! @date Apr 6, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AccessibilityAAAssignment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! pointer to amino acid for which this assignment is generated for
      util::SiPtr< const biol::AABase> m_AminoAcid;

      //! the calculated exposure
      double m_ExposureValue;

      //! the experimentally measured accessibilities
      storage::Map< AccessibilityAA::EnvironmentEnum, double> m_Accessibility;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AccessibilityAAAssignment();

      //! @brief constructor from member variables
      //! @param AMINO_ACID the amino acid that this assignment is for
      //! @param EXPOSURE_VALUE the calculated exposure
      //! @param EXPOSURE_CALCULATOR the method to use for calculating exposure of residues from structure
      //! @param ACCESSIBILITY the experimentally measured accessibilities
      AccessibilityAAAssignment
      (
        const util::SiPtr< const biol::AABase> AMINO_ACID,
        const double EXPOSURE_VALUE,
        const storage::Map< AccessibilityAA::EnvironmentEnum, double> &ACCESSIBILITY
      );

      //! @brief Clone function
      //! @return pointer to new AccessibilityAAAssignment
      AccessibilityAAAssignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the amino acid this assignment is for
      //! @return the amino acid this assignment is for
      const util::SiPtr< const biol::AABase> &GetAABase() const;

      //! @brief get the calculated exposure value of this residue
      //! @return the exposure value calculated for this residue
      const double GetExposureValue() const;

      //! @brief data access to the experimental accessibilities
      //! @return set to experimental accessibilities which is m_Accessibility
      const storage::Map< AccessibilityAA::EnvironmentEnum, double> &GetAccessibility() const;

      //! @brief provides the desired accessibility based on a given environment type
      //! @return pair of bool and Accessibility where the bool indicates if that environment type exists and
      //!         the Accessibility is the experimentally measured value
      storage::Pair< bool, double>
      GetAccessibilityByEnvironment( const AccessibilityAA::EnvironmentType &ENVIRONMENT) const;

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

    }; // class AccessibilityAAAssignment

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ACCESSIBILITY_AA_ASSIGNMENT_H_
