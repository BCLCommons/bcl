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

#ifndef BCL_RESTRAINT_HANDLER_ACCESSIBILITY_AA_H_
#define BCL_RESTRAINT_HANDLER_ACCESSIBILITY_AA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_accessibility_profile.h"
#include "bcl_restraint_handler_base.h"
#include "bcl_restraint_interface.h"
#include "assemble/bcl_assemble_aa_exposure_interface.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerAccessibilityAA
    //! @brief HandlerAccessibilityAA is for creating restraint::Accessibilities from a file.
    //!
    //! @see @link example_restraint_handler_accessibility_aa.cpp @endlink
    //! @author alexanns
    //! @date 11/16/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API HandlerAccessibilityAA :
      public HandlerBase< AccessibilityProfile>
    {
    private:

    //////////
    // data //
    //////////

      //! s_FileHeader is what should be at the first line of the restraint file
      static const std::string s_FileHeader;

      //! s_LineFormat is the format that a line should follow to be read in correctly
      static const std::string s_LineFormat;

      //! s_MinimumEntriesPerLine is the minimum number of columns a line can have and still be a real restraint
      static const size_t s_MinimumEntriesPerLine;

      //! "s_ChainIDColumn" is the column of the line where the chain id should be given
      static const size_t s_ChainIDColumn;

      //! "s_AASeqIDColumn" is the column of the line where the sequence id of the aa should be given
      static const size_t s_AASeqIDColumn;

      //! "s_FirstEnvironmentColumn" is the column of the line where the first environment should be given
      static const size_t s_FirstEnvironmentColumn;

      //! "s_FirstMeasurementColumn" is the column of the line where the first measurement should be given
      static const size_t s_FirstMeasurementColumn;

      //! the method to use for calculating exposure of residues from structure
      util::ShPtr< assemble::AAExposureInterface> m_ExposureCalculator;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      HandlerAccessibilityAA();

      //! @brief constructor taking member variables
      //! @param EXPOSURE_CALCULATOR the method to use for calculating exposure of residues from structure
      HandlerAccessibilityAA( const util::ShPtr< assemble::AAExposureInterface> &EXPOSURE_CALCULATOR);

      //! @brief virtual copy constructor
      HandlerAccessibilityAA *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief CreateRestraints is the function which creates the AccessibilityProfile from an istream
      //! @param ISTREAM is the istream from which the AccessibilityProfile will be created
      //! @return returns a AccessibilityProfile created from ISTREAM
      AccessibilityProfile ReadRestraints( std::istream &ISTREAM) const;

      //! @brief writes restraints to a stream
      //! @param OSTREAM the stream the restraints will be written to
      //! @param RESTRAINT the restraint that will be written to the stream
      //! @return stream the restraints were written to
      static std::ostream &WriteRestraints( std::ostream &OSTREAM, const AccessibilityProfile &PROFILE);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read restraint from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write restraint to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class HandlerAccessibilityAA

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_HANDLER_ACCESSIBILITY_AA_H_
