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

#ifndef BCL_FOLD_LOCATOR_UNCONNECTED_SEGMENTS_H_
#define BCL_FOLD_LOCATOR_UNCONNECTED_SEGMENTS_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_loop_parameters.h"
#include "assemble/bcl_assemble_domain_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorUnconnectedSegments
    //! @brief Locates unconnected structure segments in a protein model
    //! @detail The locator finds the segments of loop regions which don't have coordinates assigned. Segments which are
    //! beginning at the n-terminus or ending at the c-terminus are ignored. The found connecting segments are
    //! parameterized using fold::LoopParameters.
    //!
    //! @see @link example_fold_locator_unconnected_segments.cpp @endlink
    //! @author fischea
    //! @date Dec 20, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LocatorUnconnectedSegments :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorUnconnectedSegments();

      //! @brief clone function
      //! @return pointer to a new LocatorUnconnectedSegments
      LocatorUnconnectedSegments *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief locates unconnected segments in a protein model and computes a parameterization of the segment
      //! @detail two segments are considered unconnected if there is at least one residue without coordinates in
      //! between
      //! @param MODEL protein model for which to find unconnected segments
      //! @return parameterizations of unconnected segments in the given protein model
      util::ShPtrVector< LoopParameters> Locate( const assemble::ProteinModel &MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads members from a given input stream
      //! @param ISTREAM input stream to read members from
      //! @return input stream which members were read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into a given output stream
      //! @param OSTREAM output stream to write members into
      //! @param INDENT number of indentations
      //! @return output stream into which members were written
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LocatorUnconnectedSegments

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOCATOR_UNCONNECTED_SEGMENTS_H_
