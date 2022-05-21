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

#ifndef BCL_NMR_ROSETTA_RDC_HANDLER_H_
#define BCL_NMR_ROSETTA_RDC_HANDLER_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_handler_base.h"
#include "restraint/bcl_restraint_rdc.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RosettaRDCHandler
    //! @brief Handler for reading/writing Rosetta-style RDC files
    //! @details http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/file_constraints.html
    //!
    //! @see @link example_nmr_rosetta_rdc_handler.cpp @endlink
    //! @author weinerbe
    //! @date Aug 9, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RosettaRDCHandler :
      public restraint::HandlerBase< restraint::RDC>
    {

    private:

    //////////
    // data //
    //////////

      //! determines if the restraint's seq ids need to be recalculated based on their relation to the pdb ids
      size_t m_SequenceOffset;

      //! whether to adjust signs of the read-in rdcs because they are scaled correctly but values are scaled correctly,
      //! but non-N containing RDCs need to switch sign
      bool m_AdjustSigns;

      //! whether to normalize rdcs to NH values
      bool m_Normalize;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from the sequence offset
      //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
      RosettaRDCHandler( const size_t SEQUENCE_OFFSET = 0);

      //! @brief Clone function
      //! @return pointer to new RosettaRDCHandler
      RosettaRDCHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns short name for the class used over the command line
      //! @return short name for the class used over the command line
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief reads restraints from a stream
      //! @param ISTREAM the stream the restraints will be read from
      //! @return read in rdcs
      restraint::RDC ReadRestraints( std::istream &ISTREAM) const;

      //! @brief writes restraints to a stream
      //! @param OSTREAM the stream the restraints will be written to
      //! @param RESTRAINT the restraint that will be written to the stream
      //! @return stream the restraints were written to
      static std::ostream &WriteRestraints( std::ostream &OSTREAM, const restraint::RDC &RESTRAINT);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class RosettaRDCHandler

  } // namespace nmr
} // namespace bcl

#endif // BCL_NMR_ROSETTA_RDC_HANDLER_H_
