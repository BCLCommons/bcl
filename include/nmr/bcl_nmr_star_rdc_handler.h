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

#ifndef BCL_NMR_STAR_RDC_HANDLER_H_
#define BCL_NMR_STAR_RDC_HANDLER_H_

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
    //! @class StarRDCHandler
    //! @brief Handler for reading/writing Star-formatted RDCs
    //! @brief Reads in NMR-Star formatted files containing NMR restraint data.  These files are
    //!        automatically generated for NMR structures deposited to the PDB and are stored at the NMR restraints grid
    //!        operated by the BMRB.  More information regarding the file format can be found at
    //!        http://www.bmrb.wisc.edu/dictionary/3.1html/SuperGroupPage.html
    //!
    //! @see @link example_nmr_star_rdc_handler.cpp @endlink
    //! @author weinerbe
    //! @date Aug 9, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StarRDCHandler :
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

      //! @brief construct from member data / reading preferences
      //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
      StarRDCHandler
      (
        const std::string &DEFAULT_EXTENSION = ".rdc_star",
        const size_t SEQUENCE_OFFSET = 0,
        const bool &ADJUST_SIGNS = false,
        const bool &NORMALIZE = false
      );

      //! @brief Clone function
      //! @return pointer to new StarRDCHandler
      StarRDCHandler *Clone() const;

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
      //! @return stream that the restraints were read from
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

    }; // class StarRDCHandler

  } // namespace nmr

} // namespace bcl

#endif // BCL_NMR_STAR_RDC_HANDLER_H_
