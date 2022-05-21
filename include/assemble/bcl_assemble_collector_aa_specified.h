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

#ifndef BCL_ASSEMBLE_COLLECTOR_AA_SPECIFIED_H_
#define BCL_ASSEMBLE_COLLECTOR_AA_SPECIFIED_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorAASpecified
    //! @brief Class for collecting a set of residues specified by seq id and chain which are provided from a
    //!        file, list, etc.
    //! @details It can be constructed from a file which contains a list of residues and their chain. The file should be
    //!          formatted as
    //!          <chain id> <resi seq id>
    //!          with each residue on a new line. An example is
    //!          'A' 23
    //!          'D' 12
    //!
    //! @see @link example_assemble_collector_aa_specified.cpp @endlink
    //! @author alexanns
    //! @date Oct 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorAASpecified :
      public find::CollectorInterface< util::SiPtrList< const biol::AABase>, ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! the list of residues that will be collected
      storage::List< LocatorAA> m_ResidueList;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CollectorAASpecified();

      //! @brief construct from a filename which contains the list of residues to collect
      //! @param RESI_LIST_FILENAME the file which contains the list of residues
      CollectorAASpecified( const std::string &RESI_LIST_FILENAME);

      //! @brief Clone function
      //! @return pointer to new CollectorAASpecified
      CollectorAASpecified *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetResidueList gives the list of residue locators
      //! @return m_ResidueList which is the list of amino acid locators
      const storage::List< LocatorAA> &GetResidueList() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Collect the specified residues from a protein model
      //! @param MODEL the protein model from which residues will be collected
      //! @return SiPtrList of residues from "m_ResidueList" which could be collected
      virtual util::SiPtrList< const biol::AABase> Collect( const ProteinModel &MODEL) const;

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

    public:

      //! @brief ReadAALocators reads in the information from a file needed to create a list of residue locators
      //! @param FILENAME the file from which the information will be read
      //! @return list of LocatorAAs that have been created from FILENAME
      static storage::List< LocatorAA> ReadAALocators( const std::string &FILENAME);

    }; // class CollectorAASpecified

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_COLLECTOR_AA_SPECIFIED_H_ 
