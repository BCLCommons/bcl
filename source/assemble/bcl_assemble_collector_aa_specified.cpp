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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_collector_aa_specified.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "biol/bcl_biol_aa_base.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorAASpecified::s_Instance
    (
      GetObjectInstances().AddInstance( new CollectorAASpecified())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorAASpecified::CollectorAASpecified() :
      m_ResidueList()
    {
    }

    //! @brief construct from a filename which contains the list of residues to collect
    //! @param RESI_LIST_FILENAME the file which contains the list of residues
    CollectorAASpecified::CollectorAASpecified( const std::string &RESI_LIST_FILENAME) :
      m_ResidueList( ReadAALocators( RESI_LIST_FILENAME))
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorAASpecified
    CollectorAASpecified *CollectorAASpecified::Clone() const
    {
      return new CollectorAASpecified( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorAASpecified::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetResidueList gives the list of residue locators
    //! @return m_ResidueList which is the list of amino acid locators
    const storage::List< LocatorAA> &CollectorAASpecified::GetResidueList() const
    {
      return m_ResidueList;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &CollectorAASpecified::GetAlias() const
    {
      static const std::string s_name( "CollectorAASpecified");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorAASpecified::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects residues specified by sequence and chain id.");
      serializer.AddInitializer
      (
        "residues",
        "list of residues to collect",
        io::Serialization::GetAgent( &m_ResidueList)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Collect the specified residues from a protein model
    //! @param MODEL the protein model from which residues will be collected
    //! @return SiPtrList of residues from "m_ResidueList" which could be collected
    util::SiPtrList< const biol::AABase> CollectorAASpecified::Collect( const ProteinModel &MODEL) const
    {
      // create list that will hold all of the residues collected from "MODEL"
      util::SiPtrList< const biol::AABase> collected_aas;

      // iterate through "m_ResidueList" in order to collect the residues
      for
      (
        storage::List< LocatorAA>::const_iterator
          itr( m_ResidueList.Begin()), itr_end( m_ResidueList.End()); itr != itr_end; ++itr
      )
      {
        // try to locate the residue dentod by "itr" in "MODEL"
        const util::SiPtr< const biol::AABase> located_aa( itr->Locate( MODEL));

        // true if the residue was able to be located
        if( located_aa.IsDefined())
        {
          // add the residue to "collected_aas"
          collected_aas.PushBack( located_aa);
        }
      }

      // return the list of residues collected from "MODEL"
      return collected_aas;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CollectorAASpecified::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ResidueList, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CollectorAASpecified::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ResidueList, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief ReadAALocators reads in the information from a file needed to create a list of residue locators
    //! @param FILENAME the file from which the information will be read
    //! @return list of LocatorAAs that have been created from FILENAME
    storage::List< LocatorAA> CollectorAASpecified::ReadAALocators( const std::string &FILENAME)
    {
      // create list which will hold all of the residue locators
      storage::List< LocatorAA> resi_list;

      // create istream to read locators from "FILENAME"
      io::IFStream read;

      // open read and bind it to "FILENAME"
      io::File::MustOpenIFStream( read, FILENAME);

      // read in all the contents of "FILENAME" line by line
      while( read.good())
      {
        std::string current_line;
        std::getline( read, current_line);

        // empty lines indicates end
        if( util::TrimString( current_line).empty())
        {
          break;
        }

        std::stringstream str_stream( current_line);
        char chain_id;
        io::Serialize::Read( chain_id, str_stream);

        // read in the seq id from the current line
        int seq_id;
        io::Serialize::Read( seq_id, str_stream);

        // message the chain and seq ids
        BCL_MessageDbg( "chain " + util::Format()( chain_id) + " seqid " + util::Format()( seq_id));
        BCL_Assert( str_stream.eof(), "The file has wrong format, line does not have chain id and seqid: " + current_line);

        // create a locator out of "chain" and "seq_id" and add it to "resi_list"
        resi_list.PushBack( LocatorAA( chain_id, seq_id));
      }

      // return the list of residue locators
      return resi_list;
    }

  } // namespace assemble
} // namespace bcl
