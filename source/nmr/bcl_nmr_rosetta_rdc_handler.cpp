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
#include "nmr/bcl_nmr_rosetta_rdc_handler.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RosettaRDCHandler::s_Instance
    (
      util::Enumerated< restraint::HandlerBase< restraint::RDC> >::AddInstance
      (
        new RosettaRDCHandler()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the sequence offset
    //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
    RosettaRDCHandler::RosettaRDCHandler( const size_t SEQUENCE_OFFSET) :
      m_SequenceOffset( SEQUENCE_OFFSET),
      m_AdjustSigns( false),
      m_Normalize( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RosettaRDCHandler
    RosettaRDCHandler *RosettaRDCHandler::Clone() const
    {
      return new RosettaRDCHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RosettaRDCHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &RosettaRDCHandler::GetAlias() const
    {
      static const std::string s_name( "Rosetta");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads restraints from a stream
    //! @param ISTREAM the stream the restraints will be read from
    //! @return read in rdcs
    restraint::RDC RosettaRDCHandler::ReadRestraints( std::istream &ISTREAM) const
    {
      // reset the restraints
      restraint::RDC restraints;

      // store the file
      storage::Vector< std::string> line_data( util::StringLineListFromIStream( ISTREAM));

      // iterate through the line data
      for
      (
        storage::Vector< std::string>::const_iterator line_itr( line_data.Begin()), end_line_itr( line_data.End());
          line_itr != end_line_itr; ++line_itr
      )
      {
        // split each line by tabs
        storage::Vector< std::string> split_line( util::SplitString( *line_itr, "\t"));

        if( split_line.GetSize() < 5)
        {
          break;
        }

        // set chains to A
        const char chain_a( 'A');
        const char chain_b( 'A');

        // get the 1st and 2nd seq id and adjust according to the SequenceOffset
        const size_t seq_id_a( util::ConvertStringToNumericalValue< size_t>( split_line( 0)) - m_SequenceOffset);
        const size_t seq_id_b( util::ConvertStringToNumericalValue< size_t>( split_line( 2)) - m_SequenceOffset);

        // construct the locators
        const restraint::LocatorCoordinatesHydrogen locator_a( chain_a, seq_id_a, split_line( 1));
        const restraint::LocatorCoordinatesHydrogen locator_b( chain_b, seq_id_b, split_line( 3));

        // get the rdc value
        const double rdc_value( util::ConvertStringToNumericalValue< double>( split_line( 4)));

        // get the default bond length from the atom type
        const size_t bond_length( locator_a.GetAtomType()->GetBondLength( locator_b.GetAtomType()));

        // add the restraint to the restraint container
        restraints.PushBack( locator_a, locator_b, bond_length, rdc_value);
      } // end of looping through lines

      // if the normalize to NH flag is given
      if( m_Normalize)
      {
        // normalize the restraint
        restraints.NormalizetoNH();
      }
      // if the adjust sign flag is given
      else if( m_AdjustSigns)
      {
        // adjust the sign
        restraints.AdjustSigns();
      }

      // end
      return restraints;
    }

    //! @brief writes restraints to a stream
    //! @param OSTREAM the stream the restraints will be written to
    //! @param RESTRAINT the restraint that will be written to the stream
    //! @return stream the restraints were written to
    std::ostream &RosettaRDCHandler::WriteRestraints( std::ostream &OSTREAM, const restraint::RDC &RESTRAINT)
    {
      // iterate through the RDC container
      for
      (
        storage::Vector< storage::Triplet< restraint::DataPairwise, double, double> >::const_iterator
          rdc_itr( RESTRAINT.GetData().Begin()), end_itr( RESTRAINT.GetData().End());
        rdc_itr != end_itr; ++rdc_itr
      )
      {
        // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_a( rdc_itr->First().First());
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_b( rdc_itr->First().Second());

        // assert that they are defined (otherwise these are not RDC restraints)
        BCL_Assert( sp_locator_a.IsDefined() && sp_locator_b.IsDefined(), "Restraints to be written are not RDCs");

        // seqID of first Atom
        OSTREAM << sp_locator_a->GetSeqID() << "\t";
        // type of first Atom
        OSTREAM << sp_locator_a->GetAtomTypeString() << "\t";
        // seqID of second Atom
        OSTREAM << sp_locator_b->GetSeqID() << "\t";
        // type of second Atom
        OSTREAM << sp_locator_b->GetAtomTypeString() << "\t";
        // RDC value
        OSTREAM << rdc_itr->Third() << "\n";
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RosettaRDCHandler::GetSerializer() const
    {
      io::Serializer serializer( restraint::HandlerInterface::GetSerializer());
      serializer.SetClassDescription
      (
        "Reads residual dipolar coupling restraints in Rosetta format; see "
        "http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/file_constraints.html for format specification"
      );
      serializer.AddInitializer
      (
        "offset",
        "offset of pdb ids relative to the restraints ids. This may be needed when the restraints were determined on a "
        "different sequence than is folded, defines the number to be subtracted from restraint ids to get to the "
        "sequence id ",
        io::Serialization::GetAgent( &m_SequenceOffset),
        "0"
      );
      serializer.AddInitializer
      (
        "normalize",
        "Normalize (sign and magnitude) of RDCs to NH values",
        io::Serialization::GetAgent( &m_Normalize),
        "False"
      );
      serializer.AddInitializer
      (
        "adjust signs",
        "adjust signs of normalized RDCs only. Set this if the values are scaled correctly, but non-N containing RDCs "
        "need to switch sign. Note that normalize=True will automatically adjust signs as well",
        io::Serialization::GetAgent( &m_AdjustSigns),
        "False"
      );
      return serializer;
    }

  } // namespace nmr
} // namespace bcl
