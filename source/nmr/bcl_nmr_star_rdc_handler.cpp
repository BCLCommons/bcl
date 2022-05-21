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
#include "nmr/bcl_nmr_star_rdc_handler.h"

// includes from bcl - sorted alphabetically
#include "nmr/bcl_nmr_star_tags.h"
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
    const util::SiPtr< const util::ObjectInterface> StarRDCHandler::s_Instance
    (
      util::Enumerated< restraint::HandlerBase< restraint::RDC> >::AddInstance( new StarRDCHandler())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the sequence offset
    //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
    StarRDCHandler::StarRDCHandler
    (
      const std::string &DEFAULT_EXTENSION,
      const size_t SEQUENCE_OFFSET,
      const bool &ADJUST_SIGNS,
      const bool &NORMALIZE
    ) :
      restraint::HandlerBase< restraint::RDC>( DEFAULT_EXTENSION),
      m_SequenceOffset( SEQUENCE_OFFSET),
      m_AdjustSigns( ADJUST_SIGNS),
      m_Normalize( NORMALIZE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StarRDCHandler
    StarRDCHandler *StarRDCHandler::Clone() const
    {
      return new StarRDCHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarRDCHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &StarRDCHandler::GetAlias() const
    {
      static const std::string s_name( "Star");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads restraints from a stream
    //! @param ISTREAM the stream the restraints will be read from
    //! @return stream that the restraints were read from
    restraint::RDC StarRDCHandler::ReadRestraints( std::istream &ISTREAM) const
    {
      // reset the restraints
      restraint::RDC restraints;

      // read in the file.
      storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > line_data( StarTags::ReadStarFile( ISTREAM));

      // find the RDC data in the map
      storage::Map
      <
        StarTagCategory, storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      >::const_iterator map_itr_find( line_data.Find( GetStarTagCategories().e_RDCConstraint));

      // make sure the input file has been read in and has RDC data
      BCL_Assert( map_itr_find != line_data.End(), "File containing RDC restraints was not read in");

      // create vector to hold chain ids
      storage::Map< std::string, char> chain_ids;

      // iterate through the rdc data lines
      for
      (
        storage::List< std::string>::const_iterator data_itr( map_itr_find->second.Second().Begin()),
          data_itr_end( map_itr_find->second.Second().End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // split the string by ' '
        storage::Vector< std::string> split_line( util::SplitString( *data_itr, " \t"));

        // create chains
        char chain_a
        (
          StarTags::GetChainID( split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthAsymID1)), chain_ids)
        );
        char chain_b
        (
          StarTags::GetChainID( split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthAsymID2)), chain_ids)
        );

        // create seq ids
        const size_t seq_id_a
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthSeqID1))
          ) - m_SequenceOffset
        );
        const size_t seq_id_b
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthSeqID2))
          ) - m_SequenceOffset
        );

        // create rdc value
        const double rdc_value
        (
          util::ConvertStringToNumericalValue< double>
          (
            split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCVal))
          )
        );

        // create atom types
        const std::string atom_type_a_string
        (
          split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthAtomID1))
        );
        const std::string atom_type_b_string
        (
          split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthAtomID2))
        );

        // construct the locators
        const restraint::LocatorCoordinatesHydrogen locator_a( chain_a, seq_id_a, atom_type_a_string);
        const restraint::LocatorCoordinatesHydrogen locator_b( chain_b, seq_id_b, atom_type_b_string);

        // find the bond length entry
        const size_t bond_length_index( map_itr_find->second.First().Find( GetStarTags().e_RDCBondLength));
        double bond_length;

        // if it was found
        if( bond_length_index != map_itr_find->second.First().GetSize())
        {
          // set the bond length
          bond_length = util::ConvertStringToNumericalValue< double>( split_line( bond_length_index));
        }
        // if it was not found
        else
        {
          // get the default bond length from the atom type
          bond_length = locator_a.GetAtomType()->GetBondLength( locator_b.GetAtomType());
          if( !util::IsDefined( bond_length))
          {
            continue;
          }
        }

        // add the restraint
        restraints.PushBack( locator_a, locator_b, bond_length, rdc_value);
      }

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
    std::ostream &StarRDCHandler::WriteRestraints( std::ostream &OSTREAM, const restraint::RDC &RESTRAINT)
    {
      // write file header
      const std::string rdc_constraint( GetStarTagCategories().e_RDCConstraint->GetDescription());
      OSTREAM << "loop_" << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCID->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCVal->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCLowerBound->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCUpperBound->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthAsymID1->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthSeqID1->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthAtomID1->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthAsymID2->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthSeqID2->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthAtomID2->GetDescription() << '\n';

      // initialize a counter for writing out the NOE number
      int counter( 0);

      // iterate over the data
      for
      (
        storage::Vector< storage::Triplet< restraint::DataPairwise, double, double> >::const_iterator
          rdc_itr( RESTRAINT.GetData().Begin()), end_itr( RESTRAINT.GetData().End());
        rdc_itr != end_itr; ++rdc_itr, ++counter
      )
      {
        // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_a( rdc_itr->First().First());
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_b( rdc_itr->First().Second());

        // assert that they are defined (otherwise these are not RDC restraints)
        BCL_Assert( sp_locator_a.IsDefined() && sp_locator_b.IsDefined(), "Restraints to be written are not RDCs");

        // print out data to file
        OSTREAM << '\n' << counter + 1 << ' ' << util::Format().FFP( 1)( rdc_itr->Third()) << ' '; // counter and val
        OSTREAM << util::Format().FFP( 1)( rdc_itr->Third() - 1) << ' '; // lower bound
        OSTREAM << util::Format().FFP( 1)( rdc_itr->Third() + 1) << ' '; // upper bound
        OSTREAM << sp_locator_a->GetChainID() << ' '; // chain id 1
        OSTREAM << sp_locator_a->GetSeqID() << ' '; // seq id 1
        OSTREAM << sp_locator_a->GetAtomTypeString() << ' '; // atom type 1
        OSTREAM << sp_locator_b->GetChainID() << ' '; // chain id 2
        OSTREAM << sp_locator_b->GetSeqID() << ' '; // seq id 2
        OSTREAM << sp_locator_b->GetAtomTypeString(); // atom type 2
      }

      OSTREAM << "\n\tstop_\nsave_";

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer StarRDCHandler::GetSerializer() const
    {
      io::Serializer serializer( restraint::HandlerInterface::GetSerializer());
      serializer.SetClassDescription
      (
        "Reads residual dipolar coupling restraints in STAR format; see "
        "http://www.bmrb.wisc.edu/dictionary/3.1html/SuperGroupPage.html for format specification"
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
