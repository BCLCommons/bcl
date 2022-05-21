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
#include "nmr/bcl_nmr_star_noe_handler.h"

// includes from bcl - sorted alphabetically
#include "nmr/bcl_nmr_star_tags.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////
  // data //
  //////////

    //! defines the amount to be added to the NOE restraint for defining the upper limit
    const double StarNOEHandler::s_UpperLimit( 0.5);

    //! the lower limit for NOEs
    const double StarNOEHandler::s_LowerLimit( 1.8);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> StarNOEHandler::s_Instance
    (
      util::Enumerated< restraint::HandlerBase< util::ShPtrVector< restraint::AtomDistance> > >::AddInstance
      (
        new StarNOEHandler()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the sequence distance and offset
    //! @param SEQUENCE_DISTANCE size_t which is the smallest distance in sequence two residues can be if a restraint
    //!        is going to be stored
    //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
    StarNOEHandler::StarNOEHandler( const std::string &DEFAULT_EXTENSION, const size_t SEQUENCE_DISTANCE, const size_t SEQUENCE_OFFSET) :
      restraint::HandlerBase< util::ShPtrVector< restraint::AtomDistance> >( DEFAULT_EXTENSION),
      m_SequenceDistance( SEQUENCE_DISTANCE),
      m_SequenceOffset( SEQUENCE_OFFSET)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StarNOEHandler
    StarNOEHandler *StarNOEHandler::Clone() const
    {
      return new StarNOEHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarNOEHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &StarNOEHandler::GetAlias() const
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
    util::ShPtrVector< restraint::AtomDistance> StarNOEHandler::ReadRestraints( std::istream &ISTREAM) const
    {
      // reset the restraints
      util::ShPtrVector< restraint::AtomDistance> restraints;

      // read in the star file
      storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > line_data( StarTags::ReadStarFile( ISTREAM));

      // first need to find the corresponding atom pairs from the distance constraint category

      // create a variable to hold the atom pair information
      storage::Map< size_t, storage::VectorND< 2, storage::Triplet< char, int, std::string> > > atom_pairs
      (
        GetAtomPairData( line_data)
      );

      // now read in the distance information that corresponds to stored atom pairs

      // create an iterator to the dist constraint info in the map
      storage::Map
      <
        StarTagCategory, storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      >::const_iterator map_itr_find( line_data.Find( GetStarTagCategories().e_DistConstraintValue));

      // make sure the input file has been read in and has dist constraint data
      BCL_Assert
      (
        map_itr_find != line_data.End(),
        "File containing distance restraints was not read in or does not contain the _Dist_constraint_value category"
      );

      // iterate through the dist constraint data
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

        // initialize the distance, upper and lower bounds
        double distance;
        double upper_bound;
        double lower_bound;

        // if the distance is not in the correct place in the star file, define it by the upper bound
        if
        (
          split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceVal)) == "." ||
          split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceVal)) == "0.0"
        )
        {
          distance =
          util::ConvertStringToNumericalValue< double>
          (
            split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceUpperBoundVal))
          );
          upper_bound = distance + s_UpperLimit;

          // only give lower bound as standard if distance value is greater than 1.8
          if( distance > s_LowerLimit)
          {
            lower_bound = s_LowerLimit;
          }
          // else subtract the s_UpperLimit from the distance to get a lower bound
          else
          {
            lower_bound = distance - s_UpperLimit;
          }
        }

        // if the UpperBound value is missing, then calculate an upper and lower bound from the distance value
        else if
        (
          split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceUpperBoundVal)) == "."
        )
        {
          distance =
            util::ConvertStringToNumericalValue< double>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceVal))
            );
          upper_bound = distance + s_UpperLimit;
          // only give lower bound as standard if distance value is greater than 1.8
          if( distance > s_LowerLimit)
          {
            lower_bound = s_LowerLimit;
          }
          // else subtract the s_UpperLimit from the distance to get a lower bound
          else
          {
            lower_bound = distance - s_UpperLimit;
          }
        }
        else
        {
          distance =
            util::ConvertStringToNumericalValue< double>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceVal))
            );
          upper_bound =
            util::ConvertStringToNumericalValue< double>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceUpperBoundVal))
            );
          lower_bound =
            util::ConvertStringToNumericalValue< double>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceLowerBoundVal))
            );
        }

        // get the corresponding atom pair
        storage::VectorND< 2, storage::Triplet< char, int, std::string> > atom_pair
        (
          atom_pairs
          [
            util::ConvertStringToNumericalValue< size_t>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueConstraintID))
            )
          ]
        );

        const size_t seq_id_a( atom_pair.First().Second());
        const size_t seq_id_b( atom_pair.Second().Second());

        // create a new restraint
        util::ShPtr< restraint::AtomDistance> sp_restraint
        (
          new restraint::AtomDistance
          (
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen
              (
                atom_pair.First().First(),
                seq_id_a,
                atom_pair.First().Third()
              )
            ),
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen
              (
                atom_pair.Second().First(),
                seq_id_b,
                atom_pair.Second().Third()
              )
            ),
            util::ShPtr< restraint::Distance>
            (
              new restraint::Distance
              (
                distance,
                upper_bound,
                lower_bound
              )
            )
          )
        );

        // if the restraint has defined values
        if( sp_restraint->IsDefined())
        {
          // if the sequence distance is too small and the atoms are on the same chain
          if
          (
            atom_pair.First().First() == atom_pair.Second().First() &&
            math::Absolute( atom_pair.First().Second() - atom_pair.Second().Second()) < int( m_SequenceDistance)
          )
          {
            // print a message
            BCL_MessageStd
            (
              "Sequence distance is smaller than cutoff, " +
              util::Format()( m_SequenceDistance) + " : " + sp_restraint->GetIdentification()
            );
          }
          // the sequence distance is large enough
          else
          {
            // add it to the vector
            restraints.PushBack( sp_restraint);
          }

        }
        // the restraint is undefined
        else
        {
          // print a message
          BCL_MessageStd
          (
            "Star file contains undefined restraint: " + sp_restraint->GetIdentification()
          );
        }
      }

      // end
      return restraints;
    }

    //! @brief writes restraints to a stream
    //! @param OSTREAM the stream the restraints will be written to
    //! @param RESTRAINT the restraint that will be written to the stream
    //! @return stream the restraints were written to
    std::ostream &StarNOEHandler::WriteRestraints
    (
      std::ostream &OSTREAM,
      const util::ShPtrVector< restraint::AtomDistance> &RESTRAINT
    )
    {
      // write file header
      const std::string dist_constraint( GetStarTagCategories().e_DistConstraint->GetDescription());
      OSTREAM << "loop_" << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistTreeNodeMemberConstraintID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistConstraintTreeNodeMemberID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistAuthAsymID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistAuthSeqID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistAuthAtomID->GetDescription() << '\n';

      size_t number_restraints( 0);

      // continue until the number of restraints has been satisfied or there are no more left to chose from
      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator itr( RESTRAINT.Begin()),
          end_itr( RESTRAINT.End());
        itr != end_itr; ++itr
      )
      {
        // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_a( ( *itr)->GetData().First());
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_b( ( *itr)->GetData().Second());

        // assert that they are defined (otherwise these are not NOE restraints)
        BCL_Assert( sp_locator_a.IsDefined() && sp_locator_b.IsDefined(), "Restraints to be written are not NOEs");

        // make new line in output file
        OSTREAM << '\n';

        // write "Tree_node_member_constraint_ID"
        OSTREAM << number_restraints + 1 << ' ';
        // write "Constraint_tree_node_member_ID"
        OSTREAM << "1 ";
        // write Chain ID
        OSTREAM << sp_locator_a->GetChainID() << ' ';
        // write seq id
        OSTREAM << sp_locator_a->GetSeqID() << ' ';
        // write atom type
        OSTREAM << sp_locator_a->GetAtomTypeString();
        OSTREAM << '\n';
        // write "Tree_node_member_constraint_ID"
        OSTREAM << number_restraints + 1 << ' ';
        // write "Constraint_tree_node_member_ID"
        OSTREAM << "2 ";
        // write Chain ID
        OSTREAM << sp_locator_b->GetChainID() << ' ';
        // write seq id
        OSTREAM << sp_locator_b->GetSeqID() << ' ';
        // write atom type
        OSTREAM << sp_locator_b->GetAtomTypeString();

        ++number_restraints;
      }

      // write column headers
      const std::string dist_value( GetStarTagCategories().e_DistConstraintValue->GetDescription());
      OSTREAM << "\n\tstop_\n" << '\n';
      OSTREAM << "loop_" << '\n';
      OSTREAM << '\t' << dist_value << "." << GetStarTags().e_DistValueConstraintID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_value << "." << GetStarTags().e_DistValueDistanceVal->GetDescription() << '\n';
      OSTREAM << '\t' << dist_value << "." << GetStarTags().e_DistValueDistanceLowerBoundVal->GetDescription() << '\n';
      OSTREAM << '\t' << dist_value << "." << GetStarTags().e_DistValueDistanceUpperBoundVal->GetDescription() << '\n';

      // loop through the vectors of AtomDistances to add to the bottom of the star file
      size_t counter( 1);

      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator itr( RESTRAINT.Begin()),
          itr_end( RESTRAINT.End());
        itr != itr_end; ++itr, ++counter
      )
      {
        OSTREAM << '\n';
        OSTREAM << counter;
        OSTREAM << ' ' << ( *itr)->GetDistance()->GetDistance();
        OSTREAM << ' ' << ( *itr)->GetDistance()->LowerBound();
        OSTREAM << ' ' << ( *itr)->GetDistance()->UpperBound();
      }

      OSTREAM << "\n\tstop_\nsave_";

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StarNOEHandler::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SequenceDistance, ISTREAM);
      io::Serialize::Read( m_SequenceOffset, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StarNOEHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SequenceDistance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SequenceOffset, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gets atom pair data used for determining distance restraints
    //! @param LINE_DATA the data needed to get NOE values
    //! @return vector of atom pair information (chain id, seq_id, atom type string)
    storage::Map
    <
      size_t, storage::VectorND< 2, storage::Triplet< char, int, std::string> >
    > StarNOEHandler::GetAtomPairData
    (
      storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > &LINE_DATA
    ) const
    {
      // create an iterator to the dist constraint info in the map
      storage::Map
      <
        StarTagCategory, storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      >::const_iterator map_itr_find( LINE_DATA.Find( GetStarTagCategories().e_DistConstraint));

      // make sure the input file has been read in and has dist constraint data
      BCL_Assert
      (
        map_itr_find != LINE_DATA.End(),
        "File containing distance restraints was not read in or does not contain the _Dist_constraint category"
      );

      // create a variable to hold the atom pair information
      storage::Map< size_t, storage::VectorND< 2, storage::Triplet< char, int, std::string> > > atom_pairs;

      // create vector to hold chain ids
      storage::Map< std::string, char> chain_ids;

      // iterate through the dist constraint data
      for
      (
        storage::List< std::string>::const_iterator data_itr( map_itr_find->second.Second().Begin()),
          data_itr_end( map_itr_find->second.Second().End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // split the string by ' '
        storage::Vector< std::string> split_line_a( util::SplitString( *data_itr, " \t"));
        storage::Vector< std::string> split_line_b( util::SplitString( *( ++data_itr), " \t"));

        // get the current atom pair index and current atom
        const size_t current_index
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistTreeNodeMemberConstraintID))
          )
        );
        const size_t current_atom
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistConstraintTreeNodeMemberID))
          )
        );

        // create variable to hold the data
        storage::VectorND< 2, storage::Triplet< char, int, std::string> > atom_pair;

        // only read in the next line if it is a new atom (STAR format allows for multiple restraints from ambiguous
        // protons, such as 1   MET   HG2 vs 1   MET   HG3 that both may be in the restraint file, so only 1 needs to
        // be read in)
        while
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistTreeNodeMemberConstraintID))
          ) == current_index &&
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistConstraintTreeNodeMemberID))
          ) == current_atom
        )
        {
          ++data_itr;
          if( data_itr == data_itr_end)
          {
            break;
          }
          split_line_b = util::SplitString( *( data_itr), " \t");
        }

        // read in the chains
        const char chain_a
        (
          StarTags::GetChainID( split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistAuthAsymID)), chain_ids)
        );
        const char chain_b
        (
          StarTags::GetChainID( split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistAuthAsymID)), chain_ids)
        );

        // read in the seq ids
        const int seq_id_a
        (
          util::ConvertStringToNumericalValue< int>
          (
            split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistAuthSeqID))
          ) - int( m_SequenceOffset)
        );
        const int seq_id_b
        (
          util::ConvertStringToNumericalValue< int>
          (
            split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistAuthSeqID))
          ) - int( m_SequenceOffset)
        );

        BCL_Assert( util::IsDefined( m_SequenceDistance), "m_SequenceDistance is not defined");

        // read in the atom types as strings
        const std::string atom_a( split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistAuthAtomID)));
        const std::string atom_b( split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistAuthAtomID)));

        // add the atom pair to the vector
        storage::Triplet< char, int, std::string> atom_info_a( chain_a, seq_id_a, atom_a);
        storage::Triplet< char, int, std::string> atom_info_b( chain_b, seq_id_b, atom_b);
        atom_pair.First() = atom_info_a;
        atom_pair.Second() = atom_info_b;
        atom_pairs[ current_index] = atom_pair;

        // read through new lines if they are not part of a new atom pair
        while
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistTreeNodeMemberConstraintID))
          ) == current_index
        )
        {
          ++data_itr;
          if( data_itr == data_itr_end)
          {
            break;
          }
          split_line_b = util::SplitString( *( data_itr), " \t");
        }

        // found a new pair, so back up 1 so that it can be added in the next iteration
        --data_itr;
      }

      // end
      return atom_pairs;
    }

  } // namespace nmr
} // namespace bcl
