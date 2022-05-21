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
#include "nmr/bcl_nmr_rosetta_noe_handler.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"
#include "score/bcl_score_restraint_nmr_distance_interface.h"
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
    const double RosettaNOEHandler::s_UpperLimit( 0.5);

    //! the lower limit for NOEs
    const double RosettaNOEHandler::s_LowerLimit( 1.8);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RosettaNOEHandler::s_Instance
    (
      util::Enumerated< restraint::HandlerBase< util::ShPtrVector< restraint::AtomDistance> > >::AddInstance
      (
        new RosettaNOEHandler()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the sequence distance and offset
    //! @param SEQUENCE_DISTANCE size_t which is the smallest distance in sequence two residues can be if a restraint
    //!        is going to be stored
    //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
    //! @param PREFIX prefix (path) for KB potential files
    //! @param SCORE_WEIGHT score weight to use
    RosettaNOEHandler::RosettaNOEHandler
    (
      const size_t SEQUENCE_DISTANCE,
      const size_t SEQUENCE_OFFSET,
      const std::string PREFIX,
      const double &SCORE_WEIGHT
    ) :
      m_SequenceDistance( SEQUENCE_DISTANCE),
      m_SequenceOffset( SEQUENCE_OFFSET),
      m_Prefix( PREFIX),
      m_Weight( SCORE_WEIGHT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RosettaNOEHandler
    RosettaNOEHandler *RosettaNOEHandler::Clone() const
    {
      return new RosettaNOEHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RosettaNOEHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &RosettaNOEHandler::GetAlias() const
    {
      static const std::string s_name( "Rosetta");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads restraints from a stream
    //! @param ISTREAM the stream the restraints will be read from
    //! @return the read in restraints
    util::ShPtrVector< restraint::AtomDistance> RosettaNOEHandler::ReadRestraints( std::istream &ISTREAM) const
    {
      // reset the restraints
      util::ShPtrVector< restraint::AtomDistance> restraints;

      // create a place to store the file
      storage::Vector< std::string> line_data( util::StringLineListFromIStream( ISTREAM));

      // iterate through the file's data
      for
      (
        storage::Vector< std::string>::const_iterator line_itr( line_data.Begin()), end_line_itr( line_data.End());
          line_itr != end_line_itr; ++line_itr
      )
      {
        // split the line
        storage::Vector< std::string> split_line( util::SplitString( *line_itr, "\t"));

        // break if the amount of data given is too small
        if( split_line.GetSize() < 10)
        {
          break;
        }

        // set chains to A
        const char chain_a( 'A');
        const char chain_b( 'A');

        // get the 1st and 2nd seq id and adjust according to the sequence offset
        const size_t seq_id_a( util::ConvertStringToNumericalValue< size_t>( split_line( 2)) - m_SequenceOffset);
        const size_t seq_id_b( util::ConvertStringToNumericalValue< size_t>( split_line( 4)) - m_SequenceOffset);

        // get the atom type strings
        const std::string atom_type_a( split_line( 1));
        const std::string atom_type_b( split_line( 3));

        // get the lower bound
        const double lower_bound( util::ConvertStringToNumericalValue< double>( split_line( 6)));

        // get the noe distance
        const double noe_file_value( util::ConvertStringToNumericalValue< double>( split_line( 7)));

        // get the actual noe value
        const double noe_value( std::max( s_LowerLimit, noe_file_value));

        // create an upper bound
        const double upper_bound( noe_value + s_UpperLimit);

        // create the distance restraint
        util::ShPtr< restraint::AtomDistance> sp_restraint
        (
          new restraint::AtomDistance
          (
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen( chain_a, seq_id_a, atom_type_a)
            ),
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen( chain_b, seq_id_b, atom_type_b)
            ),
            util::ShPtr< restraint::Distance>( new restraint::Distance( noe_value, upper_bound, lower_bound))
          )
        );

        // only store if if the restraint has defined values
        if( sp_restraint->IsDefined())
        {
          // only store if there is one chain and the two AAs are separated by enough amino acids
          if
          (
            math::Absolute( int( seq_id_a) - int( seq_id_b)) < int( m_SequenceDistance)
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
        else
        {
          // print a message
          BCL_MessageStd
          (
            "Rosetta file contains undefined restraint: " + sp_restraint->GetIdentification()
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
    std::ostream &RosettaNOEHandler::WriteRestraints
    (
      std::ostream &OSTREAM,
      const util::ShPtrVector< restraint::AtomDistance> &RESTRAINT
    ) const
    {
      // iterate through the vector of restraints
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

        // print the type of ROSETTA restraint this is
        OSTREAM << "AtomPair\t";

        // if no prefix is given, use bounded potential
        if( m_Prefix == "")
        {
          // Give 1st Atom Type
          OSTREAM << sp_locator_a->GetAtomTypeString() << "\t";
          // Give 1st Atom SeqID
          OSTREAM << sp_locator_a->GetSeqID() << "\t";
          // Give 2nd Atom Type
          OSTREAM << sp_locator_b->GetAtomTypeString() << "\t";
          // Give 2nd Atom SeqID
          OSTREAM << sp_locator_b->GetSeqID() << "\t";
          // Required lines which tell how the restraint will be calculated (0) is the lower bound
          OSTREAM << "BOUNDED\t0\t";
          // and then the distance
          OSTREAM << ( *itr)->GetDistance()->GetDistance() << "\t";
          // tell that it's an NOE
          OSTREAM << util::Format()( m_Weight) << "\tNOE\n";
        }
        // use KB potentials
        else
        {
          // store the current atom types
          biol::AtomType atom_a( GetBaseAtomType( sp_locator_a->GetAtomTypeFromString()));
          biol::AtomType atom_b( GetBaseAtomType( sp_locator_b->GetAtomTypeFromString()));

          // get the total bonds to CB
          const size_t nr_bonds
          (
            score::RestraintNMRDistanceInterface::GetBondsFromCB( atom_a) +
            score::RestraintNMRDistanceInterface::GetBondsFromCB( atom_b)
          );

          // set the string
          std::string atom_string( "CB");
          if( atom_a == biol::GetAtomTypes().H || atom_b == biol::GetAtomTypes().H)
          {
            atom_string = "H";
          }
          if
          (
            atom_a == biol::GetAtomTypes().HA || atom_a == biol::GetAtomTypes().HA2 || atom_a == biol::GetAtomTypes().HA3 ||
            atom_b == biol::GetAtomTypes().HA || atom_b == biol::GetAtomTypes().HA2 || atom_b == biol::GetAtomTypes().HA3
          )
          {
            atom_string = "HA";
          }

          // handle HA2 and HA3
          std::string atom_a_name( atom_a.GetName());
          std::string atom_b_name( atom_b.GetName());
          if( atom_a == biol::GetAtomTypes().HA2)
          {
            atom_a_name = "1HA";
          }
          else if( atom_a == biol::GetAtomTypes().HA3)
          {
            atom_a_name = "2HA";
          }
          if( atom_b == biol::GetAtomTypes().HA2)
          {
            atom_b_name = "1HA";
          }
          else if( atom_b == biol::GetAtomTypes().HA3)
          {
            atom_b_name = "2HA";
          }

          // Give 1st Atom Type
          OSTREAM << atom_a_name << "\t";
          // Give 1st Atom SeqID
          OSTREAM << sp_locator_a->GetSeqID() << "\t";
          // Give 2nd Atom Type
          OSTREAM << atom_b_name << "\t";
          // Give 2nd Atom SeqID
          OSTREAM << sp_locator_b->GetSeqID() << "\t";
          // use bounded function if number of bonds is zero
          if( nr_bonds == 0)
          {
            // Required lines which tell how the restraint will be calculated (0) is the lower bound
            OSTREAM << "BOUNDED\t0\t";
            // and then the distance
            OSTREAM << ( *itr)->GetDistance()->GetDistance() << "\t";
            // tell that it's an NOE
            OSTREAM << util::Format()( m_Weight) << "\tNOE\n";
          }
          else
          {
            // output spline
            OSTREAM << "SPLINE\tdifference\t";
            // output filename
            OSTREAM << m_Prefix << atom_string << "_" << nr_bonds << ".potential" << '\t';
            // and then the distance
            OSTREAM << ( *itr)->GetDistance()->GetDistance() << "\t";
            OSTREAM << util::Format()( m_Weight) << "\t0.25\n";
          }
        }
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RosettaNOEHandler::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SequenceDistance, ISTREAM);
      io::Serialize::Read( m_SequenceOffset, ISTREAM);
      io::Serialize::Read( m_Prefix, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RosettaNOEHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SequenceDistance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SequenceOffset, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Prefix, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gets the base atom type (i.e. CB)
    //! @param ATOM_TYPE atom type to be converted
    //! @return base atom type
    biol::AtomType RosettaNOEHandler::GetBaseAtomType( const biol::AtomType &ATOM_TYPE)
    {
      // initialize atom type to CB
      biol::AtomType atom_type( biol::GetAtomTypes().CB);

      // if type is HA or H
      if
      (
        ATOM_TYPE == biol::GetAtomTypes().HA ||
        ATOM_TYPE == biol::GetAtomTypes().HA2 ||
        ATOM_TYPE == biol::GetAtomTypes().HA3 ||
        ATOM_TYPE == biol::GetAtomTypes().H
      )
      {
        atom_type = ATOM_TYPE;
      }

      // end
      return atom_type;
    }

  } // namespace nmr
} // namespace bcl
