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
#include "score/bcl_score_protein_model_membrane_topology.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_domain_sse_pool_overlapping.h"
#include "biol/bcl_biol_atom.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "fold/bcl_fold_setup.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelMembraneTopology::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelMembraneTopology())
    );

    //! @brief gives the flag that allows the input of the filename of the pool containing the expectd tm helices
    //! @return shptr to flag with parameter for providing the expected transmembrane helices
    const util::ShPtr< command::FlagInterface> &
    ProteinModelMembraneTopology::GetFlagExpectedTransmembraneHelicesPoolFile()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "tm_helices",
          "The pool formatted file containing the expected tm helices that the protein model will be scored against",
          command::Parameter
          (
            "pool_file", "filename of the pool formatted file containing the expected helices", "tm_helices.pool"
          )
        )
      );

      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelMembraneTopology::ProteinModelMembraneTopology() :
      m_TransmembraneDomain(),
      m_SSEPool(),
      m_Scheme( GetStaticClassName< ProteinModelMembraneTopology>())
    {
    }

    //! @brief constructor taking member variables
    //! @param LOCATOR the locator that will be used to locate the transmembrane domains of the protein model
    //! @param POOL the pool of sses identifying the transmembrane sses
    //! @param SCHEME the scheme of this score
    ProteinModelMembraneTopology::ProteinModelMembraneTopology
    (
      const util::ShPtr< find::LocatorInterface< util::ShPtr< assemble::Domain>, assemble::ProteinModel> > &LOCATOR,
      const assemble::SSEPool &POOL,
      const std::string &SCHEME
    ) :
      m_TransmembraneDomain( LOCATOR),
      m_SSEPool( POOL),
      m_Scheme( SCHEME)
    {
      BCL_Assert
      (
        !m_SSEPool.IsOverlapping(), "SSE pool has overlapping SSE definitions. "
         "Please choose which unique transmembrane segments you would like to use."
      );

      BCL_Assert
      (
        m_SSEPool.GetSSEs( biol::GetSSTypes().STRAND).GetSize() == 0, "Expected TM SSEs pool should only have helices"
      );
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelMembraneTopology
    ProteinModelMembraneTopology *ProteinModelMembraneTopology::Clone() const
    {
      return new ProteinModelMembraneTopology( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelMembraneTopology::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelMembraneTopology::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief get a more readable score scheme
    //! @return a more readable score scheme
    const std::string &ProteinModelMembraneTopology::GetReadableScheme() const
    {
      static const std::string s_readable_scheme( "MP topology");
      return s_readable_scheme;
    }

    //! @brief get score type
    //! @return score type
    ProteinModel::Type ProteinModelMembraneTopology::GetType() const
    {
      return ProteinModel::e_Structure;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize this score from the flag
    bool ProteinModelMembraneTopology::InitializeFromFlag()
    {
      if( !GetFlagExpectedTransmembraneHelicesPoolFile()->GetFlag())
      {
        return false;
      }
      io::IFStream read;
      std::string pool_filename
      (
        GetFlagExpectedTransmembraneHelicesPoolFile()->GetFirstParameter()->GetValue()
      );
      io::File::MustOpenIFStream( read, pool_filename);
      m_SSEPool.ReadSSEPool( read, *fold::GetSetup().GetEmptyModel(), 9, 3);
      io::File::CloseClearFStream( read);
      // make a domain locator
      m_TransmembraneDomain = util::ShPtr< assemble::LocatorDomainSSEPoolOverlapping>
      (
        new assemble::LocatorDomainSSEPoolOverlapping( m_SSEPool)
      );
      m_Scheme = "mp_helix_topology";
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an PROTEIN_MODEL and returning a t_ResultType object
    //! @param PROTEIN_MODEL Protein Model to be used to evaluate the function
    //! @return function value of the given argument
    double ProteinModelMembraneTopology::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // locate the transmebmrane domain
      const util::ShPtr< assemble::Domain> transmembrane_domain( m_TransmembraneDomain->Locate( PROTEIN_MODEL));
      BCL_MessageDbg
      (
        "transmembrane domain size " + util::Format()( transmembrane_domain->GetSSEs().GetSize())
      );

      // to hold the tm sses coming from the sse pool and the associated sses from the located tm domain
      storage::Map
      <
        util::SiPtr< const assemble::SSE>, //< expected tm segment
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>, //< model tm sses
        assemble::SSELessThan
      > tm_assignments( AssociateModelTMSSEsWithExpectedTMSegments( *transmembrane_domain, m_SSEPool));

      double score( 0);

      // check n-termini
      {
        // to hold the side of the membrane the nterminal ca atom is on for each transmembrane sse in the model
        // and the chain the sse is in
        storage::Map< size_t, storage::Vector< storage::Pair< double, char> > > membrane_side_coordinates
        (
          GetExpectedMembraneSideAndModelSSETerminusCoordinates( tm_assignments, m_SSEPool.GetSSEs(), true)
        );

        // get the score it is not normalized; penalize ntermini that should be on the same side but are not
        score += ScoreTopology( membrane_side_coordinates);

        BCL_MessageDbg( "n-termini ScoreTopology score is " + util::Format()( score));

        // penalize ntermini that are on the same side but shouldn't be
        score += ScoreTopologyOpposingSides( membrane_side_coordinates);

        BCL_MessageDbg( "n-termini ScoreTopologyOpposingSides score is " + util::Format()( score));
      }

      // check c-termini
      {
        // to hold the side of the membrane the cterminal ca atom is on for each transmembrane sse in the model
        // and the chain the sse is in
        storage::Map< size_t, storage::Vector< storage::Pair< double, char> > > membrane_side_coordinates
        (
          GetExpectedMembraneSideAndModelSSETerminusCoordinates( tm_assignments, m_SSEPool.GetSSEs(), false)
        );

        score += ScoreTopology( membrane_side_coordinates);

        BCL_MessageDbg( "c-termini ScoreTopology score is " + util::Format()( score));

        // penalize ntermini that are on the same side but shouldn't be
        score += ScoreTopologyOpposingSides( membrane_side_coordinates);

        BCL_MessageDbg( "c-termini ScoreTopologyOpposingSides score is " + util::Format()( score));
      }

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelMembraneTopology::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_TransmembraneDomain, ISTREAM);
      io::Serialize::Read( m_SSEPool, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelMembraneTopology::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_TransmembraneDomain, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSEPool, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief associate the expected tm sses coming from the sse pool with the sses from the located tm domain
    //! @param TRANSMEMBRANE_DOMAIN
    //! @param EXPECTED_TM_SEGMENTS
    //! @return map with the expected tm sse as the key and the set of associated model tm sses as the value.
    //!         there can be multiple model sses associated with an expected tm sse since in the model the whole tm
    //!         segment might be broken
    storage::Map
    <
      util::SiPtr< const assemble::SSE>,
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>,
      assemble::SSELessThan
    > ProteinModelMembraneTopology::AssociateModelTMSSEsWithExpectedTMSegments
    (
      const assemble::Domain &TRANSMEMBRANE_DOMAIN, const assemble::SSEPool &EXPECTED_TM_SEGMENTS
    )
    {
      storage::Map
      <
        util::SiPtr< const assemble::SSE>,
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>,
        assemble::SSELessThan
      > tm_assignments;

      // iterate through the sses in the domain
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          domain_itr( TRANSMEMBRANE_DOMAIN.GetData().Begin()), domain_itr_end( TRANSMEMBRANE_DOMAIN.GetData().End());
        domain_itr != domain_itr_end;
        ++domain_itr
      )
      {
        // reference on current sse
        const assemble::SSE &current_sse( **domain_itr);

        BCL_MessageDbg( "checking domain sse " + current_sse.GetIdentification());

        // get the best matching sse out of the sse pool
        const util::SiPtr< const assemble::SSE> best_tm_match
        (
          EXPECTED_TM_SEGMENTS.FindBestMatchFromPool( current_sse, 10).First()
        );

        // true if no match could be found
        if( !best_tm_match.IsDefined())
        {
          BCL_MessageDbg( "no match could be found for sse " + current_sse.GetIdentification());
          // go to next sse
          continue;
        }

        // assign the current sse to the best tm match from the sse pool
        tm_assignments[ best_tm_match].Insert( *domain_itr);
        BCL_MessageDbg
        (
          "associated " + ( *domain_itr)->GetIdentification() + " with expected tm sse "
          + best_tm_match->GetIdentification()
        );
      }

      return tm_assignments;
    }

    //! @brief For a side of the membrane, gathers the coordinates of the n-termini of SSEs that should be on that side
    //!        The side is implicitly defined by the fact that all expected TM sses are provided, and every other SSE
    //!        should have its n-terminus on the same side
    //! @param EXPECTED_AND_MODEL_TM_SSES the assignment of an expected TM sse with corresponding sses from the model
    //! @param ALL_EXPECTED_TM_SSES the list of all expected SSEs
    //! @return map which has a size_t representing a side of the membrane, and the corresponding coordinates from
    //!         model sses assigned to expected tm sses whose n-termini should be on the same side of the membrane
    //!         The char is the chain id of the sse, because if the sses aren't in the same chain they should not be
    //!         compared together
    storage::Map< size_t, storage::Vector< storage::Pair< double, char> > >
    ProteinModelMembraneTopology::GetExpectedMembraneSideAndModelSSETerminusCoordinates
    (
      const storage::Map
      <
        util::SiPtr< const assemble::SSE>,
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>,
        assemble::SSELessThan
      > &EXPECTED_AND_MODEL_TM_SSES,
      const util::SiPtrVector< const assemble::SSE> &ALL_EXPECTED_TM_SSES,
      const bool N_TERMINUS
    )
    {
      // to hold the side of the membrane the nterminal ca atom is on for each transmembrane sse in the model
      storage::Map< size_t, storage::Vector< storage::Pair< double, char> > > membrane_side_coordinates;

      // iterate through the sse pool of tm segments
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          pool_tm_sse_itr( ALL_EXPECTED_TM_SSES.Begin()), pool_tm_sse_itr_end( ALL_EXPECTED_TM_SSES.End());
        pool_tm_sse_itr != pool_tm_sse_itr_end;
        ++pool_tm_sse_itr
      )
      {
        storage::Map
        <
          util::SiPtr< const assemble::SSE>,
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>,
          assemble::SSELessThan
        >::const_iterator assignment_itr( EXPECTED_AND_MODEL_TM_SSES.Find( *pool_tm_sse_itr));

        // true if there are no sses in the model corresponding to the current tm sse
        if( assignment_itr == EXPECTED_AND_MODEL_TM_SSES.End())
        {
          continue;
        }

        // build up set of sses from the fragments overlapping with the tm sse to ensure they are ordered by sequence
        const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan> &model_tm_segments
        (
          assignment_itr->second
        );

        // the z-coordinate of the terminus of the sse of interest
        const assemble::SSE &model_tm_segment( N_TERMINUS ? **model_tm_segments.Begin() : **( --model_tm_segments.End()));
        const double z_coordinate
        (
          N_TERMINUS ? model_tm_segment.GetFirstAA()->GetCA().GetCoordinates().Z() :
                       model_tm_segment.GetLastAA()->GetCA().GetCoordinates().Z()
        );

        BCL_MessageDbg( "current expected TM SSE is " + ( *pool_tm_sse_itr)->GetIdentification());
        BCL_MessageDbg
        (
          "current membrane side is tm sse # "
          + util::Format()( pool_tm_sse_itr - ALL_EXPECTED_TM_SSES.Begin()) + " = membrane side "
          + util::Format()( ( pool_tm_sse_itr - ALL_EXPECTED_TM_SSES.Begin()) % 2)
        );
        BCL_MessageDbg
        (
          "current model segment is "
          + model_tm_segment.GetIdentification()
        );
        BCL_MessageDbg
        (
          "current model nterminal coordinate is "
          + util::Format()( z_coordinate)
        );

        // ntermini should be on alternating sides of the membrane - arbitrarily assign even indices to one side, and
        // odd to the other. on either side, z coordinate of nterminal ca atom of each sse should have the same sign
        membrane_side_coordinates[ ( pool_tm_sse_itr - ALL_EXPECTED_TM_SSES.Begin()) % 2].PushBack
        (
          storage::Pair< double, char>( z_coordinate, model_tm_segment.GetChainID())
        );
      }

      return membrane_side_coordinates;
    }

    //! @brief scores the agreement of the actual arrangement of sses with the expected arrangement
    //! @param TOPOLOGY map that implicitly holds the expected and model topologies. The size_t indicates a side of
    //!        the membrane, the vector of doubles is the coordinates of the nterminus of sses that should be on that
    //!        side of the membrane.
    //! @return double which is the score of the agreement of the model sse topology with expected topology
    double ProteinModelMembraneTopology::ScoreTopology
    (
      const storage::Map< size_t, storage::Vector< storage::Pair< double, char> > > &TOPOLGY
    )
    {
      double score( 0);

      // count number of times z coordinates don't match sign even though they should be on the same side of membrane
      for
      (
        storage::Map< size_t, storage::Vector< storage::Pair< double, char> > >::const_iterator
          membrane_side_itr( TOPOLGY.Begin()),
          membrane_side_itr_end( TOPOLGY.End());
        membrane_side_itr != membrane_side_itr_end;
        ++membrane_side_itr
      )
      {
        // get the coordinates for this side
        const storage::Vector< storage::Pair< double, char> > &coordinates( membrane_side_itr->second);

        // iterate through the coordinates
        for
        (
          storage::Vector< storage::Pair< double, char> >::const_iterator
            coordinate_itr( coordinates.Begin()), coordinate_itr_end( coordinates.End());
          coordinate_itr != coordinate_itr_end;
          ++coordinate_itr
        )
        {
          // iterate through the coordinates
          for
          (
            storage::Vector< storage::Pair< double, char> >::const_iterator coordinate_itr_b( coordinate_itr + 1);
            coordinate_itr_b != coordinate_itr_end;
            ++coordinate_itr_b
          )
          {
            // true if the coordiantes are on opposite sides of the membrane i.e. z coordinates have opposite signs
            // and they are in the same chain
            if
            (
              ( coordinate_itr_b->First() / coordinate_itr->First()) < 0 &&
              coordinate_itr_b->Second() == coordinate_itr->Second()
            )
            {
              ++score;
            }
          }
        }
      }

      return score;
    }

    //! @brief scores the agreement of the actual arrangement of sses with the expected arrangement
    //!        Penalizes occurances of ntermini being on the same side of the membrane although they shouldn't be,
    //!        since their arrangement would force a loop to go through the membrane
    //! @param TOPOLOGY map that implicitly holds the expected and model topologies. The size_t indicates a side of
    //!        the membrane, the vector of doubles is the coordinates of the nterminus of sses that should be on that
    //!        side of the membrane.
    //! @return double which is the score of the agreement of the model sse topology with expected topology
    double ProteinModelMembraneTopology::ScoreTopologyOpposingSides
    (
      const storage::Map< size_t, storage::Vector< storage::Pair< double, char> > > &TOPOLGY
    )
    {
      double score( 0);

      // count number of times z coordinates match sign even though they shouldn't be on the same side of membrane
      for
      (
        storage::Map< size_t, storage::Vector< storage::Pair< double, char> > >::const_iterator
          membrane_side_itr( TOPOLGY.Begin()),
          membrane_side_itr_end( TOPOLGY.End());
        membrane_side_itr != membrane_side_itr_end;
        ++membrane_side_itr
      )
      {
        // get the coordinates for this side
        const storage::Vector< storage::Pair< double, char> > &coordinates( membrane_side_itr->second);

        // iterate through membrane side again
        for
        (
          storage::Map< size_t, storage::Vector< storage::Pair< double, char> > >::const_iterator
            membrane_side_itr_b
            (
              ++storage::Map< size_t, storage::Vector< storage::Pair< double, char> > >::const_iterator( membrane_side_itr)
            );
            membrane_side_itr_b != membrane_side_itr_end;
          ++membrane_side_itr_b
        )
        {
          // get the coordinates for this side
          const storage::Vector< storage::Pair< double, char> > &coordinates_b( membrane_side_itr_b->second);

          // iterate through the coordinates
          for
          (
            storage::Vector< storage::Pair< double, char> >::const_iterator
              coordinate_itr( coordinates.Begin()), coordinate_itr_end( coordinates.End());
            coordinate_itr != coordinate_itr_end;
            ++coordinate_itr
          )
          {
            // iterate through the coordinates
            for
            (
              storage::Vector< storage::Pair< double, char> >::const_iterator
                coordinate_itr_b( coordinates_b.Begin()), coordinate_itr_b_end( coordinates_b.End());
              coordinate_itr_b != coordinate_itr_b_end;
              ++coordinate_itr_b
            )
            {
              // true if the coordinates are on same sides of the membrane i.e. z coordinates have same signs
              // and they have the same chain id
              if
              (
                ( coordinate_itr_b->First() / coordinate_itr->First()) > 0 &&
                coordinate_itr_b->Second() == coordinate_itr->Second()
              )
              {
                ++score;
              }
            }
          }
        }
      }

      return score;
    }

  } // namespace score
} // namespace bcl
