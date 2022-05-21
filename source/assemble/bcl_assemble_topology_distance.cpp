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
#include "assemble/bcl_assemble_topology_distance.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> TopologyDistance::s_Instance
    (
      GetObjectInstances().AddInstance( new TopologyDistance())
    );

    //! @brief returns default distance cutoff
    //! @return default distance cutoff
    double TopologyDistance::GetDefaultDistanceCutoff()
    {
      return 5.0;
    }

    //! @brief returns default angle cutoff
    //! @return default angle cutoff
    double TopologyDistance::GetDefaultAngleCutoff()
    {
      return math::g_Pi / 2.0;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    TopologyDistance::TopologyDistance() :
      m_DistanceCutoff( GetDefaultDistanceCutoff()),
      m_AngleCutoff( GetDefaultAngleCutoff())
    {
    }

    //! @brief constructor from a distance cutoff and angle cutoff
    //! @param DISTANCE_CUTOFF distance cutoff
    //! @param ANGLE_CUTOFF angle cutoff
    TopologyDistance::TopologyDistance( const double DISTANCE_CUTOFF, const double ANGLE_CUTOFF) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_AngleCutoff( ANGLE_CUTOFF)
    {
    }

    //! @brief Clone function
    //! @return pointer to new TopologyDistance
    TopologyDistance *TopologyDistance::Clone() const
    {
      return new TopologyDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &TopologyDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns distance cutoff
    //! @return distance cutoff
    double TopologyDistance::GetDistanceCutoff() const
    {
      return m_DistanceCutoff;
    }

    //! @brief returns angle cutoff
    //! @return angle cutoff
    double TopologyDistance::GetAngleCutoff() const
    {
      return m_AngleCutoff;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates the distance between two given protein models
    //! @param MODEL ProteinModel for which distance to given TEMPLATE_MODEL will be calculated
    //! @param TEMPLATE_MODEL ProteinModel that will be used as a template
    double TopologyDistance::operator()( const ProteinModel &MODEL, const ProteinModel &TEMPLATE_MODEL) const
    {

      // static variable to hold missing SSE penalty
      static const double s_penalty_missing_penalty( 2.0);

      // calculate maximum penalty
      const double maximum_penalty
      (
        ( TEMPLATE_MODEL.GetNumberSSE( biol::GetSSTypes().HELIX) + TEMPLATE_MODEL.GetNumberSSE( biol::GetSSTypes().STRAND))
        * s_penalty_missing_penalty
      );

      // initialize distance
      double model_topology_distance( 0.0);

      // iterate over chains in model
      for
      (
        util::ShPtrVector< Chain>::const_iterator temp_chain_itr( TEMPLATE_MODEL.GetChains().Begin()),
          temp_chain_itr_end( TEMPLATE_MODEL.GetChains().End());
        temp_chain_itr != temp_chain_itr_end; ++temp_chain_itr
      )
      {
        // make a reference to this chain
        const Chain &temp_chain( **temp_chain_itr);

        // get the corresponding chain from the model
        const util::ShPtr< Chain> sp_chain( MODEL.GetChain( temp_chain.GetChainID()));

        // if the chain do not exist
        if( !sp_chain.IsDefined())
        {
          // add penalty to score
          BCL_MessageDbg( "Following Chain is missing: " + util::Format()( temp_chain.GetChainID()));
          model_topology_distance +=
            s_penalty_missing_penalty *
            ( temp_chain.GetNumberSSE( biol::GetSSTypes().HELIX) + temp_chain.GetNumberSSE( biol::GetSSTypes().STRAND));
          continue;
        }

        // get the distance between these two chains
        model_topology_distance += operator()( *sp_chain, temp_chain);
      }

      // normalize by the maximum penalty
      model_topology_distance = 100 * ( model_topology_distance / maximum_penalty);

      // end
      return model_topology_distance;
    }

    //! @brief calculates the distance between two given chains
    //! @param CHAIN Chain for which distance to given TEMPLATE_CHAIN will be calculated
    //! @param TEMPLATE_CHAIN Chain that will be used as a template
    double TopologyDistance::operator()( const Chain &CHAIN, const Chain &TEMPLATE_CHAIN) const
    {
      // initialize static penalties
      static const double penalty_missing_sse( 2.0);

      // initialize distance
      double chain_topology_distance( 0.0);

      // otherwise iterate over SSEs in the template model
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          temp_sse_itr( TEMPLATE_CHAIN.GetData().Begin()), temp_sse_itr_end( TEMPLATE_CHAIN.GetData().End());
        temp_sse_itr != temp_sse_itr_end; ++temp_sse_itr
      )
      {
        // skip if not a helix or a strand
        if( !( *temp_sse_itr)->GetType()->IsStructured())
        {
          continue;
        }

        // find this SSE in the model
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr
        (
          std::find_if( CHAIN.GetData().Begin(), CHAIN.GetData().End(), SSECompare( **temp_sse_itr))
        );

        // if such an SSE is not found
        if( sse_itr == CHAIN.GetData().End())
        {
          BCL_MessageDbg( "Following SSE is missing: " + ( *temp_sse_itr)->GetIdentification());
          chain_topology_distance += penalty_missing_sse;
          continue;
        }

        // calculate the topology distance for these two SSEs and sum it up
        chain_topology_distance += operator()( **sse_itr, **temp_sse_itr);
      }

      // end
      return chain_topology_distance;
    }

    //! @brief calculates the distance between two given SSEs
    //! @param SS_ELEMENT for which distance to given TEMPLATE_SS_ELEMENT will be calculated
    //! @param TEMPLATE_SS_ELEMENT SSE that will be used as a template
    double TopologyDistance::operator()( const SSE &SS_ELEMENT, const SSE &TEMPLATE_SS_ELEMENT) const
    {
      // initialize static penalties
      static const double penalty_dislocated_sse( 1.0);
      static const double penalty_flipped_sse( 0.5);

      // initialize distance
      double sse_topology_distance( 0.0);

      // calculate the distance between centers
      const double sse_pair_distance( linal::Distance( TEMPLATE_SS_ELEMENT.GetCenter(), SS_ELEMENT.GetCenter()));

      // if distance is more than the threshold
      if( sse_pair_distance > m_DistanceCutoff)
      {
        BCL_MessageDbg( "Following SSE is dislocated: " + SS_ELEMENT.GetIdentification());
        // add the dislocated penalty
        sse_topology_distance += penalty_dislocated_sse;
      }
      // if the distance is less than cutoff meaning SSE is in the same spot relatively
      else
      {
        // check if the SSE was flipped
        // calculate the angle deviation between SSEs
        const double sse_pair_angle
        (
          linal::ProjAngle
          (
            TEMPLATE_SS_ELEMENT.GetAxis( coord::GetAxes().e_Z), SS_ELEMENT.GetAxis( coord::GetAxes().e_Z)
          )
        );

        // if the angle deviation is larger than the angle cutoff
        if( sse_pair_angle > m_AngleCutoff)
        {
          BCL_MessageDbg( "Following SSE is flipped: " + SS_ELEMENT.GetIdentification());
          // add the flipped SSE threshold
          sse_topology_distance += penalty_flipped_sse;
        }
      }

      BCL_MessageDbg
      (
        "topology distance for SSE " + SS_ELEMENT.GetIdentification() + " is : " + util::Format()( sse_topology_distance)
      );

      // end
      return sse_topology_distance;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TopologyDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceCutoff, ISTREAM);
      io::Serialize::Read( m_AngleCutoff, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &TopologyDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AngleCutoff, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble
} // namespace bcl
