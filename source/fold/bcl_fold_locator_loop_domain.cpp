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
#include "fold/bcl_fold_locator_loop_domain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_loop_domain_c_to_n.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorLoopDomain::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorLoopDomain())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorLoopDomain::LocatorLoopDomain() :
      m_LoopSegments(),
      m_NToCSequenceDirection()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param SEGMENT_LOCATORS the loop segment locators that will be used to find the domain loop segments
    //! @param N_TO_C true if the loop domain goes from n to c - false if LocatorLoopDomainCToN should be used instead
    LocatorLoopDomain::LocatorLoopDomain
    (
      const storage::List< LocatorLoopSegment> &SEGMENT_LOCATORS,
      const bool N_TO_C
    ) :
      m_LoopSegments( SEGMENT_LOCATORS),
      m_NToCSequenceDirection( N_TO_C)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LoopSegment
    LocatorLoopDomain *LocatorLoopDomain::Clone() const
    {
      return new LocatorLoopDomain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LocatorLoopDomain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetIdentification
    //! @return GetIdentification
    const std::string LocatorLoopDomain::GetIdentification() const
    {
      std::string identification( "LocatorLoopDomain identification\n");

      for
      (
        storage::List< LocatorLoopSegment>::const_iterator itr( m_LoopSegments.Begin()), itr_end( m_LoopSegments.End());
          itr != itr_end; ++itr
      )
      {
        identification += "loop segment " + itr->GetIdentification() + '\n';
      }

      identification += " n to c direction ? " + util::Format()( m_NToCSequenceDirection);

      return identification;
    }

    //! @brief GetLoopSegments gives the list of locators for the loop segments that make up the loop domain
    //! @return the list of locators for the loop segments that make up the loop domain
    const storage::List< LocatorLoopSegment> &LocatorLoopDomain::GetLoopSegments() const
    {
      return m_LoopSegments;
    }

    //! @brief give vector of all the locator sses that correspond to this loop domain
    //! @return vector of all the locator sses that correspond to this loop domain
    util::SiPtrVector< const assemble::LocatorSSE> LocatorLoopDomain::GetLocatorSSEs() const
    {
      // to hold the locator sses of this loop domain
      util::SiPtrVector< const assemble::LocatorSSE> locator_sses;

      // iterate through the list of locator loop segments
      for
      (
          storage::List< LocatorLoopSegment>::const_iterator
            segment_itr( m_LoopSegments.Begin()), segment_itr_end( m_LoopSegments.End());
          segment_itr != segment_itr_end;
          ++segment_itr
      )
      {
        locator_sses.PushBack( segment_itr->GetLocatorSSE());
      }

      return locator_sses;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate locates the desired LoopDomain in a domain
    //! @param SSE_DOMAIN the domain in which the LoopDomain will be located
    //! @return the desired loop domain as located in SSE_DOMAIN
    util::ShPtr< LoopDomain> LocatorLoopDomain::Locate( const assemble::DomainInterface &SSE_DOMAIN) const
    {
      // create list which will hold all of the located loop segments from "PROTEIN_MODEL"
      storage::List< LoopSegment> located_loop_segments;

      // iterate through "m_LoopSegments" in order to fill "located_loop_segments"
      for
      (
        storage::List< LocatorLoopSegment>::const_iterator
          segment_itr( m_LoopSegments.Begin()), segment_itr_end( m_LoopSegments.End());
        segment_itr != segment_itr_end;
        ++segment_itr
      )
      {
        // create a LoopSegment from ProteinModel using the LocatorLoopSegment currently denoted by "segment_itr"
        const LoopSegment loop_segment( segment_itr->Locate( SSE_DOMAIN));

        // add "loop_segment" to "located_loop_segments"
        located_loop_segments.PushBack( loop_segment);
      }

      // create LoopDomain from "located_loop_segments", "m_PseudoResidue", and "located_sse"
      util::ShPtr< LoopDomain> located_domain( new LoopDomain( located_loop_segments));

      if( m_NToCSequenceDirection)
      {
        located_domain = util::ShPtr< LoopDomain>
        (
          new LoopDomain( located_loop_segments)
        );
      }
      else
      {
        located_domain = util::ShPtr< LoopDomain>
        (
          new LoopDomainCToN( located_loop_segments)
        );
      }

      // return the created loop domain
      return located_domain;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorLoopDomain::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_LoopSegments, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorLoopDomain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_LoopSegments, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates the sum of the square of the distances between target and moving points
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param LOCATOR_LOOP_DOMAIN Loop domain locator
    //! @return calculates the sum of the square of the distances between target and moving points
    double LocatorLoopDomain::CalculateSquareDistanceSum
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const LocatorLoopDomain &LOCATOR_LOOP_DOMAIN
    )
    {
      util::ShPtr< LoopDomain> sp_loop_domain( LOCATOR_LOOP_DOMAIN.Locate( PROTEIN_MODEL));
      return CalculateSquareDistanceSum( sp_loop_domain->TargetAndMovingPointsForCCD( PROTEIN_MODEL));
    }

    //! @brief calculates the root mean square distance between target and moving points
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param LOCATOR_LOOP_DOMAIN Loop domain locator
    //! @return the RMS of the points to be superimposed
    double LocatorLoopDomain::CalculateRMS
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const LocatorLoopDomain &LOCATOR_LOOP_DOMAIN
    )
    {
      util::ShPtr< LoopDomain> sp_loop_domain( LOCATOR_LOOP_DOMAIN.Locate( PROTEIN_MODEL));
      const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> points( sp_loop_domain->TargetAndMovingPointsForCCD( PROTEIN_MODEL));

      // sum of square distances
      double rms( CalculateSquareDistanceSum( points));

      // rmsd
      rms /= points.GetSize();
      rms = math::Sqrt( rms);

      // end
      return rms;
    }

    //! @brief calculates the sum of the square of the distances between target and moving points
    //! @param TARGET_AND_MOVING_POINTS List of target and moving points
    double LocatorLoopDomain::CalculateSquareDistanceSum
    (
      const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> &TARGET_AND_MOVING_POINTS
    )
    {
      // create double which will hold the sum of the square of the distances between the target and moving points
      double distance_sum( 0);

      // iterate through the target and moving points in order to calculate the square distance between each pair
      for
      (
        storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair>::const_iterator
          points_itr( TARGET_AND_MOVING_POINTS.Begin()), points_itr_end( TARGET_AND_MOVING_POINTS.End());
        points_itr != points_itr_end;
        ++points_itr
      )
      {
        BCL_MessageDbg
        (
          "current GetTargetPoint is " + util::Format()( points_itr->GetTargetPoint())
        );
        BCL_MessageDbg
        (
          "current GetMovingPoint is " + util::Format()( points_itr->GetMovingPoint())
        );
        // add the current square distance to "starting_sum_distance"
        const double sqr_dist
        (
          linal::SquareDistance( points_itr->GetTargetPoint(), points_itr->GetMovingPoint())
        );

        distance_sum += sqr_dist;

        BCL_MessageDbg( "current square distance is " + util::Format()( sqr_dist));
      }

      // message the distance sum
      BCL_MessageDbg( "total sum square distance is " + util::Format()( distance_sum));

      // return the sum of the square of the distances between the target and moving points
      return distance_sum;
    }

    //! @brief determines if a loop domain is closed or not
    //! @param DOMAIN_LOCATOR the domain that will be checked for closure
    //! @param MODEL the protein model needed to see if loop is closed or not
    //! @param RMSD_CLOSURE_THRESHOLD the threshold in rmsd Angstrom for considering the loop closed
    //! @return bool true if the loop domain is closed - false otherwise
    bool LocatorLoopDomain::IsClosed
    (
      const LocatorLoopDomain &DOMAIN_LOCATOR,
      const assemble::ProteinModel &MODEL,
      const double RMSD_CLOSURE_THRESHOLD
    )
    {
      // create double with the current sum of the square of the distances between the target and moving points
      const double rmsd( CalculateRMS( MODEL, DOMAIN_LOCATOR));

      // return false if the distance sum more than the threshold
      if( rmsd > RMSD_CLOSURE_THRESHOLD)
      {
        return false;
      }

      // otherwise return true indicating the loop is closed
      return true;
    }

  } // namespace fold
} // namespace bcl
