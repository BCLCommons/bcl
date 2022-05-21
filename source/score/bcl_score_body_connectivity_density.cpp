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
#include "score/bcl_score_body_connectivity_density.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BodyConnectivityDensity::s_Instance
    (
      GetObjectInstances().AddInstance( new BodyConnectivityDensity())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &BodyConnectivityDensity::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "connectivity");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BodyConnectivityDensity::BodyConnectivityDensity() :
      m_Scores()
    {
    }

    //! @brief constructor taking a storage::Map and  ShPtr< restraint::Body>
    //! @param CONNECTIVITY_SCORES the storage::Map which will be "m_Scores"
    //! @param BODY_RESTRAINT util::ShPtr< restraint::Body> which will be "m_BodyRestraint"
    //! @param SCHEME
    BodyConnectivityDensity::BodyConnectivityDensity
    (
      const storage::Map< density::Connectivity, double, density::Connectivity::LessThan> &
        CONNECTIVITY_SCORES,
      const util::ShPtr< restraint::Body> &BODY_RESTRAINT,
      const std::string &SCHEME
    ) :
      m_Scores( CONNECTIVITY_SCORES),
      m_BodyRestraint( BODY_RESTRAINT),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking a restraint::Body and a Density Map
    //!        This constructor calculates the scores that go into "m_Scores" and fills "m_Scores"
    //! @param BODY_RESTRAINT contains the information about which bodies densities will be calculated between and
    //!        which will be "m_BodyRestraint"
    //! @param DENSITY_MAP experimental density map which is used to calculate the minimal density between each of
    //!        the bodies in "BODY_RESTRAINT"
    //! @param SCHEME
    BodyConnectivityDensity::BodyConnectivityDensity
    (
      const util::ShPtr< restraint::Body> &BODY_RESTRAINT,
      const util::ShPtr< density::Map> &DENSITY_MAP,
      const std::string &SCHEME
    ) :
      m_Scores( InitializeScores( BODY_RESTRAINT, DENSITY_MAP)),
      m_BodyRestraint( BODY_RESTRAINT),
      m_Scheme( SCHEME)
    {
      for
      (
        storage::Map< density::Connectivity, double, density::Connectivity::LessThan>::const_iterator
          list_itr( m_Scores.Begin()), list_itr_end( m_Scores.End());
        list_itr != list_itr_end; ++list_itr
      )
      {
        BCL_MessageDbg
        (
          "score found in map: " + util::Format()( list_itr->second) +
            " , distance: " + util::Format()( list_itr->first.GetDistance()) +
            " , connectivity: " + util::Format()( list_itr->first.GetConnectivity())
        );
      }

    }

    //! @brief virtual copy constructor
    BodyConnectivityDensity *BodyConnectivityDensity::Clone() const
    {
      return new BodyConnectivityDensity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BodyConnectivityDensity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &BodyConnectivityDensity::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetScores gives a const reference to "m_Scores"
    //! @return a const reference to "m_Scores"
    const storage::Map< density::Connectivity, double, density::Connectivity::LessThan> &
    BodyConnectivityDensity::GetScores() const
    {
      return m_Scores;
    }

    //! @brief GetBodyRestraint returns a const reference to "m_BodyRestraint"
    //! @return a const reference to "m_BodyRestraint"
    const util::ShPtr< restraint::Body> &BodyConnectivityDensity::GetBodyRestraint() const
    {
      return m_BodyRestraint;
    }

  ///////////////
  // operators //
  ///////////////

    // use GetNeighborSSE function to get neighboring SSEs!!!!!!!!!!!!

    //! @brief operator() which takes an Assignment and returns the density connectivity score.
    //!        The assignment contains all the body restraints and the corresponding protein model bodies.
    //!        For each pair of restraint bodies which are related by sequence (i.e. would be connected by a loop)
    //!        according to the protein model bodies which fill them,
    //!        the two points between which the density needs to be calculated are determined. The score for the
    //!        density between the two points is then determined.
    //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the connectivity
    //! @return return a double which is the score of the agreement of the ProteinModel with the connectivity
    double BodyConnectivityDensity::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // create double "connectivity_score" and initialize with zero
      double connectivity_score( 0.0);

      // iterate over chains of "PROTEIN_MODEL"
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if the number of SSEs in current chain is smaller than two, skip considering this chain
        if( ( *chain_itr)->GetData().GetSize() < 2)
        {
          BCL_MessageVrb( "Skip chain for BodyConnectivityDensity score, SSE count < 2");
          continue;
        }

        // Iterate over sse neighbors in the Set of SSEs which means they are also neighbors in sequence
        // i.e. sse1 and sse2; sse2 and sse3; sse3 and sse4... So that the connection points between the two SSEs
        // adjacent in sequence can be gotten and the density connectivity score between those two connection points
        // can be determined. No discrimination is taken to automatically exclude
        // connections between sses which are, in the protein model, adjacent in sequence, but, in space, very
        // far apart. The scores are left to take care of this.
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_iter_a( ( *chain_itr)->GetData().Begin()),
            sse_iter_b( ++( *chain_itr)->GetData().Begin()),
            sse_itr_end( ( *chain_itr)->GetData().End());
          sse_iter_b != sse_itr_end;
          ++sse_iter_a, ++sse_iter_b
        )
        {
          // create SiPtr to SSE "sse_a" initialize with the SSE denoted by "sse_iter_a"
          const util::SiPtr< const assemble::SSE> sse_a( **sse_iter_a);

          // create SiPtr to SSE "sse_b" initialize with the SSE denoted by "sse_iter_b"
          const util::SiPtr< const assemble::SSE> sse_b( **sse_iter_b);

          // in agreement with the implementation in the old BCL, only calculate scores for neighboring sses that
          // are not further than 10 aa in sequence apart
          if( sse_b->GetFirstAA()->GetSeqID() - sse_a->GetLastAA()->GetSeqID() > 10)
          {
            BCL_MessageVrb( "Skip this pair of SSEs, distance > 10AA");
            continue;
          }
          // create const VectorND "connection_points" initialize with the connection points between
          // the SSEs denoted by "sse_iter_a" and "sse_iter_b"
          const storage::Pair< storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> >, storage::VectorND< 2, bool> >
          restraint_bodies_and_orientation
          (
            GetRestraintBodiesandOrientations
            (
              storage::VectorND< 2, util::SiPtr< const assemble::SSE> >( sse_a, sse_b)
            )
          );

          if( !restraint_bodies_and_orientation.First().First().IsDefined() || !restraint_bodies_and_orientation.First().Second().IsDefined())
          {
            BCL_MessageVrb( "Skip this pair of SSEs, at least one does not occupy a body");
            continue;
          }

          if( restraint_bodies_and_orientation.First().First() == restraint_bodies_and_orientation.First().Second())
          {
            BCL_MessageVrb( "Skip this pair of SSEs, both occupy the same body");
            continue;
          }

          // create const double "current_score" and initialize with the score of the connectivity between
          // the points specified by restraint_bodies_and_orientation
          const double current_score( CalculateScore( restraint_bodies_and_orientation));

          // add "current_score" to "connectivity_score"
          connectivity_score += current_score;
        }
      }

      // return "connectivity_score" which is the total score of all the connectivities in "PROTEIN_MODEL"
      return connectivity_score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BodyConnectivityDensity::Read( std::istream &ISTREAM)
    {
      // read in members
      io::Serialize::Read( m_Scores, ISTREAM);
      io::Serialize::Read( m_BodyRestraint, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &BodyConnectivityDensity::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scores, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BodyRestraint, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief InitializeScores initializes "m_Scores" upon construction
    //! @param BODY_RESTRAINT contains the information about which bodies densities will be calculated between and
    //!        which will be "m_BodyRestraint"
    //! @param DENSITY_MAP experimental density map which is used to calculate the minimal density between each of
    //!        the bodies in "BODY_RESTRAINT"
    //! @return storage::Map which will be "m_Scores"
    storage::Map< density::Connectivity, double, density::Connectivity::LessThan>
    BodyConnectivityDensity::InitializeScores
    (
      const util::ShPtr< restraint::Body> &BODY_RESTRAINT,
      const util::ShPtr< density::Map> &DENSITY_MAP
    ) const
    {
      const util::ShPtrList< assemble::SSEGeometryInterface> body_list
      (
        BODY_RESTRAINT->GetBody()->Begin(), BODY_RESTRAINT->GetBody()->End()
      );

      //initialize storage list to hold all the DensityConnectivity objects
      storage::List< density::Connectivity> all_connectivities
      (
        density::Connectivity::DetermineConnectivities( *DENSITY_MAP, body_list)
      );

      // create storage::Map "connection_intensities" to hold the raw intensities between all the bodies of
      // "BODY_RESTRAINT" by calculating the scores for all connectivities
      storage::Map< density::Connectivity, double, density::Connectivity::LessThan> connection_intensities
      (
        ConvertIntensitiesToScores( all_connectivities)
      );

      //  return the resulting map of all connectivities and scores in density map
      return connection_intensities;
    }

    //! @brief ConvertIntensitiesToScores takes the raw intensities determined between restraint::Bodies and
    //!        converts them into a Z-score
    //! @param CONNECTIVITIES a storage::List which has all the DensityConnectivity objects to be converted to scores.
    //! @return returns a storage::Map which has the Z-scores and will be "m_Scores"
    storage::Map< density::Connectivity, double, density::Connectivity::LessThan>
    BodyConnectivityDensity::ConvertIntensitiesToScores
    (
      const storage::List< density::Connectivity> &CONNECTIVITIES
    ) const
    {
      // initialize empty RunningAverageSD< double> object for sd and mean calculation
      math::RunningAverageSD< double> mean_sd;
      // calculate the sd and mean by iterating over all connectivity intensities
      for
      (
        storage::List< density::Connectivity>::const_iterator
          intensity_iter( CONNECTIVITIES.Begin()), intensity_iter_end( CONNECTIVITIES.End());
        intensity_iter != intensity_iter_end;
        ++intensity_iter
      )
      {
        // consider the current intensity
        mean_sd += intensity_iter->GetConnectivity();
      }

      // create double "mean" and initialize with the mean of the intensities in "intensities"
      const double mean( mean_sd.GetAverage());

      // create double std_dev and initialize with the std deviation of the intensities in "intensities"
      const double std_dev( mean_sd.GetStandardDeviation());

      BCL_MessageCrt( "mean: " + util::Format()( mean) + " , stddev: " + util::Format()( std_dev));

      // create storage::Map "scores" and initialize with "CONNECTIVITIES"
      storage::Map< density::Connectivity, double, density::Connectivity::LessThan> scores;

      // iterate through "CONNECTIVITIES" to convert the raw intensities into z-scores
      for
      (
        storage::List< density::Connectivity>::const_iterator
          score_iter( CONNECTIVITIES.Begin()), score_iter_end( CONNECTIVITIES.End());
        score_iter != score_iter_end;
        ++score_iter
      )
      {
        const double current_intensity( score_iter->GetConnectivity());
        const double score( ( mean - current_intensity) / std_dev);

        // assert when insert into map fails (because the same element is already in the map)
        BCL_Assert
        (
          scores.Insert( std::pair< density::Connectivity, double>( *score_iter, score)).second,
          "pair of connectivity object and z-score could not be inserted into map\nBody A bool" +
            util::Format()( score_iter->GetOrientationA()) +
            "Body B bool" + util::Format()( score_iter->GetOrientationB()) +
            "connectivity" + util::Format()( score_iter->GetConnectivity()) +
            "distance" + util::Format()( score_iter->GetDistance()) +
            "score:" + util::Format()( score)
        );

        const storage::Map< density::Connectivity, double, density::Connectivity::LessThan>::const_iterator
        connectivity_itr( scores.Find( *score_iter));
        BCL_Assert
        (
          connectivity_itr != scores.End(),
          "could not find particular DensityConnectivity object:" + util::Format()( *score_iter)
        );
      }

      // return "scores" which has the z-scores
      return scores;
    }

    //! @brief GetRestraintBodiesandOrientations determines the two restraint bodies (i.e. density rods) occupied by
    //!        two neighboring SSEs in the protein model and at which ends of the restraint bodies the protein model
    //!        connects them (in the form of two bools)
    //! @param SSES are the two SSEs for which the two points of connection are desired
    //! @return returns a storage::Pair of vector of two bodies and vector of two bools
    storage::Pair< storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> >, storage::VectorND< 2, bool> >
    BodyConnectivityDensity::GetRestraintBodiesandOrientations
    (
      const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > &SSES
    ) const
    {
      // obtain the bodies of the density rods that sses are occupying
      const util::ShPtr< assemble::SSEGeometryInterface> body_restraint_first( m_BodyRestraint->GetOccupiedBody( *SSES.First()));
      const util::ShPtr< assemble::SSEGeometryInterface> body_restraint_second( m_BodyRestraint->GetOccupiedBody( *SSES.Second()));

      if( !body_restraint_first.IsDefined() || !body_restraint_second.IsDefined())
      {
        return storage::Pair< storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> >, storage::VectorND< 2, bool> >();
      }

      bool orientation_first, orientation_second;

      // true if the first SSE comes before the second SSE in sequence
      if( assemble::SSELessThan().operator()( *SSES.First(), *SSES.Second()))
      {
        orientation_first = DetermineSSEOrientationInBodyRestraint( SSES.First()->EndOfZ(), *body_restraint_first);
        orientation_second = DetermineSSEOrientationInBodyRestraint( SSES.Second()->BeginOfZ(), *body_restraint_second);

        // return the end of the first SSE and the beginning of the second SSE
        return storage::Pair< storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> >, storage::VectorND< 2, bool> >
        (
          storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> > ( body_restraint_first, body_restraint_second),
          storage::VectorND< 2, bool>( orientation_first, orientation_second)
        );
      }

      // the first SSE comes after the second SSE in sequence so return the beginning of the first SSE and the end
      // of the second SSE
      orientation_first = DetermineSSEOrientationInBodyRestraint( SSES.First()->BeginOfZ(), *body_restraint_first);
      orientation_second = DetermineSSEOrientationInBodyRestraint( SSES.Second()->EndOfZ(), *body_restraint_second);

      // return the beginning of the first SSE and the end of the second SSE
      return storage::Pair< storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> >, storage::VectorND< 2, bool> >
      (
        storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> > ( body_restraint_first, body_restraint_second),
        storage::VectorND< 2, bool>( orientation_first, orientation_second)
      );
    }

    //! @brief DetermineSSEOrientationInBodyRestraint is used to match a coordinate with one of the end points of
    //!        a body by checking which end the coordinate is closer to. This is used to determine which
    //!        direction an SSE is pointing in a restraining body, for example.
    //! @param SSE_POINT is the coordinate for which its position in the restraining body is desired
    //! @param BODY the body for which the placement position of SSE_POINT is desired
    //! @return returns a bool to indicate to which of the ends of "BODY" "SSE_POINT" is closest (true for beginning
    //!         of BODY, false for end of BODY)
    bool BodyConnectivityDensity::DetermineSSEOrientationInBodyRestraint
    (
      const linal::Vector3D &SSE_POINT,
      const coord::GeometryInterface &BODY
    ) const
    {
      // get the end of "BODY"
      const linal::Vector3D body_end( BODY.GetGeometries().LastElement()->GetCenter());

      // get the beginning of "BODY"
      const linal::Vector3D body_begin( BODY.GetGeometries().FirstElement()->GetCenter());

      // true if "SSE_POINT" is closer to "body_end" than "body_begin", then return "false"
      // else return "true" because "SSE_POINT" is closer to "body_begin" than "body_end"
      return !( linal::SquareDistance( SSE_POINT, body_end) < linal::SquareDistance( SSE_POINT, body_begin));
    }

    //! @brief GetScore returns the calculated Z-score connectivity score for connecting two bodies (at specific ends)
    //! @param CONNECTIVITY_INFORMATION a storage::Pair of vector of two bodies and vector of two bools (specifying
    //!        which two bodies are to be connected and at which ends of the bodies)
    //! @return returns a double which is the Z-score associated with connecting CONNECTIVITY_INFORMATION
    double BodyConnectivityDensity::CalculateScore
    (
      const storage::Pair< storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> >, storage::VectorND< 2, bool> >
      &CONNECTIVITY_INFORMATION
    ) const
    {
      // initialize distance cutoff beyond which the connection will not be scored
      const double distance_cutoff( 10.0);

      // construct half filled DensityConnectivity object (in one specific order)
      density::Connectivity temp_connectivity
      (
        CONNECTIVITY_INFORMATION.First().First(),
        CONNECTIVITY_INFORMATION.Second().First(),
        CONNECTIVITY_INFORMATION.First().Second(),
        CONNECTIVITY_INFORMATION.Second().Second()
      );

      const storage::Map< density::Connectivity, double, density::Connectivity::LessThan>::const_iterator
      connectivity_itr( m_Scores.Find( temp_connectivity));

      // if that particular temporary density connectivity object could not be found in the map
      if( connectivity_itr == m_Scores.End())
      {
        // try constructing the swapped object and finding this in the map
        density::Connectivity temp_connectivity_swapped
        (
          CONNECTIVITY_INFORMATION.First().Second(),
          CONNECTIVITY_INFORMATION.Second().Second(),
          CONNECTIVITY_INFORMATION.First().First(),
          CONNECTIVITY_INFORMATION.Second().First()
        );

        const storage::Map< density::Connectivity, double, density::Connectivity::LessThan>::const_iterator
        connectivity_itr_swapped( m_Scores.Find( temp_connectivity_swapped));

        // if the swapped object also isn't found in the map, then something is seriously wrong
        BCL_Assert
        (
          connectivity_itr_swapped != m_Scores.End(),
          "could not find particular DensityConnectivity object:" + util::Format()( temp_connectivity_swapped)
        );

        // return the score of the identified connectivity object if the two points are within a certain distance
        if( connectivity_itr_swapped->first.GetDistance() <= distance_cutoff)
        {
          return connectivity_itr_swapped->second;
        }
        // if too far away, just return 0.0
        return double( 0.0);
      }

      // return the score of the identified connectivity object if the two points are within a certain distance
      if( connectivity_itr->first.GetDistance() <= distance_cutoff)
      {
        // return the score of the identified connectivity object
        return connectivity_itr->second;
      }
      // if too far away, just return 0.0
      return double( 0.0);
    }

  } // namespace score
} // namespace bcl

