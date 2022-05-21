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
#include "score/bcl_score_data_set_pairwise_coordinate_exclusion.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseCoordinateExclusion::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseCoordinateExclusion())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseCoordinateExclusion::GetDefaultScheme()
    {
      static const std::string s_scheme( "exclusion");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseCoordinateExclusion::DataSetPairwiseCoordinateExclusion( const std::string &SCHEME) :
      m_Radius(),
      m_DistanceMap(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param EXCLUSION_RADIUS data with atoms inside this radius of any exclusion coordinates will be penalized
    //! @param READ the stream the coordinates will be read from
    //! @param X_COORD_COLUMN column in istream that has the x-coordinate - columns start at 0
    //! @param Y_COORD_COLUMN column in istream that has the y-coordinate - columns start at 0
    //! @param Z_COORD_COLUMN column in istream that has the z-coordinate - columns start at 0
    //! @param ALL_POSSIBLE_DATA_POINTS the set of data points that should be subjected to this filter
    //! @param ENSEMBLE ensemble for which coordinates will be checked to make sure they aren't near exclusion coords
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseCoordinateExclusion::DataSetPairwiseCoordinateExclusion
    (
      const double &EXCLUSION_RADIUS,
      std::istream &READ,
      const size_t X_COORD_COLUMN,
      const size_t Y_COORD_COLUMN,
      const size_t Z_COORD_COLUMN,
      const storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > &ALL_POSSIBLE_DATA_POINTS,
      const assemble::ProteinEnsemble &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_Radius( EXCLUSION_RADIUS),
      m_DistanceMap(),
      m_Scheme( SCHEME)
    {
      m_DistanceMap =
        FillDistanceMap( READ, X_COORD_COLUMN, Y_COORD_COLUMN, Z_COORD_COLUMN, ALL_POSSIBLE_DATA_POINTS, ENSEMBLE);
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseCoordinateExclusion
    DataSetPairwiseCoordinateExclusion *DataSetPairwiseCoordinateExclusion::Clone() const
    {
      return new DataSetPairwiseCoordinateExclusion( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseCoordinateExclusion::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseCoordinateExclusion::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score of a data set
    //! @param DATA data set to be scored
    //! @return the score of the current data set
    double DataSetPairwiseCoordinateExclusion::operator()( const restraint::DataSetPairwise &DATA) const
    {
      double score( 0);

      // iterate through the data set
      for
      (
        restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // calculate exclusion scores for each data point
        const double score_first( CalculateExclusionScore( data_itr->First(), m_DistanceMap, m_Radius));
        const double score_second( CalculateExclusionScore( data_itr->Second(), m_DistanceMap, m_Radius));

        // true if the first data point is defined
        if( util::IsDefined( score_first))
        {
          score += score_first;
        }

        // true if the second data point is defined
        if( util::IsDefined( score_second))
        {
          score += score_second;
        }
      }

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseCoordinateExclusion::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Radius, ISTREAM);
      io::Serialize::Read( m_DistanceMap, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseCoordinateExclusion::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Radius, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceMap, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief fills the distance map of undesirable coordinates by calculating the distance between a set of data
    //!        points and undesirable coordinates that are read in from an istream.
    //! @param READ the stream from which undesirable coordinates will be read
    //! @param X_COORD_COLUMN the column the x coordinate are in - starts at 0
    //! @param Y_COORD_COLUMN the column the y coordinate are in - starts at 0
    //! @param Z_COORD_COLUMN the column the z coordinate are in - starts at 0
    //! @param ALL_POSSIBLE_DATA_POINTS the data points that should be subjected to this score
    //! @param ENSEMBLE the protein models for which distances will be calculated
    //! @return distance map of data point to vector of mean sd info calculated from ENSEMBLE from that data point to
    //!         each excluded coordinate
    DataSetPairwiseCoordinateExclusion::DistanceMap DataSetPairwiseCoordinateExclusion::FillDistanceMap
    (
      std::istream &READ,
      const size_t X_COORD_COLUMN,
      const size_t Y_COORD_COLUMN,
      const size_t Z_COORD_COLUMN,
      const storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > &ALL_POSSIBLE_DATA_POINTS,
      const assemble::ProteinEnsemble &ENSEMBLE
    )
    {
      // get the exclusion coordinate from the stream
      const storage::List< linal::Vector3D> exclusion_coordinates
      (
        ReadCoordinates( READ, X_COORD_COLUMN, Y_COORD_COLUMN, Z_COORD_COLUMN)
      );

      // locate the coordinates of data points for each model
      storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, //< data points
        storage::Vector< linal::Vector3D>, //< coodinates for each model
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > coordinate_map( LocateCoordinates( ALL_POSSIBLE_DATA_POINTS, ENSEMBLE));

      return CalculateAverageDistanceToExcludedCoordinates( coordinate_map, exclusion_coordinates);
    }

    //! @brief gives list of coordinates read from an istream
    //! @param READ the stream from which undesirable coordinates will be read
    //! @param X_COORD_COLUMN the column the x coordinate are in - starts at 0
    //! @param Y_COORD_COLUMN the column the y coordinate are in - starts at 0
    //! @param Z_COORD_COLUMN the column the z coordinate are in - starts at 0
    //! @return list of Vector3D which are the coordinates read from a file
    storage::List< linal::Vector3D> DataSetPairwiseCoordinateExclusion::ReadCoordinates
    (
      std::istream &READ,
      const size_t X_COORD_COLUMN,
      const size_t Y_COORD_COLUMN,
      const size_t Z_COORD_COLUMN
    )
    {
      // get the coord list from the file
      const storage::Vector< storage::Vector< std::string> >
        coord_lines( util::SplittedStringLineListFromIStream( READ));

      storage::List< linal::Vector3D> coordinates;

      // iterate through the lines to extract coordinates and build the coordinate list
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr( coord_lines.Begin()), itr_end( coord_lines.End());
        itr != itr_end; ++itr
      )
      {
        // get x, y, z coordinates
        const double x_coord( util::ConvertStringToNumericalValue< double>( itr->operator()( X_COORD_COLUMN)));
        const double y_coord( util::ConvertStringToNumericalValue< double>( itr->operator()( Y_COORD_COLUMN)));
        const double z_coord( util::ConvertStringToNumericalValue< double>( itr->operator()( Z_COORD_COLUMN)));

        // add vector to coordinates
        const linal::Vector3D coord( x_coord, y_coord, z_coord);
        coordinates.PushBack( coord);
      }

      return coordinates;
    }

    //! @brief locates coordinates for a set of locators for an ensemble of proteins
    //! @param ALL_POSSIBLE_DATA_POINTS the data points whose coordinates will be located
    //! @param ENSEMBLE the protein models from which coordinates will be found
    //! @return map of locator to its associated vector or coordinates, one coordinate for each model in the ensemble
    storage::Map
    <
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
      storage::Vector< linal::Vector3D>,
      assemble::LocatorAtomCoordinatesInterface::PtrLessThan
    > DataSetPairwiseCoordinateExclusion::LocateCoordinates
    (
      const storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > &ALL_POSSIBLE_DATA_POINTS,
      const assemble::ProteinEnsemble &ENSEMBLE
    )
    {
      // map to hold the located coordinates of the ensemble
      storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        storage::Vector< linal::Vector3D>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > coord_map;

      // iterate through the data points
      for
      (
        storage::Set
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator data_itr( ALL_POSSIBLE_DATA_POINTS.Begin()), data_itr_end( ALL_POSSIBLE_DATA_POINTS.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // add the current locator to the distance map
        std::pair< storage::Map
          <
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
            storage::Vector< linal::Vector3D>,
            assemble::LocatorAtomCoordinatesInterface::PtrLessThan
          >::iterator, bool> insert_status
        (
          coord_map.Insert
          (
            std::pair< util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Vector< linal::Vector3D> >
            (
              *data_itr, storage::Vector< linal::Vector3D>()
            )
          )
        );

        // make sure data point can be inserted
        BCL_Assert( insert_status.second, "could not insert " + ( *data_itr)->GetIdentification());

        // iterate through the ensemble to calculate the coordinates and fill coordinate map
        for
        (
          assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
          ensemble_itr != ensemble_itr_end;
          ++ensemble_itr
        )
        {
          // locate coordinate
          const linal::Vector3D current_coords( ( *data_itr)->Locate( **ensemble_itr));

          // true if coordinate could not be found
          if( !current_coords.IsDefined())
          {
            // go to next protein model
            continue;
          }

          // add the current coords to the coordinate map
          insert_status.first->second.PushBack( current_coords);
        }
      }

      return coord_map;
    }

    //! @brief calculates the mean and sd info for data points distances to exclusion coordinates
    //! @param COORDINATE_MAP map of locator to associated vector of coords; 1 coord for each model in the ensemble
    //! @param EXCLUSION_COORDINATES list of coordinates which will be used to calculate distances to them
    //! @return distance map of data point to vector of mean sd info calculated from ENSEMBLE from that data point to
    //!         each excluded coordinate
    DataSetPairwiseCoordinateExclusion::DistanceMap
    DataSetPairwiseCoordinateExclusion::CalculateAverageDistanceToExcludedCoordinates
    (
      const storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        storage::Vector< linal::Vector3D>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > &COORDINATE_MAP,
      const storage::List< linal::Vector3D> &EXCLUSION_COORDINATES
    )
    {
      DistanceMap distance_map;

      // iterate through the coordinate map
      for
      (
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
          storage::Vector< linal::Vector3D>,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator locator_itr( COORDINATE_MAP.Begin()), locator_itr_end( COORDINATE_MAP.End());
        locator_itr != locator_itr_end;
        ++locator_itr
      )
      {
        // the current vector of coordinates for the current atom
        const storage::Vector< linal::Vector3D> &coords( locator_itr->second);

        if( coords.IsEmpty())
        {
          continue;
        }

        // add the current locator to the distance map
        std::pair< storage::Map
          <
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
            storage::Vector< math::RunningAverageSD< double> >,
            assemble::LocatorAtomCoordinatesInterface::PtrLessThan
          >::iterator, bool> insert_status
        (
          distance_map.Insert
          (
            std::pair
            <
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Vector< math::RunningAverageSD< double> >
            >
            (
              locator_itr->first, storage::Vector< math::RunningAverageSD< double> >()
            )
          )
        );

        // make sure the data point could be inserted
        BCL_Assert( insert_status.second, "could not insert " + locator_itr->first->GetIdentification());

        // iterate through the excluded coordinates
        for
        (
          storage::List< linal::Vector3D>::const_iterator
            excluded_itr( EXCLUSION_COORDINATES.Begin()), excluded_itr_end( EXCLUSION_COORDINATES.End());
          excluded_itr != excluded_itr_end;
          ++excluded_itr
        )
        {
          const linal::Vector3D &exclusion_coord( *excluded_itr);

          // true if the excluded coordinate is undefined
          if( !exclusion_coord.IsDefined())
          {
            continue;
          }

          // to keep track of ensemble average distance to current exclusion coordinate
          math::RunningAverageSD< double> mean_sd;

          // iterate through the coordinates of the current locator of the ensemble to calculate average distance to
          // excluded coordinate
          for
          (
            storage::Vector< linal::Vector3D>::const_iterator coord_itr( coords.Begin()), coord_itr_end( coords.End());
            coord_itr != coord_itr_end;
            ++coord_itr
          )
          {
            const linal::Vector3D &coord( *coord_itr);

            // true if the coordinate is undefined
            if( !coord.IsDefined())
            {
              continue;
            }

            // calculate current distance between data point coord and exclusion coord
            const double current_distance( linal::Distance( coord, exclusion_coord));

            // true if the current distance is defined
            if( util::IsDefined( current_distance))
            {
              // add the distance to the considered distances
              mean_sd += current_distance;
            }
          } //< iterate through coordinates of current locator

          // add the average distance to the distance map
          insert_status.first->second.PushBack( mean_sd);
        } //< iterate through excluded coordinates
      } //< iterate through the data point locators

      return distance_map;
    }

    //! @brief for a given data point, calculates the score it has according to its distance from exclusion points
    //! @param LOCATOR the data point whose score will be calculated
    //! @param DISTANCE_MAP map of data point to vector of mean sd info calculated from ENSEMBLE from that data point
    //!        to each excluded coordinate
    //! @param RADIUS if either component of data pair is within this radius the data pair will be scored unfavorably
    //! @return double which is the score of LOCATOR according to its distances from exclusion coordinates
    double DataSetPairwiseCoordinateExclusion::CalculateExclusionScore
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR,
      const DistanceMap &DISTANCE_MAP,
      const double RADIUS
    )
    {
      // try to find the datapoint in the distance map
      DistanceMap::const_iterator itr_first( DISTANCE_MAP.Find( LOCATOR));

      // true if the data point is not found in the exposure map
      if( itr_first == DISTANCE_MAP.End())
      {
        return util::GetUndefinedDouble();
      }

      // initialize score to zero
      double score( 0);

      // iterate through the distance statistics of the locator for all coordinate exclusions
      for
      (
        storage::Vector< math::RunningAverageSD< double> >::const_iterator
          stats_itr( itr_first->second.Begin()), stats_itr_end( itr_first->second.End());
        stats_itr != stats_itr_end;
        ++stats_itr
      )
      {
        const double distance_mean( stats_itr->GetAverage());

        // distance_mean should be greater than m_Radius - give penalty if it is not above the cutoff
        const double current_score( distance_mean > RADIUS ? 0 : RADIUS - distance_mean);

        // add current score to score
        score += current_score;
      }

      return score;
    }

  } // namespace score
  
} // namespace bcl
