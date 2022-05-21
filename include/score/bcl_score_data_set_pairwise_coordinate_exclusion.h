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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_COORDINATE_EXCLUSION_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_COORDINATE_EXCLUSION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwiseCoordinateExclusion
    //! @brief Scores a DataSetPairwise to select for data pairs that are away from undesirable coordinates
    //! @details Selects for DataSetPairwise that have data pairs whose located coordinates are outside of an input
    //!          radius around input coordinates that are undesirable to be around. Undesirable coordinates are file
    //!          input and the columns that contain the x, y, z coordinate can be defined. Column numbering starts at 0.
    //!
    //! @see @link example_score_data_set_pairwise_coordinate_exclusion.cpp @endlink
    //! @author alexanns
    //! @date May 9, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseCoordinateExclusion :
      public math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double>
    {

    private:

      // type def for a map that contains for LocatorAtomCoordinatesInterface objects the mean and sd of how close
      // close that coordinate comes to any of possible multiple undesirable coordinates. The mean and sd comes from
      // the possibility of an ensemble and therefore an ensmble of distances to each undesirable coordinate
      typedef storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        // mean and sd of distance to an exclusion coordinate and vector holds this info for each exclusion coordinate
        storage::Vector< math::RunningAverageSD< double> >,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > DistanceMap;

    //////////
    // data //
    //////////

      //! if either component of a data pair is within this radius the data pair will be scored unfavorably
      double m_Radius;

      //! map of distances to input exclusion coordinates
      DistanceMap m_DistanceMap;

      //! the scheme of this mutate
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking optional scheme
      //! @param SCHEME the scheme for this scoring function
      explicit DataSetPairwiseCoordinateExclusion( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking members
      //! @param EXCLUSION_RADIUS data with atoms inside this radius of any exclusion coordinates will be penalized
      //! @param READ the stream the coordinates will be read from
      //! @param X_COORD_COLUMN column in istream that has the x-coordinate - columns start at 0
      //! @param Y_COORD_COLUMN column in istream that has the y-coordinate - columns start at 0
      //! @param Z_COORD_COLUMN column in istream that has the z-coordinate - columns start at 0
      //! @param ALL_POSSIBLE_DATA_POINTS the set of data points that should be subjected to this filter
      //! @param ENSEMBLE ensemble for which coordinates will be checked to make sure they aren't near exclusion coords
      //! @param SCHEME the scheme for this scoring function
      DataSetPairwiseCoordinateExclusion
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
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseCoordinateExclusion
      DataSetPairwiseCoordinateExclusion *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the score of a data set
      //! @param DATA data set to be scored
      //! @return the score of the current data set
      double operator()( const restraint::DataSetPairwise &DATA) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

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
      static DataSetPairwiseCoordinateExclusion::DistanceMap FillDistanceMap
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
      );

      //! @brief gives list of coordinates read from an istream
      //! @param READ the stream from which undesirable coordinates will be read
      //! @param X_COORD_COLUMN the column the x coordinate are in - starts at 0
      //! @param Y_COORD_COLUMN the column the y coordinate are in - starts at 0
      //! @param Z_COORD_COLUMN the column the z coordinate are in - starts at 0
      //! @return list of Vector3D which are the coordinates read from a file
      static storage::List< linal::Vector3D> ReadCoordinates
      (
        std::istream &READ,
        const size_t X_COORD_COLUMN,
        const size_t Y_COORD_COLUMN,
        const size_t Z_COORD_COLUMN
      );

      //! @brief locates coordinates for a set of locators for an ensemble of proteins
      //! @param ALL_POSSIBLE_DATA_POINTS the data points whose coordinates will be located
      //! @param ENSEMBLE the protein models from which coordinates will be found
      //! @return map of locator to its associated vector or coordinates, one coordinate for each model in the ensemble
      static storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        storage::Vector< linal::Vector3D>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > LocateCoordinates
      (
        const storage::Set
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        > &ALL_POSSIBLE_DATA_POINTS,
        const assemble::ProteinEnsemble &ENSEMBLE
      );

      //! @brief calculates the mean and sd info for data points distances to exclusion coordinates
      //! @param COORDINATE_MAP map of locator to associated vector of coords; 1 coord for each model in the ensemble
      //! @param EXCLUSION_COORDINATES list of coordinates which will be used to calculate distances to them
      //! @return distance map of data point to vector of mean sd info calculated from ENSEMBLE from that data point to
      //!         each excluded coordinate
      static DistanceMap CalculateAverageDistanceToExcludedCoordinates
      (
        const storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
          storage::Vector< linal::Vector3D>,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        > &COORDINATE_MAP,
        const storage::List< linal::Vector3D> &EXCLUSION_COORDINATES
      );

      //! @brief for a given data point, calculates the score it has according to its distance from exclusion points
      //! @param LOCATOR the data point whose score will be calculated
      //! @param DISTANCE_MAP map of data point to vector of mean sd info calculated from ENSEMBLE from that data point
      //!        to each excluded coordinate
      //! @param RADIUS if either component of data pair is within this radius the data pair will be scored unfavorably
      //! @return double which is the score of LOCATOR according to its distances from exclusion coordinates
      static double CalculateExclusionScore
      (
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR,
        const DistanceMap &DISTANCE_MAP,
        const double RADIUS
      );

    }; // class DataSetPairwiseCoordinateExclusion

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DATA_SET_PAIRWISE_COORDINATE_EXCLUSION_H_ 
