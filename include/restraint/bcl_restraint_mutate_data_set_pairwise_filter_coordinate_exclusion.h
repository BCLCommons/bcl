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

#ifndef BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_COORDINATE_EXCLUSION_H_
#define BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_COORDINATE_EXCLUSION_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_running_average_sd.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDataSetPairwiseFilterCoordinateExclusion
    //! @brief Mutates a DataSetPairwise by filtering out data pairs that are too close to undesirable coordinates
    //! @details Data pairs are excluded from a given data set if either of the data points comes too close to any of
    //!          undesired coordinates. Coordinates are input from a file. The data points
    //!          ( i.e. LocatorAtomCoordinatesInterfaces) that are subjected to this filter are user definable.
    //!
    //! @see @link example_restraint_mutate_data_set_pairwise_filter_coordinate_exclusion.cpp @endlink
    //! @author alexanns
    //! @date May 11, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDataSetPairwiseFilterCoordinateExclusion :
      public math::MutateInterface< DataSetPairwise>
    {

    private:

      //! type def for a map that contains for LocatorAtomCoordinatesInterface objects the mean and sd of how close
      //! close that coordinate comes to any of possible multiple undesirable coordinates. The mean and sd comes from
      //! the possibility of an ensemble and therefore an ensmble of distances to each undesirable coordinate
      typedef storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        // RunningAverageSD< double> is stats for ensemble average distance to an exclusion coordinate and
        // the vector holds this info for each exclusion coordinate
        storage::Vector< math::RunningAverageSD< double> >,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > DistanceMap;

    //////////
    // data //
    //////////

      //! if either component of a data pair is within this radius the data pair will be removed
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
      explicit MutateDataSetPairwiseFilterCoordinateExclusion( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking members
      //! @param EXCLUSION_RADIUS data with atoms inside this radius of any exclusion coordinates will be penalized
      //! @param READ the stream the coordinates will be read from
      //! @param X_COORD_COLUMN column in istream that has the x-coordinate - columns start at 0
      //! @param Y_COORD_COLUMN column in istream that has the y-coordinate - columns start at 0
      //! @param Z_COORD_COLUMN column in istream that has the z-coordinate - columns start at 0
      //! @param ALL_POSSIBLE_DATA_POINTS the set of data points that should be subjected to this filter
      //! @param ENSEMBLE ensemble for which coordinates will be checked to make sure they aren't near exclusion coords
      //! @param SCHEME the scheme for this scoring function
      MutateDataSetPairwiseFilterCoordinateExclusion
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
      //! @return pointer to new MutateDataSetPairwiseFilterCoordinateExclusion
      MutateDataSetPairwiseFilterCoordinateExclusion *Clone() const;

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

      //! @brief virtual operator taking an DATA and returning a mutated object of DATA
      //! @param DATA DATA of interest that will be mutated
      //! @return MutateResult that results from mutating to the argument DATA
      math::MutateResult< DataSetPairwise> operator()( const DataSetPairwise &DATA) const;

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

    private:

    }; // class MutateDataSetPairwiseFilterCoordinateExclusion

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_COORDINATE_EXCLUSION_H_ 
