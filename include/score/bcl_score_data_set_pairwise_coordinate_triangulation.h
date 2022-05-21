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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_COORDINATE_TRIANGULATION_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_COORDINATE_TRIANGULATION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwiseCoordinateTriangulation
    //! @brief Scores a DataSetPairwise to select for data points that fall far away from other data points
    //! @details An ensemble provides the coordinate positions of the data points. If the average distance between
    //!          two data points is closure than the provided threshold, the data pair containing one of them is
    //!          unfavorable.
    //!
    //! @see @link example_score_data_set_pairwise_coordinate_triangulation.cpp @endlink
    //! @author alexanns
    //! @date May 10, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseCoordinateTriangulation :
      public math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double>
    {

    private:

    //////////
    // data //
    //////////

      //! the distance points need to be away from one another to not be unfavorable
      double m_RadiusCutoff;

      //! the ensemble from which the distances will be calculated
      util::ShPtr< assemble::ProteinEnsemble> m_Ensemble;

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
      explicit DataSetPairwiseCoordinateTriangulation( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking members
      //! @param RADIUS_CUTOFF coordinates are considered far enough apart when above this distance
      //! @param ENSEMBLE ensemble for which exposures will be calculated
      //! @param SCHEME the scheme for this scoring function
      DataSetPairwiseCoordinateTriangulation
      (
        const double &RADIUS_CUTOFF,
        const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseCoordinateTriangulation
      DataSetPairwiseCoordinateTriangulation *Clone() const;

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

      //! @brief determines if a locator can be used or if one is already close in space to it so don't use it
      //! @param USED_LOCATORS locators that already exist in a data set
      //! @param ENSEMBLE structures that will be used to determine average spatial proximity
      //! @param LOCATOR locator which will be checked to see if it can be used or not based on already used locators
      //! @param DISTANCE_THRESHOLD the minimum spatial distance allowed between LOCATOR and any other used locator
      static bool UsableLocator
      (
        storage::Set
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        > &USED_LOCATORS,
        const assemble::ProteinEnsemble &ENSEMBLE,
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR,
        const double DISTANCE_THRESHOLD
      );

    }; // class DataSetPairwiseCoordinateTriangulation

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DATA_SET_PAIRWISE_COORDINATE_TRIANGULATION_H_ 
