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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_EUCLIDIAN_DISTANCE_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_EUCLIDIAN_DISTANCE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwiseEuclidianDistance
    //! @brief Scores DataSetPairwise to select for data pairs such that the two data points are within a distance range
    //! @details The data points within a datapair should be within the provided distance range of one another in
    //!          order to be considered favorable. The average distance between the two data points is calculated from
    //!          a provided ensemble.
    //!
    //! @see @link example_score_data_set_pairwise_euclidian_distance.cpp @endlink
    //! @author alexanns
    //! @date May 8, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseEuclidianDistance :
      public math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double>
    {

    private:

    //////////
    // data //
    //////////

      //! the distance range that the two data points in data pair should be within
      math::Range< double> m_DistanceRange;

      //! the ensemble from which distances will be calculated
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
      explicit DataSetPairwiseEuclidianDistance( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking optional scheme
      //! @param RANGE desired range of euclidian separation between data points
      //! @param ENSEMBLE ensemble distances will be calculated from
      //! @param SCHEME the scheme for this scoring function
      DataSetPairwiseEuclidianDistance
      (
        const math::Range< double> &RANGE,
        const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief constructor taking optional scheme
      //! @param RANGE desired range of euclidian separation between data points
      //! @param MODEL protein model distances will be calculated from
      //! @param SCHEME the scheme for this scoring function
      DataSetPairwiseEuclidianDistance
      (
        const math::Range< double> &RANGE,
        const assemble::ProteinModel &MODEL,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseEuclidianDistance
      DataSetPairwiseEuclidianDistance *Clone() const;

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

    private:

    }; // class DataSetPairwiseEuclidianDistance

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DATA_SET_PAIRWISE_EUCLIDIAN_DISTANCE_H_ 
