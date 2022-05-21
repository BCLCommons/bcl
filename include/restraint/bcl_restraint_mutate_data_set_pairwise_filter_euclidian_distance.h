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

#ifndef BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_EUCLIDIAN_DISTANCE_H_
#define BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_EUCLIDIAN_DISTANCE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDataSetPairwiseFilterEuclidianDistance
    //! @brief Mutates a DataSetPairwise by removing data pairs that don't have distances within a given range
    //! @details If the distance indicated by a data pair is longer than a desired range the data pair will be removed.
    //!
    //! @see @link example_restraint_mutate_data_set_pairwise_filter_euclidian_distance.cpp @endlink
    //! @author alexanns
    //! @date May 11, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDataSetPairwiseFilterEuclidianDistance :
      public math::MutateInterface< DataSetPairwise>
    {

    private:

    //////////
    // data //
    //////////

      //! the distance range data pairs must be within
      math::Range< double> m_DistanceRange;

      //! ensemble of models distances will be calculated from
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
      explicit MutateDataSetPairwiseFilterEuclidianDistance( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking optional scheme
      //! @param RANGE desired range of euclidian separation between data points
      //! @param ENSEMBLE ensemble distances will be calculated from
      //! @param SCHEME the scheme for this scoring function
      MutateDataSetPairwiseFilterEuclidianDistance
      (
        const math::Range< double> &RANGE,
        const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new MutateDataSetPairwiseFilterEuclidianDistance
      MutateDataSetPairwiseFilterEuclidianDistance *Clone() const;

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

      //! @brief virtual operator taking an DataSetPairwise and returning a mutated object
      //! @param DATA DataSetPairwise of interest that will be mutated
      //! @return MutateResult that results from mutating to the argument DATA_SET
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

    }; // class MutateDataSetPairwiseFilterEuclidianDistance

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_EUCLIDIAN_DISTANCE_H_ 
