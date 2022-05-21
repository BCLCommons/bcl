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

#ifndef BCL_RESTRAINT_SAXS_DATA_REDUCTION_H_
#define BCL_RESTRAINT_SAXS_DATA_REDUCTION_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_scattering_data.h"
#include "bcl_restraint_sas_scattering_point.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_flag_interface.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SaxsDataReduction
    //! @brief Reduces a large data set to a smaller set using Shannon sampling and Noisy Data Reduction
    //! @details see Accurate Assessment of mass, models and resolution by small-angle scattering, Robert Rambo
    //!          John Tainer, Nature 2013
    //!
    //! @see @link example_restraint_saxs_data_reduction.cpp @endlink
    //! @author putnamdk
    //! @date Aug 1, 2013
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SaxsDataReduction :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    /////////////////
    // operations  //
    /////////////////

      //! @brief compute the shannon numbers for saxs profile generation from experimental data
      //! @param DMAX - Maximum dimension of protein inferred from GNOM
      //! @param QMAX - Maximum momentum transfer value
      //! @return number of bins necessary to represent given protein model
      static size_t ComputeShannonNumber( const double &DMAX, const double &QMAX);

      //! @brief samples shannon bins n times ( n typically is 1000 - 3000) to produce noise free scattering profile
      //! @param PROTEIN_MODEL - protein model that contains the coordinates for current model
      //! @param NUMBER_OF_ITERATIONS - the number of iterations to compute ( n typically is 1000 - 3000)
      //! @param ORIGINAL_DATA - pointer to the preprocessed experimental data
      //! @param DENSITY_DATA - pointer to the transformed experimental data to P(r) domain
      static util::ShPtr< SasScatteringData> SasSignalRecovery
      (
        assemble::ProteinModel &PROTEIN_MODEL,
        size_t NUMBER_OF_ITERATIONS,
        const util::ShPtr< SasScatteringData> ORIGINAL_DATA,
        const double DMAX
      );

      //! @brief splits scatting profile into n (shannon) bins.  Selects datapoint to represent the curve from
      //!        each bin that has the least error associated with it.
      //! @param ORIGINAL_DATA - pointer to the preprocessed experimental data
      //! @param DMAX - maximum dimension of the protein
      static util::ShPtr< SasScatteringData> SasSignalRecoveryEstimate
      (
        const util::ShPtr< SasScatteringData> ORIGINAL_DATA,
        const double DMAX
      );

      //! @brief randomly select a point inside the given boundary
      //! @param LEFT - the left boundary of the interval
      //! @param RIGHT - the right boundary of the interval
      //! @return the randomly selected value inside the interval [ LEFT, RIGHT]
      static int SelectIndex( int LEFT, int RIGHT);

      //! @brief randomly select a point inside the given boundary
      //! @param ORIGINAL_DATA - pointer to the experimental data
      //! @param LEFT - the left boundary of the interval
      //! @param RIGHT - the right boundary of the interval
      //! @return the randomly selected value inside the interval [ LEFT, RIGHT]
      static int SelectIndexMinError( const util::ShPtr< SasScatteringData> ORIGINAL_DATA, int LEFT, int RIGHT);

    }; // class SaxsDataReduction

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAXS_DATA_REDUCTION_H_
