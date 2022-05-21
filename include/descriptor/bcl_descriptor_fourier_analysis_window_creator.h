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

#ifndef BCL_DESCRIPTOR_FOURIER_ANALYSIS_WINDOW_CREATOR_H_
#define BCL_DESCRIPTOR_FOURIER_ANALYSIS_WINDOW_CREATOR_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_window_weighting_interface.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FourierAnalysisWindowCreator
    //! @brief distributes cache-able fourier-series coefficients for descriptor classes
    //! Does not derive from Object interface because this is a pure singleton class used by internal algorithms and for
    //! caching
    //!
    //! @see @link example_descriptor_power_spectrum.cpp @endlink
    //! @author mendenjl
    //! @date Mar 11, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FourierAnalysisWindowCreator
    {

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FourierAnalysisWindowCreator();

      //! @param copy constructor; private and undefined because this is a singleton class
      FourierAnalysisWindowCreator( const FourierAnalysisWindowCreator &);

    public:

    //////////
    // data //
    //////////

      typedef util::OwnPtr< storage::Pair< linal::Matrix< float>, storage::Vector< float> > > WaveletsFrequenciesPtr;

      //! @brief set up the window matrix
      //! @param ALIGNMENT alignment to use for the matrix
      //! @param SIZE size of the window to use
      //! @param TARGET_PERIODS periods of interest of interest
      //! @param MAX_PERIODS maximum # of periods to insert into the signal; also controls windowing
      //! @param WEIGHTING object that creates the weights vector
      //! @param MAX_OVERLAP the maximum overlap to allow for any two adjacent signals
      //!        Values >= 1 ensure that all frequencies are included, even if they are not orthogonal
      //! @param CACHEABLE true to allow caching the matrix and retrieval of the matrix from the cache, if it is available
      //! @return wavelets matrix and actual periods represented by each pair of rows in the matrix
      static WaveletsFrequenciesPtr CreateWindowMatrix
      (
        const WindowAlignmentEnum &ALIGNMENT,
        const size_t &SIZE,
        const storage::Vector< float> &TARGET_PERIODS,
        const size_t &MAX_PERIODS,
        const util::Implementation< WindowWeightingInterface> &WEIGHTING,
        const float &MAX_OVERLAP = float( 1.0),
        const bool &CACHEABLE = true
      );

      //! @brief get unweighted sin/cosine vectors
      //! @param ALIGNMENT alignment to use; determines meaning of each index in the vector
      //! @param SIZE radius of the window to use
      //! @param PERIOD period desired
      //! @return VectorND2 to linal::VectorConstRef to Cosine (first) and Sine (second) unweighted components
      static storage::Pair< linal::VectorConstReference< float>, linal::VectorConstReference< float> >
        GetUnweightedCosSinCoefficients
        (
          const WindowAlignmentEnum &ALIGNMENT,
          const size_t &SIZE,
          const float &PERIOD
        );

    }; // class FourierAnalysisWindowCreator

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_FOURIER_ANALYSIS_WINDOW_CREATOR_H_
