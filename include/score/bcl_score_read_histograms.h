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

#ifndef BCL_SCORE_READ_HISTOGRAMS_H_
#define BCL_SCORE_READ_HISTOGRAMS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ReadHistograms
    //! @brief TODO: add a brief comment
    //! @details TODO: add an general comment to this class
    //!
    //! @see @link example_score_read_histograms.cpp @endlink
    //! @author woetzen
    //! @date Feb 15, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ReadHistograms
    {

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor should not be called
      ReadHistograms()
      {
      }

    public:

    ////////////////
    // operations //
    ////////////////

      //! read histograms from stream, containing 20 histograms for each amino acids in blocks for each membrane region
      static
      storage::Vector< storage::Vector< math::Histogram> >
      ReadMembraneDependentEnvironmentHistograms
      (
        std::istream &ISTREAM
      );

      //! read histograms from stream, containing 20 histograms for each amino acid
      static
      storage::Vector< math::Histogram>
      ReadEnvironmentHistograms
      (
        std::istream &ISTREAM
      );

    }; // ReadHistograms

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_READ_HISTOGRAMS_H_
