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

#ifndef BCL_SCORE_RESIDUAL_DIPOLAR_COUPLING_HISTOGRAM_H_
#define BCL_SCORE_RESIDUAL_DIPOLAR_COUPLING_HISTOGRAM_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "nmr/bcl_nmr_rdc_container.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ResidualDipolarCouplingHistogram
    //! @brief This class scores RDCs between experimental and theoretical values using a global
    //! histogram approach
    //! @details Clore, G. M. et. al. J. of Magnetic Resonance 133, 216-221 (1998)
    //!
    //! @see @link example_score_residual_dipolar_coupling_histogram.cpp @endlink
    //! @author akinlr, weinerbe
    //! @date Jun 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ResidualDipolarCouplingHistogram :
      public math::FunctionInterfaceSerializable< nmr::RDCContainer, double>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ResidualDipolarCouplingHistogram
      ResidualDipolarCouplingHistogram *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return GetDefaultScheme();
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() for calculating the agreement between theoretical and experimental RDCs
      //! @param RDC_CONTAINER contains both the theoretical and experimental RDCs
      //! @return double which is the agreement between the theoretical and experimental RDCs
      double operator()( const nmr::RDCContainer &RDC_CONTAINER) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ResidualDipolarCouplingHistogram

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESIDUAL_DIPOLAR_COUPLING_HISTOGRAM_H_
