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

#ifndef BCL_SCORE_RESIDUAL_DIPOLAR_COUPLING_Q_VALUE_H_
#define BCL_SCORE_RESIDUAL_DIPOLAR_COUPLING_Q_VALUE_H_

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
    //! @class ResidualDipolarCouplingQValue
    //! @brief calculates the "q-value" agreement between theoretical and experimental RDCs
    //!        See Meiler et al. Journal of Biomolecular Nmr, 2000. The DipoCoup paper
    //!
    //! @see @link example_score_residual_dipolar_coupling_q_value.cpp @endlink
    //! @author alexanns, weinerbe
    //! @date May 31, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ResidualDipolarCouplingQValue :
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
      //! @return pointer to new ResidualDipolarCouplingQValue
      ResidualDipolarCouplingQValue *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
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
      //! @param RDC_RESTRAINTS contains the experimental RDCs and calculated RDCs
      //! @return double which is the agreement between the theoretical and experimental RDCs
      double operator()( const nmr::RDCContainer &RDC_RESTRAINTS) const;

    //////////////////////
    // input and output //
    //////////////////////

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param RDC_RESTRAINTS contains the experimental RDCs and calculated RDCs
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const nmr::RDCContainer &RDC_RESTRAINTS,
        std::ostream &OSTREAM
      ) const;

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ResidualDipolarCouplingQValue

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESIDUAL_DIPOLAR_COUPLING_Q_VALUE_H_
