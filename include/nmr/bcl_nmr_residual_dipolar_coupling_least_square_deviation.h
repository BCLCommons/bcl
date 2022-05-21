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

#ifndef BCL_NMR_RESIDUAL_DIPOLAR_COUPLING_LEAST_SQUARE_DEVIATION_H_
#define BCL_NMR_RESIDUAL_DIPOLAR_COUPLING_LEAST_SQUARE_DEVIATION_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_nmr_rdc_container.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_rdc_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ResidualDipolarCouplingLeastSquareDeviation
    //! @brief for calculating theoretical RDCs from coordinates of a structure and experimental RDCs
    //! @details The restraints are assumed to be correctly normalized prior
    //!        See Meiler et al. Journal of Biomolecular Nmr, 2000. The DipoCoup paper
    //!
    //! @see @link example_nmr_residual_dipolar_coupling_least_square_deviation.cpp @endlink
    //! @author alexanns, weinerbe
    //! @date June 1, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ResidualDipolarCouplingLeastSquareDeviation :
      public math::FunctionInterfaceSerializable< restraint::RDCAssignment, RDCContainer>
    {

    private:

    //////////
    // data //
    //////////

      //! const size_t "number_independent_parameters" to hold the number of parameters that need to be solved for
      static const size_t s_NumberIndependentParameters;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ResidualDipolarCouplingLeastSquareDeviation();

      //! @brief Clone function
      //! @return pointer to new ResidualDipolarCoupling
      ResidualDipolarCouplingLeastSquareDeviation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() taking a list of RDC assignments and returning a ResidualDipolarCouplingContainer
      //! The ResidualDipolarCouplingContainer will contain theoretical RDCs which have been calculated by doing a
      //! least squares deviation fitting to the experimental RDCs and the ResidualDipolarCouplingContainer will
      //! also contain the experimental RDCs
      //! @param RESTRAINTS list of RDC assignments
      //! @return experimental and calculated RDCs
      RDCContainer operator()( const restraint::RDCAssignment &RESTRAINTS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ResidualDipolarCouplingLeastSquareDeviation

  } // namespace nmr
} // namespace bcl

#endif // BCL_NMR_RESIDUAL_DIPOLAR_COUPLING_LEAST_SQUARE_DEVIATION_H_ 
