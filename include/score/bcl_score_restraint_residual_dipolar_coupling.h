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

#ifndef BCL_SCORE_RESTRAINT_RESIDUAL_DIPOLAR_COUPLING_H_
#define BCL_SCORE_RESTRAINT_RESIDUAL_DIPOLAR_COUPLING_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "nmr/bcl_nmr.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintResidualDipolarCoupling
    //! @brief is for scoring a protein model's agreement with residual dipolar coupling data
    //!        This is the upper class for scoring residual dipolar coupling data
    //!        It generates the assignments from a given protein model and a member variable list of Restraints
    //!        It then uses a member variable to calculate theoretical RDCs from the assignments
    //!        Then it uses a member variable to calculate the agreement between the theoretical and experimental rdcs
    //!
    //! @see @link example_score_restraint_residual_dipolar_coupling.cpp @endlink
    //! @author alexanns, weinerbe
    //! @date Apr 12, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintResidualDipolarCoupling :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! ShPtr to a function interface which takes a list of assigned RDCs and calculates theoretical RDCs
      //! the calculated and experimental RDCs are put into a ResidualDipolarCouplingContainer
      util::Implementation< math::FunctionInterfaceSerializable< restraint::RDCAssignment, nmr::RDCContainer> > m_CalculatorOfRDCs;

      //! ShPtr to a function interface which takes a ResidualDipolarCouplingContainer and returns an agreement value
      //! since the ResidualDipolarCouplingContainer contains the experimental and calculated RDCs
      util::Implementation< math::FunctionInterfaceSerializable< nmr::RDCContainer, double> > m_RDCAgreementCalculator;

      //! RDC restraints
      util::ShPtr< restraint::RDC> m_RDCs;

      //! scheme to be used
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintResidualDipolarCoupling();

      //! @brief construct from member variables
      //! @param RDCS experimental rdcs
      //! @param RDC_FROM_STRUCTURE_CALCULATOR function interface for calculating theoretical RDCs
      //! @param RDC_AGREEMENT_CALCULATOR function interface for calculating an RDC agreement value
      //! @param SCHEME Scheme to be used
      RestraintResidualDipolarCoupling
      (
        const util::ShPtr< restraint::RDC> &RDCS,
        const math::FunctionInterfaceSerializable< restraint::RDCAssignment, nmr::RDCContainer> &RDC_FROM_STRUCTURE_CALCULATOR,
        const math::FunctionInterfaceSerializable< nmr::RDCContainer, double> &RDC_AGREEMENT_CALCULATOR,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new RestraintResidualDipolarCoupling
      RestraintResidualDipolarCoupling *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() which takes an ProteinModel for calculating its agreement with the RDC data
      //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the RDC data
      //! @return return a double which is the score of the agreement of the ProteinModel with the RDC data
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

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

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the RDC data
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        std::ostream &OSTREAM
      ) const;

    }; // class RestraintResidualDipolarCoupling

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_RESIDUAL_DIPOLAR_COUPLING_H_
