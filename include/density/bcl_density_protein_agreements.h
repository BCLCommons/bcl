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

#ifndef BCL_DENSITY_PROTEIN_AGREEMENTS_H_
#define BCL_DENSITY_PROTEIN_AGREEMENTS_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_density_protein_agreement_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinAgreements
    //! @brief class that enumerates implementations of ProteinAgreementInterface
    //! @details A convenience function will create a ShPtr to an implementation setting all the parameters
    //!
    //! @see @link example_density_protein_agreements.cpp @endlink
    //! @author woetzen
    //! @date Aug 8, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinAgreements :
      public util::Enumerate< util::ShPtr< ProteinAgreementInterface>, ProteinAgreements>
    {
      friend class util::Enumerate< util::ShPtr< ProteinAgreementInterface>, ProteinAgreements>;

    public:

    //////////
    // data //
    //////////

      const ProteinAgreement e_CCC;                //!< cross correlation coefficient
      const ProteinAgreement e_CCCScaled;          //!< cross correlation coefficient scaled by protein size
      const ProteinAgreement e_LikelihoodCa;       //!< Likelihood CA atom
      const ProteinAgreement e_LikelihoodCb;       //!< Likelihood CB atom
      const ProteinAgreement e_LikelihoodBackBone; //!< Likelihood Backbone atoms

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinAgreements();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add a new protein agreement instance
      //! @param SP_PROTEIN_AGREEMENT shptr to instance of a ProteinAgreementInterface derived class
      //! @return the new enum constructed from the agreement
      ProteinAgreement AddAgreement( const util::ShPtr< ProteinAgreementInterface> &SP_PROTEIN_AGREEMENT);

      //! @brief create ProteinAgreement from enum, density map and resolution
      //! @param PROTEIN_AGREEMENT enum of the agreement to be used
      //! @param SP_DENSITY density map to be used
      //! @param RESOLUTION resolution in Angstrom to be used
      //! @return ShPtr to density::ProteinAgreement
      util::ShPtr< ProteinAgreementInterface> CreateProteinAgreement
      (
        const ProteinAgreement &PROTEIN_AGREEMENT,
        const Simulator &SIMULATOR,
        const util::SiPtr< const Map> &SP_DENSITY,
        const double RESOLUTION
      ) const;

    }; // class ProteinAgreements

    //! @brief construct on access function for all ProteinAgreements
    //! @return reference to only instances of ProteinAgreements
    BCL_API
    ProteinAgreements &GetProteinAgreements();

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< density::ProteinAgreementInterface>, density::ProteinAgreements>;

  } // namespace util
} // namespace bcl

#endif // BCL_DENSITY_PROTEIN_AGREEMENTS_H_ 
