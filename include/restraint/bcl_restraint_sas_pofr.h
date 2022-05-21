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

#ifndef BCL_RESTRAINT_SAS_POFR_H_
#define BCL_RESTRAINT_SAS_POFR_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_experimental_and_calculated_density.h"
#include "bcl_restraint_sas_pofr_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasPofR
    //! @brief Compute the density distribution function from a protein model
    //! @details
    //!
    //! @see @link example_restraint_saxs_pofr.cpp @endlink
    //! @author putnamdk
    //! @date July 1, 2015
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasPofR :
      public SasPofRInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SasDebye
      SasPofR *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Get the positions for all the atoms in the protein model
      //! @param MODEL the protein model of interest
      //! @return a vector containing positions
      storage::Vector< linal::Vector3D> GetAtoms( const assemble::ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief overloaded () operator to calculate Intensity from Q
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return returns Intensity for given Q for both the experimental and calculated data
      SasExperimentalAndCalculatedDensity operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class SasPofR

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_POFR_H_
