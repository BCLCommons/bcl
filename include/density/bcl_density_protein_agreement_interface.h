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

#ifndef BCL_DENSITY_PROTEIN_AGREEMENT_INTERFACE_H_
#define BCL_DENSITY_PROTEIN_AGREEMENT_INTERFACE_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "score/bcl_score_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinAgreementInterface
    //! @brief Interface class for scoring the agreement of a protein model with a density map
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Aug 8, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinAgreementInterface :
      public score::ProteinModel
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ProteinAgreementInterface
      virtual ProteinAgreementInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief access to the density simulator
      //! @return the density simulator used
      virtual const util::ShPtr< SimulateInterface> &GetSimulator() const = 0;

      //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
      //! @param SP_SIMULATOR ShPtr to SimulatInterface
      virtual void SetSimulator( const util::ShPtr< SimulateInterface> &SP_SIMULATOR) = 0;

      //! @brief access to the density used for agreement calculation
      //! @return SiPtr to the density
      virtual const util::SiPtr< const Map> &GetDensity() const = 0;

      //! @brief set the density used for agreement calculation
      //! @param SP_DENSITY SiPtr to the density map
      virtual void SetDensityMap( const util::SiPtr< const Map> &SP_DENSITY) = 0;

    }; // class ProteinAgreementInterface

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_PROTEIN_AGREEMENT_INTERFACE_H_ 
