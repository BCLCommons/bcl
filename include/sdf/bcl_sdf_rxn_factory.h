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

#ifndef BCL_SDF_RXN_FACTORY_H_
#define BCL_SDF_RXN_FACTORY_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sdf_mdl_handler.h"
#include "bcl_sdf_molecule_reading_pref.h"
#include "bcl_sdf_rxn_handler.h"
#include "chemistry/bcl_chemistry_reaction_complete.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RXNFactory
    //! @brief create a fragment object that contains all the, configuration and conformation
    //! @details This class functionality for creating fragment entitiy that has, configuration and
    //! conformation information from file sources like SDF or CSD files
    //!
    //! @see @link example_sdf_rxn_factory.cpp @endlink
    //! @author geanesar, combss, mendenjl
    //! @date Jun 02, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RXNFactory
    {

    public:

      //! @brief construct a ReactionComplete from a rxn file
      //! @param HANDLER the handler to use
      static chemistry::ReactionComplete MakeReactionComplete( const RXNHandler &HANDLER);

      //! @brief construct a ReactionStructure from a molfile
      static chemistry::ReactionStructure MakeReactionStructure( const CTabHandler &HANDLER);

    };
  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_RXN_FACTORY_H_
