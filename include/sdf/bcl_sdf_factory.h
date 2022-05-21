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

#ifndef BCL_SDF_FACTORY_H_
#define BCL_SDF_FACTORY_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sdf_mdl_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Factory
    //! @brief creates molecule objects, either for a specific layer (e.g. constitution, configuration, conformation)
    //! or all layers (complete)
    //!
    //! @see @link example_sdf_factory.cpp @endlink
    //! @author mendenjl
    //! @date Mar 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Factory
    {

    public:

      //! @brief read complete molecule from an mdl file
      //! @param MDL_HANDLER handler that has the information
      //! @return a molecule complete
      static chemistry::MoleculeComplete             MakeMolecule( const MdlHandler &HANDLER);

      //! @brief read molecule conformation from MDL file
      //! @param MDL_HANDLER handler that has the conformation information
      //! @return a molecule conformation shared
      static chemistry::MolecularConformationShared MakeConformation( const MdlHandler &HANDLER);

      //! @brief read molecule configuration from MDL file
      //! @param MDL_HANDLER handler that has the configuration information
      //! @return a molecule configuration shared
      static chemistry::MolecularConfigurationShared MakeConfiguration( const MdlHandler &HANDLER);

      //! @brief read molecule constitution from MDL file
      //! @param MDL_HANDLER handler that has the constitution information
      //! @return a molecule constitution shared
      static chemistry::MolecularConstitutionShared MakeConstitution( const MdlHandler &HANDLER);

    };
  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_FACTORY_H_
