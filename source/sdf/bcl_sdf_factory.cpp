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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "sdf/bcl_sdf_factory.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_molecular_configuration_shared.h"
#include "chemistry/bcl_chemistry_molecular_conformation_shared.h"
#include "chemistry/bcl_chemistry_molecule_complete.h"
#include "sdf/bcl_sdf_fragment_factory.h"

namespace bcl
{
  namespace sdf
  {

    //! @brief read complete molecule from an mdl file
    //! @param MDL_HANDLER handler that has the information
    //! @return a molecule complete
    chemistry::MoleculeComplete Factory::MakeMolecule( const MdlHandler &HANDLER)
    {
      return FragmentFactory::MakeFragment( HANDLER, e_Saturate);
    }

    //! @brief read molecule conformation from MDL file
    //! @param MDL_HANDLER handler that has the conformation information
    //! @return a molecule conformation shared
    chemistry::MolecularConformationShared Factory::MakeConformation( const MdlHandler &HANDLER)
    {
      return chemistry::MolecularConformationShared( MakeMolecule( HANDLER));
    }

    //! @brief read molecule configuration from MDL file
    //! @param MDL_HANDLER handler that has the configuration information
    //! @return a molecule configuration shared
    chemistry::MolecularConfigurationShared Factory::MakeConfiguration( const MdlHandler &HANDLER)
    {
      return chemistry::MolecularConfigurationShared( MakeMolecule( HANDLER));
    }

    //! @brief read molecule constitution from MDL file
    //! @param MDL_HANDLER handler that has the constitution information
    //! @return a molecule constitution shared
    chemistry::MolecularConstitutionShared Factory::MakeConstitution( const MdlHandler &HANDLER)
    {
      return chemistry::MolecularConstitutionShared( MakeMolecule( HANDLER));
    }

  } // namespace sdf
} // namespace bcl

