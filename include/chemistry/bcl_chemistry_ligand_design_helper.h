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

#ifndef BCL_CHEMISTRY_LIGAND_DESIGN_HELPER_H_
#define BCL_CHEMISTRY_LIGAND_DESIGN_HELPER_H_

// include the namespace header
#include "bcl_defines.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include <string>

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LigandDesignHelper
    //! @brief Functionality needed by the (java-based) BCL::LigandDesign tool
    //!
    //! @see @link example_chemistry_ligand_design_helper.cpp @endlink
    //! @author mendenjl
    //! @date Mar 02, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LigandDesignHelper
    {
    public:

      //! @brief get the comma separated conformation-independent properties
      //! @return the comma-seperated conformation-independent properties
      static const std::string &GetAvailableProperties();

      //! @brief get the comma separated predictions list
      //! @note prediction is defined here as anything that requires a conformation
      //! @return the comma-seperated predictions list
      static const std::string &GetAvailablePredictions();

      //! @brief set the directory of the jar file. This is used to locate the bcl-related files
      //! @param JAR_FILE_DIRECTORY the directory containing the jar file; should also contain a bcl folder
      static void SetJarDirectory( const std::string &JAR_FILE_DIR);

      //! @brief calculate conformation-independent properties.
      //! @param MOLECULE molecule of interest, encoded as a string in SDF format
      //! @param DESCRIPTORS descriptors (comma seperated) to calculate for the molecule
      //! @return DESCRIPTORS given as values in a comma separated list.
      //!         Elements in vector-valued descriptors separated by spaces
      //!         Note that if the atom types are undefined, returns empty string
      static std::string CalculateProperties
      (
        const std::string &MOLECULE,
        const std::string &PROPERTIES
      );

      //! @brief processes the molecule; producing a 3d-conformation,
      //! @param MOLECULE molecule of interest, encoded as a string in SDF format
      //! @param DESCRIPTORS descriptors (comma seperated) to calculate for the molecule
      //! @return molecule in sdf formation with 3D conformation and DESCRIPTORS as misc properties
      //!         Note that if the atom types are undefined, an empty string is returned
      static std::string ProcessMolecule
      (
        const std::string &MOLECULE,
        const std::string &DESCRIPTORS
      );
    };
  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_LIGAND_DESIGN_HELPER_H_
