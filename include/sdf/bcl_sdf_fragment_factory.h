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

#ifndef BCL_SDF_FRAGMENT_FACTORY_H_
#define BCL_SDF_FRAGMENT_FACTORY_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sdf_mdl_handler.h"
#include "bcl_sdf_molecule_reading_pref.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentFactory
    //! @brief create a fragment object that contains all the, configuration and conformation
    //! @details This class functionality for creating fragment entitiy that has, configuration and
    //! conformation information from file sources like SDF or CSD files
    //!
    //! @see @link example_sdf_fragment_factory.cpp @endlink
    //! @author mendenjl
    //! @date Mar 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentFactory
    {

    public:

      //! @brief make complete fragment from an mdl file
      //! @param MDL_HANDLER handler that has the information
      //! @return a fragment complete
      static chemistry::FragmentComplete
        MakeFragment
        (
          const CTabHandler &HANDLER,
          const HydrogenHandlingPref &H_PREF = e_Maintain,
          const NeutralizationPref &NEUTRALIZATION_PREF = e_None,
          const std::string &MOLECULE_ID = std::string()
        );

      //! @brief make complete fragment from an mdl file
      //! @param MDL_HANDLER handler that has the information
      //! @return a fragment complete
      static chemistry::FragmentComplete
        MakeFragment
        (
          const MolfileHandler &HANDLER,
          const HydrogenHandlingPref &H_PREF = e_Maintain,
          const NeutralizationPref &NEUTRALIZATION_PREF = e_None,
          const std::string &MOLECULE_ID = std::string()
        );

      //! @brief make complete fragment from an mdl file
      //! @param MDL_HANDLER handler that has the information
      //! @return a fragment complete
      static chemistry::FragmentComplete
        MakeFragment
        (
          const MdlHandler &HANDLER,
          const HydrogenHandlingPref &H_PREF = e_Maintain,
          const NeutralizationPref &NEUTRALIZATION_PREF = e_None,
          const std::string &MOLECULE_ID = std::string()
        );

      //! @brief make fragment conformation from MDL file
      //! @param MDL_HANDLER handler that has the conformation information
      //! @return a fragment conformation shared
      static chemistry::FragmentConformationShared
        MakeConformation( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF = e_Maintain);

      //! @brief make fragment configuration from MDL file
      //! @param MDL_HANDLER handler that has the configuration information
      //! @return a fragment configuration shared
      static chemistry::FragmentConfigurationShared
        MakeConfiguration( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF = e_Maintain);

      //! @brief make fragment constitution from MDL file
      //! @param MDL_HANDLER handler that has the constitution information
      //! @return a fragment constitution shared
      static chemistry::FragmentConstitutionShared
        MakeConstitution( const MdlHandler &HANDLER, const HydrogenHandlingPref &H_PREF = e_Maintain);

    };
  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_FRAGMENT_FACTORY_H_
