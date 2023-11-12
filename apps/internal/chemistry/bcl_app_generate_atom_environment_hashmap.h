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

#ifndef BCL_APP_GENERATE_ATOM_ENVIRONMENT_HASHMAP_H_
#define BCL_APP_GENERATE_ATOM_ENVIRONMENT_HASHMAP_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "storage/bcl_storage_hash_map.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_running_average_sd.h"
#include "storage/bcl_storage_vector.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateAtomEnvironmentHashmap
    //! @brief Application for generating a hashmap reference file for atom environments
    //!
    //! @author brownbp1, mendenjl
    //! @date April 12, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateAtomEnvironmentHashmap :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! Output filename base
      util::ShPtr< command::FlagInterface> m_OutputFilenameBase;

      //! count element_type-element_type-bond_type statistics (normalized)
      util::ShPtr< command::FlagInterface> m_WriteHashMap;

      //! bond radius for counts; valid 2-4
      util::ShPtr< command::FlagInterface> m_BondRadius;

      //! map of atom environment counts
      mutable storage::Vector< storage::HashMap< std::string, size_t> > m_HashMap;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      GenerateAtomEnvironmentHashmap();

    public:

      // instantiate enumerator for GenerateAtomEnvironmentHashmap class
      static const ApplicationType GenerateAtomEnvironmentHashmap_Instance;

      //! @brief Clone function
      //! @return pointer to new GenerateAtomEnvironmentHashmap
      GenerateAtomEnvironmentHashmap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief generate counts of element-element bonds for statistics
      //! @param ENSEMBLE the fragment ensemble for which to calculate the counts
      void CountAtomEnvironments( const chemistry::FragmentEnsemble &ENSEMBLE) const;

    }; // GenerateAtomEnvironmentHashmap

  } // namespace app
} // namespace bcl

#endif // BCL_APP_GENERATE_ATOM_ENVIRONMENT_HASHMAP_H_
