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

#ifndef BCL_APP_MAP_PARAMS_H_
#include "sched/bcl_sched_mutex.h"

#define BCL_APP_MAP_PARAMS_H_

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"

// app header
#include "app/bcl_app.h"
#include "app/bcl_app_interface.h"

// include headers from the bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MapParams
    //! @brief Maps Rosetta atom names to Amber prep files; intended that this app will evolve to handle multiple
    //! params type conversions; not terribly generic at the moment
    //!
    //! @author brownbp1
    //! @date 10/16/2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class MapParams :
        public InterfaceRelease
        {
    private:

    //////////
    // data //
    //////////

      //! flag that specifies the input Rosetta SDF
      util::ShPtr< command::FlagInterface> m_RosettaSDF;

      //! flag that specifies the MDL property containing Rosetta atom names
      util::ShPtr< command::FlagInterface> m_RosettaAtomNamesMDL;

      //! flag that specifies the input AmberTools SDF
      util::ShPtr< command::FlagInterface> m_AmberToolsSDF;

      //! flag that specifies the input AmberTools PREPI
      util::ShPtr< command::FlagInterface> m_AmberToolsPrepi;

      //! flag that specifies the input AmberTools MC
      util::ShPtr< command::FlagInterface> m_AmberToolsMC;

      //! Scheme used for comparing whether two atoms are equivalent
      util::ShPtr< command::FlagInterface> m_AtomComparisonType;

      //! Scheme used for comparing whether two bonds are equivalent
      util::ShPtr< command::FlagInterface> m_BondComparisonType;

      //! flag that specifies the output filename
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

      //! rosetta molecule
      mutable chemistry::FragmentEnsemble m_RosettaMol;

      //! ambertools molecule
      mutable chemistry::FragmentEnsemble m_AmberToolsMol;

      //! rosetta atom names
      mutable storage::Vector< std::string> m_RosettaSDFAtomNames;

      //! ambertools atom names
      mutable storage::Vector< std::string> m_AmberToolsSDFAtomNames;

      //! ambertools prepi atom names
      mutable storage::Vector< std::string> m_AmberToolsPrepiAtomNames;

      //! ambertools mc omit atom names
      mutable storage::Vector< std::string> m_AmberToolsOmitAtomNames;

      //! map between ambertools sdf molecule indices and ambertools prepi names
      mutable storage::Map< size_t, std::string> m_AmberToolsMap;

      //! mutex for i/o thread safety
      mutable sched::Mutex m_Mutex;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief Default constructor
      MapParams();

    public:

      //! @brief Clone function
      //! @return pointer to new MapParams
      MapParams *Clone() const
      {
        return new MapParams( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

    //////////////////////
    //   operations     //
    //////////////////////

      //! @brief initializes the command object for this application
      //! @return a ShPtr to a Command containing all of this applications parameters
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief reads in the Rosetta SDF
      void ReadRosettaSDF() const;

      //! @brief reads in the AmberTools SDF
      void ReadAmberToolsSDF() const;

      //! @brief reads in the AmberTools PREPI
      void ReadAmberToolsAtomNames() const;

      //! @brief reads in the AmberTools MC
      void ReadAmberToolsOmitAtomNames() const;

      //! @brief gets rosetta atom names from SDF
      void GetRosettaSDFAtomNames() const;

      //! @brief gets ambertools atom names from SDF
      void GetAmberToolsSDFAtomNames() const;

      //! @brief maps the rosetta molecule to the ambertools molecule
      //! @param ROSETTA_MOL the input rosetta molecule
      //! @param AMBERTOOLS_MOL the input ambertools molecule
      //! @return the common subgraph isomorphism between the two
      graph::CommonSubgraphIsomorphism< size_t, size_t> MapRosettaToAmberToolsSDF
      (
        const chemistry::FragmentComplete &ROSETTA_MOL,
        const chemistry::FragmentComplete &AMBERTOOLS_MOL
      ) const;

      //! @brief map the rosetta atom names to a new prepi file
      storage::Vector< std::string> MapRosettaToPrepi() const;

      //! @brief gets the indices of the AmberTools molecule that
      //! correspond to the atom names in the AmberTools PREPI
      void MapAmberToolsSDFPrepi() const;

    //////////////////////
    //    operators     //
    //////////////////////

      //! @brief the Main function
      //! @return 0 for success
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

      // Static instance of MapParams
      static const ApplicationType MapParams_Instance;

        }; // class MapParams

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MAP_PARAMS_H_
