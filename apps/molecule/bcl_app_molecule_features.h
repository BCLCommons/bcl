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

#ifndef BCL_APP_MOLECULE_FEATURES_H_
#define BCL_APP_MOLECULE_FEATURES_H_
// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "app/bcl_app.h"
#include "app/bcl_app_interface.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"

// include headers from the bcl - sorted alphabetically

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeFeatures
    //! @brief Used for visualization of input sensitivity on given parts of a molecule
    //!
    //! @author geanesar, brownbp1
    //! @date 08/31/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class MoleculeFeatures :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! flag sets the filename for the molecules to be analyzed
      util::ShPtr< command::FlagInterface> m_MoleculesFlag;

      //! flag sets the filename for the reference molecule(s)
      util::ShPtr< command::FlagInterface> m_ReferenceFlag;

      //! whether to consider mutually matched atom pairs instead of doing MCS search
      util::ShPtr< command::FlagInterface> m_MutuallyMatchingAtomsFlag;

      //! flag that takes the models to use for scoring
      util::ShPtr< command::FlagInterface> m_ScorerFlag;

      //! flag that allows the user to specify element types to perturb
      util::ShPtr< command::FlagInterface> m_ElementPerturberFlag;

      //! flag that specifies the output filename
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! whether to include hydrogen atoms in analysis
      util::ShPtr< command::FlagInterface> m_IgnoreHFlag;

      //! whether to normalize naive perturbation analysis by number of atoms in fragment
      util::ShPtr< command::FlagInterface> m_NormalizeFlag;

      //! whether to normalize final values by the number of molecules in the original ensemble
      util::ShPtr< command::FlagInterface> m_StatisticFlag;

      //! whether to remove disconnected atoms of tiny substructures when doing certain perturbations
      util::ShPtr< command::FlagInterface> m_SplitLargestFlag;

      //! set the perturbation type
      util::ShPtr< command::FlagInterface> m_PerturbTypeFlag;

      //! set the atom comparison type
      util::ShPtr< command::FlagInterface> m_AtomTypeFlag;

      //! set the bond comparison type
      util::ShPtr< command::FlagInterface> m_BondTypeFlag;

      //! set the min and max PyMol spectrum visualization values
      util::ShPtr< command::FlagInterface> m_ColorFlag;

      //! set Rigorous vs. Naive distribution of atom contributions to score
      util::ShPtr< command::FlagInterface> m_ScoreDistributionTypeFlag;

      //! set the PyMol color scheme for output visualization script
      util::ShPtr< command::FlagInterface> m_ColorSpectrumFlag;

      //! whether to output intermediate perturbations to SDF files
      util::ShPtr< command::FlagInterface> m_OutputIntermediatesFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief Default constructor
      MoleculeFeatures();

    public:

      //! @brief Clone function
      //! @return pointer to new MoleculeFeatures
      MoleculeFeatures *Clone() const
      {
        return new MoleculeFeatures( *this);
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

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief initializes the command object for this application
      //! @return a ShPtr to a Command containing all of this applications parameters
      util::ShPtr< command::Command> InitializeCommand() const;

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

      // Static instance of MoleculeFeatures
      static const ApplicationType MoleculeFeatures_Instance;

    }; // class MoleculeFeatures

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_FEATURES_H_
