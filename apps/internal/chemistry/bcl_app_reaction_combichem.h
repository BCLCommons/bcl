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

#ifndef BCL_APP_REACTION_COMBICHEM_H_
#define BCL_APP_REACTION_COMBICHEM_H_

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_reaction_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter_check_file_existence.h"
namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ReactionCombichem
    //! @brief Reacts fragments together to create larger molecules
    //! @details Provided one or more SDF files, reacts together fragments in a pairwise manner enumeratively for a specified
    //! number of iterations. Scores the resulting molecules with a specified score function
    //!
    //! @author brownbp1
    //! @date April 06, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ReactionCombichem :
      public InterfaceRelease
    {
    public:

      // static member of this class
      static const ApplicationType ReactionCombichem_Instance;

    private:

    //////////
    // data //
    //////////

      //! molecules to be linked
      util::ShPtr< command::FlagInterface> m_StartingMoleculesFlag;

      //! output SDF file
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! the protein-ligand interface scoring function to use for pose refinement
      util::ShPtr< command::FlagInterface> m_LigandLocalDockFlag;

      //! the directory containing the reaction files
      util::ShPtr< command::FlagInterface> m_ReactionsFlag;

      //! the molecules to be reacted with the input molecules
      util::ShPtr< command::FlagInterface> m_ReagentsFlag;

      //! the number of rounds of combichem to perform
      util::ShPtr< command::FlagInterface> m_NRoundsFlag;

      //! limit the number of reactions that can be performed each round to 1
      util::ShPtr< command::FlagInterface> m_LimitOneRxnPerRoundFlag;

      //! save all products from all rounds for final output
      util::ShPtr< command::FlagInterface> m_SaveAllRoundsFlag;

      //! the MDL property name specifying the protein pocket filename
      //! note that this must match the model
      util::ShPtr< command::FlagInterface> m_MDLStringFlag;

      //! the drug-likeness metric to use
      util::ShPtr< command::FlagInterface> m_DrugLikenessTypeFlag;

      //! flag indication of first molecule to load from ensemble A
      util::ShPtr< command::FlagInterface> m_Start;

      //! flag indication of last molecule to load from ensemble A
      util::ShPtr< command::FlagInterface> m_MaxMols;

      //! flag to enable corina conformer generation (if installed)
      util::ShPtr< command::FlagInterface> m_CorinaFlag;

      //! The reaction operation class, used for modifying structures
      mutable chemistry::FragmentReact m_ReactOp;

      //! The reaction operation class, used for modifying structures
      mutable chemistry::ReactionSearch m_ReactionSearch;

      //! the MDL property name specifying the protein pocket filename
      mutable std::string m_MDLString;

      //! get the filename of the protein pocket
      mutable std::string m_PocketFilename;

      //! the path to the reactions directory
      mutable std::string m_ReactionsDirectory;

      //! molecule ensembles
      mutable storage::Vector< chemistry::FragmentComplete> m_StartEnsemble;
      mutable storage::Vector< chemistry::FragmentComplete> m_Products;
      mutable chemistry::FragmentEnsemble m_ReactantEnsemble;

      //! reactions
      mutable chemistry::ReactionEnsemble m_ReactionEnsemble;
      mutable storage::Vector< chemistry::ReactionComplete> m_Reactions;
      mutable storage::Map< size_t, chemistry::ReactionComplete> m_UnusedReactions;

      //! descriptor for pose scoring
      mutable descriptor::CheminfoProperty m_LigandDockScorer;

      static chemistry::ConformationComparisonPsiField s_Aligner;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      ReactionCombichem();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      ReactionCombichem *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

    /////////////////
    // operations  //
    /////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // helper functions //
    /////////////////////

      //! @brief remove whitespace (via isspace) from a string
      //! @param STR the string to remove whitespace from
      //! @return STR without any whitespace
      std::string RemoveWhitespace( const std::string &STR) const;

      //! @brief wrapper that calls the FragmentMapConformer::Clean function to generate a legitimate 3D conformer and fix bond lengths
      //! @param FRAGMENT the molecule that needs to be fixed
      //! @param REFERENCE scaffold molecule for substructure alignment reference
      //! @return pointer to cleaned molecule
      util::ShPtr< chemistry::FragmentComplete> CallCleaner
      (
        const chemistry::FragmentComplete &FRAGMENT,
        const chemistry::FragmentComplete &REFERENCE
      ) const;

      //! @brief wrapper that calls ResolveClashes and OptimizePose with a combined global/local conf ensemble for each mol
      //! @param ENSEMBLE the molecule that needs to be fixed
      //! @param SCORER the property to be used to score the interface
      //! @param REFERENCE the scaffold for pre-alignment
      //! @return pose-optimized ensemble
      chemistry::FragmentEnsemble LigandLocalDock
      (
        const chemistry::FragmentEnsemble &ENSEMBLE,
        const descriptor::CheminfoProperty &SCORER,
        const chemistry::FragmentComplete &REFERENCE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // ReactionCombichem

  } // namespace app
} // namespace bcl

#endif // BCL_APP_REACTION_COMBICHEM_H_
