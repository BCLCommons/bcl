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

#ifndef BCL_APP_LINK_FRAGMENTS_H_
#define BCL_APP_LINK_FRAGMENTS_H_

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
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
    //! @class LinkFragments
    //! @brief Links fragments together to create larger molecules
    //! @details Provided one or more SDF files, links together fragments in a pairwise manner either enumeratively or at specified atom indices.
    //! Linking occurs either directly, or with the specified linker type and number of linker repeats
    //!
    //! @author brownbp1
    //! @date March 20, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LinkFragments :
      public InterfaceRelease
    {
    public:

      // static member of this class
      static const ApplicationType LinkFragments_Instance;

    private:

    //////////
    // data //
    //////////

      //! molecules to be linked
      util::ShPtr< command::ParameterInterface> m_InputMoleculesAFlag;

      //! molecules to be linked
      util::ShPtr< command::ParameterInterface> m_InputMoleculesBFlag;

      //! atom indices in molecules from m_InputMoleculesA
      util::ShPtr< command::FlagInterface> m_LinkIndicesAFlag;

      //! //! atom indices in molecules from InputMoleculesA
      util::ShPtr< command::FlagInterface> m_LinkIndicesBFlag;

      //! output SDF file
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! directly link the fragments at provided indices
      util::ShPtr< command::FlagInterface> m_DirectLinkFlag;

      //! link the fragments via a single element
      util::ShPtr< command::FlagInterface> m_SingleElementLinkFlag;

      //! link the fragments via an alkyl chain
      util::ShPtr< command::FlagInterface> m_AlkylLinkFlag;

      //! link the fragments via methoxy repeats
      util::ShPtr< command::FlagInterface> m_MethoxyLinkFlag;

      //! link the fragments via ethoxy repeats
      util::ShPtr< command::FlagInterface> m_EthoxyLinkFlag;

      //! link the fragments via amide repeats
      util::ShPtr< command::FlagInterface> m_AmideLinkFlag;

      //! the protein-ligand interface scoring function to use for pose refinement
      util::ShPtr< command::FlagInterface> m_LigandLocalDockFlag;

      //! the MDL property name specifying the protein pocket filename
      //! note that this must match the model
      util::ShPtr< command::FlagInterface> m_MDLStringFlag;

      //! the drug-likeness metric to use
      util::ShPtr< command::FlagInterface> m_DrugLikenessTypeFlag;

      //! flag indication of first molecule to load from ensemble A
      util::ShPtr< command::FlagInterface> m_StartA;

      //! flag indication of first molecule to load from ensemble B
      util::ShPtr< command::FlagInterface> m_StartB;

      //! flag indication of last molecule to load from ensemble A
      util::ShPtr< command::FlagInterface> m_MaxMolsA;

      //! flag indication of last molecule to load from ensemble B
      util::ShPtr< command::FlagInterface> m_MaxMolsB;

      //! the MDL property name specifying the protein pocket filename
      mutable std::string m_MDLString;

      //! get the filename of the protein pocket
      mutable std::string m_PocketFilename;

      //! link point atom indices for each ensemble
      mutable storage::Vector< size_t> m_LinkAtomIndicesA;
      mutable storage::Vector< size_t> m_LinkAtomIndicesB;

      // member data for threading and efficiency - modeled after molecule:Compare
      //! # pairs that have been considered
      mutable size_t m_NumberPairsConsidered;

      //! # pairs that will be considered
      mutable size_t m_NumberPairsToConsider;

      //! molecule ensembles
      mutable storage::Vector< chemistry::FragmentComplete> m_EnsembleA;
      mutable storage::Vector< chemistry::FragmentComplete> m_EnsembleB;
      mutable storage::Vector< chemistry::FragmentEnsemble> m_LinkedMolEnsembles;

      //! ensemble sizes
      mutable size_t m_EnsembleASize;
      mutable size_t m_EnsembleBSize;

      //! to avoid repeating fragment combinations
      mutable bool m_IdenticalEnsembles;

      //! descriptor for pose scoring
      mutable descriptor::CheminfoProperty m_LigandDockScorer;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      LinkFragments();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      LinkFragments *Clone() const;

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

      //! @brief directly link two fragments
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDICES_A link indices in first molecule
      //! @param LINK_INDICES_B link indices in second molecule
      //! @return the newly generated molecules
      chemistry::FragmentEnsemble DirectLink
      (
        const chemistry::FragmentComplete &FRAGMENT_A,
        const chemistry::FragmentComplete &FRAGMENT_B,
        const storage::Vector< size_t> &LINK_INDICES_A,
        const storage::Vector< size_t> &LINK_INDICES_B
      ) const;

      //! @brief link two fragments via a single element
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDICES_A link indices in first molecule
      //! @param LINK_INDICES_B link indices in second molecule
      //! @param ELEMENT_TYPE element type serving as a link
      //! @return the newly generated molecules
      chemistry::FragmentEnsemble SingleElementLink
      (
        const chemistry::FragmentComplete &FRAGMENT_A,
        const chemistry::FragmentComplete &FRAGMENT_B,
        const storage::Vector< size_t> &LINK_INDICES_A,
        const storage::Vector< size_t> &LINK_INDICES_B,
        const storage::Vector< std::string> &ELEMENT_TYPE
      ) const;

      //! @brief link two fragments with an alkyl chain
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDICES_A link indices in first molecule
      //! @param LINK_INDICES_B link indices in second molecule
      //! @param REPEATS the number of linker repeats
      //! @return the newly generated molecules
      chemistry::FragmentEnsemble AlkylLink
      (
        const chemistry::FragmentComplete &FRAGMENT_A,
        const chemistry::FragmentComplete &FRAGMENT_B,
        const storage::Vector< size_t> &LINK_INDICES_A,
        const storage::Vector< size_t> &LINK_INDICES_B,
        const size_t &REPEATS
      ) const;

      //! @brief link two fragments with a methoxy repeat
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDICES_A link indices in first molecule
      //! @param LINK_INDICES_B link indices in second molecule
      //! @param REPEATS the number of linker repeats
      //! @return the newly generated molecules
      chemistry::FragmentEnsemble MethoxyLink
      (
        const chemistry::FragmentComplete &FRAGMENT_A,
        const chemistry::FragmentComplete &FRAGMENT_B,
        const storage::Vector< size_t> &LINK_INDICES_A,
        const storage::Vector< size_t> &LINK_INDICES_B,
        const size_t &REPEATS
      ) const;

      //! @brief link two fragments with an ethoxy repeat
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDICES_A link indices in first molecule
      //! @param LINK_INDICES_B link indices in second molecule
      //! @param REPEATS the number of linker repeats
      //! @return the newly generated molecules
      chemistry::FragmentEnsemble EthoxyLink
      (
        const chemistry::FragmentComplete &FRAGMENT_A,
        const chemistry::FragmentComplete &FRAGMENT_B,
        const storage::Vector< size_t> &LINK_INDICES_A,
        const storage::Vector< size_t> &LINK_INDICES_B,
        const size_t &REPEATS
      ) const;

      //! @brief link two fragments with an amide repeat
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDICES_A link indices in first molecule
      //! @param LINK_INDICES_B link indices in second molecule
      //! @param REPEATS the number of linker repeats
      //! @return the newly generated molecules
      chemistry::FragmentEnsemble AmideLink
      (
        const chemistry::FragmentComplete &FRAGMENT_A,
        const chemistry::FragmentComplete &FRAGMENT_B,
        const storage::Vector< size_t> &LINK_INDICES_A,
        const storage::Vector< size_t> &LINK_INDICES_B,
        const size_t &REPEATS
      ) const;

      //! @brief link two fragments with a urea repeat
      //! @param FRAGMENT_A first molecule
      //! @param FRAGMENT_B second molecule
      //! @param LINK_INDICES_A link indices in first molecule
      //! @param LINK_INDICES_B link indices in second molecule
      //! @param REPEATS the number of linker repeats
      //! @return the newly generated molecules
      chemistry::FragmentEnsemble UreaLink
      (
        const chemistry::FragmentComplete &FRAGMENT_A,
        const chemistry::FragmentComplete &FRAGMENT_B,
        const storage::Vector< size_t> &LINK_INDICES_A,
        const storage::Vector< size_t> &LINK_INDICES_B,
        const size_t &REPEATS
      ) const;

      //! @brief link two fragments with an amide repeat
      //! @param CONNECTION whether to connect the amide to molecule A via C or N or both
      //! @param REPEATS the number of linker repeats
      //! @return the newly generated molecules
      chemistry::FragmentEnsemble GenerateAmideLinker
      (
        const std::string &CONNECTION,
        const size_t &REPEATS
      ) const;

      //! @brief remove a hydrogen atom from a target atom
      //! @param FRAGMENT the molecule of interest
      //! @param ATOM_INDEX the index of the atom in the molecule of interest
      //! @return the new molecule and index of the desired atom
      storage::Pair< chemistry::FragmentComplete, size_t> OpenValence
      (
        const chemistry::FragmentComplete &FRAGMENT,
        const size_t &ATOM_INDEX
      ) const;

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
      //! @return pose-optimized ensemble
      chemistry::FragmentEnsemble LigandLocalDock
      (
        const chemistry::FragmentEnsemble &ENSEMBLE,
        const descriptor::CheminfoProperty &SCORER
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

    }; // LinkFragments

  } // namespace app
} // namespace bcl

#endif // BCL_APP_LINK_FRAGMENTS_H_
