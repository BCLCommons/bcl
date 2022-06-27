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
#include "chemistry/bcl_chemistry_fragment_mutate_react.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateReact::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateReact())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateReact::FragmentMutateReact() :
        m_LigandBased( false),
        m_CorrectGeometry( false),
        m_CorrectNonReferenceRingGeometry( false),
        m_AdditionalAdjacentAtoms( size_t( 0))
    {
      // important to get options from base class and initialize reaction search
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief constructor with fragment react object
    FragmentMutateReact::FragmentMutateReact
    (
      const FragmentReact &REACT
    ) :
      FragmentReact( REACT),
      m_LigandBased( false),
      m_CorrectGeometry( false),
      m_CorrectNonReferenceRingGeometry( false),
      m_AdditionalAdjacentAtoms( size_t( 0))
    {
      // important to get options from base class and initialize reaction search
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief clone constructor
    FragmentMutateReact *FragmentMutateReact::Clone() const
    {
      return new FragmentMutateReact( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateReact::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateReact::GetAlias() const
    {
      static const std::string s_name( "React");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an SmallMolecule and returning a grown SmallMolecule
    //! @param FRAGMENT small molecule of interest
    //! @return Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateReact::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "FragmentReact!");

      // partially prevent sneaky access to reagent-only reactions
      if( !FRAGMENT.GetSize())
      {
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // add starting molecule to reagents
      m_Reagents.PushBack( FRAGMENT);

      // must run initialize to make the reaction structure tree
      m_ReactionSearch.Initialize();

      // collect valid reactions
      m_Reactions = ReactionEnsemble(
        storage::List< ReactionComplete>
        (
          m_ReactionSearch.GetReactions()->Begin(),
          m_ReactionSearch.GetReactions()->End()
        )
      );

      // Reset the reaction search object with the new molecules
      m_ReactionSearch.Reset();
      m_ReactionSearch = ReactionSearch
          (
            m_Reagents,
            m_Reactions
          );

      // re-initialize
      m_ReactionSearch.Initialize();

      // provide user messages
      if( !m_LigandBased)
      {
        BCL_MessageStd( "Setting pose-dependent options");
        BCL_MessageStd( "Ligand-based: " + util::Format()( m_LigandBased ? "true" : "false"));
        BCL_MessageStd( "Fix bad geometry: " + util::Format()( m_CorrectGeometry ? "true" : "false"));
        BCL_MessageStd("Fix bad ring geometry: " + util::Format()( m_CorrectNonReferenceRingGeometry ? "true" : "false"));
        BCL_MessageStd("Extend adjacent atoms: " + util::Format()( m_AdditionalAdjacentAtoms));
      }

      // try a few times
      for( size_t i( 0); i < m_NumberMaxAttempts; ++i)
      {
        // perform reaction
        auto products( ReactRandom( FRAGMENT));
        if( !products.Second().GetSize())
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // for cleaning and optimizing the new molecule
        AtomVector< AtomComplete> new_atom_vector( products.Second().GetMolecules().FirstElement().GetAtomVector());
        FragmentMapConformer cleaner
        (
          m_DrugLikenessType,
          m_MDL,
          FRAGMENT.GetMDLProperty( m_MDL),
          m_PropertyScorer,
          m_ResolveClashes,
          m_BFactors,
          m_Corina,
          storage::Vector< size_t>(),
          false,
          m_CorrectGeometry,
          m_AdditionalAdjacentAtoms
        );

        // remove hydrogen atoms so ease burden on the isomorphism search during cleaning
        static HydrogensHandler hydrogens_handler;
        hydrogens_handler.Remove( new_atom_vector);

        // also known as the "fuck it, just give me a conformer" option
        if( m_LigandBased)
        {
          // remove any restrictions to sampling
          FragmentComplete unclean_mol( new_atom_vector, "");
          unclean_mol.RemoveProperty( "SampleByParts");

          // finalize
          util::ShPtr< FragmentComplete> new_mol_ptr
          (
            m_ScaffoldFragment.GetSize() ?
                cleaner.Clean( unclean_mol.GetAtomVector(), m_ScaffoldFragment, m_DrugLikenessType) :
                cleaner.Clean( unclean_mol.GetAtomVector(), FRAGMENT, m_DrugLikenessType)
          );
          return math::MutateResult< FragmentComplete>( new_mol_ptr, *this);
        }

        // clean atoms only; 3D conformer was adjusted in ReactionWorker
        AtomVector< AtomComplete> clean_mol( cleaner.CleanAtoms( new_atom_vector, m_DrugLikenessType));
        FragmentComplete clean_labeled_mol( clean_mol, "");
        this->SetPropertiesFromOther( clean_labeled_mol, products.Second().GetMolecules().FirstElement());
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>( new FragmentComplete( clean_labeled_mol)), *this);
      }
      // if we never get it right then just return null
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // helper functions //
  //////////////////////

     io::Serializer FragmentMutateReact::GetSerializer() const
     {
       // initialize with base class serializers
       io::Serializer parameters( FragmentMutateInterface::GetSerializer());
       parameters.Merge( FragmentReact::GetSerializer());

       // add class description and parameters
       parameters.SetClassDescription
       (
         "Reacts a starting target fragment with one or more additional reagents"
       );

       parameters.AddInitializer
       (
         "reagents",
         "reagents that can be used to react with the input fragment according "
         "to reaction schemes available in 'reactions_directory'",
         io::Serialization::GetAgent( &m_ReagentsFilename),
         ""
       );

       parameters.AddInitializer
       (
         "reactions_directory",
         "path to the directory containing reaction (.rxn) files of interest",
         io::Serialization::GetAgent( &m_ReactionsDirectory),
         RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "functional_reactions"
       );

       parameters.AddInitializer
       (
         "allowed_reactant_positions",
         "reactant position indices (0-indexed) indicating which parts of the "
         "reaction(s) the target molecule is allowed to match; i.e., if you only "
         "want the target fragment to be considered a potential candidate for the "
         "second position of a reaction, pass '1'; default is to test all positions",
         io::Serialization::GetAgent( &m_TargetReactantPositionsStr),
         ""
       );

       parameters.AddInitializer
       (
         "ligand_based",
         "overrides 3D conformer settings to just produce an arbitrary conformer without preserving spatial information; "
         "by default all of the mutates derived from FragmentMutateInterface attempt to preserve as much coordinate "
         "information as possible while still yielding high quality 3D conformers. In ligand-based drug discovery we "
         "frequently do not care about real space coordinates. Because preserving the real space information is "
         "considerably more expensive than generating a de novo conformer with complex reaction-based mutates and "
         "the resulting pose is less likely to be biologically relevant (because usually with reactions we have "
         "multiple smaller fragments that we start with and frequently lose atoms along the way), it may be desirable to "
         "simply generate a random 3D conformer.",
         io::Serialization::GetAgent( &m_LigandBased),
         "false"
       );

       parameters.AddInitializer
       (
         "fix_geometry",
         "pose-dependent; "
         "if 3D conformer matters, fix atoms with bad geometry even if they are in reference structure",
         io::Serialization::GetAgent( &m_CorrectGeometry),
         "false"
       );

       parameters.AddInitializer
       (
         "fix_ring_geometry",
         "pose-dependent; "
         "if 3D conformer matters, add all ring atoms from non-reference scaffolds to mobile selection",
         io::Serialization::GetAgent( &m_CorrectNonReferenceRingGeometry),
         "false"
       );

       parameters.AddInitializer
       (
         "extend_adjacent_atoms",
         "pose-dependent; "
         "include adjacent atoms out this many bonds from any perturbed atom when generating a new 3D conformer",
         io::Serialization::GetAgent( &m_AdditionalAdjacentAtoms),
         "0"
       );

       return parameters;
     }

     //! @brief Set the members of this property from the given LABEL
     //! @param LABEL the label to parse
     //! @param ERROR_STREAM the stream to write errors to
     bool FragmentMutateReact::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
     {
       // static initialization check
       if( command::CommandState::IsInStaticInitialization())
       {
         return true;
       }

       // call RISH function of the base class
       if
       (
           !FragmentMutateInterface::ReadInitializerSuccessHook( LABEL, ERROR_STREAM) ||
           !FragmentReact::ReadInitializerSuccessHook( LABEL, ERROR_STREAM)
       )
       {
         return false;
       }

       // read in reagents
       if( m_ReagentsFilename.size())
       {
         s_Mutex.Lock();
         io::IFStream input;
         io::File::MustOpenIFStream( input, m_ReagentsFilename);
         m_Reagents.ReadMoreFromMdl( input);
         io::File::CloseClearFStream( input);
         s_Mutex.Unlock();
       }

       // we must have reagents
       BCL_Assert( m_Reagents.GetSize(), "No reagents! Exiting...");

       // we need some reactions
       if( m_Reactions.GetSize())
       {
         // get the default reactions directory if nothing is specified
         if( !m_ReactionsDirectory.size())
         {
           m_ReactionsDirectory = RotamerLibraryFile::GetRotamerFinder().FindFile( "") + ( "functional_reactions");
           BCL_MessageStd( "No reagents specified. Defaulting to reactions in " + util::Format()( m_ReactionsDirectory));
         }
       }

       // this is the easiest way to read in the reactions given the current reaction framework
       m_ReactionSearch = ReactionSearch
       (
         m_ReagentsFilename,
         m_ReactionsDirectory
       );

       // read in reactant position indices
       if( m_TargetReactantPositionsStr.size())
       {
         m_TargetReactantPositions.Reset();
         m_TargetReactantPositions = util::SplitStringToNumerical< size_t>( m_TargetReactantPositionsStr);
       }

       // set ligand- or structure-based conformer reqs
       m_ReactionWorker.SetProductConformerArbitrary( m_LigandBased);

       // for pose-dependent conformers set some details
       if( !m_LigandBased)
       {
         // set pose-dependent options
         m_ReactionWorker.SetCorrectGeometry( m_CorrectGeometry);
         m_ReactionWorker.SetCorrectNonReferenceRingGeometry( m_CorrectNonReferenceRingGeometry);
         m_ReactionWorker.SetAdditionalAdjacentAtoms( m_AdditionalAdjacentAtoms);
       }

       // done
       return true;
     }

  } // namespace chemistry
} // namespace bcl
