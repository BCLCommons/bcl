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
#include "chemistry/bcl_chemistry_configurational_bond_types.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_stochastic_pose_optimizer.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_default.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_directory_entry.h"
#include "sched/bcl_sched_job_interface.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "bcl_app_link_fragments.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality.h"
#include "assemble/bcl_assemble_voxel_grid_atom.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extensions_file_existence.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {

    // Static instance initialization
    const ApplicationType LinkFragments::LinkFragments_Instance
    (
      GetAppGroups().AddAppToGroup( new LinkFragments(), GetAppGroups().e_ChemInfo)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    LinkFragments::LinkFragments() :
      m_InputMoleculesAFlag
      (
        new command::Parameter
        (
          "input_filenames_a",
          "if provided without input_filenames_b, pairwise combination of all fragments in the provided file will occur",
          command::ParameterCheckFileExistence()
        )
      ),
      m_InputMoleculesBFlag
      (
        new command::Parameter
        (
          "input_filenames_b",
          "fragments to be connected to fragments in input_filenames_a",
          command::ParameterCheckFileExistence()
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_filename",
          "linked molecules",
          command::Parameter
          (
            "output_filename",
            "SDF containing linked molecules",
            "linked_molecules.sdf.gz"
          )
        )
      ),
      m_LinkIndicesAFlag
      (
        new command::FlagDynamic
        (
          "link_atoms_a",
          "atoms indices from molecules in input_filenames_a at which to link fragments (0-indexed)",
          command::Parameter
          (
            "link_atoms",
            "if not atoms are provided then the default is to consider "
            "all atoms in the maximum common substructure of the molecules in Ensemble A,"
            "or if the ensemble only contains a single molecule every heavy atom",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
            ""
          )
        )
      ),
      m_LinkIndicesBFlag
      (
        new command::FlagDynamic
        (
          "link_atoms_b",
          "atoms indices from molecules in input_filenames_b at which to link fragments (0-indexed)",
          command::Parameter
          (
            "link_atoms",
            "if not atoms are provided then the default is to consider "
            "all atoms in the maximum common substructure of the molecules in Ensemble B,"
            "or if the ensemble only contains a single molecule every heavy atom",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
            ""
          )
        )
      ),
      m_DirectLinkFlag
      (
        new command::FlagStatic
        (
          "direct_link",
          "connect fragments directly"
        )
      ),
      m_SingleElementLinkFlag
      (
        new command::FlagDynamic
        (
          "single_element_linker",
          "connect fragments with alkyl linkers",
          command::Parameter
          (
            "element_type",
            "the element to insert as a linker",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create
              (
                "C", "O", "N", "S"
              )
            ),
            "C"
          )
        )
      ),
      m_AlkylLinkFlag
      (
        new command::FlagStatic
        (
          "alkyl_linker",
          "connect fragments with alkyl linkers",
          command::Parameter
          (
            "repeats",
            "number of saturated carbon repeats",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "1"
          )
        )
      ),
      m_MethoxyLinkFlag
      (
        new command::FlagStatic
        (
          "methoxy_linker",
          "connect fragments with alkyl linkers",
          command::Parameter
          (
            "repeats",
            "number of methoxy repeats",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "1"
          )
        )
      ),
      m_EthoxyLinkFlag
      (
        new command::FlagStatic
        (
          "ethoxy_linker",
          "connect fragments with alkyl linkers",
          command::Parameter
          (
            "repeats",
            "number of ethoxy repeats",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "1"
          )
        )
      ),
      m_AmideLinkFlag
      (
        new command::FlagStatic
        (
          "amide_linker",
          "connect fragments with alkyl linkers",
          command::Parameter
          (
            "repeats",
            "number of amide repeats",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "1"
          )
        )
      ),
      m_LigandLocalDockFlag
      (
        new command::FlagDynamic
        (
          "ligand_dock_score",
          "the protein-ligand interface scoring function to use for pose refinement",
          command::Parameter
          (
            "function",
            "the scoring function implementation to use"
          ),
          0,
          1
        )
      ),
      m_MDLStringFlag
      (
        new command::FlagDynamic
        (
          "MDL_property",
          "the MDL property specifying the protein binding pocket used for pose-dependent scoring",
          command::Parameter
          (
            "property_name",
            "the name of the property"
          ),
          0,
          1
        )
      ),
      m_DrugLikenessTypeFlag
      (
        new command::FlagStatic
        (
          "druglikeness_type", "the type of druglikeness filter to use to determine when a molecule is skipped by the Monte Carlo algorithm",
          command::Parameter
          (
            "type", "the type of druglikenes to use",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "IsConstitutionDruglike", "IsConstitutionDruglikeAndHitlike", "None")),
            "IsConstitutionDruglike"
          )
        )
      ),
      m_StartA
      (
        new command::FlagStatic
        (
          "ensemble_a_start",
          "flag for indicating which mol from ensemble a to load in first",
          command::Parameter
          (
            "index",
            "index of ensemble a molecules to start with",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            "0"
          )
        )
      ),
      m_StartB
      (
        new command::FlagStatic
        (
          "ensemble_b_start",
          "flag for indicating which mol from ensemble a to load in first",
          command::Parameter
          (
            "index",
            "index of ensemble b molecules to start with",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            "0"
          )
        )
      ),
      m_MaxMolsA
      (
        new command::FlagStatic
        (
          "ensemble_a_max",
          "flag for indicating maximum number of molecules to take from ensemble a",
          command::Parameter
          (
            "max",
            "max number of molecules from ensemble a",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            util::Format()( math::GetHighestBoundedValue< size_t>())
          )
        )
      ),
      m_MaxMolsB
      (
        new command::FlagStatic
        (
          "ensemble_b_max",
          "flag for indicating maximum number of molecules to take from ensemble b",
          command::Parameter
          (
            "max",
            "max number of molecules from ensemble b",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            util::Format()( math::GetHighestBoundedValue< size_t>())
          )
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    LinkFragments *LinkFragments::Clone() const
    {
      return new LinkFragments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LinkFragments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string LinkFragments::GetDescription() const
    {
      return "connect molecule fragments with various linkers";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &LinkFragments::GetReadMe() const
    {
      static std::string s_read_me =
        "LinkFragments connects molecule fragments at specified atom indices "
        "via any of the following options: "
        "direct linkage, single element (N, O, C) linkage, alkyl linkage of specified length, "
        "and more.";
      return s_read_me;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> LinkFragments::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // molecule i/o
      sp_cmd->AddParameter( m_InputMoleculesAFlag);
      sp_cmd->AddParameter( m_InputMoleculesBFlag);
      sp_cmd->AddFlag( m_OutputFilenameFlag);

      // specify link indices
      sp_cmd->AddFlag( m_LinkIndicesAFlag);
      sp_cmd->AddFlag( m_LinkIndicesBFlag);

      // linker types
      sp_cmd->AddFlag( m_DirectLinkFlag);
      sp_cmd->AddFlag( m_SingleElementLinkFlag);
      sp_cmd->AddFlag( m_AlkylLinkFlag);
      sp_cmd->AddFlag( m_MethoxyLinkFlag);
      sp_cmd->AddFlag( m_EthoxyLinkFlag);
      sp_cmd->AddFlag( m_AmideLinkFlag);

      // protein-ligand interface scorer for pose refinement
      sp_cmd->AddFlag( m_LigandLocalDockFlag);
      sp_cmd->AddFlag( m_MDLStringFlag);

      // drug-likeness filter
      sp_cmd->AddFlag( m_DrugLikenessTypeFlag);

      // add defaults
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled command object
      return sp_cmd;
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief the Main function
    //! @return error code - 0 for success
    int LinkFragments::Main() const
    {

      // read in ensemble_a
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_InputMoleculesAFlag->GetValue());

      // get the starting molecule for the ensemble_a input
      math::Range< size_t> ens_a_load_rng( size_t( 0), math::GetHighestBoundedValue< size_t>());
      ens_a_load_rng.SetMin( util::ConvertStringToNumericalValue< size_t>( m_StartA->GetFirstParameter()->GetValue()));

      // set max number of molecules to load for ensemble_a
      const size_t n_to_load_a( util::ConvertStringToNumericalValue< size_t>( m_MaxMolsA->GetFirstParameter()->GetValue()));
      ens_a_load_rng.SetMax
      (
        math::GetHighestBoundedValue< size_t>() - n_to_load_a > ens_a_load_rng.GetMin() ?
            ens_a_load_rng.GetMin() + n_to_load_a - 1 :
            math::GetHighestBoundedValue< size_t>()
      );

      // get command-line preference for hydrogen atom handling
      chemistry::FragmentEnsemble ensemble_a( input, sdf::GetCommandLineHydrogensPref(), ens_a_load_rng);
      if( m_MDLStringFlag->GetFlag())
      {
        m_MDLString = m_MDLStringFlag->GetFirstParameter()->GetValue();
        m_PocketFilename = ensemble_a.GetMolecules().Begin()->GetStoredPropertiesNonConst().GetMDLProperty( m_MDLString);
      }
      if( m_LigandLocalDockFlag->GetFlag())
      {
        BCL_Assert
        (
          m_PocketFilename.size(),
          "Either no MDL property was specified for the protein binding pocket filename, or "
          "the specified MDL property did not contain a filename. Check your input SDF and try again."
        );
      }

      // close input stream
      io::File::CloseClearFStream( input);

      // finalize ensemble_a
      m_EnsembleA = storage::Vector< chemistry::FragmentComplete>( ensemble_a.Begin(), ensemble_a.End());
      m_EnsembleASize = m_EnsembleA.GetSize();

      // test whether only one ensemble was given or whether the same filename was given for both ensembles
      m_IdenticalEnsembles =
        !m_InputMoleculesBFlag->GetWasSetInCommandLine()
        ||
        (
          ( io::DirectoryEntry( m_InputMoleculesBFlag->GetValue()).GetFullName() == io::DirectoryEntry( m_InputMoleculesAFlag->GetValue()).GetFullName())
          && m_StartA->GetFirstParameter()->GetValue() == m_StartB->GetFirstParameter()->GetValue()
          && m_MaxMolsA->GetFirstParameter()->GetValue() == m_MaxMolsB->GetFirstParameter()->GetValue()
        );

      // initialize ensemble_b
      chemistry::FragmentEnsemble ensemble_b;

      // if the ensembles are not identical then we need to read in ensemble_b the same way we read in ensemble_a
      if( !m_IdenticalEnsembles)
      {
        // read in ensemble_b
        io::File::MustOpenIFStream( input, m_InputMoleculesBFlag->GetValue());

        // get the starting molecule for the ensemble_b input
        math::Range< size_t> ens_b_load_rng( size_t( 0), math::GetHighestBoundedValue< size_t>());
        ens_b_load_rng.SetMin( util::ConvertStringToNumericalValue< size_t>( m_StartB->GetFirstParameter()->GetValue()));

        // set max number of molecules to load for ensemble_b
        const size_t n_to_load_b( util::ConvertStringToNumericalValue< size_t>( m_MaxMolsB->GetFirstParameter()->GetValue()));
        ens_b_load_rng.SetMax
        (
          math::GetHighestBoundedValue< size_t>() - n_to_load_b > ens_b_load_rng.GetMin() ?
              ens_b_load_rng.GetMin() + n_to_load_b - 1 :
              math::GetHighestBoundedValue< size_t>()
        );

        // get command-line preference for hydrogen atom handling
        ensemble_b.ReadMoreFromMdl( input, sdf::GetCommandLineHydrogensPref(), ens_b_load_rng);

        // close input stream
        io::File::CloseClearFStream( input);

        // finalize ensemble_b
        m_EnsembleB = storage::Vector< chemistry::FragmentComplete>( ensemble_b.Begin(), ensemble_b.End());
        m_EnsembleBSize = m_EnsembleB.GetSize();

        // set # of pairs to consider
        m_NumberPairsToConsider = m_EnsembleASize * m_EnsembleBSize;
      }

      // if the ensembles are identical then set our number of pairs to consider to the correct number
      else
      {
        m_EnsembleB = m_EnsembleA;
        m_EnsembleBSize = m_EnsembleASize;

        // set # of pairs to consider
        m_NumberPairsToConsider = m_EnsembleASize * ( m_EnsembleASize + 1) / 2;
      }

      // if the files are provided then they cannot be empty
      if( m_EnsembleASize == 0 || m_EnsembleBSize == 0)
      {
        BCL_MessageCrt( "One of the provided files did not contain any molecules; no output will be given");
        return 0;
      }

      // initialize atom indices for each ensemble (assumes that the indices are common within an ensemble)
      // A
      if( m_LinkIndicesAFlag->GetFlag())
      {
        // convert the functionalization points to numeric values
        storage::Vector< std::string> fxnl_pts_str_a( m_LinkIndicesAFlag->GetStringList());

        // output linkable indices to terminal
        std::string link_atoms_a;
        for
        (
            auto itr( fxnl_pts_str_a.Begin()), itr_end( fxnl_pts_str_a.End());
            itr != itr_end;
            ++itr
        )
        {
          link_atoms_a.append(util::Format()( *itr));
          link_atoms_a.append(",");
        }
        BCL_MessageStd( "Linkable atom indices in molecule ensemble A: " + util::Format()( link_atoms_a));

        storage::Set< size_t> fxnl_pts_set_a;

        // A
        for( size_t i( 0), l( fxnl_pts_str_a.GetSize()); i < l; ++i)
        {
          size_t point;
          if( !util::TryConvertFromString( point, fxnl_pts_str_a( i), util::GetLogger()))
          {
            BCL_MessageStd( "Could not parse \"" + fxnl_pts_str_a( i) + "\" as a number");
            continue;
          }
//          if( point < m_EnsembleASize)
//          {
            fxnl_pts_set_a.Insert( point);
//          }
//          else
//          {
//            BCL_MessageStd
//            (
//              "Warning: specified point \"" + util::Format()( point) + "\""
//              " has an index greater than the number of atoms in the molecule, not using this point"
//            );
//          }
        }
        m_LinkAtomIndicesA = storage::Vector< size_t>( fxnl_pts_set_a.Begin(), fxnl_pts_set_a.End());
      }

      // B
      if( m_LinkIndicesBFlag->GetFlag())
      {
        // convert the functionalization points to numeric values
        storage::Vector< std::string> fxnl_pts_str_b( m_LinkIndicesBFlag->GetStringList());

        // output linkable indices to terminal
        std::string link_atoms_b;
        for
        (
            auto itr( fxnl_pts_str_b.Begin()), itr_end( fxnl_pts_str_b.End());
            itr != itr_end;
            ++itr
        )
        {
          link_atoms_b.append(util::Format()( *itr));
          link_atoms_b.append(",");
        }
        BCL_MessageStd( "Linkable atom indices in molecule ensemble B: " + util::Format()( link_atoms_b));

        storage::Set< size_t> fxnl_pts_set_b;

        // B
        for( size_t i( 0), l( fxnl_pts_str_b.GetSize()); i < l; ++i)
        {
          size_t point;
          if( !util::TryConvertFromString( point, fxnl_pts_str_b( i), util::GetLogger()))
          {
            BCL_MessageStd( "Could not parse \"" + fxnl_pts_str_b( i) + "\" as a number");
            continue;
          }
//          if( point < m_EnsembleBSize)
//          {
            fxnl_pts_set_b.Insert( point);
//          }
//          else
//          {
//            BCL_MessageStd
//            (
//              "Warning: specified point \"" + util::Format()( point) + "\""
//              " has an index greater than the number of atoms in the molecule, not using this point"
//            );
//          }
        }
        m_LinkAtomIndicesB = storage::Vector< size_t>( fxnl_pts_set_b.Begin(), fxnl_pts_set_b.End());
      }

      // initialize number of pairs of molecules that need to be connected
      m_NumberPairsConsidered = size_t( 0);

      // link fragments of non-identical input files
        // go over each molecule in ensemble A
        for
        (
            auto ens_a_itr( m_EnsembleA.Begin()), ens_a_itr_end( m_EnsembleA.End());
            ens_a_itr != ens_a_itr_end;
            ++ens_a_itr
        )
        {
          // if not previously specified, set link atom indices for current molecule
          if( !m_LinkIndicesAFlag->GetFlag())
          {
            m_LinkAtomIndicesA.Reset();
            for
            (
                auto a_itr( ens_a_itr->GetAtomVector().Begin()), a_itr_end( ens_a_itr->GetAtomVector().End());
                a_itr != a_itr_end;
                ++a_itr
            )
            {
              if( a_itr->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
              {
                continue;
              }
              m_LinkAtomIndicesA.PushBack( ens_a_itr->GetAtomVector().GetAtomIndex( *a_itr));
            }
          }

          // go over each molecule in ensemble B
          for
          (
              auto ens_b_itr( m_EnsembleB.Begin()), ens_b_itr_end( m_EnsembleB.End());
              ens_b_itr != ens_b_itr_end;
              ++ens_b_itr
          )
          {
            // if not previously specified, set link atom indices for current molecule
            if( !m_LinkIndicesBFlag->GetFlag())
            {
              m_LinkAtomIndicesB.Reset();
              for
              (
                  auto b_itr( ens_b_itr->GetAtomVector().Begin()), b_itr_end( ens_b_itr->GetAtomVector().End());
                  b_itr != b_itr_end;
                  ++b_itr
              )
              {
                if( b_itr->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
                {
                  continue;
                }
                m_LinkAtomIndicesB.PushBack( ens_b_itr->GetAtomVector().GetAtomIndex( *b_itr));
              }
            }

            // direct links
            if( m_DirectLinkFlag->GetFlag())
            {
              m_LinkedMolEnsembles.PushBack
              (
                chemistry::FragmentEnsemble
                (
                  DirectLink
                  (
                    *ens_a_itr,
                    *ens_b_itr,
                    m_LinkAtomIndicesA,
                    m_LinkAtomIndicesB
                  )
                )
              );
            } // end direct link
            // single element links
            if( m_SingleElementLinkFlag->GetFlag())
            {
              m_LinkedMolEnsembles.PushBack
              (
                chemistry::FragmentEnsemble
                (
                  SingleElementLink
                  (
                    *ens_a_itr,
                    *ens_b_itr,
                    m_LinkAtomIndicesA,
                    m_LinkAtomIndicesB,
                    m_SingleElementLinkFlag->GetStringList() // element type
                  )
                )
              );
            } // end single element link
            if( m_AmideLinkFlag->GetFlag())
            {
              m_LinkedMolEnsembles.PushBack
              (
                chemistry::FragmentEnsemble
                (
                  AmideLink
                  (
                    *ens_a_itr,
                    *ens_b_itr,
                    m_LinkAtomIndicesA,
                    m_LinkAtomIndicesB,
                    m_AmideLinkFlag->GetFirstParameter()->GetNumericalValue< size_t>()
                  )
                )
              );
            } // end amide link
          } // end ensemble B
        } // end ensemble A

      // output molecules and done
      io::OFStream output;
      io::File::MustOpenOFStream( output, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
      for
      (
          auto final_ens_itr( m_LinkedMolEnsembles.Begin()), final_ens_itr_end( m_LinkedMolEnsembles.End());
          final_ens_itr != final_ens_itr_end;
          ++final_ens_itr
      )
      {
        //refine docked pose if relevant
        chemistry::FragmentEnsemble mols;
        if( m_LigandLocalDockFlag->GetFlag())
        {
          m_LigandDockScorer = m_LigandLocalDockFlag->GetFirstParameter()->GetValue();
          mols = chemistry::FragmentEnsemble( LigandLocalDock( *final_ens_itr, m_LigandDockScorer));
        }
        else
        {
          mols = *final_ens_itr;
        }
        for
        (
            auto mol_itr( mols.Begin()), mol_itr_end( mols.End());
            mol_itr != mol_itr_end;
            ++mol_itr
        )
        {
          BCL_MessageStd( "Writing output!");
          mol_itr->WriteMDL( output);
        }
      }
      io::File::CloseClearFStream( output);
      return 0;
    }

  //////////////////////
  // helper functions //
  /////////////////////

     //! @brief directly link two fragments
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDICES_A link indices in first molecule
     //! @param LINK_INDICES_B link indices in second molecule
     //! @return the newly generated molecules
     chemistry::FragmentEnsemble LinkFragments::DirectLink
     (
       const chemistry::FragmentComplete &FRAGMENT_A,
       const chemistry::FragmentComplete &FRAGMENT_B,
       const storage::Vector< size_t> &LINK_INDICES_A,
       const storage::Vector< size_t> &LINK_INDICES_B
     ) const
     {
       BCL_MessageCrt( "Begin directly linking fragments!");

       // initialize output ensemble
       chemistry::FragmentEnsemble linked_mols;

       // get atom vectors
       chemistry::AtomVector< chemistry::AtomComplete>
       atom_v_a( FRAGMENT_A.GetAtomVector()), atom_v_b( FRAGMENT_B.GetAtomVector());

       // connect atoms from fragment a
       for
       (
           auto a_index_itr( LINK_INDICES_A.Begin()), a_index_itr_end( LINK_INDICES_A.End());
           a_index_itr != a_index_itr_end;
           ++a_index_itr
       )
       {
         // cannot connect hydrogen atoms to anything
         if( atom_v_a( *a_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
         {
           BCL_MessageStd( "Whoops! Specified atom index in A is a hydrogen atom! Ignoring index...");
           continue;
         }

         // connect atoms from from fragment b
         for
         (
             auto b_index_itr( LINK_INDICES_B.Begin()), b_index_itr_end( LINK_INDICES_B.End());
             b_index_itr != b_index_itr_end;
             ++b_index_itr
         )
         {
           // cannot connect hydrogen atoms to anything
           if( atom_v_b( *b_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
           {
             BCL_MessageStd( "Whoops! Specified atom index in B is a hydrogen atom! Ignoring index...");
             continue;
           }

           // open valences
           storage::Pair< chemistry::FragmentComplete, size_t> pair_a( OpenValence( FRAGMENT_A, *a_index_itr));
           storage::Pair< chemistry::FragmentComplete, size_t> pair_b( OpenValence( FRAGMENT_B, *b_index_itr));

           // link fragments
           storage::Pair< bool, chemistry::FragmentComplete> new_fragment
           (
             chemistry::MergeFragmentComplete::MergeFragments
             (
               pair_a.First(),
               pair_b.First(),
               chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
               storage::Pair< size_t, size_t>( pair_a.Second(), pair_b.Second())
             )
           );

           // clean the new molecule
           BCL_MessageStd( "Clean the molecule");
           util::ShPtr< chemistry::FragmentComplete> clean_frag( CallCleaner( new_fragment.Second(), FRAGMENT_A));

           // add new molecule to the ensemble
           if( new_fragment.First() && clean_frag.IsDefined())
           {
             linked_mols.PushBack( *clean_frag);
           }
         }
       }
       BCL_MessageCrt( "Done directly linking fragments!");
       return linked_mols;
     }

     //! @brief link two fragments via a single element
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDICES_A link indices in first molecule
     //! @param LINK_INDICES_B link indices in second molecule
     //! @param ELEMENT_TYPE element type serving as a link
     //! @return the newly generated molecules
     chemistry::FragmentEnsemble LinkFragments::SingleElementLink
     (
       const chemistry::FragmentComplete &FRAGMENT_A,
       const chemistry::FragmentComplete &FRAGMENT_B,
       const storage::Vector< size_t> &LINK_INDICES_A,
       const storage::Vector< size_t> &LINK_INDICES_B,
       const storage::Vector< std::string> &ELEMENT_TYPE
     ) const
     {
       // initialize output ensemble
       chemistry::FragmentEnsemble linked_mols;

       // loop over the element type string
       for
       (
           auto string_itr( ELEMENT_TYPE.Begin()), string_itr_end( ELEMENT_TYPE.End());
           string_itr != string_itr_end;
           ++string_itr
       )
       {
         BCL_MessageCrt( "Begin linking fragments with single element: " + util::Format()( *string_itr));

         // get link atom
         sdf::AtomInfo link_atom;

         // carbon
         if( *string_itr == "C")
         {
           link_atom = sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_NonChiral);
         }
         // oxygen
         else if( *string_itr == "O")
         {
           link_atom = sdf::AtomInfo( chemistry::GetAtomTypes().O_Te2Te2TeTe, chemistry::e_NonChiral);
         }
         // nitrogen
         else if( *string_itr == "N")
         {
           link_atom = sdf::AtomInfo( chemistry::GetAtomTypes().N_Te2Te2TeTe, chemistry::e_NonChiral);
         }
         // sulfur
         else if( *string_itr == "S")
         {
           link_atom = sdf::AtomInfo( chemistry::GetAtomTypes().S_Te2Te2TeTe, chemistry::e_NonChiral);
         }
         // this is not needed because of the constructor parameter check, but oh well
         else
         {
           BCL_MessageStd( "Invalid element type specified for link atom");
           return chemistry::FragmentEnsemble();
         }

         // connect atoms from fragment a
         for
         (
             auto a_index_itr( LINK_INDICES_A.Begin()), a_index_itr_end( LINK_INDICES_A.End());
             a_index_itr != a_index_itr_end;
             ++a_index_itr
         )
         {

           // get atom vector
           chemistry::AtomVector< chemistry::AtomComplete> atom_v_a( FRAGMENT_A.GetAtomVector());

           // cannot connect hydrogen atoms to anything
           if( atom_v_a( *a_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
           {
             BCL_MessageStd( "Whoops! Specified atom index in A is a hydrogen atom! Ignoring index...");
             continue;
           }

//           // skip obviously bad heavy atom connections; this would get processed with the drug-likeness filter,
//           // but I don't even want it to get that far. It is a waste of time.
//           if
//           (
//               link_atom.GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Sulfur &&
//               atom_v_a( *a_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Sulfur
//           )
//           {
//             continue;
//           }
//           else if
//           (
//               link_atom.GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Oxygen &&
//               atom_v_a( *a_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Oxygen
//           )
//           {
//             continue;
//           }
//           else if
//           (
//               link_atom.GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Nitrogen &&
//               atom_v_a( *a_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Nitrogen
//           )
//           {
//             continue;
//           }
//           else if
//           (
//               (link_atom.GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Nitrogen &&
//               atom_v_a( *a_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Oxygen) ||
//               (link_atom.GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Oxygen &&
//               atom_v_a( *a_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Nitrogen)
//           )
//           {
//             continue;
//           }

           // add link atom to fragment A
           storage::Vector< sdf::AtomInfo> atominfo_a( atom_v_a.GetAtomInfo());
           storage::Vector< sdf::BondInfo> bondinfo( atom_v_a.GetBondInfo());
           atominfo_a.PushBack( link_atom);
           chemistry::AtomVector< chemistry::AtomComplete> atom_v_a_mod( atominfo_a, bondinfo);

           // connect new link atom to fragment A with a single bond
           bondinfo.PushBack
           (
             sdf::BondInfo
             (
               *a_index_itr,
               atom_v_a_mod.GetSize() - 1, // 0-index correction for last atom
               chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond
             )
           );

           // re-make this atom vector now that we have the extra bond
           atom_v_a_mod = chemistry::AtomVector< chemistry::AtomComplete>( atominfo_a, bondinfo);
           chemistry::FragmentComplete atom_v_a_mod_frag( atom_v_a_mod, "");

           // connect atoms from from fragment b to the link atom
           for
           (
               auto b_index_itr( LINK_INDICES_B.Begin()), b_index_itr_end( LINK_INDICES_B.End());
               b_index_itr != b_index_itr_end;
               ++b_index_itr
           )
           {

             // get atom vector
             chemistry::AtomVector< chemistry::AtomComplete> atom_v_b( FRAGMENT_B.GetAtomVector());

             // cannot connect hydrogen atoms to anything
             if( atom_v_b( *b_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
             {
               BCL_MessageStd( "Whoops! Specified atom index in B is a hydrogen atom! Ignoring index...");
               continue;
             }

             // open valences
             storage::Pair< chemistry::FragmentComplete, size_t> pair_a( OpenValence( atom_v_a_mod_frag, atom_v_a_mod_frag.GetSize() - 1));
             storage::Pair< chemistry::FragmentComplete, size_t> pair_b( OpenValence( FRAGMENT_B, *b_index_itr));

             // link fragment A to the link atom
             storage::Pair< bool, chemistry::FragmentComplete> new_fragment
             (
               chemistry::MergeFragmentComplete::MergeFragments
               (
                 pair_a.First(),
                 pair_b.First(),
                 chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
                 storage::Pair< size_t, size_t>( pair_a.Second(), pair_b.Second())
               )
             );

             // clean the new molecule
             util::ShPtr< chemistry::FragmentComplete> clean_frag( CallCleaner( new_fragment.Second(), FRAGMENT_A));

             // add new molecule to the ensemble
             if( new_fragment.First() && clean_frag.IsDefined())
             {
               linked_mols.PushBack( *clean_frag);
             }
           }
         }
         BCL_MessageCrt( "Done linking fragments with single element: " + util::Format()( *string_itr));
       }
       // return final molecules
       return linked_mols;
     }

     //! @brief link two fragments with an alkyl chain
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDICES_A link indices in first molecule
     //! @param LINK_INDICES_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     chemistry::FragmentEnsemble LinkFragments::AlkylLink
     (
       const chemistry::FragmentComplete &FRAGMENT_A,
       const chemistry::FragmentComplete &FRAGMENT_B,
       const storage::Vector< size_t> &LINK_INDICES_A,
       const storage::Vector< size_t> &LINK_INDICES_B,
       const size_t &REPEATS
     ) const
     {

       return chemistry::FragmentEnsemble();
     }

     //! @brief link two fragments with a methoxy repeat
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDICES_A link indices in first molecule
     //! @param LINK_INDICES_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     chemistry::FragmentEnsemble LinkFragments::MethoxyLink
     (
       const chemistry::FragmentComplete &FRAGMENT_A,
       const chemistry::FragmentComplete &FRAGMENT_B,
       const storage::Vector< size_t> &LINK_INDICES_A,
       const storage::Vector< size_t> &LINK_INDICES_B,
       const size_t &REPEATS
     ) const
     {

       return chemistry::FragmentEnsemble();
     }

     //! @brief link two fragments with an ethoxy repeat
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDICES_A link indices in first molecule
     //! @param LINK_INDICES_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     chemistry::FragmentEnsemble LinkFragments::EthoxyLink
     (
       const chemistry::FragmentComplete &FRAGMENT_A,
       const chemistry::FragmentComplete &FRAGMENT_B,
       const storage::Vector< size_t> &LINK_INDICES_A,
       const storage::Vector< size_t> &LINK_INDICES_B,
       const size_t &REPEATS
     ) const
     {

       return chemistry::FragmentEnsemble();
     }

     //! @brief link two fragments with an amide repeat
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDICES_A link indices in first molecule
     //! @param LINK_INDICES_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     chemistry::FragmentEnsemble LinkFragments::AmideLink
     (
       const chemistry::FragmentComplete &FRAGMENT_A,
       const chemistry::FragmentComplete &FRAGMENT_B,
       const storage::Vector< size_t> &LINK_INDICES_A,
       const storage::Vector< size_t> &LINK_INDICES_B,
       const size_t &REPEATS
     ) const
     {
       // initialize output ensemble
       chemistry::FragmentEnsemble linked_mols;

       BCL_MessageCrt( "Begin linking fragments with an amide repeat of size " + util::Format()( REPEATS));
       BCL_MessageStd( "A");

       // get atom vectors
       chemistry::AtomVector< chemistry::AtomComplete>
       atom_v_a( FRAGMENT_A.GetAtomVector()), atom_v_b( FRAGMENT_B.GetAtomVector());

       // connect atoms from fragment a
       for
       (
           auto a_index_itr( LINK_INDICES_A.Begin()), a_index_itr_end( LINK_INDICES_A.End());
           a_index_itr != a_index_itr_end;
           ++a_index_itr
       )
       {
         // cannot connect hydrogen atoms to anything
         if( atom_v_a( *a_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
         {
           BCL_MessageStd( "Whoops! Specified atom index in A is a hydrogen atom! Ignoring index...");
           continue;
         }

         // connect atoms from from fragment b
         for
         (
             auto b_index_itr( LINK_INDICES_B.Begin()), b_index_itr_end( LINK_INDICES_B.End());
             b_index_itr != b_index_itr_end;
             ++b_index_itr
         )
         {
           // cannot connect hydrogen atoms to anything
           if( atom_v_b( *b_index_itr).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
           {
             BCL_MessageStd( "Whoops! Specified atom index in B is a hydrogen atom! Ignoring index...");
             continue;
           }

           // open valences
           BCL_MessageStd( "B");
           storage::Pair< chemistry::FragmentComplete, size_t> pair_a( OpenValence( FRAGMENT_A, *a_index_itr));
           storage::Pair< chemistry::FragmentComplete, size_t> pair_b( OpenValence( FRAGMENT_B, *b_index_itr));

           // retrieve amide
           chemistry::FragmentEnsemble amide( GenerateAmideLinker( "N", REPEATS));

           // get indices for connection to molecule
           size_t n_index( util::GetUndefined< size_t>()), c_index( util::GetUndefined< size_t>());
           size_t n_undef_index( util::GetUndefined< size_t>()), c_undef_index( util::GetUndefined< size_t>());
           for
           (
               auto atom_itr( amide.GetMolecules().FirstElement().GetAtomVector().Begin()),
               atom_itr_end( amide.GetMolecules().FirstElement().GetAtomVector().End());
               atom_itr != atom_itr_end;
               ++atom_itr
           )
           {
             // get the indices for the attachment atoms (they are bonded to undefined atoms)
             for
             (
                 auto bond_itr( atom_itr->GetBonds().Begin()), bond_itr_end( atom_itr->GetBonds().End());
                 bond_itr != bond_itr_end;
                 ++bond_itr
             )
             {
               // go over bonds of each atom
               if( bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Undefined)
               {
                 if( atom_itr->GetElementType() == chemistry::GetElementTypes().e_Nitrogen)
                 {
                   n_index = amide.GetMolecules().FirstElement().GetAtomVector().GetAtomIndex( *atom_itr);
                   n_undef_index = amide.GetMolecules().FirstElement().GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom());
                 }
                 else if( atom_itr->GetElementType() == chemistry::GetElementTypes().e_Carbon)
                 {
                   c_index = amide.GetMolecules().FirstElement().GetAtomVector().GetAtomIndex( *atom_itr);
                   c_undef_index = amide.GetMolecules().FirstElement().GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom());
                 }
               } // end bond_itr
             } // end atom_itr
           }

           BCL_MessageStd( "C");
           // make two fragments for attachment to molecule A
           // one fragment will be attached via nitrogen
           // the other via carbon
           // remove the necessary undefined atoms
           storage::Vector< size_t> n_vec_indices, c_vec_indices;
           chemistry::AtomVector< chemistry::AtomComplete>
               n_vec( amide.GetMolecules().FirstElement().GetAtomVector()),
               c_vec( amide.GetMolecules().FirstElement().GetAtomVector());
           for( size_t n_defined( 0); n_defined < n_vec.GetSize(); ++n_defined)
           {
             if( n_vec( n_defined).GetElementType() != chemistry::GetElementTypes().e_Undefined)
             {
               n_vec_indices.PushBack( n_defined);
             }
           }
           for( size_t c_defined( 0); c_defined < c_vec.GetSize(); ++c_defined)
           {
             if( c_vec( c_defined).GetElementType() != chemistry::GetElementTypes().e_Undefined)
             {
               c_vec_indices.PushBack( c_defined);
             }
           }

           // make corrections to attachment indices
           if( n_undef_index < n_index)
           {
             n_index -= size_t( 1);
           }
           if( c_undef_index < c_index)
           {
             c_index -= size_t( 1);
           }

           // build new amide fragments without the undefined atom at the attachment site
           n_vec.Reorder( n_vec_indices);
           chemistry::FragmentComplete new_n_amide( n_vec, "");
           c_vec.Reorder( c_vec_indices);
           chemistry::FragmentComplete new_c_amide( c_vec, "");

           BCL_MessageStd( "D");
           // make an n-amide link to molecule A
           storage::Pair< bool, chemistry::FragmentComplete> n_fragment
           (
             chemistry::MergeFragmentComplete::MergeFragments
             (
               pair_a.First(),
               new_n_amide,
               chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
               storage::Pair< size_t, size_t>( pair_a.Second(), n_index)
             )
           );
           // make a c-amide link to molecule A
           storage::Pair< bool, chemistry::FragmentComplete> c_fragment
           (
             chemistry::MergeFragmentComplete::MergeFragments
             (
               pair_a.First(),
               new_c_amide,
               chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
               storage::Pair< size_t, size_t>( pair_a.Second(), c_index)
             )
           );

           BCL_MessageStd( "E");
           // now we want to find the index of the other attachment site for each new fragment
           // this will be the atom attached to the undefined atom in the amide linker
           size_t n_frag_connection_index( util::GetUndefined< size_t>()), c_frag_connection_index( util::GetUndefined< size_t>());
           size_t n_frag_undefined_index( util::GetUndefined< size_t>()), c_frag_undefined_index( util::GetUndefined< size_t>());
           BCL_Assert( n_fragment.First(), "Failed N");
           if( n_fragment.First())
           {
             // each atom in fragment
             for
             (
                 auto atom_itr( n_fragment.Second().GetAtomVector().Begin()), atom_itr_end( n_fragment.Second().GetAtomVector().End());
                 atom_itr != atom_itr_end;
                 ++atom_itr
             )
             {
               // go over each bond
               for
               (
                   auto bond_itr( atom_itr->GetBonds().Begin()), bond_itr_end( atom_itr->GetBonds().End());
                   bond_itr != bond_itr_end;
                   ++bond_itr
               )
               {
                 if( bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Undefined)
                 {
                   n_frag_connection_index = n_fragment.Second().GetAtomVector().GetAtomIndex( *atom_itr);
                   n_frag_undefined_index = n_fragment.Second().GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom());
                   BCL_Debug( n_frag_connection_index);
                   BCL_Debug( n_frag_undefined_index);
                   break;
                 }
               }
             }
           }
           BCL_Assert( c_fragment.First(), "Failed C");
           if( c_fragment.First())
           {
             // each atom in fragment
             for
             (
                 auto atom_itr( c_fragment.Second().GetAtomVector().Begin()), atom_itr_end( c_fragment.Second().GetAtomVector().End());
                 atom_itr != atom_itr_end;
                 ++atom_itr
             )
             {
               // go over each bond
               for
               (
                   auto bond_itr( atom_itr->GetBonds().Begin()), bond_itr_end( atom_itr->GetBonds().End());
                   bond_itr != bond_itr_end;
                   ++bond_itr
               )
               {
                 if( bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Undefined)
                 {
                   c_frag_connection_index = c_fragment.Second().GetAtomVector().GetAtomIndex( *atom_itr);
                   c_frag_undefined_index = c_fragment.Second().GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom());
                   break;
                 }
               }
             }
           }

           BCL_MessageStd( "F");
           // remove the necessary undefined atoms from each fragments
           n_vec_indices.Reset();
           c_vec_indices.Reset();
           chemistry::AtomVector< chemistry::AtomComplete> n_frag_vec( n_fragment.Second().GetAtomVector()), c_frag_vec( c_fragment.Second().GetAtomVector());
           for( size_t n_defined( 0); n_defined < n_frag_vec.GetSize(); ++n_defined)
           {
             if( n_frag_vec( n_defined).GetElementType() != chemistry::GetElementTypes().e_Undefined)
             {
               n_vec_indices.PushBack( n_defined);
             }
           }
           for( size_t c_defined( 0); c_defined < c_frag_vec.GetSize(); ++c_defined)
           {
             if( c_frag_vec( c_defined).GetElementType() != chemistry::GetElementTypes().e_Undefined)
             {
               c_vec_indices.PushBack( c_defined);
             }
           }

           // make corrections to attachment indices
           if( n_frag_undefined_index < n_frag_connection_index)
           {
             n_frag_connection_index -= size_t( 1);
           }
           if( c_frag_undefined_index < c_frag_connection_index)
           {
             c_frag_connection_index -= size_t( 1);
           }

           // build new amide fragments without the undefined atom at the attachment site
           n_frag_vec.Reorder( n_vec_indices);
           chemistry::FragmentComplete new_n_frag( n_frag_vec, "");
           c_frag_vec.Reorder( c_vec_indices);
           chemistry::FragmentComplete new_c_frag( c_frag_vec, "");

           BCL_MessageStd( "G");
           // link fragments
           BCL_Debug( n_frag_connection_index);
           BCL_Debug( pair_b.Second());
           storage::Pair< bool, chemistry::FragmentComplete> n_fully_connected
           (
             chemistry::MergeFragmentComplete::MergeFragments
             (
               new_n_frag,
               pair_b.First(),
               chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
               storage::Pair< size_t, size_t>( n_frag_connection_index, pair_b.Second())
             )
           );
           BCL_Debug( c_frag_connection_index);
           BCL_Debug( pair_b.Second());
           storage::Pair< bool, chemistry::FragmentComplete> c_fully_connected
           (
             chemistry::MergeFragmentComplete::MergeFragments
             (
               new_c_frag,
               pair_b.First(),
               chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
               storage::Pair< size_t, size_t>( c_frag_connection_index, pair_b.Second())
             )
           );

           // clean the new molecules
           BCL_MessageStd( "H");
           util::ShPtr< chemistry::FragmentComplete> n_clean_frag( CallCleaner( n_fully_connected.Second(), FRAGMENT_A));
           util::ShPtr< chemistry::FragmentComplete> c_clean_frag( CallCleaner( c_fully_connected.Second(), FRAGMENT_A));

           // add new molecules to the ensemble
           if( n_fully_connected.First() && n_clean_frag.IsDefined())
           {
             linked_mols.PushBack( *n_clean_frag);
           }
           if( c_fully_connected.First() && c_clean_frag.IsDefined())
           {
             linked_mols.PushBack( *c_clean_frag);
           }
         }
       }

       // return final molecules
       BCL_MessageCrt( "End linking fragments with an amide repeat of size " + util::Format()( REPEATS));
       return linked_mols;
     }

     //! @brief link two fragments with an amide repeat
     //! @param CONNECTION whether to connect the amide to molecule A via C or N or both
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     chemistry::FragmentEnsemble LinkFragments::GenerateAmideLinker
     (
       const std::string &CONNECTION,
       const size_t &REPEATS
     ) const
     {
       // Begin
       chemistry::FragmentEnsemble amide;

       // Read in amide linker file
       io::IFStream file;
       io::File::MustOpenIFStream
       (
         file,
         chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "combichem_linkers/amide.sdf.gz"
       );
       amide.ReadMoreFromMdl( file, sdf::e_Maintain);

       io::File::CloseClearFStream( file);

       // amide linker repeat length of 1
       if( REPEATS == size_t( 1))
       {
         return amide;
       }
       else
       {
         return amide;
       }
       return chemistry::FragmentEnsemble();
     }

     //! @brief remove a hydrogen atom from a target atom
     //! @param FRAGMENT the molecule of interest
     //! @param ATOM_INDEX the index of the atom in the molecule of interest
     //! @return the new molecule and index of the desired atom
     storage::Pair< chemistry::FragmentComplete, size_t> LinkFragments::OpenValence
     (
       const chemistry::FragmentComplete &FRAGMENT,
       const size_t &ATOM_INDEX
     ) const
     {
       // find a hydrogen atom attached to specified atom index
       size_t h_index( util::GetUndefined< size_t>()), new_atom_index( ATOM_INDEX);
       for
       (
           auto bond_itr( FRAGMENT.GetAtomVector()( ATOM_INDEX).GetBonds().Begin()),
           bond_itr_end( FRAGMENT.GetAtomVector()( ATOM_INDEX).GetBonds().End());
           bond_itr != bond_itr_end;
           ++bond_itr
       )
       {
         if( bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
         {
           h_index = FRAGMENT.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom());
           break;
         }
       }

       // if not hydrogen atoms then return input
       if( !util::IsDefined( h_index))
       {
         return storage::Pair< chemistry::FragmentComplete, size_t>( FRAGMENT, ATOM_INDEX);
       }

       // adjust the reference atom index as needed
       if( h_index < ATOM_INDEX)
       {
         new_atom_index -= size_t( 1);
       }

       // generate new atom indices excluding the hydrogen atom
       storage::Vector< size_t> keep_indices;
       for
       (
           size_t i( 0); i < FRAGMENT.GetSize(); ++i
       )
       {
         if( i == h_index)
         {
           continue;
         }
         keep_indices.PushBack( i);
       }

       // make new molecule without hydrogen atom
       chemistry::AtomVector< chemistry::AtomComplete> new_atom_v( FRAGMENT.GetAtomVector());
       new_atom_v.Reorder( keep_indices);
       chemistry::FragmentComplete new_frag( new_atom_v, "");

       // return new molecule with the updated atom index
       return storage::Pair< chemistry::FragmentComplete, size_t>( new_frag, new_atom_index);
     }

     //! @brief wrapper that calls the FragmentMapConformer::Clean function to generate a legitimate 3D conformer and fix bond lengths
     //! @param FRAGMENT the molecule that needs to be fixed
     //! @param REFERENCE scaffold molecule for substructure alignment reference
     //! @return pointer to cleaned molecule
     util::ShPtr< chemistry::FragmentComplete> LinkFragments::CallCleaner
     (
       const chemistry::FragmentComplete &FRAGMENT,
       const chemistry::FragmentComplete &REFERENCE
     ) const
     {
       // hold on to this
       static chemistry::HydrogensHandler hdyrogens_handler;
       static chemistry::FragmentMapConformer cleaner;

       // remove hydrogen atoms
       chemistry::AtomVector< chemistry::AtomComplete> new_frag_v( FRAGMENT.GetAtomVector());
       hdyrogens_handler.Remove( new_frag_v);

       // call cleaner
       util::ShPtr< chemistry::FragmentComplete> clean_frag
       (
         cleaner.Clean
         (
           new_frag_v, REFERENCE, m_DrugLikenessTypeFlag->GetFirstParameter()->GetValue()
         )
       );
       return clean_frag;
     }

     //! @brief wrapper that calls ResolveClashes and OptimizePose with a combined global/local conf ensemble for each mol
     //! @param ENSEMBLE the molecule that needs to be fixed
     //! @param SCORER the property to be used to score the interface
     //! @return pointer to cleaned molecule
     chemistry::FragmentEnsemble LinkFragments::LigandLocalDock
     (
       const chemistry::FragmentEnsemble &ENSEMBLE,
       const descriptor::CheminfoProperty &SCORER
     ) const
     {
       // so we can add the properties
       chemistry::FragmentEnsemble ensemble( ENSEMBLE), best_confs;

       // we only want to make once
       static chemistry::RotamerLibraryFile rotamer_library_file;
       static chemistry::SampleConformations sample_confs
       (
         rotamer_library_file,
         "",
         0.0,   // tolerance
         100,    // number of conformations
         2000,   // number of iterations
         false,   // change chirality
         0.0,    // random dihedral change weight
         true,   // generate 3d
         0.5     // clash tolerance
       );
       static chemistry::FragmentMapConformer conf_mapper
       (
         m_DrugLikenessTypeFlag->GetFirstParameter()->GetValue(),
         m_MDLString,
         m_PocketFilename,
         SCORER,
         true,
         storage::Vector< float>(),
         false
       );
       static storage::Vector< float> bfactors( conf_mapper.GetBFactors());

       // Sample conformers of each molecule
       for
       (
           auto mol_itr( ensemble.Begin()), mol_itr_end( ensemble.End());
           mol_itr != mol_itr_end;
           ++mol_itr
       )
       {
         // store MDL for interface scoring and generate confs
         mol_itr->GetStoredPropertiesNonConst().SetMDLProperty( conf_mapper.GetMDL(), conf_mapper.GetPocketFilename());
         chemistry::FragmentEnsemble confs( sample_confs( *mol_itr).First());

         // optimize pose
         util::ShPtr< chemistry::FragmentComplete> opti_pose;
         if( confs.GetSize())
         {
           BCL_MessageStd( "Getting best scoring non-clashing conformer!");
           static chemistry::FragmentStochasticPoseOptimizer pose_optimizer
           (
             SCORER,
             bfactors,
             conf_mapper.GetMDL(),
             conf_mapper.GetPocketFilename(),
             size_t( 20),
             size_t( 100),
             size_t( 100),
             float( 1.0),
             opti::e_LargerIsBetter,
             float( 5.0)
           );
           opti_pose = pose_optimizer.StochasticPoseOptimization( confs);
         }
         else
         {
           opti_pose = util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete( *mol_itr));
           continue;
         }

         // save best pose
         if( opti_pose->GetSize())
         {
           best_confs.PushBack( *opti_pose);
         }
       }
       return best_confs;
     }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LinkFragments::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &LinkFragments::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
