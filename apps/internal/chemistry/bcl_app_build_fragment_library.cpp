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
#include "bcl_app_build_fragment_library.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_molecule.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_template_instantiations.h"

namespace bcl
{
  namespace app
  {

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    BuildFragmentLibrary::BuildFragmentLibrary() :
      m_OutputFileFlag
      (
        new command::FlagStatic
        (
          "output",
          "filename to output all fragments",
          command::Parameter
          (
            "output",
            "filename to output sdf of all component fragments"
          )
        )
      ),
      m_MaxNumberofBonds
      (
        new command::FlagStatic
        (
          "max_bonds",
          "maximum number of non-ring bonds (which will broken) that a molecule can have for it to be fragmented",
          command::Parameter
          (
            "max_bonds",
            "maximum number of non-ring bonds (which will broken) that a molecule can have for it to be fragmented",
            command::ParameterCheckRanged< size_t>(),
            "100"
          )
        )
      ),
      m_MaxRotatableBonds
      (
        new command::FlagStatic
        (
          "max_rot",
          "maximum number of rotatable bonds that a fragment should have so that it is written to output",
          command::Parameter
          (
            "max_rot",
            "maximum number of rotatable bonds that a fragment should have so that it is written to output",
            command::ParameterCheckRanged< size_t>(),
            "11"
          )
        )
      ),
      m_FindConstitutions
      (
        new command::FlagStatic
        (
          "constitutions",
          "flag to find all constitutions rather than configurations"
        )
      ),
      m_IgnoreUncommonRings
      (
        new command::FlagStatic
        (
          "ignore_uncommon_rings",
          "If set, ignore all but the 242 most common ring systems (which are present in at least 50 structures of the COD database)"
        )
      )
    {
    }

    //! @brief copy constructor; only copy flags for applications
    BuildFragmentLibrary::BuildFragmentLibrary( const BuildFragmentLibrary &PARENT) :
        m_OutputFileFlag( PARENT.m_OutputFileFlag),
        m_MaxNumberofBonds( PARENT.m_MaxNumberofBonds),
        m_MaxRotatableBonds( PARENT.m_MaxRotatableBonds),
        m_FindConstitutions( PARENT.m_FindConstitutions),
        m_IgnoreUncommonRings( PARENT.m_IgnoreUncommonRings)
    {
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> BuildFragmentLibrary::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // ensembles containing the molecules to be edited
      chemistry::FragmentFeed::AddFlags( *sp_cmd, "molecules");

      // insert all the flags and params
      sp_cmd->AddFlag( m_OutputFileFlag);
      sp_cmd->AddFlag( m_MaxNumberofBonds);
      sp_cmd->AddFlag( m_MaxRotatableBonds);
      sp_cmd->AddFlag( m_FindConstitutions);
      sp_cmd->AddFlag( m_IgnoreUncommonRings);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    BuildFragmentLibrary *BuildFragmentLibrary::Clone() const
    {
      return new BuildFragmentLibrary( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &BuildFragmentLibrary::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string BuildFragmentLibrary::GetDescription() const
    {
      return "BuildFragmentLibrary does what it suggests";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &BuildFragmentLibrary::GetReadMe() const
    {
      static std::string s_read_me
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::BuildFragmentLibrary, terms of use, "
        "appropriate citation, installation procedures, BCL::BuildFragmentLibrary execution, "
        "technical support, and future research directions.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::BuildFragmentLibrary?\n"
        "BCL::BuildFragmentLibrary is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::BuildFragmentLibrary is a utility that fragments "
        "molecules into fragments by breaking all the non-ring bonds."
        "\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::BuildFragmentLibrary.\n"
        "\n"
        ""
        ""
        ""
        ""
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::BuildFragmentLibrary.\n"
        "Running BCL::BuildFragmentLibrary requires an sdf file containing molecules that need to be fragmented.\n"
        "\n"
        "2) Run BCL::BuildFragmentLibrary to fragment molecules\n"
        "\n"
        "At a command prompt, navigate to the location of your BCL::BuildFragmentLibrary executable program. The syntax for"
        "running the application looks like the following"
        "\n"
        "bcl.exe BuildFragmentLibrary -molecules_filename <filename.sdf> -output <filename>"
        "\n\nFLAGS:\n\n"
        "-molecules_filenames <filename> -> file containing ensemble of molecules that need to be fragmented\n"
        "-output <filename> -> file to which fragments will be written out\n"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe BuildFragmentLibrary -help\n"
        "\n"
        "For more general information about the product, type bcl.exe BuildFragmentLibrary -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::BuildFragmentLibrary.\n"
        "BCL::BuildFragmentLibrary is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );
      return s_read_me;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &BuildFragmentLibrary::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::BuildFragmentLibrary is a utility for fragmenting molecules"
        "Features of BCL::BuildFragmentLibrary\n"
        "<ul>"
        "  <li>Fragment molecules by breaking non ring bonds as\n"
        "  </li>"
        "  <li>Compressed molecule files (bz2, gzip) are supported</li>"
        "</ul>\n\n"
        "!build_fragment_library.png!"
      );

      return s_web_text;
    }

    namespace
    {
      //! @brief remove bonds that are simple terminal bonds that hang off a ring. The positions of such atoms is trivial to
      //!        reconstruct
      //! @param MOL molecule to remove the unininteresting terminal bonds from rings from
      //! @return molecule without any terminal bonds off of rings
      chemistry::FragmentComplete RemoveTerminalBondsOfRings( const chemistry::FragmentComplete &FRAGMENT)
      {
        typename chemistry::ConformationGraphConverter::t_AtomGraph atom_graph
        (
          chemistry::ConformationGraphConverter::CreateGraphWithAtoms( FRAGMENT)
        );
        storage::Vector< size_t> indices_to_keep;
        indices_to_keep.AllocateMemory( FRAGMENT.GetSize());
        size_t i( 0);
        for( auto itr( FRAGMENT.GetAtomsIterator()); itr.NotAtEnd(); ++itr, ++i)
        {
          if
          (
            itr->GetBonds().GetSize() != size_t( 1)
          )
          {
            indices_to_keep.PushBack( i);
          }
          else if
          (
            itr->GetBonds().Begin()->GetTargetAtom().CountNonValenceBondsWithProperty
            (
              chemistry::ConfigurationalBondTypeData::e_IsInRing,
              size_t( 1)
            ) >= size_t( 2)
          )
          {
            if
            (
              itr->GetBonds().Begin()->GetTargetAtom().CountNonValenceBondsWithProperty
              (
                chemistry::ConfigurationalBondTypeData::e_IsAromatic,
                size_t( 1)
              ) < size_t( 2)
              && itr->GetAtomType() == chemistry::GetAtomTypes().H_S
            )
            {
              // bonds off of non-aromatic rings. Only keep the hydrogens for exact substitution pattern matching, which
              // can influence the most common rotamer of the ring
              indices_to_keep.PushBack( i);
            }
          }
          else if
          (
            itr->GetAtomType() != chemistry::GetAtomTypes().H_S
            || itr->GetBonds().Begin()->GetTargetAtom().GetBonds().GetSize()
               > itr->GetBonds().Begin()->GetTargetAtom().GetNumberCovalentlyBoundHydrogens() + size_t( 1)
          )
          {
            // skip terminal hydrogens to reduce computational complexity
            indices_to_keep.PushBack( i);
          }
        }
        return indices_to_keep.GetSize() == FRAGMENT.GetSize()
               ? FRAGMENT
               : chemistry::FragmentComplete( chemistry::ConformationGraphConverter::CreateAtomsFromGraph( atom_graph.GetSubgraph( indices_to_keep)), FRAGMENT.GetName());
      }

      //! @get a list of all the rings that have been seen at least 100 times
      util::ShPtr< chemistry::ConstitutionSet> GetCommonRingGraphs()
      {
        io::IFStream input;
        const std::string filename( chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "rings_mincnt50.sdf.gz");
        io::File::MustOpenIFStream( input, filename);
        chemistry:: FragmentEnsemble common_rings( input, sdf::e_Remove);
        io::File::CloseClearFStream( input);
        storage::Vector< graph::ConstGraph< size_t, size_t> > common_ring_graphs;
        common_ring_graphs.AllocateMemory( common_rings.GetSize());
        util::ShPtr< chemistry::ConstitutionSet> set_sp( new chemistry::ConstitutionSet);
        for( auto itr( common_rings.Begin()), itr_end( common_rings.End()); itr != itr_end; ++itr)
        {
          set_sp->Insert( chemistry::FragmentConstitutionShared( *itr));
        }
        return set_sp;
      }

      //! @brief returns list of ring or chain fragments
      //! @param MOLECULE molecule of interest
      //! @param MOLECULE_GRAPH graph of molecule of interest
      //! @return list of ring or chain fragments
      chemistry::FragmentEnsemble RemoveUncommonRings( const chemistry::ConformationInterface &MOLECULE)
      {
        static const util::ShPtr< chemistry::ConstitutionSet> s_common_rings( GetCommonRingGraphs());

        typename chemistry::ConformationGraphConverter::t_AtomGraph atom_graph
        (
          chemistry::ConformationGraphConverter::CreateGraphWithAtoms( MOLECULE)
        );

        typename chemistry::ConformationGraphConverter::t_AtomGraph atom_graph_copy_for_rings( atom_graph);

        // declare a list for storing atom indices of rings/chain
        chemistry::FragmentSplitRings fsr( true);
        storage::List< storage::Vector< size_t> > ring_indices( fsr.GetComponentVertices( MOLECULE, atom_graph_copy_for_rings));
        chemistry::FragmentEnsemble all_ring_systems( fsr.ConvertComponentsIntoEnsemble( MOLECULE, ring_indices, atom_graph_copy_for_rings));
        storage::Set< size_t> atoms_in_uncommon_rings;
        auto itr_macrocycle_indices( ring_indices.Begin());
        for
        (
          auto itr_all_rings( all_ring_systems.Begin()), itr_all_rings_end( all_ring_systems.End());
          itr_all_rings != itr_all_rings_end;
          ++itr_all_rings, ++itr_macrocycle_indices
        )
        {
          if( s_common_rings->Find( chemistry::FragmentConstitutionShared( *itr_all_rings)) == s_common_rings->End())
          {
            // not a common ring, remove edges between the atoms from the graph
            for
            (
              auto itr_a( itr_macrocycle_indices->Begin()), itr_a_end( itr_macrocycle_indices->End());
              itr_a != itr_a_end;
              ++itr_a
            )
            {
              for
              (
                auto itr_b( itr_macrocycle_indices->Begin()), itr_b_end( itr_macrocycle_indices->End());
                itr_b != itr_a;
                ++itr_b
              )
              {
                atom_graph.RemoveEdge( *itr_a, *itr_b);
              }
            }
          }
        }
        return fsr.ConvertComponentsIntoEnsemble( MOLECULE, graph::Connectivity::GetComponents( atom_graph), atom_graph);
      }
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int BuildFragmentLibrary::Main() const
    {
      // get a feed to get molecules that need to be fragmented
      chemistry::FragmentFeed feed( "molecules", sdf::e_Saturate);

      // set high tolerance so that only one conformation is obtained for each configuration
      util::ShPtr< chemistry::ConfigurationSet> configuration_set( new chemistry::ConfigurationSet());
      util::ShPtr< chemistry::ConstitutionSet> constitution_set( new chemistry::ConstitutionSet());

      // get parameters
      const size_t max_bonds( m_MaxNumberofBonds->GetFirstParameter()->GetNumericalValue< size_t>());
      const size_t max_rot( m_MaxRotatableBonds->GetFirstParameter()->GetNumericalValue< size_t>());

      // create object that fragments molecules
      chemistry::FragmentMolecule fragmenter;
      if( m_FindConstitutions->GetFlag())
      {
        fragmenter = chemistry::FragmentMolecule( constitution_set, max_rot);
      }
      else
      {
        fragmenter = chemistry::FragmentMolecule( configuration_set, max_rot);
      }

      math::RunningAverage< double> ave_seconds;
      const util::Time start_time( util::Time::GetCurrent());
      util::Time last_mol_start_time( util::Time::GetCurrent());
      for( size_t mol_number( 0); feed.NotAtEnd(); ++feed, ++mol_number)
      {
        // remove terminal bonds off rings
        chemistry::FragmentComplete newmol( RemoveTerminalBondsOfRings( *feed));
        chemistry::FragmentEnsemble ensemble
        (
          m_IgnoreUncommonRings->GetFlag()
          ? RemoveUncommonRings( newmol)
          : chemistry::FragmentEnsemble( storage::List< chemistry::FragmentComplete>( size_t( 1), newmol))
        );
        const size_t ens_size( ensemble.GetSize());
        size_t frag_number( 0);
        for( auto itr_ens( ensemble.Begin()), itr_ens_end( ensemble.End()); itr_ens != itr_ens_end; ++itr_ens, ++frag_number)
        {
          // get number of breakable bonds and if it is more than the specified maximum number of bonds, then skip the
          // molecule
          size_t breakable_bonds
          (
            itr_ens->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, 0)
          );

          if( itr_ens->GetNumberAtoms() > 15 && breakable_bonds > max_bonds)
          {
            util::GetLogger().LogStatus
            (
              " molecule # " + util::Format()( mol_number) + " with "
              + util::Format()( breakable_bonds) + " breakable bonds, skipping."
            );
            continue;
          }
          util::Time end_time( util::Time::GetCurrent());
          const double last_mol_time( ( end_time - last_mol_start_time).GetSecondsFractional());
          util::GetLogger().LogStatus
          (
            " molecule # " + util::Format()( mol_number)
            + " fragment # " + util::Format()( frag_number) + "/" + util::Format()( ens_size)
            + " with "
            + util::Format()( breakable_bonds) + " breakable bonds. Last mol took "
            + util::Format()( last_mol_time) + " seconds. Ave seconds per molecule: "
            + util::Format()( ( end_time - start_time).GetSecondsFractional() / ( mol_number + 0.5))
          );
          last_mol_start_time = end_time;

          // call fragment molecule function
          fragmenter( *itr_ens);
        }
      }

      // open the output file for fragments that are generated by fragmenting molecules
      io::OFStream output;
      io::File::MustOpenOFStream( output, m_OutputFileFlag->GetFirstParameter()->GetValue());

      // maximum number of chain dihedral angles for fragments to be written out
      static const size_t s_max_dihedrals( 4);

      sdf::MdlHandler::GetAddAtomMdlLineFlag()->SetFlag();
      if( !m_FindConstitutions->GetFlag())
      {
        for
        (
          util::ShPtrList< chemistry::FragmentConfigurationShared>::const_iterator
            itr_config( fragmenter.GetConfigurations()->Begin()),
            itr_config_end( fragmenter.GetConfigurations()->End());
          itr_config != itr_config_end;
          ++itr_config
        )
        {
          chemistry::PriorityDihedralAngles pdi;
          chemistry::FragmentComplete comp( **itr_config);
          auto edges( pdi.GetDihedralEdges( comp));
          size_t n_chain( 0);
          for( auto itr_edges( edges.Begin()), itr_edges_end( edges.End()); itr_edges != itr_edges_end; ++itr_edges)
          {
            if( !itr_edges->GetEdgeData()->IsBondInRing() && ++n_chain > s_max_dihedrals)
            {
              break;
            }
          }
          if( n_chain > s_max_dihedrals)
          {
            continue;
          }
          // store counts for fragments
          ( *itr_config)->WriteMDL( output);
        }
      }
      else
      {
        for
        (
          util::ShPtrList< chemistry::FragmentConstitutionShared>::const_iterator
            itr_config( fragmenter.GetConstitutions()->Begin()),
            itr_config_end( fragmenter.GetConstitutions()->End());
          itr_config != itr_config_end;
          ++itr_config
        )
        {
          chemistry::PriorityDihedralAngles pdi;
          chemistry::FragmentComplete comp( **itr_config);
          auto edges( pdi.GetDihedralEdges( comp));
          size_t n_chain( 0);
          for( auto itr_edges( edges.Begin()), itr_edges_end( edges.End()); itr_edges != itr_edges_end; ++itr_edges)
          {
            if( !itr_edges->GetEdgeData()->IsBondInRing() && ++n_chain > s_max_dihedrals)
            {
              break;
            }
          }
          if( n_chain > s_max_dihedrals)
          {
            continue;
          }
          // store counts for fragments
          ( *itr_config)->WriteMDL( output);
        }
      }

      // now that application has run close the output stream
      io::File::CloseClearFStream( output);

      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    const ApplicationType BuildFragmentLibrary::BuildFragmentLibrary_Instance
    (
      GetAppGroups().AddAppToGroup( new BuildFragmentLibrary(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl

