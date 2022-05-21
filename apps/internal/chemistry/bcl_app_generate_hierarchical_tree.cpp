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
#include "bcl_app_generate_hierarchical_tree.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_probability_score.h"
#include "chemistry/bcl_chemistry_rotamer_library_interface.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_isomorphism.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"

namespace bcl
{
  namespace app
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    GenerateHierarchicalTree::GenerateHierarchicalTree() :
      m_LibraryFlag
      (
        new command::FlagStatic
        (
          "library",
          "file that contains the library of molecules for which one desires hierarchical tree. Ignored if -prune is given",
          command::Parameter
          (
            "library",
            "file that contains the library of molecules for which one desires hierarchical tree",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      ),
      m_FormatFlag
      (
        new command::FlagStatic
        (
          "format",
          "format in which the hierarchical tree is desired",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< chemistry::RotamerLibraryInterface>()),
            "File(prefix=rotlib)"
          )
        )
      ),
      m_PruneFlag
      (
        new command::FlagStatic
        (
          "prune",
          "option to prune an existing library",
          command::Parameter
          (
            "existing library to prune",
            "",
            command::ParameterCheckSerializable( util::Implementation< chemistry::RotamerLibraryInterface>()),
            "File(prefix=rotlib)"
          )
        )
      )
    {
    }

      //! copy constructor
      GenerateHierarchicalTree::GenerateHierarchicalTree( const GenerateHierarchicalTree &APP) :
        m_LibraryFlag( APP.m_LibraryFlag),
        m_FormatFlag( APP.m_FormatFlag),
        m_PruneFlag( APP.m_PruneFlag)
      {
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> GenerateHierarchicalTree::InitializeCommand() const
      {
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // hydrogen preferences
        sdf::AddMoleculeIOPrefFlags( *sp_cmd);

        // ensemble containing the molecules to be aligned to scaffold
        sp_cmd->AddFlag( m_LibraryFlag);
        // ensemble that will be written out
        sp_cmd->AddFlag( m_FormatFlag);
        // ensemble that will be written out
        sp_cmd->AddFlag( m_PruneFlag);

        // add default bcl parameters
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      GenerateHierarchicalTree *GenerateHierarchicalTree::Clone() const
      {
        return new GenerateHierarchicalTree( *this);
      }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GenerateHierarchicalTree::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string GenerateHierarchicalTree::GetDescription() const
    {
      return "CreateRotamerLibraryTreeFormat creates a hierarchical tree for an ensemble of molecules where parent node"
          "is immediate substructure of child node";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &GenerateHierarchicalTree::GetReadMe() const
    {
      static std::string s_read_me
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::CreateRotamerLibraryTreeFormat, terms of use, "
        "appropriate citation, installation procedures, BCL::CreateRotamerLibraryTreeFormat execution, "
        "technical support, and future research directions.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::CreateRotamerLibraryTreeFormat?\n"
        "BCL::CreateRotamerLibraryTreeFormat is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::CreateRotamerLibraryTreeFormat is a utility that "
        "hierarchical tree for an ensemble of molecules where parent node is immediate substructure of child node."
        "\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::CreateRotamerLibraryTreeFormat.\n"
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
        "VI. RUNNING BCL::CreateRotamerLibraryTreeFormat.\n"
        "Running BCL::CreateRotamerLibraryTreeFormat requires an sdf file containing the ensemble of molecules for which hierarchical tree is desired.\n"
        "\n"
        "2) Run BCL::CreateRotamerLibraryTreeFormat generate hierarchical tree format for an ensemble of molecules\n"
        "\n"
        "At a command prompt, navigate to the location of your BCL::CreateRotamerLibraryTreeFormat executable program. The syntax for"
        "running the application looks like the following"
        "\n"
        "bcl.exe CreateRotamerLibraryTreeFormat -rotamer_library <filename> -rotamer_library_format <file or db>"
        "\n\nFLAGS:\n\n"
        "-library <filename> -> file containing ensemble of molecules\n"
        "-format <file or db> -> specify flat file format or database as output. See help for more information \n"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe CreateRotamerLibraryTreeFormat -help\n"
        "\n"
        "For more general information about the product, type bcl.exe CreateRotamerLibraryTreeFormat -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::CreateRotamerLibraryTreeFormat.\n"
        "BCL::CreateRotamerLibraryTreeFormat is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );
      return s_read_me;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &GenerateHierarchicalTree::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::CreateRotamerLibraryTreeFormat is a utility for generating hierarchical tree of molecules "
        "Features of BCL::CreateRotamerLibraryTreeFormat\n"
        "<ul>"
        "  <li>Parent nodes are immediate substructures of child nodes.\n"
        "  </li>"
        "  <li>Output in flat file or database.\n"
        "  </li>"
        "  <li>Compressed molecule files (bz2, gzip) are supported</li>"
        "</ul>\n\n"
        "!create_tree_format.png!"
      );

      return s_web_text;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int GenerateHierarchicalTree::Main() const
    {
      //read in scaffold sdf and ensemble sdf
      if( !m_PruneFlag->GetFlag())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_LibraryFlag->GetFirstParameter()->GetValue());
        chemistry::FragmentEnsemble rotamer_library( input, sdf::e_Maintain);
        io::File::CloseClearFStream( input);
//        io::OFStream output;
//        io::File::MustOpenOFStream( output, "pruned_library.sdf.gz");
//        size_t mol_index( 0);
//        size_t n_kept( 0), n_removed( 0), n_removed_scoring( 0);
//        for
//        (
//          auto itr_rot_lib( rotamer_library.GetMolecules().Begin()), itr_rot_lib_end( rotamer_library.GetMolecules().End());
//          itr_rot_lib != itr_rot_lib_end;
//          ++mol_index
//        )
//        {
//          chemistry::PriorityDihedralAngles pdi;
//          auto edges( pdi.GetDihedralEdges( *itr_rot_lib));
//          size_t n_ring( 0), n_chain( 0);
//          for( auto itr_edges( edges.Begin()), itr_edges_end( edges.End()); itr_edges != itr_edges_end; ++itr_edges)
//          {
//            if( itr_edges->GetEdgeData()->IsBondInRing())
//            {
//              ++n_ring;
//            }
//            else if( ++n_chain > 4)
//            {
//              break;
//            }
//          }
//          if( n_chain > size_t( 4))
//          {
//            ++n_removed;
//            util::GetLogger().LogStatus
//            (
//              "Would delete: " + util::Format()( mol_index)
//              + " " + util::Format()( n_removed) + " / " + util::Format()( n_removed + n_kept) + " due to more than 4 dihedral chains"
//            );
//            auto itr_rot_lib_prev( itr_rot_lib);
//            ++itr_rot_lib;
//            rotamer_library.GetMolecules().Remove( itr_rot_lib_prev);
//            continue;
//          }
//
//          chemistry::ConformationGraphConverter graph_maker
//          (
//            chemistry::ConformationGraphConverter::e_AtomType,
//            chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
//          );
//          graph::SubgraphIsomorphism< size_t, size_t> si;
//          auto graph_ab( graph_maker( *itr_rot_lib));
//          si.SetGraphExternalOwnership( graph_ab);
//          si.SetSubgraphExternalOwnership( graph_ab);
//          si.FindIsomorphism();
//          util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism> smfi
//          (
//            size_t( 1),
//            util::ShPtr< chemistry::SmallMoleculeFragmentIsomorphism>
//            (
//              new chemistry::SmallMoleculeFragmentIsomorphism( *itr_rot_lib, *itr_rot_lib, si)
//            )
//          );
//          auto vec_data
//          (
//            chemistry::SmallMoleculeFragmentMapping().MapFragmentIsomorphisms
//            (
//              *itr_rot_lib, smfi, storage::Set< size_t>()
//            )
//          );
//          if( vec_data.IsEmpty())
//          {
//            auto itr_rot_lib_prev( itr_rot_lib);
//            ++itr_rot_lib;
//            ++n_removed;
//            rotamer_library.GetMolecules().Remove( itr_rot_lib_prev);
//          }
//          else
//          {
//            ++n_kept;
//            itr_rot_lib->WriteMDL( output);
//            ++itr_rot_lib;
//          }
//        }
//        io::File::CloseClearFStream( output);

        m_Format = m_FormatFlag->GetFirstParameter()->GetValue();
        m_Format->Create( rotamer_library);
      }
      else
      {
        util::Implementation< chemistry::RotamerLibraryInterface> original_lib( m_PruneFlag->GetFirstParameter()->GetValue());
        auto configuration_mapping( original_lib->RetrieveConfigurationMapping());
        size_t constitution_index( 0);
        chemistry::FragmentEnsemble ensemble_final;
        size_t n_kept( 0), n_removed( 0), n_removed_scoring( 0);
        io::OFStream output;
        io::File::MustOpenOFStream( output, "pruned_library.sdf.gz");
        for
        (
          auto itr_confm( configuration_mapping.Begin()), itr_confm_end( configuration_mapping.End());
          itr_confm != itr_confm_end;
          ++itr_confm, ++constitution_index
        )
        {
          auto ensembl( original_lib->RetrieveAssociatedConfigurations( *itr_confm));
          auto itr_set_index( itr_confm->Begin());
          for( auto itr_ensem( ensembl.Begin()), itr_ensem_end( ensembl.End()); itr_ensem != itr_ensem_end; ++itr_ensem, ++itr_set_index)
          {
            chemistry::PriorityDihedralAngles pdi;
            auto edges( pdi.GetDihedralEdges( *itr_ensem));
            size_t n_ring( 0), n_chain( 0);
            for( auto itr_edges( edges.Begin()), itr_edges_end( edges.End()); itr_edges != itr_edges_end; ++itr_edges)
            {
              if( itr_edges->GetEdgeData()->IsBondInRing())
              {
                ++n_ring;
              }
              else if( ++n_chain > 4)
              {
                break;
              }
            }
            if( n_chain > size_t( 4))
            {
              ++n_removed;
              util::GetLogger().LogStatus
              (
                "Would delete: " + util::Format()( *itr_set_index)
                + " " + util::Format()( n_removed) + " / " + util::Format()( n_removed + n_kept)
                + " due to more than 4 dihedral chains"
              );
              continue;
            }

            chemistry::ConformationGraphConverter graph_maker
            (
              chemistry::ConformationGraphConverter::e_AtomType,
              chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
            );
            graph::SubgraphIsomorphism< size_t, size_t> si;
            auto graph_ab( graph_maker( *itr_ensem));
            si.SetGraphExternalOwnership( graph_ab);
            si.SetSubgraphExternalOwnership( graph_ab);
            si.FindIsomorphism();
            util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism> smfi
            (
              size_t( 1),
              util::ShPtr< chemistry::SmallMoleculeFragmentIsomorphism>
              (
                new chemistry::SmallMoleculeFragmentIsomorphism( *itr_ensem, *itr_ensem, si, false)
              )
            );
            auto vec_data
            (
              chemistry::SmallMoleculeFragmentMapping().MapFragmentIsomorphisms
              (
                *itr_ensem, smfi, storage::Set< size_t>()
              )
            );
            if
            (
              vec_data.IsEmpty()
            )
            {
              ++n_removed;
              util::GetLogger().LogStatus
              (
                "Would delete: " + util::Format()( *itr_set_index)
                + " " + util::Format()( n_removed) + " / " + util::Format()( n_removed + n_kept)
              );
            }
            else
            {
              ++n_kept;
              util::GetLogger().LogStatus
              (
                "Keeping: " + util::Format()( *itr_set_index)
                + " " + util::Format()( n_removed) + " / " + util::Format()( n_removed + n_kept) + " "
                + util::Format()( n_kept) + " / " + util::Format()( n_removed + n_kept) + " retained"
              );
              if( vec_data( 0)->GetFragment().GetNumberRingBonds() && vec_data( 0)->GetFragment().GetNumberDihedralChainBonds())
              {
                for( size_t rot_n( 1), n_rot( vec_data( 0)->GetFragment().GetRotamerNumbers()); rot_n <= n_rot; ++rot_n)
                {
                  itr_ensem->GetStoredPropertiesNonConst().RemoveProperty( "Rotamer" + util::Format()( rot_n) + "Coordinates");
                }
              }
              itr_ensem->SetName( "");
              itr_ensem->RemoveProperty( "ContainsRings");
              itr_ensem->RemoveProperty( "BinSize");
              ensemble_final.PushBack( *itr_ensem);
              itr_ensem->WriteMDL(output);
            }
          }
        }

        io::File::CloseClearFStream( output);

        m_Format = m_FormatFlag->GetFirstParameter()->GetValue();
        m_Format->Create( ensemble_final);
      }

      // end
      return 0;
    }

    const ApplicationType GenerateHierarchicalTree::GenerateHierarchicalTree_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateHierarchicalTree(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
