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

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "app/bcl_app_groups.h"
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_rotamer_library_interface.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_subgraph.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_pair.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ExtractFocusedLibraryDesignFragments
    //! @brief Extracts fragments added to a core scaffold during FocusedLibraryDesign...actually its a bit more general than that
    //!
    //! @author geanesar
    //! @date Sept 25 2014
    //!
    //! TODO merge into AlignToScaffold as functionality has overlap
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class ExtractFocusedLibraryDesignFragments :
      public Interface
    {

    public:

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

      util::ShPtr< command::FlagInterface> m_EnsembleFilenameFlag;

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      ExtractFocusedLibraryDesignFragments();

      ExtractFocusedLibraryDesignFragments *Clone() const
      {
        return new ExtractFocusedLibraryDesignFragments( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_InputFilenamesFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_EnsembleFilenameFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      int Main() const
      {

        // Read in the molecules
        chemistry::FragmentFeed scaffold_feed( m_InputFilenamesFlag->GetStringList(), sdf::e_Remove);

        chemistry::FragmentEnsemble scaffolds;

        // Read in scaffolds
        for( ; scaffold_feed.NotAtEnd(); ++scaffold_feed)
        {
          scaffolds.PushBack( *scaffold_feed);
        }

        size_t n_scaffs( scaffolds.GetSize());
        storage::Vector< graph::ConstGraph< size_t, size_t> > scaffold_graphs( n_scaffs);

        chemistry::ConformationGraphConverter graph_maker
        (
          chemistry::ConformationGraphConverter::e_AtomType,
          chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic
        );

        size_t pos( 0);
        for
        (
          chemistry::FragmentEnsemble::const_iterator itr_scaff( scaffolds.Begin()), itr_scaff_end( scaffolds.End());
          itr_scaff != itr_scaff_end;
          ++itr_scaff, ++pos
        )
        {
          scaffold_graphs( pos) = graph_maker( *itr_scaff);
        }

        chemistry::FragmentEnsemble fragments;

        chemistry::FragmentFeed mol_feed( m_EnsembleFilenameFlag->GetStringList(), sdf::e_Saturate);
        graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
        // Check each molecule for instances of the scaffolds
        io::OFStream out1;
        io::File::MustOpenOFStream( out1, "matched_mols.sdf");
        io::File::CloseClearFStream( out1);
        for( ; mol_feed.NotAtEnd(); ++mol_feed)
        {
          // Make graphs
          const graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( *mol_feed));
          const chemistry::ConformationGraphConverter::t_AtomGraph mol_atom_graph( graph_maker.CreateGraphWithAtoms( *mol_feed));
          util::OwnPtr< const graph::ConstGraph< size_t, size_t> > mol_graph_ptr( &mol_graph, false);

          // Checking for scaffolds requires an isomorphism.  Set the supergraph to the molecule we are checking
          isomorphism.SetGraph( mol_graph);

          // Check the molecule for occurrences of each scaffold
          for( size_t scaff_no( 0); scaff_no < n_scaffs; ++scaff_no)
          {
            isomorphism.SetSubgraph( scaffold_graphs( scaff_no));

            // If the scaffold was found, extract its fragments
            if( isomorphism.FindIsomorphism())
            {
              io::File::MustOpenOFStream( out1, "matched_mols.sdf", std::ios::app);
              mol_feed->WriteMDL( out1);
              io::File::CloseClearFStream( out1);
              // Get the isomorphism
              storage::Vector< size_t> mol_iso( isomorphism.GetIsomorphism());

              // Make a subgraph from that isomorphism
              graph::Subgraph< size_t, size_t> scaff_subgraph( mol_graph_ptr, mol_iso);

              // Copy the graph
              graph::ConstGraph< size_t, size_t> mol_graph_mod( mol_graph);

              // Remove edges between the scaffold and the fragments
              storage::List< storage::Pair< size_t, size_t> > edges( scaff_subgraph.GetAdjacentEdgeIndices());
              for
              (
                storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_edge( edges.Begin()), itr_edge_end( edges.End());
                itr_edge != itr_edge_end;
                ++itr_edge
              )
              {
                mol_graph_mod.RemoveEdge( itr_edge->First(), itr_edge->Second());
              }

              // Sort the isomorphism so that we can determine if one of the fragments is the scaffold
              mol_iso.Sort( std::less< size_t>());
              storage::List< storage::Vector< size_t> > frags( graph::Connectivity::GetComponents( mol_graph_mod));
              for
              (
                storage::List< storage::Vector< size_t> >::iterator itr_frag( frags.Begin()), itr_frag_end( frags.End());
                itr_frag != itr_frag_end;
                ++itr_frag
              )
              {
                itr_frag->Sort( std::less< size_t>());
                if( *itr_frag == mol_iso)
                {
                  continue;
                }
                chemistry::AtomVector< chemistry::AtomComplete> frag_atom_vector
                (
                  graph_maker.CreateAtomsFromGraph( mol_atom_graph.GetSubgraph( *itr_frag))
                );
                fragments.PushBack( chemistry::FragmentComplete( frag_atom_vector, ""));
              }
            }
          }
        }

        // Write out the discovered fragments
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
        fragments.WriteMDL( output);
        io::File::CloseClearFStream( output);

        return 0;
      } // Main()

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

    private:

      static const ApplicationType ExtractFocusedLibraryDesignFragments_Instance;

    }; // ExtractFocusedLibraryDesignFragments

      //! @brief standard constructor
    ExtractFocusedLibraryDesignFragments::ExtractFocusedLibraryDesignFragments() :
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output", "flag selecting the output file name",
          command::Parameter
          (
            "filename", "filename for output sdf"
          )
        )
      ),
      m_InputFilenamesFlag
      (
        new command::FlagDynamic
        (
          "scaffold",
          "filename of an sdf file containing the scaffolds to use",
          command::Parameter
          (
            "file containing the scaffold",
            "name of files containing derived structures",
            command::ParameterCheckFileExistence()
          ),
          1,
          21
        )
      ),
      m_EnsembleFilenameFlag
      (
        new command::FlagDynamic
        (
          "ensemble",
          "filenames for the ensemble of molecules to fragment",
          command::Parameter
          (
            "filenames of derived structure sdfs",
            "name of files containing derived structures",
            command::ParameterCheckFileExistence()
          ),
          1,
          21
        )
      )
    {
    }

    const ApplicationType ExtractFocusedLibraryDesignFragments::ExtractFocusedLibraryDesignFragments_Instance
    (
      GetAppGroups().AddAppToGroup( new ExtractFocusedLibraryDesignFragments(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
