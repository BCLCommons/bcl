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
    //! @class MakeFocusedLibraryDesignScaffold
    //! @brief Application for identifying core scaffold from substituted molecules
    //!
    //! @author geanesar
    //! @date Sept 25 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class MakeFocusedLibraryDesignScaffold :
      public Interface
    {

    public:

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

      util::ShPtr< command::FlagInterface> m_DatasetFilenameFlag;

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      MakeFocusedLibraryDesignScaffold();

      MakeFocusedLibraryDesignScaffold *Clone() const
      {
        return new MakeFocusedLibraryDesignScaffold( *this);
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
        sp_cmd->AddFlag( m_DatasetFilenameFlag);

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
        storage::Vector< graph::ConstGraph< size_t, size_t> > scaffold_graphs;
        storage::Vector< storage::Vector< size_t> > scaffold_saturate_atoms;

        chemistry::ConformationGraphConverter graph_maker
        (
          chemistry::ConformationGraphConverter::e_AtomType,
          chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic
        );

        for( ; scaffold_feed.NotAtEnd(); ++scaffold_feed)
        {
          const chemistry::FragmentComplete &mol( *scaffold_feed);
          scaffolds.PushBack( mol);
          scaffold_graphs.PushBack( graph_maker( mol));
          scaffold_saturate_atoms.PushBack( storage::Vector< size_t>( mol.GetNumberAtoms(), size_t( 1)));
        }

        chemistry::FragmentFeed mol_feed( m_DatasetFilenameFlag->GetStringList(), sdf::e_Remove);
        for( ; mol_feed.NotAtEnd(); ++mol_feed)
        {
          const chemistry::FragmentComplete &mol( *mol_feed);
          graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( mol));
          util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_graph_ptr( &mol_graph, false);

          // See if it matches any scaffold
          chemistry::FragmentEnsemble::const_iterator itr_scaff( scaffolds.Begin());
          for( size_t scaff_no( 0), num_scaffs( scaffold_graphs.GetSize()); scaff_no < num_scaffs; ++scaff_no, ++itr_scaff)
          {
            storage::Vector< size_t> &scaff_atoms( scaffold_saturate_atoms( scaff_no));
            graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
            isomorphism.SetGraph( mol_graph);
            isomorphism.SetSubgraph( scaffold_graphs( scaff_no));
            if( !isomorphism.FindIsomorphism())
            {
              continue;
            }
            storage::Vector< size_t> iso( isomorphism.GetIsomorphism());
            graph::Subgraph< size_t, size_t> subgraph( mol_graph_ptr, iso);
            storage::List< storage::Pair< size_t, size_t> > adj_edges( subgraph.GetAdjacentEdgeIndices());
            for
            (
              storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_edge( adj_edges.Begin()), itr_edge_end( adj_edges.End());
              itr_edge != itr_edge_end;
              ++itr_edge
            )
            {
              scaff_atoms( iso.Find( itr_edge->First())) = 0;
            }
          }
        }

        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
        size_t scaff_no( 0);
        for
        (
          chemistry::FragmentEnsemble::const_iterator itr_scaff( scaffolds.Begin()), itr_scaff_end( scaffolds.End());
          itr_scaff != itr_scaff_end;
          ++itr_scaff, ++scaff_no
        )
        {
          chemistry::AtomVector< chemistry::AtomComplete> atom_vector( itr_scaff->GetAtomVector());
          chemistry::HydrogensHandler::SaturatePartial( atom_vector, scaffold_saturate_atoms( scaff_no));
          chemistry::FragmentComplete new_frag( atom_vector, "");
          new_frag.WriteMDL( output);
        }
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

      static const ApplicationType MakeFocusedLibraryDesignScaffold_Instance;

    }; // MakeFocusedLibraryDesignScaffold

      //! @brief standard constructor
    MakeFocusedLibraryDesignScaffold::MakeFocusedLibraryDesignScaffold() :
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
      m_DatasetFilenameFlag
      (
        new command::FlagDynamic
        (
          "dataset",
          "filenames for dataset molecules",
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

    const ApplicationType MakeFocusedLibraryDesignScaffold::MakeFocusedLibraryDesignScaffold_Instance
    (
      GetAppGroups().AddAppToGroup( new MakeFocusedLibraryDesignScaffold(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
