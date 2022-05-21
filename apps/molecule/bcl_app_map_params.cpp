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

// app header
#include "app/bcl_app_apps.h"

// include header for this application
#include "molecule/bcl_app_map_params.h"

// include headers from the bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_types.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_configurational_bond_types.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_element_types.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_default.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "sdf/bcl_sdf_bond_info.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_string_functions.h"

namespace bcl
{
  namespace app
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    MapParams::MapParams() :
      m_RosettaSDF
      (
        new command::FlagStatic
        (
          "rosetta_sdf",
          "input rosetta sdf containing rosetta params atom names in MDL properties",
          command::Parameter
          (
            "",
            ""
          )
        )
      ),
      m_RosettaAtomNamesMDL
      (
        new command::FlagStatic
        (
          "rosetta_mdl",
          "MDL property name for atom names in 'rosetta_sdf'",
          command::Parameter
          (
            "rosetta_mdl",
            "",
            command::ParameterCheckDefault(),
            "Atom Names"
          )
        )
      ),
      m_AmberToolsSDF
      (
        new command::FlagStatic
        (
          "ambertools_sdf",
          "input ambertools sdf corresponding to the 'ambertools_prepi'",
          command::Parameter
          (
            "ambertools_sdf",
            ""
          )
        )
      ),
      m_AmberToolsPrepi
      (
        new command::FlagStatic
        (
          "ambertools_prepi",
          "input ambertools prepi corresponding to the 'ambertools_sdf'",
          command::Parameter
          (
            "ambertools_prepi",
            ""
          )
        )
      ),
      m_AmberToolsMC
      (
        new command::FlagStatic
        (
          "ambertools_mc",
          "input ambertools mc if one was used to generate the prepi file",
          command::Parameter
          (
            "ambertools_mc",
            "",
            command::ParameterCheckDefault(),
            ""
          )
        )
      ),
      m_AtomComparisonType
      (
        new command::FlagStatic
        (
          "atom_type",
          "atom type comparison; ",
          command::Parameter
          (
            "atom_comparison_type",
            "string value of the atom type info for substructure matching calculation",
            command::ParameterCheckSerializable
            (
              chemistry::ConformationGraphConverter::AtomComparisonTypeEnum()
            ),
            "ElementType"
          )
        )
      ),
      m_BondComparisonType
      (
        new command::FlagStatic
        (
          "bond_type",
          "bond type comparison; ",
          command::Parameter
          (
            "bond_comparison_type",
            "string value of the bond type info for substructure matching calculation",
            command::ParameterCheckSerializable
            (
              chemistry::ConfigurationalBondTypeData::DataEnum()
            ),
            "Identity"
          )
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix",
          "prefix for output files",
          command::Parameter
          (
            "output_prefix",
            "",
            "AmberToRosettaPrepi"
          )
        )
      )
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MapParams::GetDescription() const
    {
      return "MapParams is a helper application to reduce manual labor in NCAA design. "
          "It maps atom names from a Rosetta-derived SDF to an AmberTools prepi file.";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MapParams::GetReadMe() const
    {
      static std::string s_read_me =
          "MapParams takes an input SDF from Rosetta (must contain MDL property mapping "
          "atom indices to atom names in Rosetta params file), an input SDF from AmberTools, "
          "and an input PREPI file from AmberTools corresponding to the AmberTools SDF. MapParams "
          "returns a new PREPI file with the original AmberTools atom names replaced by the "
          "Rosetta params atom names.";
      return s_read_me;
    }

  //////////////////////
  //   operations     //
  //////////////////////

    //! @brief initializes the command object for this application
    //! @return a ShPtr to a Command containing all of this applications parameters
    util::ShPtr< command::Command> MapParams::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // primary flags
      sp_cmd->AddFlag( m_RosettaSDF);
      sp_cmd->AddFlag( m_RosettaAtomNamesMDL);
      sp_cmd->AddFlag( m_AmberToolsSDF);
      sp_cmd->AddFlag( m_AmberToolsPrepi);
      sp_cmd->AddFlag( m_AmberToolsMC);
      sp_cmd->AddFlag( m_AtomComparisonType);
      sp_cmd->AddFlag( m_BondComparisonType);
      sp_cmd->AddFlag( m_OutputPrefixFlag);

      // default flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

    //! @brief reads in the AmberTools SDF
    void MapParams::ReadRosettaSDF() const
    {
      // load the Rosetta SDF
      m_Mutex.Lock();
      io::IFStream input_sdf;
      io::File::MustOpenIFStream( input_sdf, m_RosettaSDF->GetFirstParameter()->GetValue());
      m_RosettaMol.ReadMoreFromMdl( input_sdf);
      io::File::CloseClearFStream( input_sdf);
      m_Mutex.Unlock();
    }

    //! @brief reads in the Rosetta SDF
    void MapParams::ReadAmberToolsSDF() const
    {
      // load the AmberTools SDF
      m_Mutex.Lock();
      io::IFStream input_sdf;
      io::File::MustOpenIFStream( input_sdf, m_AmberToolsSDF->GetFirstParameter()->GetValue());
      m_AmberToolsMol.ReadMoreFromMdl( input_sdf);
      io::File::CloseClearFStream( input_sdf);
      m_Mutex.Unlock();
    }

    void MapParams::ReadAmberToolsAtomNames() const
    {
      // load the AmberTools PREPI
      m_Mutex.Lock();
      io::IFStream input_prepi;
      io::File::MustOpenIFStream( input_prepi, m_AmberToolsPrepi->GetFirstParameter()->GetValue());

      // read line by line
      size_t line_index( 0);
      for( std::string line; std::getline( input_prepi, line); ++line_index)
      {
        // skip the first 11 lines
        if( line_index < 10)
        {
          continue;
        }

        // split the atom names
        storage::Vector< std::string> split_line( util::SplitString( line, " "));
        if( !split_line.GetSize())
        {
          break;
        }
        std::string atom_name( util::Format().W( 2).L()( split_line( 1)));
        m_AmberToolsPrepiAtomNames.PushBack( atom_name);
      }
      io::File::CloseClearFStream( input_prepi);
      m_Mutex.Unlock();
    }

    void MapParams::ReadAmberToolsOmitAtomNames() const
    {
      // load the AmberTools PREPI
      m_Mutex.Lock();
      io::IFStream input_mc;
      io::File::MustOpenIFStream( input_mc, m_AmberToolsMC->GetFirstParameter()->GetValue());

      // read line by line
      size_t line_index( 0);
      for( std::string line; std::getline( input_mc, line); ++line_index)
      {
        // skip lines that are not for omit atom names
        if( std::string( line.substr( 0, 4)).compare( "OMIT") != 0)
        {
          continue;
        }

        // split the atom names
        storage::Vector< std::string> split_line( util::SplitString( line, " "));
        if( !split_line.GetSize())
        {
          break;
        }
        std::string atom_name( util::Format().W( 2).L()( split_line( 1)));
        m_AmberToolsOmitAtomNames.PushBack( atom_name);
      }
      io::File::CloseClearFStream( input_mc);
      m_Mutex.Unlock();
    }

    //! @brief gets rosetta atom names from SDF
    void MapParams::GetRosettaSDFAtomNames() const
    {
      // assumes we have an MDL property on the SDF
      chemistry::FragmentComplete &mol( m_RosettaMol.GetMolecules().FirstElement());
      std::string atom_names
      (
        mol.GetMDLProperty
        (
          m_RosettaAtomNamesMDL->GetFirstParameter()->GetValue()
        )
      );

      // split the atom names
      storage::Vector< std::string> split_atom_names( util::SplitString( atom_names, " (),"));

      // remove the indices
      m_RosettaSDFAtomNames.Reset();
      for( size_t i( 1), sz( split_atom_names.GetSize()); i < sz; ++++i)
      {
        m_RosettaSDFAtomNames.PushBack( util::Format().W( 2).L()( split_atom_names( i)));
      }
    }

    //! @brief gets ambertools atom names from SDF
    void MapParams::GetAmberToolsSDFAtomNames() const
    {
      // running map of elements to number of times element appears in SDF
      storage::Map< chemistry::ElementType, size_t> element_counts;
      storage::Vector< storage::Pair< chemistry::AtomType, size_t>> element_names;

      // loop over all atoms and get counts
      chemistry::FragmentComplete &mol( m_AmberToolsMol.GetMolecules().FirstElement());
      for
      (
          auto atom_itr( mol.GetAtomVector().Begin()), atom_itr_end( mol.GetAtomVector().End());
          atom_itr != atom_itr_end;
          ++atom_itr
      )
      {
        // insert the element into the map
        if( !element_counts.Has( atom_itr->GetElementType()))
        {
          element_counts.InsertElement
          (
            std::make_pair
            (
              atom_itr->GetElementType(),
              1
            )
          );
        }
        else
        {
          element_counts.Find( atom_itr->GetElementType())->second += 1;
        }

        // add the element name
        element_names.PushBack
        (
          storage::Pair< chemistry::AtomType, size_t>
          (
            std::make_pair
            (
              atom_itr->GetAtomType(),
              element_counts.Find( atom_itr->GetElementType())->second
            )
          )
        );
      }

      // write final atom names
      m_AmberToolsSDFAtomNames.Reset();
      for
      (
          auto itr( element_names.Begin()), itr_end( element_names.End());
          itr != itr_end;
          ++itr
      )
      {
        // get the two letter code
        std::string name( util::Format().W( 2).L()( itr->First()->GetTwoLetterCode()));

        // get element name; strip second character if BCNS or second character is X
        if
        (
            itr->First()->GetElementType() == chemistry::GetElementTypes().e_Boron ||
            itr->First()->GetElementType() == chemistry::GetElementTypes().e_Carbon ||
            itr->First()->GetElementType() == chemistry::GetElementTypes().e_Nitrogen ||
            itr->First()->GetElementType() == chemistry::GetElementTypes().e_Sulfur ||
            std::string( 1, name[ 1]).compare( "X") == 0
        )
        {
          name.pop_back();
        }
        m_AmberToolsSDFAtomNames.PushBack( name + util::Format().W( 2).L()( itr->Second()));
      }

      // if we pass an ambertools mc file used to help generate the prepi, exclude omit atoms
      if( m_AmberToolsOmitAtomNames.GetSize())
      {
        // loop over our SDF atom names
        size_t i( 0);
        storage::Vector< std::string> pruned_names( m_AmberToolsSDFAtomNames);
        auto tracker( pruned_names.Begin());
        for
        (
            auto sdf_atoms_itr( m_AmberToolsSDFAtomNames.Begin()),
            sdf_atoms_itr_end( m_AmberToolsSDFAtomNames.End());
            sdf_atoms_itr != sdf_atoms_itr_end;
            ++sdf_atoms_itr, ++i, ++tracker
        )
        {
          // check to see if that atom name is in the mc file list
          if( m_AmberToolsOmitAtomNames.Find( *sdf_atoms_itr) < m_AmberToolsOmitAtomNames.GetSize())
          {
            // remove it from the tracker
            pruned_names.RemoveElement( tracker);
          }
        }

        // now the final atom names are from the pruned list
        m_AmberToolsSDFAtomNames = pruned_names;
      }
    }

    //! @brief maps the rosetta molecule to the ambertools molecule
    //! @param ROSETTA_MOL the input rosetta molecule
    //! @param AMBERTOOLS_MOL the input ambertools molecule
    //! @return the common subgraph isomorphism between the two
    graph::CommonSubgraphIsomorphism< size_t, size_t> MapParams::MapRosettaToAmberToolsSDF
    (
      const chemistry::FragmentComplete &ROSETTA_MOL,
      const chemistry::FragmentComplete &AMBERTOOLS_MOL
    ) const
    {
      // initialize iso with solution type
      graph::CommonSubgraphIsomorphism< size_t, size_t> csi
      (
        graph::CommonSubgraphIsomorphismBase::e_GreedyUnconnected
      );

      // generate graphs of each molecule
      chemistry::ConformationGraphConverter graph_converter
      (
        chemistry::ConformationGraphConverter::AtomComparisonTypeEnum( m_AtomComparisonType->GetFirstParameter()->GetValue()),
        chemistry::ConfigurationalBondTypeData::DataEnum( m_BondComparisonType->GetFirstParameter()->GetValue()),
        false
      );
      graph::ConstGraph< size_t, size_t> rosetta_graph( graph_converter( ROSETTA_MOL));
      graph::ConstGraph< size_t, size_t> ambertools_graph( graph_converter( AMBERTOOLS_MOL));

      // set the graphs for the iso search
      csi.SetGraphA
      (
        util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &rosetta_graph, false)
      );
      csi.SetGraphB
      (
        util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &ambertools_graph, false)
      );

      // find the iso between the two graphs
      csi.FindIsomorphism( csi.EstimateUpperBounds());
      return csi;
    }

    //! @brief map the rosetta atom names to a new prepi file
    storage::Vector< std::string> MapParams::MapRosettaToPrepi() const
    {
      // find isomorphism between rosetta and ambertools molecules
      graph::CommonSubgraphIsomorphism< size_t, size_t> isomorphism
      (
        MapRosettaToAmberToolsSDF
        (
          m_RosettaMol.GetMolecules().FirstElement(),
          m_AmberToolsMol.GetMolecules().FirstElement()
        )
      );

      // debug
      for
      (
          auto iso_itr( isomorphism.GetIsomorphism().Begin()), iso_itr_end( isomorphism.GetIsomorphism().End());
          iso_itr != iso_itr_end;
          ++iso_itr
      )
      {
        std::string atom_name( m_AmberToolsPrepiAtomNames( iso_itr->second));
        std::string rosetta_name( m_RosettaSDFAtomNames( iso_itr->first));
        BCL_MessageStd( "Mapping " + util::Format().W( 2).L()( atom_name) + " --> " + util::Format().W( 2).L()( rosetta_name));
      }

      // load the AmberTools PREPI
      m_Mutex.Lock();
      io::IFStream input_prepi;
      io::File::MustOpenIFStream( input_prepi, m_AmberToolsPrepi->GetFirstParameter()->GetValue());

      // for each rosetta atom, find the corresponding atom in ambertools
      size_t line_index( 0);
      std::string line, original_line, rosetta_prepi;
      for( ; std::getline( input_prepi, line); ++line_index)
      {
        original_line = line;
        bool found( false);
        for
        (
            auto iso_itr( isomorphism.GetIsomorphism().Begin()), iso_itr_end( isomorphism.GetIsomorphism().End());
            iso_itr != iso_itr_end;
            ++iso_itr
        )
        {
          std::string atom_name( m_AmberToolsPrepiAtomNames( iso_itr->second));
          std::string rosetta_name( m_RosettaSDFAtomNames( iso_itr->first));

          // find the ambertools prepi atom name to replace
          std::size_t pos( line.find( atom_name));

          // replace with the rosetta atom name
          if( pos != std::string::npos)
          {
            rosetta_name.length() > atom_name.length() ?
                line.replace( pos, rosetta_name.length(), rosetta_name) :
                line.replace( pos, atom_name.length(), rosetta_name);
            found = true;
          }
        }
        // if we never replace the line just use the original
        if( !found)
        {
          std::cout << original_line << std::endl;
          rosetta_prepi.append( line);
          rosetta_prepi.append( "\n");
        }
        else
        {
          std::cout << line << std::endl;
          rosetta_prepi.append( line);
          rosetta_prepi.append( "\n");
        }
      }
      io::File::CloseClearFStream( input_prepi);
      m_Mutex.Unlock();

      // write output
      m_Mutex.Lock();
      io::OFStream output_prepi;
      io::File::MustOpenOFStream( output_prepi, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".rosetta.prepi");
      output_prepi << rosetta_prepi;
      io::File::CloseClearFStream( output_prepi);
      m_Mutex.Unlock();

      return storage::Vector< std::string>();
    }

    //! @brief gets the indices of the AmberTools molecule that
    //! correspond to the atom names in the AmberTools PREPI
    void MapParams::MapAmberToolsSDFPrepi() const
    {
      // track which names we have matched by removing from this vector
      storage::Vector< std::string> leftover_names( m_AmberToolsSDFAtomNames);

      // loop over our SDF atom names
      size_t i( 0);
      auto tracker( leftover_names.Begin());
      for
      (
          auto sdf_atoms_itr( m_AmberToolsSDFAtomNames.Begin()),
          sdf_atoms_itr_end( m_AmberToolsSDFAtomNames.End());
          sdf_atoms_itr != sdf_atoms_itr_end;
          ++sdf_atoms_itr, ++i, ++tracker
      )
      {
        // check to see if that atom name is in the prepi file
        if( m_AmberToolsPrepiAtomNames.Find( *sdf_atoms_itr) < m_AmberToolsPrepiAtomNames.GetSize())
        {
          // remove it from the tracker
          leftover_names.RemoveElement( tracker);

          // add value to map
          m_AmberToolsMap.Insert
          (
            std::make_pair
            (
              i,
              m_AmberToolsPrepiAtomNames( m_AmberToolsPrepiAtomNames.Find( *sdf_atoms_itr))
            )
          );
        }
      }
    }

  //////////////////////
  //    operators     //
  //////////////////////

    //! @brief the Main function
    //! @return 0 for success
    int MapParams::Main() const
    {
      // get all member data set correctly
      ReadRosettaSDF();
      if( m_AmberToolsMC->GetFlag())
      {
        ReadAmberToolsOmitAtomNames();
      }
      ReadAmberToolsSDF();
      ReadAmberToolsAtomNames();

      // compute iso between rosetta and ambertools molecules
      GetRosettaSDFAtomNames();
      GetAmberToolsSDFAtomNames();
      MapAmberToolsSDFPrepi();

      // now that we have all of the pieces, generate a new prepi file with Rosetta atom names
      MapRosettaToPrepi();

      return 0;
    } // Main

  //////////////////////
  // helper functions //
  //////////////////////

    // Construct the static instance of this application, and add it to the ChemInfo group
    const ApplicationType MapParams::MapParams_Instance
    (
      GetAppGroups().AddAppToGroup( new MapParams(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
