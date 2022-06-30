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
#include "chemistry/bcl_chemistry_atom_environment.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "bcl_app_generate_atom_environment_hashmap.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_file.h"
#include "iterate/bcl_iterate_generic.h"
#include "math/bcl_math_linear_least_squares.h"
#include "model/bcl_model_retrieve_dataset_subset.h"
#include "sdf/bcl_sdf_bond_info.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    const ApplicationType GenerateAtomEnvironmentHashmap::GenerateAtomEnvironmentHashmap_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateAtomEnvironmentHashmap(), GetAppGroups().e_InternalChem)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    GenerateAtomEnvironmentHashmap::GenerateAtomEnvironmentHashmap() :
      m_OutputFilenameBase
      (
        new command::FlagStatic
        (
          "output",
          "base name for output of histograms",
          command::Parameter( "output", "base name for output of histograms")
        )
      ),
      m_WriteHashMap
      (
        new command::FlagStatic
        (
          "write_atom_environment_map",
          "count the number of times each atom environment occurs"
        )
      ),
      m_BondRadius
      (
        new command::FlagStatic
        (
          "bond_radius",
          "bond radius for atom environment counts; valid range 2 - 4",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckRanged< size_t>( 2, 4),
            "2"
          )
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new GenerateAtomEnvironmentHashmap
    GenerateAtomEnvironmentHashmap *GenerateAtomEnvironmentHashmap::Clone() const
    {
      return new GenerateAtomEnvironmentHashmap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GenerateAtomEnvironmentHashmap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> GenerateAtomEnvironmentHashmap::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // add flags for input
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // Output filename base
      sp_cmd->AddFlag( m_OutputFilenameBase);

      // whether to collect element-element-bond statistics
      sp_cmd->AddFlag( m_WriteHashMap);

      // get bond radius
      sp_cmd->AddFlag( m_BondRadius);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and read-me are useful
      command::GetAppDefaultFlags().AddRequiredCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &GenerateAtomEnvironmentHashmap::GetReadMe() const
    {
      static std::string s_read_me =
        "GenerateAtomEnvironmentHashmap generates atom environments from input molecular ensembles "
        "and stores them in a hashmap for reference later. Practical applications include restricting"
        "designed molecular structures to the fragments in the map and the ability to build local atom"
        "environment-based statistical potentials for drug design and discovery.";
      return s_read_me;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int GenerateAtomEnvironmentHashmap::Main() const
    {
      io::OFStream output;

      // iterate over the fragments from whatever feeds were given
      chemistry::FragmentEnsemble molecules_from_feed;
      for( chemistry::FragmentFeed feed; feed.NotAtEnd(); ++feed)
      {
        // read in molecules
        chemistry::FragmentComplete frag( *feed);

        // Skip molecules with non-gasteiger atom types
        if( frag.HasNonGasteigerAtomTypes())
        {
          BCL_MessageStd
          (
            "Molecule #" + util::Format()( feed.GetPosition())
            + " has atoms that cannot be described with gasteiger atom types; skipping"
          );
          continue;
        }
        molecules_from_feed.PushBack( frag);
      }

      // Collect stats
      if( m_WriteHashMap->GetFlag())
      {
        CountAtomEnvironments( molecules_from_feed);
      }

      // end
      return 0;
    }

    //! @brief generate counts of element-element bonds for statistics
    //! @param FRAG the fragment for which to calculate the counts
    void GenerateAtomEnvironmentHashmap::CountAtomEnvironments( const chemistry::FragmentEnsemble &ENSEMBLE) const
    {
      // create the atom environment object that will convert our atom environments to strings
      chemistry::AtomEnvironment generate_atom_environment;

      // reset hashmap before we start
      m_HashMap = storage::Vector< storage::HashMap< std::string, size_t> >( size_t( 3));

      // loop over all molecules in ensemble
      for
      (
          auto mol_itr( ENSEMBLE.Begin()), mol_itr_end( ENSEMBLE.End());
          mol_itr != mol_itr_end;
          ++mol_itr
      )
      {
        // iterate over all atoms in molecule
        for
        (
            auto atom_itr( mol_itr->GetAtomVector().Begin()), atom_itr_end( mol_itr->GetAtomVector().End());
            atom_itr != atom_itr_end;
            ++atom_itr
        )
        {
          // convert atom environment to string
          std::string atom_env_str_two, atom_env_str_three, atom_env_str_four;
          if( m_BondRadius->GetFirstParameter()->GetNumericalValue< size_t>() == size_t( 2))
          {
              atom_env_str_two = std::string( generate_atom_environment.MakeAtomEnvironmentStringTwo( *atom_itr));
          }
          if( m_BondRadius->GetFirstParameter()->GetNumericalValue< size_t>() == size_t( 3))
          {
              atom_env_str_three = std::string( generate_atom_environment.MakeAtomEnvironmentStringThree( *atom_itr));
          }
          if( m_BondRadius->GetFirstParameter()->GetNumericalValue< size_t>() == size_t( 4))
          {
              atom_env_str_four = std::string( generate_atom_environment.MakeAtomEnvironmentStringFour( *atom_itr));
          }

          // add it to hashmap
          auto insert_itr_two( m_HashMap( 0).Insert( std::make_pair( atom_env_str_two, size_t( 1))));
          auto insert_itr_three( m_HashMap( 1).Insert( std::make_pair( atom_env_str_three, size_t( 1))));
          auto insert_itr_four( m_HashMap( 2).Insert( std::make_pair( atom_env_str_four, size_t( 1))));

          // if we failed to insert it (i.e. because we already have one) then increment the count
          if( !insert_itr_two.second)
          {
            insert_itr_two.first->second += 1;
          }
          if( !insert_itr_three.second)
          {
            insert_itr_three.first->second += 1;
          }
          if( !insert_itr_four.second)
          {
            insert_itr_four.first->second += 1;
          }
        }
      }

      // output counts
      if( m_OutputFilenameBase->GetFlag() && m_BondRadius->GetFirstParameter()->GetNumericalValue< size_t>() == size_t( 2))
      {
        io::OFStream raw_counts_out;
        io::File::MustOpenOFStream( raw_counts_out, m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".atom_environment_hashmap.two_bonds.txt");
        raw_counts_out << m_HashMap( 0) << '\n';
        io::File::CloseClearFStream( raw_counts_out);
      }
      if( m_OutputFilenameBase->GetFlag() && m_BondRadius->GetFirstParameter()->GetNumericalValue< size_t>() == size_t( 3))
      {
        io::OFStream raw_counts_out;
        io::File::MustOpenOFStream( raw_counts_out, m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".atom_environment_hashmap.three_bonds.txt");
        raw_counts_out << m_HashMap( 1) << '\n';
        io::File::CloseClearFStream( raw_counts_out);
      }
      if( m_OutputFilenameBase->GetFlag() && m_BondRadius->GetFirstParameter()->GetNumericalValue< size_t>() == size_t( 4))
      {
        io::OFStream raw_counts_out;
        io::File::MustOpenOFStream( raw_counts_out, m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".atom_environment_hashmap.four_bonds.txt");
        raw_counts_out << m_HashMap( 2) << '\n';
        io::File::CloseClearFStream( raw_counts_out);
      }
    }

  } // namespace app
} // namespace bcl

