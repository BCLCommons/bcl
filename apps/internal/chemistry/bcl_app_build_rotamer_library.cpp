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
#include "bcl_app_build_rotamer_library.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_graph_marker.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_rotamer_cluster_center.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "chemistry/bcl_chemistry_valence_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_matrix3x3.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "quality/bcl_quality_rmsd.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "storage/bcl_storage_triplet.h"

namespace bcl
{
  namespace app
  {

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    BuildRotamerLibrary::BuildRotamerLibrary() :
      m_StructureDatabaseFlag
      (
        new command::FlagStatic
        (
          "structure_database",
          "sdf file containing the molecules from which to find fragments",
          command::Parameter
          (
            "structure_database",
            "an sdf input file",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_DihedralFragmentFlag
      (
        new command::FlagDynamic
        (
          "dihedral_fragments",
          "get fragments for dihedral bonds",
          command::Parameter
          (
            "dihedral_fragments",
            "get fragments for dihedral bonds"
          ),
          0,
          0
        )
      ),
      m_BondAngleFragmentFlag
      (
        new command::FlagStatic
        (
          "bond_angles",
          "get fragments for bond angles"
        )
      ),
      m_IdentifyNewConformersFlag
      (
        new command::FlagDynamic
        (
          "add_rotamers",
          "identify new conformers in conformation ensemble that are not present in native ensemble",
          command::Parameter
          (
            "add_rotamers",
            "identify new conformers in conformation ensemble that are not present in native ensemble"
          ),
          0,
          1
        )
      ),
      m_BinSizeFlag
      (
        new command::FlagStatic
        (
          "bin_size",
          "Bin size to be used for determining conformations using dihedralbins, usually 15, 30 or 60",
          command::Parameter
          (
            "bin_size",
            "Bin size to be used for determining conformations using dihedralbins, usually 15, 30 or 60",
            command::ParameterCheckRanged< double>( 0.0, 360.0),
            "30.0"
          )
        )
      ),
      m_OutputFileFlag
      (
        new command::FlagStatic
        (
          "output",
          "filename to output sdf containing conformers",
          command::Parameter
          (
            "output_filename",
            "filename to output sdf containing conformers"
          )
        )
      ),
      m_MinFragmentCount
      (
        new command::FlagStatic
        (
          "min_fragment_count",
          "minimum count that a fragment should have to be written out",
          command::Parameter
          (
            "min_fragment_count",
            "minimum count that a fragment should have to be written out",
            "3"
          )
        )
      ),
      m_MaxFragmentCount
      (
        new command::FlagStatic
        (
          "max_fragment_count",
          "maximum count that a fragment should have to be written out since clustering will be limiting factor",
          command::Parameter
          (
            "max_fragment_count",
            "maximum count that a fragment should have to be written out since clustering will be limiting factor",
            "1000"
          )
        )
      ),
      m_MinRotamerCount
      (
        new command::FlagStatic
        (
          "min_rotamer_count",
          "minimum count that a rotamer should have to be accepted",
          command::Parameter
          (
            "min_rotamer_count",
            "minimum count that a rotamer should have to be accepted",
            "3"
          )
        )
      ),
      m_CountFactor
      (
        new command::FlagDynamic
        (
          "count_weight",
          "weight to be given to rotamers seen in the conformer library",
          command::Parameter
          (
            "count_weight",
            "weight to be given to rotamers seen in the conformer library",
            "1.0"
          ),
          0,
          1
        )
      ),
      m_GetIndividualRotamers
      (
        new command::FlagDynamic
        (
          "get_individual_rotamers",
          "if individual rotamers are desired then provide the flag with filename",
          command::Parameter
          (
            "get_individual_rotamers",
            "if individual rotamers are desired then provide the flag with filename rotamers"
          ),
          0,
          1
        )
      ),
      m_AddSimulatedCounts
      (
        new command::FlagStatic
        (
          "add_simulated_counts",
          "if set, add simulated rotamer counts as well for use in computing energies"
        )
      )
    {
    }

    //! copy constructor
    BuildRotamerLibrary::BuildRotamerLibrary( const BuildRotamerLibrary &APP) :
      m_StructureDatabaseFlag( APP.m_StructureDatabaseFlag),
      m_DihedralFragmentFlag( APP.m_DihedralFragmentFlag),
      m_BondAngleFragmentFlag( APP.m_BondAngleFragmentFlag),
      m_IdentifyNewConformersFlag( APP.m_IdentifyNewConformersFlag),
      m_BinSizeFlag( APP.m_BinSizeFlag),
      m_OutputFileFlag( APP.m_OutputFileFlag),
      m_MinFragmentCount( APP.m_MinFragmentCount),
      m_MaxFragmentCount( APP.m_MaxFragmentCount),
      m_MinRotamerCount( APP.m_MinRotamerCount),
      m_CountFactor( APP.m_CountFactor),
      m_GetIndividualRotamers( APP.m_GetIndividualRotamers),
      m_AddSimulatedCounts( APP.m_AddSimulatedCounts)
    {
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> BuildRotamerLibrary::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // add command line options to add mdl lines
      sdf::MdlHandler::AddAtomMdlLineFlag( *sp_cmd);

      // ensembles containing the molecules to be edited
      chemistry::FragmentFeed::AddFlags( *sp_cmd, "fragments");

      // insert all the flags and params
      sp_cmd->AddFlag( m_StructureDatabaseFlag);
      sp_cmd->AddFlag( m_DihedralFragmentFlag);
      sp_cmd->AddFlag( m_BondAngleFragmentFlag);
      sp_cmd->AddFlag( m_IdentifyNewConformersFlag);
      sp_cmd->AddFlag( m_BinSizeFlag);
      sp_cmd->AddFlag( m_OutputFileFlag);
      sp_cmd->AddFlag( m_MinFragmentCount);
      sp_cmd->AddFlag( m_MaxFragmentCount);
      sp_cmd->AddFlag( m_MinRotamerCount);
      sp_cmd->AddFlag( m_CountFactor);
      sp_cmd->AddFlag( m_GetIndividualRotamers);
      sp_cmd->AddFlag( m_AddSimulatedCounts);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    BuildRotamerLibrary *BuildRotamerLibrary::Clone() const
    {
      return new BuildRotamerLibrary( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &BuildRotamerLibrary::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string BuildRotamerLibrary::GetDescription() const
    {
      return "BuildRotamerLibrary does what it does";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &BuildRotamerLibrary::GetReadMe() const
    {
      static std::string s_read_me
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::BuildRotamerLibrary, terms of use, "
        "appropriate citation, installation procedures, BCL::BuildRotamerLibrary execution, "
        "technical support, and future research directions.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::BuildRotamerLibrary?\n"
        "BCL::BuildRotamerLibrary is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons. BCL::BuildRotamerLibrary is a utility that creates a "
        "rotamer library for given fragments using a database of interest."
        "\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::BuildRotamerLibrary.\n"
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
        "VI. RUNNING BCL::BuildRotamerLibrary.\n"
        "Running BCL::BuildRotamerLibrary requires sdf files containing fragments and structure database.\n"
        "\n"
        "2) Run BCL::BuildRotamerLibrary to generate rotamer library for fragments of interest\n"
        "\n"
        "At a command prompt, navigate to the location of your BCL::BuildRotamerLibrary executable program. The syntax for"
        "running the application looks like the following"
        "\n"
        "bcl.exe BuildRotamerLibrary -structure_database <filename> -fragments <filename>-output <filename>"
        "\n\nFLAGS:\n\n"
        "-structure_database <filename> -> file containing structure database from which rotamer statistics will be derived \n"
        "-fragments <filename> -> file containing fragments for which rotamer library needs to be generated \n"
        "-output <filename> -> file to which the rotamer library will be written out\n"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe BuildRotamerLibrary -help\n"
        "\n"
        "For more general information about the product, type bcl.exe BuildRotamerLibrary -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::BuildRotamerLibrary.\n"
        "BCL::BuildRotamerLibrary is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );
      return s_read_me;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &BuildRotamerLibrary::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::BuildRotamerLibrary is a utility for generating rotamer library for fragments using a database of interest"
        "Features of BCL::BuildRotamerLibrary\n"
        "<ul>"
        "  <li>Generates rotamer library for fragments\n"
        "  </li>"
        "  <li>Compressed molecule files (bz2, gzip) are supported</li>"
        "</ul>\n\n"
        "!build_rotamer_library.png!"
      );

      return s_web_text;
    }

    //! @brief determine rings and bonds contained in the fragment
    //! @param FRAGMENT molecule whose rings are to be determined
    //! @return a list of list of indices of atoms. Inner list has atom indices contained in a single ring.
    storage::List< storage::List< size_t> > BuildRotamerLibrary::DetermineRings
    (
      const chemistry::FragmentComplete &FRAGMENT
    ) const
    {
      // get list of vertices of fragment that are in a ring and contained in the molecule of interest
      storage::List< storage::Vector< size_t> > ring_vertices
      (
        chemistry::SmallMoleculeFragmentMapping::GetRingVertices( FRAGMENT)
      );

      storage::Vector< storage::VectorND< 4, size_t> > fragment_priority( chemistry::PriorityDihedralAngles()( FRAGMENT).Second());

      storage::List< storage::List< size_t> > ring_bond_indices;
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator
            itr_vert( ring_vertices.Begin()), itr_vert_end( ring_vertices.End());
        itr_vert != itr_vert_end;
        ++itr_vert
      )
      {
        size_t bond_index( 0);
        storage::List< size_t> current_ring_bonds;
        for
        (
          storage::Vector< storage::VectorND< 4, size_t> >::const_iterator
            itr_priority_bonds( fragment_priority.Begin()), itr_priority_bonds_end( fragment_priority.End());
          itr_priority_bonds != itr_priority_bonds_end;
          ++itr_priority_bonds, ++bond_index
        )
        {
          // if ring vertices contains center bond atoms then the bond is contained in the current ring
          if
          (
              chemistry::RotamerDihedralBondData::IfFullIntersection
              (
                *itr_vert,
                storage::Vector< size_t>::Create( itr_priority_bonds->Second(), itr_priority_bonds->Third())
              )
          )
          {
            current_ring_bonds.PushBack( bond_index);
          }
        }
        ring_bond_indices.PushBack( current_ring_bonds);
      }
      return ring_bond_indices;
    }

    //! @brief determine if rings of molecule contain multiple ring rotamers
    //! @param RING_BONDS a list of list of indices of atoms. Inner list has atom indices contained in a single ring
    //! @param BIN_VECTOR rotamer bins for the corresponding rings
    //! @return true if any ring has multiple rotamers else false
    bool BuildRotamerLibrary::IfContainsRingConformations
    (
      const storage::List< storage::List< size_t> > &RING_BONDS,
      const storage::List< storage::Vector< size_t> > &BIN_VECTOR
    ) const
    {
      // go through the ring bonds and from each rotamer get the dihedral bin for the bond
      for
      (
        storage::List< storage::List< size_t> >::const_iterator itr( RING_BONDS.Begin()), itr_end( RING_BONDS.End());
        itr != itr_end;
        ++itr
      )
      {
        // store the only binning seen for this ring so far
        storage::Vector< size_t> binning;
        for
        (
          storage::List< storage::Vector< size_t> >::const_iterator
            itr_rotamer_bins( BIN_VECTOR.Begin()), itr_rotamer_bins_end( BIN_VECTOR.End());
          itr_rotamer_bins != itr_rotamer_bins_end;
          ++itr_rotamer_bins
        )
        {
          storage::Vector< size_t> ring_bins;
          const storage::Vector< size_t> &current_rotamer( *itr_rotamer_bins);
          for
          (
            storage::List< size_t>::const_iterator itr_inner( itr->Begin()), itr_inner_end( itr->End());
            itr_inner != itr_inner_end;
            ++itr_inner
          )
          {
            ring_bins.PushBack( current_rotamer( *itr_inner));
          }
          if( binning.IsEmpty())
          {
            binning = ring_bins;
          }
          else if( !( binning == ring_bins))
          {
            return true;
          }
        }
      }
      return false;
    }

    //! @brief finds all the conformers of a fragment inside an ensemble
    //! @param FRAGMENT the fragment to find conformers of
    //! @param FRAGMENT_INDEX index of the fragment in the fragment ensemble
    //! @return void
    void BuildRotamerLibrary::FindConformers
    (
      const chemistry::FragmentComplete &FRAGMENT,
      const size_t &FRAGMENT_INDEX
    ) const
    {
      // make graphs of the small molecule
      // this graph will be used to compute the actual isomorphism
      // Vertices are colored by atom type
      // Edges are colored by bond order (including aromatic)
      const graph::ConstGraph< size_t, size_t> scaffold_simple_graph
      (
        m_FragmentSimpleGraphs( FRAGMENT_INDEX)
      );

      // the isomorphism will then be applied to extract the appropriate scaffold from this graph, which
      // is colored by pointers to the actual atoms, and size_t's for the bond types
      // bond types can't be used directly here due to limitations of math::Matrix, which ConstGraph uses to hold the
      // edges
      graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> scaffold_graph
      (
        chemistry::ConformationGraphConverter::CreateGraphWithAtoms( FRAGMENT)
      );

      // store all the different conformations of the FRAGMENT in conformers
      chemistry::FragmentEnsemble conformers;

      // object that will compute the fragment isomorphism
      graph::SubgraphIsomorphism< size_t, size_t> simple_csi;

      // after the isomorphism is complete, we will convert it back into this type of graph
      graph::SubgraphIsomorphism< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> csi_with_atoms;

      // set the smaller graph of the substructure to the scaffold graph
      simple_csi.SetSubgraphExternalOwnership( scaffold_simple_graph);
      csi_with_atoms.SetSubgraphExternalOwnership( scaffold_graph);

      // create instance of class conformerclustercenter for the molecule of interest
      chemistry::RotamerClusterCenter cluster_center( FRAGMENT, m_BinSize);

      bool seen_enough_times( false);
      storage::List< size_t> super_fragment_mols;

      size_t minimum_molecules_seen( m_MinFragmentCount->GetFirstParameter()->GetNumericalValue< size_t>());
      // now go over all the molecules provided in the ensemble to find isomorphism between the molecule of interest and
      // any of the molecules from ensemble
      for
      (
        size_t ensemble_number( 0), ensemble_size( m_EnsembleGraphs.GetSize());
        ensemble_number < ensemble_size;
        ++ensemble_number
      )
      {
        if
        (
          scaffold_simple_graph.GetSize() > m_MoleculeSimpleGraphs( ensemble_number).GetSize()
          || scaffold_simple_graph.NumEdges() > m_MoleculeSimpleGraphs( ensemble_number).NumEdges()
        )
        {
          continue;
        }
        simple_csi.SetGraphExternalOwnership( m_MoleculeSimpleGraphs( ensemble_number));

        if( !seen_enough_times)
        {
          if( simple_csi.FindIsomorphism())
          {
            super_fragment_mols.PushBack( ensemble_number);
            if( super_fragment_mols.GetSize() >= minimum_molecules_seen - size_t( 1))
            {
              seen_enough_times = true;

              for
              (
                storage::List< size_t>::const_iterator
                  itr_list( super_fragment_mols.Begin()), itr_list_end( super_fragment_mols.End());
                itr_list != itr_list_end;
                ++itr_list
              )
              {
                simple_csi.SetGraphExternalOwnership( m_MoleculeSimpleGraphs( *itr_list));
                // check if there are enough edges and vertices of the right type in the ensemble's graph
                if( simple_csi.FindAllIsomorphisms())
                {
                  if( util::GetMessenger().GetCurrentMessageLevel() >= util::Message::e_Verbose)
                  {
                    util::GetLogger().LogStatus
                    (
                      "Found: " + util::Format()( simple_csi.GetIsomorphisms().GetSize())
                      + " isomorphisms for molecule: " + util::Format()( *itr_list)
                    );
                  }
                  csi_with_atoms.SetGraphExternalOwnership( m_EnsembleGraphs( *itr_list));
                  const storage::Vector< storage::Vector< size_t> > &simple_isomorphisms( simple_csi.GetIsomorphisms());
                  csi_with_atoms.SetIsomorphisms( simple_isomorphisms);
                  cluster_center( csi_with_atoms, simple_isomorphisms, FRAGMENT.GetName());
                }
              }
            }
          }
        }
        // check if there are enough edges and vertices of the right type in the ensemble's graph
        else if( simple_csi.FindAllIsomorphisms())
        {
          if( util::GetMessenger().GetCurrentMessageLevel() >= util::Message::e_Verbose)
          {
            util::GetLogger().LogStatus
            (
              "Found: " + util::Format()( simple_csi.GetIsomorphisms().GetSize())
              + " isomorphisms for molecule: " + util::Format()( ensemble_number)
            );
          }
          csi_with_atoms.SetGraphExternalOwnership( m_EnsembleGraphs( ensemble_number));
          const storage::Vector< storage::Vector< size_t> > &simple_isomorphisms( simple_csi.GetIsomorphisms());
          csi_with_atoms.SetIsomorphisms( simple_isomorphisms);
          cluster_center( csi_with_atoms, simple_isomorphisms, FRAGMENT.GetName());
        }
      }
      if( cluster_center.GetClusters().IsEmpty())
      {
        if( m_IdentifyNewConformersFlag->GetFlag())
        {
          FRAGMENT.WriteMDL( m_Output);
        }
        return;
      }

      size_t number_chain_dihedrals( 0), number_ring_dihedrals( 0);
      chemistry::PriorityDihedralAngles pda;
      auto dihedrals( pda.GetDihedralEdges( FRAGMENT));
      for( auto itr_dih( dihedrals.Begin()), itr_dih_end( dihedrals.End()); itr_dih != itr_dih_end; ++itr_dih)
      {
        itr_dih->GetEdgeData()->IsBondInRing() ? ++number_ring_dihedrals : ++number_chain_dihedrals;
      }

      bool contains_rings( false);
      const storage::List< storage::List< size_t> > ring_bonds
      (
        !number_chain_dihedrals
        ? DetermineRings( FRAGMENT)
        : storage::List< storage::List< size_t> >()
      );

      if( !ring_bonds.IsEmpty())
      {
        contains_rings = true;
      }

      size_t max_instance( cluster_center.GetMaxInstance());
      size_t minimum_rotamer_count( m_MinRotamerCount->GetFirstParameter()->GetNumericalValue< size_t>());
      if( max_instance < double( minimum_rotamer_count))
      {
        if( m_IdentifyNewConformersFlag->GetFlag())
        {
          FRAGMENT.WriteMDL( m_Output);
        }
        return;
      }
      math::RunningAverageSD< double> rotamer_stats( cluster_center.GetStats());
      double rotamer_sum( rotamer_stats.GetAverage() * rotamer_stats.GetWeight());
      if( rotamer_sum < minimum_molecules_seen)
      {
        if( m_IdentifyNewConformersFlag->GetFlag())
        {
          FRAGMENT.WriteMDL( m_Output);
        }
        return;
      }
      if( number_ring_dihedrals || m_GetIndividualRotamers->GetFlag())
      {
        // once all instances have been collected determine rotamer centeres
        cluster_center.CalculateClusterCenters
        (
          contains_rings && !number_chain_dihedrals
          ? m_MaxFragmentCount->GetFirstParameter()->GetNumericalValue< size_t>()
          : size_t( 10)
        );
        cluster_center.PruneNonPlanarAromaticRings();
      }
      if( minimum_rotamer_count)
      {
        cluster_center.PruneRotamersByInstances( minimum_rotamer_count);
      }
      chemistry::FragmentComplete molecule_first( cluster_center.GetAverageStucture( m_MaxFragmentCount->GetFirstParameter()->GetNumericalValue< size_t>()));

      if( m_IdentifyNewConformersFlag->GetFlag())
      {
        AddRotamers( molecule_first, FRAGMENT_INDEX, cluster_center, ring_bonds);
        return;
      }
      if( m_AddSimulatedCounts->GetFlag() && number_chain_dihedrals)
      {
        chemistry::FragmentEnsemble simulated( m_SampleConformations( molecule_first).First());
        simple_csi.SetGraphExternalOwnership( scaffold_simple_graph);
        simple_csi.FindAllIsomorphisms();
        BCL_MessageStd
        (
          "Adding counts for molecule with " + util::Format()( simple_csi.GetIsomorphisms().GetSize())
          + " internal isomorphisms"
        );
        size_t sim_number( 0), nsim( simulated.GetSize());
        for( auto itrsim( simulated.Begin()), itrsim_end( simulated.End()); itrsim != itrsim_end; ++itrsim, ++sim_number)
        {
          graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> simulated_graph
          (
            chemistry::ConformationGraphConverter::CreateGraphWithAtoms( *itrsim)
          );
          if( sim_number < 100 || ( sim_number % 100) == 0)
          {
            util::GetLogger().LogStatus( "On simulated molecule: " + util::Format()( sim_number) + " / " + util::Format()( nsim));
          }
          csi_with_atoms.SetGraphExternalOwnership( simulated_graph);
          csi_with_atoms.SetIsomorphisms( simple_csi.GetIsomorphisms());
          cluster_center( csi_with_atoms, simple_csi.GetIsomorphisms(), FRAGMENT.GetName(), true);
        }
      }
      OutputRotamerLibrary
      (
        molecule_first,
        FRAGMENT_INDEX,
        cluster_center,
        ring_bonds
      );
    }

    //! @brief updates rotamer information for given fragment given a set of rotamers in the form of RotamerClusterCenter
    //! @param FRAGMENT the fragment whose rotamers need to be updated
    //! @param FRAGMENT_INDEX index of the fragment in the fragment ensemble
    //! @param CLUSTER_CENTER contains rotamer information that needs to be added to the fragment of interest
    //! @param RING_BONDS ring bonds of fragment.
    //! @return void
    void BuildRotamerLibrary::AddRotamers
    (
      const chemistry::FragmentComplete &FRAGMENT,
      const size_t FRAGMENT_INDEX,
      const chemistry::RotamerClusterCenter &CLUSTER_CENTER,
      const storage::List< storage::List< size_t> > &RING_BONDS
    ) const
    {
      chemistry::FragmentComplete molecule_first( FRAGMENT);

      // conformation tolerance crieria
      chemistry::ConformationComparisonByDihedralBins conformation_comparison( m_BinSize);

      size_t rotamer_index_a( 0);
      double rotamer_counts( 0);
      storage::List< linal::Vector< int> > fragment_rotamers_a;
      bool contains_ring_conformations( false);
      if( FRAGMENT.IsPropertyStored( "NRotamers"))
      {
        // get number of rotamers of the fragment in question
        rotamer_index_a = util::ConvertStringToNumericalValue< size_t>( FRAGMENT.GetStoredProperties().GetMDLProperty( "NRotamers"));

        rotamer_counts = FRAGMENT.GetStoredProperties().GetMDLPropertyAsVector( "ScaffoldCount").First() / rotamer_index_a;

        if( m_CountFactor->GetFlag())
        {
          rotamer_counts = rotamer_counts * util::ConvertStringToNumericalValue< double>( m_CountFactor->GetFirstParameter()->GetValue());
        }

        // for each rotamer get rotamer counts, average angle, rotamer bin signature, coordinates
        for( size_t rotamer_count( 1); rotamer_count <= rotamer_index_a; ++rotamer_count)
        {
          // rotamer count of a particular rotamer
          fragment_rotamers_a.PushBack
          (
            FRAGMENT.GetStoredProperties().GetMDLPropertyAsVector
            (
              "Rotamer" + util::Format()( rotamer_count) + "Bins"
            )
          );
        }
      }

      // now go over each rotamer of molecule and add data about all rotamers to non-const version of molecule
      size_t rotamer_index( rotamer_index_a + 1);
      double scaffold_count( 0);
      size_t common_rotamers( 0);

      // iterate through each cluster and create rotamer library from the information
      for
      (
        storage::List< chemistry::RotamerEnsemble>::const_iterator
          itr_cluster( CLUSTER_CENTER.GetClusters().Begin()), itr_cluster_end( CLUSTER_CENTER.GetClusters().End());
        itr_cluster != itr_cluster_end;
        ++itr_cluster, ++rotamer_index
      )
      {
        // get bin string, average and standard deviation for rotamer
        const linal::Vector< double> &average( itr_cluster->GetAverage());

        bool conformer_exists( false);
        // get bin vector for each dihedral angle of fragment
        const linal::Vector< int> bin_vector( conformation_comparison.DetermineDihedralKeys( average, true, true));

        for
        (
          storage::List< linal::Vector< int> >::const_iterator
            itr_a( fragment_rotamers_a.Begin()), itr_a_end( fragment_rotamers_a.End());
          itr_a != itr_a_end;
          ++itr_a
        )
        {
          if( bin_vector == *itr_a)
          {
            ++common_rotamers;
            conformer_exists = true;
            break;
          }
        }

        if( conformer_exists)
        {
          --rotamer_index;
          continue;
        }

        bool has_ring_conformations( HasRingConformations( *itr_cluster, fragment_rotamers_a));
        // if the rotamer from given library doesnt have ring conformatin but a ring conformation is found in this structure
        // database, there is inconsistency. Hence dont add this rotamer
        if( !contains_ring_conformations && has_ring_conformations)
        {
          --rotamer_index;
          continue;
        }

        // get rotamer index string
        const std::string rotamer_st( "Rotamer" + util::Format()( rotamer_index));

        // sum rotamer counts to get total count for fragment of interest
        scaffold_count += rotamer_counts;

        // store rotamer count and coordinates for current rotamer
        molecule_first.StoreProperty( rotamer_st + "Count", storage::Vector< double>( 1, double( rotamer_counts)));

        if( contains_ring_conformations)
        {
          // get coordinates for current rotamer in a string
          std::stringstream atom_vector;
          const storage::Vector< sdf::AtomInfo> atom_info_cur( itr_cluster->GetClusterCenter().GetAtomInfo());
          for
          (
              storage::Vector< sdf::AtomInfo>::const_iterator
              itr_info( atom_info_cur.Begin()), itr_info_end( atom_info_cur.End());
              itr_info != itr_info_end;
              ++itr_info
          )
          {
            atom_vector << itr_info->ToMdlAtomLine() << "\n";
          }
          molecule_first.StoreProperty( rotamer_st + "Coordinates", atom_vector.str());
        }

        std::stringstream average_string;
        std::stringstream bin_string;

        // store bin signature, average angle and std dev of angles
        size_t dihedral_index( 0);
        for
        (
          auto itr_index( bin_vector.Begin()), itr_index_end( bin_vector.End());
          itr_index != itr_index_end;
          ++itr_index, ++dihedral_index
        )
        {
          // format the string for representing bins, average and standard deviation
          bin_string << util::Format().W( 8).L()( *itr_index);
          average_string << util::Format().FFP( 2).W( 8).L()( average( dihedral_index));
        }
        // store bin_string, average dihedral angle and standard deviation for current rotamer
        molecule_first.StoreProperty( rotamer_st + "Bins", bin_string.str());
        molecule_first.StoreProperty( rotamer_st + "Ave", average_string.str());
      }
      if( common_rotamers)
      {
        BCL_MessageStd
        (
          "Unexpected: " + util::Format()( common_rotamers)
          + " out of " + util::Format()( CLUSTER_CENTER.GetClusters().GetSize())
          + " average rotamers mapped to same rotamer"
        );
      }

      if( FRAGMENT.IsPropertyStored( "NRotamers"))
      {
        // add number of different rotamers to the non-const version of molecule
        molecule_first.RemoveProperty( "NRotamers");
        size_t prev_counts( FRAGMENT.GetStoredProperties().GetMDLPropertyAsVector( "ScaffoldCount").First());
        molecule_first.RemoveProperty( "ScaffoldCount");
        molecule_first.StoreProperty( "NRotamers", storage::Vector< double>( 1, double( rotamer_index - 1)));
        molecule_first.StoreProperty( "ScaffoldCount", storage::Vector< double>( 1, double( prev_counts + scaffold_count)));
        molecule_first.WriteMDL( m_Output);
      }
    }

    //! @brief determine wheter rotamers have different ring conformations
    //! @param FRAGMENT rotamers for a fragment of interest
    //! @param FRAGMENT_INDEX rotamers that have been determined for the fragment of interest
    //! @return true if different ring conformations exist, otherwise false
    bool BuildRotamerLibrary::HasRingConformations
    (
      const chemistry::RotamerEnsemble &FRAGMENT,
      const storage::List< linal::Vector< int> > &FRAGMENT_ROTAMERS
    ) const
    {
      // get list of vertices of fragment that are in a ring and contained in the molecule of interest
      storage::List< storage::Vector< size_t> > ring_vertices
      (
        chemistry::SmallMoleculeFragmentMapping::GetRingVertices( FRAGMENT.GetRotamers().FirstElement())
      );

      storage::Vector< graph::UndirectedEdge< chemistry::ConfigurationalBondType> > priority_edges
      (
        chemistry::PriorityDihedralAngles::GetDihedralEdges( FRAGMENT.GetRotamers().FirstElement())
      );

      storage::List< storage::List< size_t> > ring_bonds;

      // iterate through the rings contained in the fragment of interest
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator
          itr_ring( ring_vertices.Begin()), itr_ring_end( ring_vertices.End());
        itr_ring != itr_ring_end;
        ++itr_ring
      )
      {
        // get vertices contained in the current ring
        const storage::Vector< size_t> &cur_ring( *itr_ring);
        storage::List< size_t> cur_ring_bonds;

        // go through the center bonds of the molecule that the fragment contains
        size_t bond_index( 0);
        for
        (
          storage::Vector< graph::UndirectedEdge< chemistry::ConfigurationalBondType> >::const_iterator
            itr_center_bonds( priority_edges.Begin()), itr_center_bonds_end( priority_edges.End());
          itr_center_bonds != itr_center_bonds_end;
          ++itr_center_bonds, ++bond_index
        )
        {
          // if ring vertices contains center bond atoms then the bond is contained in the current ring
          if
          (
            chemistry::RotamerDihedralBondData::IfFullIntersection
            (
              cur_ring,
              storage::Vector< size_t>::Create( itr_center_bonds->GetIndexLow(), itr_center_bonds->GetIndexHigh())
            )
          )
          {
            cur_ring_bonds.PushBack( bond_index);
          }
        }
        ring_bonds.PushBack( cur_ring_bonds);
      }

      // get bin vector for each dihedral angle of fragment
      const storage::Vector< int> bin_vector
      (
        chemistry::ConformationComparisonByDihedralBins( m_BinSize).DetermineDihedralKeys( FRAGMENT.GetAverage(), true, true)
      );

      // go through the ring bonds and from each rotamer get the dihedral bin for the bond
      for
      (
        storage::List< storage::List< size_t> >::const_iterator itr( ring_bonds.Begin()), itr_end( ring_bonds.End());
        itr != itr_end;
        ++itr
      )
      {
        // store the only binning seen for this ring so far
        storage::Vector< int> binning;
        for
        (
          storage::List< size_t>::const_iterator itr_inner( itr->Begin()), itr_inner_end( itr->End());
          itr_inner != itr_inner_end;
          ++itr_inner
        )
        {
          binning.PushBack( bin_vector( *itr_inner));
        }

        for
        (
          storage::List< linal::Vector< int> >::const_iterator
            itr_rotamer_bins( FRAGMENT_ROTAMERS.Begin()),
            itr_rotamer_bins_end( FRAGMENT_ROTAMERS.End());
          itr_rotamer_bins != itr_rotamer_bins_end;
          ++itr_rotamer_bins
        )
        {
          storage::Vector< int> ring_bins;
          const linal::Vector< int> &current_rotamer( *itr_rotamer_bins);
          for
          (
            storage::List< size_t>::const_iterator itr_inner( itr->Begin()), itr_inner_end( itr->End());
            itr_inner != itr_inner_end;
            ++itr_inner
          )
          {
            ring_bins.PushBack( current_rotamer( *itr_inner));
          }

          if( !( binning == ring_bins))
          {
            return true;
          }
        }
      }
      return false;
    }

    //! @brief finds all the conformers of a fragment inside an ensemble
    //! @param FRAGMENT the fragment to find conformers of
    //! @param FRAGMENT_INDEX index of the fragment in the fragment ensemble
    //! @param CLUSTER_CENTER object containing all rotamers for the fragment of interest
    //! @param RING_BONDS true if fragment contains ring, false if not.
    void BuildRotamerLibrary::OutputRotamerLibrary
    (
      const chemistry::FragmentComplete &FRAGMENT,
      const size_t FRAGMENT_INDEX,
      const chemistry::RotamerClusterCenter &CLUSTER_CENTER,
      const storage::List< storage::List< size_t> > &RING_BONDS
    ) const
    {
      chemistry::FragmentComplete molecule_output( FRAGMENT);
      storage::Vector< storage::VectorND< 4, size_t> > priority_dihedral( CLUSTER_CENTER.GetPriority()( FRAGMENT).Second());

      // conformation tolerance crieria
      chemistry::ConformationComparisonByDihedralBins conformation_comparison( m_BinSize);

      // now go over each rotamer of molecule and add data about all rotamers to non-const version of molecule
      size_t rotamer_index( 1);
      double scaffold_count( 0);

      storage::List< storage::Vector< size_t> > all_rotamer_bins;

      double count_weight( 1.0);
      if( m_CountFactor->GetFlag())
      {
        count_weight = util::ConvertStringToNumericalValue< double>( m_CountFactor->GetFirstParameter()->GetValue());
      }

      io::OFStream individual_rotamers;
      if( m_GetIndividualRotamers->GetFlag())
      {
        // open the output file so that we can write rotamer library as we go
        io::File::MustOpenOFStream( individual_rotamers, m_GetIndividualRotamers->GetFirstParameter()->GetValue() + "_" + util::Format()( FRAGMENT_INDEX) + ".sdf");
      }
      storage::Vector< double> rotamer_counts_vec;
      rotamer_counts_vec.AllocateMemory( CLUSTER_CENTER.GetClusters().GetSize());
      storage::List< double> average_vec;
      storage::List< double> std_dev_vec;
      storage::List< double> bin_vec;
      storage::Vector< double> simulated_vec;
      // iterate through each cluster and create rotamer library from the information
      for
      (
        storage::List< chemistry::RotamerEnsemble>::const_iterator
          itr_cluster( CLUSTER_CENTER.GetClusters().Begin()), itr_cluster_end( CLUSTER_CENTER.GetClusters().End());
        itr_cluster != itr_cluster_end;
        ++itr_cluster, ++rotamer_index
      )
      {
        // get rotamer count and skip adding rotamer if count is less than minimum required
        double rotamer_counts( itr_cluster->GetWeights());

        rotamer_counts = double( rotamer_counts * count_weight);

        const std::string rotamer_st( "Rotamer" + util::Format()( rotamer_index));

        // sum rotamer counts to get total count for fragment of interest
        scaffold_count += rotamer_counts;
        rotamer_counts_vec.PushBack( rotamer_counts);

        // get bin string, average and standard deviation for rotamer
        const linal::Vector< double> &average( itr_cluster->GetAverage());
        const storage::Vector< size_t> bin_vector( conformation_comparison.DetermineDihedralKeys( average, true, true));

        average_vec.InsertElements( average_vec.End(), average.Begin(), average.End());
        bin_vec.InsertElements( bin_vec.End(), bin_vector.Begin(), bin_vector.End());
        simulated_vec.PushBack( itr_cluster->GetSimulatedInstances());

        // get bin vector for each dihedral angle of fragment

        if( !RING_BONDS.IsEmpty())
        {
          storage::Vector< sdf::AtomInfo> cur_atom_info( itr_cluster->GetClusterCenter().GetAtomInfo());
          std::stringstream atom_vector;

          for
          (
            storage::Vector< sdf::AtomInfo>::const_iterator
            itr_info( cur_atom_info.Begin()), itr_info_end( cur_atom_info.End());
            itr_info != itr_info_end;
            ++itr_info
          )
          {
            atom_vector << itr_info->ToMdlAtomLine() << "\n";
          }

          molecule_output.StoreProperty( rotamer_st + "Coordinates", atom_vector.str());

          all_rotamer_bins.PushBack( bin_vector);
        }

        // store bin_string, average dihedral angle and standard deviation for current rotamer
        if( m_GetIndividualRotamers->GetFlag())
        {
          chemistry::FragmentComplete cur_rot( itr_cluster->GetClusterCenter());

          cur_rot.StoreProperty( rotamer_st + "Count", storage::Vector< double>( 1, double( rotamer_counts)));
          cur_rot.StoreProperty( rotamer_st + "Bins", bin_vector);
          cur_rot.StoreProperty( rotamer_st + "Ave", average);
          if( m_AddSimulatedCounts->GetFlag())
          {
            cur_rot.StoreProperty( rotamer_st + "CountsSimulated", util::Format()( itr_cluster->GetSimulatedInstances()));
          }

          if( rotamer_counts > 4)
          {
            cur_rot.WriteMDL( individual_rotamers);
          }
        }
      }
      molecule_output.StoreProperty( "RotamerCounts", rotamer_counts_vec);
      molecule_output.StoreProperty( "Bins", linal::Vector< size_t>( bin_vec.Begin(), bin_vec.End()));
      molecule_output.StoreProperty( "Average", linal::Vector< size_t>( average_vec.Begin(), average_vec.End()));
      if( m_AddSimulatedCounts->GetFlag())
      {
        molecule_output.StoreProperty( "RotamerCountsSimulated", simulated_vec);
        molecule_output.StoreProperty( "NovelRotamersSimulated", util::Format()( CLUSTER_CENTER.GetSimulatedWeightUnseenRotamers()));
      }
      if( m_GetIndividualRotamers->GetFlag())
      {
        // open the output file so that we can write rotamer library as we go
        io::File::CloseClearFStream( individual_rotamers);
      }

      // if no rotamer has been seen more than minimum number of times, then do not add fragment to library
//      if( rotamer_index == 1)
//      {
//        return;
//      }

      if( !RING_BONDS.IsEmpty())
      {
        bool contains_ring_conformations( IfContainsRingConformations( RING_BONDS, all_rotamer_bins));
        if( !contains_ring_conformations)
        {
          for( size_t count( 1), number_of_rotamers( rotamer_index); count != number_of_rotamers; ++count)
          {
            // get rotamer index string
            const std::string rotamer_st( "Rotamer" + util::Format()( count));

            molecule_output.RemoveProperty( rotamer_st + "Coordinates");
          }
        }
      }

      molecule_output.StoreProperty( "ScaffoldCount", storage::Vector< size_t>( 1, size_t( scaffold_count)));
      molecule_output.WriteMDL( m_Output);
    }

    //! @brief convert dihedral bond information into a string
    //! @param TYPE dihedral bond information
    template< unsigned int t_N>
    std::string BuildRotamerLibrary::ConvertAtomBondTypeIntoString( const storage::VectorND< t_N, size_t> &TYPE) const
    {
      std::stringstream fragment_information;
      fragment_information << TYPE( 0) << ' ';
      for( unsigned int i( 1); i < t_N; i += 2)
      {
        fragment_information << TYPE( i) << ' ';
        fragment_information << TYPE( i + 1) << ' ';
      }
      return fragment_information.str();
    }

    //! @brief convert dihedral information into molecule
    //! @TYPE dihedral bond information
    //! @return molecule from dihedral bond
    chemistry::FragmentComplete BuildRotamerLibrary::ConvertDihedralToMolecule( const storage::VectorND< 7, size_t> &TYPE) const
    {
      storage::Vector< sdf::AtomInfo> atom_info( size_t( 4));
      storage::Vector< sdf::BondInfo> bond_info( size_t( 3));

      // iterate over atoms
      for( size_t atom_id( 0), number_atoms( 4); atom_id < number_atoms; ++atom_id)
      {
        sdf::AtomInfo current_info( chemistry::AtomType( TYPE( atom_id * 2)), chemistry::e_NonChiral);

        // set the atom info up properly
        atom_info( atom_id) = current_info;

        if( atom_id < size_t( 3))
        {
          bond_info( atom_id) = sdf::BondInfo( atom_id, atom_id + 1, chemistry::ConfigurationalBondType( TYPE( atom_id * 2 + 1))->WithoutIsometry());
        }
      }
      chemistry::FragmentComplete molecule( chemistry::AtomVector< chemistry::AtomComplete>( atom_info, bond_info), "");
      return molecule;
    }

    namespace
    {
      //! @brief check whether the given vector nd refers to a non-gasteiger atom type
      //! @param TYPE vector containing alternating atom/bond types, starting with an atom type
      //! @return true if all the atom types are proper gasteiger types
      bool ContainsOnlyGasteigerAtomTypes( const storage::VectorND< 7, size_t> &TYPE)
      {
        // determine the last atom type that is a gasteiger type
        static size_t n_gasteiger_types( chemistry::AtomTypes::GetNumberGasteigerTypes());

        return
          TYPE( 0) < n_gasteiger_types
          && TYPE( 2) < n_gasteiger_types
          && TYPE( 4) < n_gasteiger_types
          && TYPE( 6) < n_gasteiger_types;
      }

    }

    //! @brief get dihedral bonds as fragments from the given molecule
    //! @param MOLECULE the small molecule from which dihedral information will be obtained
    void BuildRotamerLibrary::GetDihedralFragments( const chemistry::ConformationInterface &MOLECULE) const
    {
      storage::Map< storage::VectorND< 7, size_t>, storage::Vector< double> >
        dihedrals( MOLECULE.GetDihedralAnglesByType());

      // add the small molecule's dihedral angle info to the histograms
      std::stringstream vertex_edge_data;
      vertex_edge_data << size_t( 4) << ' ';
      vertex_edge_data << size_t( 3);
      storage::Vector< std::string> data_vector( util::SplitString( vertex_edge_data.str(), " "));
      const size_t h_atom_type( chemistry::GetAtomTypes().H_S.GetIndex());
      for
      (
        storage::Map< storage::VectorND< 7, size_t>, storage::Vector< double> >::const_iterator
          itr( dihedrals.Begin()), itr_end( dihedrals.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !ContainsOnlyGasteigerAtomTypes( itr->first))
        {
          continue;
        }
        if( chemistry::ConfigurationalBondType( itr->first( 3))->IsBondInRing())
        {
          continue;
        }
        if
        (
          itr->first( 0) != h_atom_type
          && itr->first( 6) != h_atom_type
          && !chemistry::ConfigurationalBondType( itr->first( 1))->IsBondInRing()
          && !chemistry::ConfigurationalBondType( itr->first( 5))->IsBondInRing()
        )
        {
          continue;
        }
        std::string angle_string( ConvertAtomBondTypeIntoString( itr->first));
        if( !m_DihedralMap.Has( angle_string))
        {
          m_DihedralMap[ angle_string] = ConvertDihedralToMolecule( itr->first);
        }
      }
    }

    //! @brief determine rings and bonds contained in the fragment
    void BuildRotamerLibrary::FindDihedralFragments() const
    {
      // read in ensemble, don't add any h because they slow down the fragment search and are unnecessary
      // If the user has added h to their molecules, then remove hydrogens
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_StructureDatabaseFlag->GetFirstParameter()->GetValue());
      chemistry::FragmentEnsemble molecules( input, sdf::e_Saturate);
      io::File::CloseClearFStream( input);
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( molecules.Begin()), itr_end( molecules.End());
        itr != itr_end;
        ++itr
      )
      {
        // add the small molecules' information to all the histograms
        // Skip molecules with bad 3d coordinates
        GetDihedralFragments( *itr);
      }
      storage::Vector< chemistry::FragmentComplete> dihedral_fragments( m_DihedralMap.GetMappedValues());

      m_FragmentSimpleGraphs.Reset();
      m_MoleculeSimpleGraphs.AllocateMemory( dihedral_fragments.GetSize());
      m_MoleculeSimpleGraphs.Reset();
      m_MoleculeSimpleGraphs.AllocateMemory( molecules.GetSize());
      m_EnsembleGraphs.Reset();
      m_EnsembleGraphs.AllocateMemory( molecules.GetSize());

      // create conformation graph maker object
      const chemistry::ConformationGraphConverter graph_maker
      (
        chemistry::ConformationGraphConverter::e_AtomType,
        chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
      );

      // store atom graphs for molecules
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( molecules.Begin()), itr_end( molecules.End());
        itr != itr_end;
        ++itr
      )
      {
        m_MoleculeSimpleGraphs.PushBack( graph_maker( *itr));

        // create graphs of molecules in structure db
        m_EnsembleGraphs.PushBack
        (
          chemistry::ConformationGraphConverter::CreateGraphWithAtoms( *itr)
        );
      }

      size_t scaffold_number( 0);
      // go through all molecules
      for
      (
        storage::Vector< chemistry::FragmentComplete>::const_iterator
          itr( dihedral_fragments.Begin()), itr_end( dihedral_fragments.End());
        itr != itr_end;
        ++itr, ++scaffold_number
      )
      {
        m_FragmentSimpleGraphs.PushBack( chemistry::ConstitutionGraphConverter()( chemistry::FragmentConstitutionShared( *itr)));
        BCL_MessageStd
        (
          "Molecule number : " + util::Format()( scaffold_number) + " / " + util::Format()( dihedral_fragments.GetSize())
          + " " + itr->GetAtomTypesString()
        );
        // call the main function that does all the work
        FindConformers( *itr, scaffold_number);
      }
      io::File::CloseClearFStream( m_Output);
    }

    namespace
    {
      storage::Vector< util::SiPtrVector< const linal::Vector3D> > GetHybridOrbitalCoords()
      {
        static storage::Vector< storage::Vector< linal::Vector3D> > ret( 4);
        storage::Vector< util::SiPtrVector< const linal::Vector3D> > refs( 4);
        for( size_t v( 1); v < size_t( 4); ++v)
        {
          chemistry::HybridOrbitalType hot( v);
          auto vec_vecs( chemistry::ValenceHandler::GetIdealizedGeometry( hot, 1.0));
          ret( v) = vec_vecs;
          refs( v).Resize( v + 1);
          for( size_t i( 0); i <= v; ++i)
          {
            refs( v)( i) = ret( v)( i);
          }
        }
        return refs;
      }

      void InsertMatrixIntoBondAngleMap
      (
        const linal::Matrix< double> &MAT,
        const linal::Vector< double> &LENGTHS,
        storage::List< storage::Pair< linal::Matrix< double>, math::RunningAverage< linal::Vector< double> > > > &LIST
      )
      {
        static const double s_tolerance( 0.025);
        double closest( 100000.0);
        auto itr_closest( LIST.Begin());
        int parity( 0);
        bool is_special_ring( MAT( 0, 0) < 1.0 && MAT.GetNumberRows() < size_t( 3));
        bool check_parity( MAT.GetNumberRows() > size_t( 2));
        int parity_row( 2), pos_start( 3);
        if( check_parity)
        {
          parity = MAT( 2, 2) > 0.0
                   ? 1
                   : MAT( 2, 2) < 0.0
                     ? -1
                     : 0;
        }
        else if( is_special_ring)
        {
          check_parity = true;
          parity = MAT( 0, 2) > 0.0
                   ? 1
                   : MAT( 0, 2) < 0.0
                     ? -1
                     : 0;
          parity_row = 0;
          pos_start = 0;
        }
        const double tolerance
        (
          math::Sqrt( MAT.GetNumberRows() - ( is_special_ring ? 0.0 : 1.0)) * s_tolerance
        );
        bool found( false);
        for( auto itr_list( LIST.Begin()), itr_list_end( LIST.End()); itr_list != itr_list_end; ++itr_list)
        {
          if( check_parity)
          {
            int list_parity( itr_list->First()( parity_row, 2) > 0.0 ? 1 : itr_list->First()( parity_row, 2) < 0.0 ? -1 : 0);
            if( list_parity != parity)
            {
              continue;
            }
          }
          double cosine( linal::Distance( itr_list->First().AsVector().Slice( pos_start), MAT.AsVector().Slice( pos_start)));
          if( cosine < tolerance)
          {
            closest = cosine;
            itr_closest = itr_list;
            found = true;
            break;
          }
        }
        if( found)
        {
          itr_closest->Second() += LENGTHS;

          double divisor( 1.0 / double( itr_closest->Second().GetWeight()));
          const double *itr_m( MAT.Begin());
          for( double *itr_a( itr_closest->First().Begin()), *itr_a_end( itr_closest->First().End()); itr_a != itr_a_end; ++itr_a, ++itr_m)
          {
            const double sa( *itr_a < 0.0 ? -math::Sqr( *itr_a) : math::Sqr( *itr_a));
            const double sm( *itr_m < 0.0 ? -math::Sqr( *itr_m) : math::Sqr( *itr_m));
            double d( sa + divisor * ( sm - sa));
            *itr_a = d < 0.0 ? -math::Sqrt( -d) : math::Sqrt( d);
          }
          //itr_closest->First() += ( MAT - itr_closest->First()) / double( itr_closest->Second());
        }
        else
        {
          LIST.PushBack();
          LIST.LastElement().First() = MAT;
          LIST.LastElement().Second() += LENGTHS;
        }
      }

      void FilterList( storage::List< storage::Pair< linal::Matrix< double>, math::RunningAverage< linal::Vector< double> > > > &LIST)
      {
        size_t largest( 0);
        for( auto itr_list( LIST.Begin()), itr_list_end( LIST.End()); itr_list != itr_list_end; ++itr_list)
        {
          largest = std::max( largest, size_t( itr_list->Second().GetWeight()));
        }
        if( largest > 120)
        {
          size_t cutoff( std::max( largest / size_t( 40), size_t( 3)));
          if( cutoff)
          {
            for( auto itr_list( LIST.Begin()), itr_list_end( LIST.End()); itr_list != itr_list_end;)
            {
              if( itr_list->Second().GetWeight() < cutoff)
              {
                auto itr_list_del( itr_list);
                ++itr_list;
                LIST.Remove( itr_list_del);
              }
              else
              {
                ++itr_list;
              }
            }
          }
        }
      }

    }

    //! @brief get dihedral bonds as fragments from the given molecule
    //! @param MOLECULE the small molecule from which dihedral information will be obtained
    void BuildRotamerLibrary::GetBondAngleFragments( const chemistry::FragmentComplete &MOLECULE) const
    {
      // get the canonical coordinates for the hybrid orbitals
      static storage::Vector< util::SiPtrVector< const linal::Vector3D> > s_canonical_coords( GetHybridOrbitalCoords());
      // get the origin
      static const linal::Vector3D origin( 0, 0, 0);

      if( MOLECULE.HasNonGasteigerAtomTypes())
      {
        return;
      }

      // Key <- atom-type <> Vector< bond type, atom type>
      // Value <- Matrix with unit-vector coordinates of all atoms after the first.
      // The first atom is always moved to 1 0 0, second atom is moved such that it is 0 in
      for( auto itr_atm( MOLECULE.GetAtomsIterator()); itr_atm.NotAtEnd(); ++itr_atm)
      {
        if
        (
          itr_atm->GetBonds().GetSize() <= size_t( 1)
          || itr_atm->GetBonds().GetSize() > size_t( 4)
          || !itr_atm->GetAtomType()->GetHybridOrbitalType().IsDefined()
          || itr_atm->GetAtomType()->GetHybridOrbitalType() == chemistry::GetHybridOrbitalTypes().e_Unhybridized
          || itr_atm->GetBonds().GetSize() == itr_atm->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
        )
        {
          continue;
        }
        for( size_t j( 0), three( 3); j < three; ++j)
        {
          chemistry::RotamerLibraryInterface::t_BondAngleMapKey key
          (
            chemistry::ConformationGraphConverter::AtomComparisonType( j),
            itr_atm->GetAtomType().GetIndex(),
            storage::Vector< storage::Pair< size_t, size_t> >( itr_atm->GetBonds().GetSize())
          );

          storage::Vector< storage::Triplet< size_t, size_t, size_t> > neighbors;
          for
          (
            auto itr_bond( itr_atm->GetBonds().Begin()), itr_bond_end( itr_atm->GetBonds().End());
            itr_bond != itr_bond_end;
            ++itr_bond
          )
          {
            neighbors.PushBack();
            neighbors.LastElement().First() =
              itr_bond->GetBondType()->GetBondData
              (
                chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
              );
            neighbors.LastElement().Second() =
              chemistry::ConformationGraphConverter::ConvertAtomTypeData
              (
                itr_bond->GetTargetAtom().GetAtomType(),
                key.First()
              );
            neighbors.LastElement().Third() = MOLECULE.GetAtomIndex( itr_bond->GetTargetAtom());
          }
          neighbors.Sort( std::less< storage::Triplet< size_t, size_t, size_t> >());
          neighbors = storage::Vector< storage::Triplet< size_t, size_t, size_t> >( neighbors.ReverseBegin(), neighbors.ReverseEnd());

          linal::Matrix< double> coords( neighbors.GetSize(), size_t( 3));

          for( size_t i( 0), nbond( neighbors.GetSize()); i < nbond; ++i)
          {
            key.Third()( i).First() = neighbors( i).First();
            key.Third()( i).Second() = neighbors( i).Second();
          }

          coords( 0, 0) = 1.0;
          coords( 0, 1) = coords( 0, 2) = coords( 1, 2) = 0.0;
          linal::Vector< double> bond_lengths( itr_atm->GetBonds().GetSize());
          for( size_t i( 0), nb( itr_atm->GetBonds().GetSize()); i < nb; ++i)
          {
            bond_lengths( i) = linal::Distance( itr_atm->GetPosition(), MOLECULE.GetAtomVector()( neighbors( i).Third()).GetPosition());
          }
          if( itr_atm->GetBonds().GetSize() == size_t( 2))
          {
            double radians
            (
              linal::ProjAngle
              (
                itr_atm->GetPosition(),
                MOLECULE.GetAtomVector()( neighbors( 0).Third()).GetPosition(),
                MOLECULE.GetAtomVector()( neighbors( 1).Third()).GetPosition()
              )
            );
            coords( 1, 0) = std::cos( radians);
            coords( 1, 1) = std::sin( radians);
            InsertMatrixIntoBondAngleMap( coords, bond_lengths, m_BondAngles[ key]);
            continue;
          }

          storage::Vector< linal::Vector3D> tmp_coords( neighbors.GetSize());
          for( size_t i( 0), nbond( neighbors.GetSize()); i < nbond; ++i)
          {
            tmp_coords( i) =
              linal::UnitVector
              (
                itr_atm->GetPosition(),
                MOLECULE.GetAtomVector()( neighbors( i).Third()).GetPosition()
              );
          }

          size_t n_in_ring
          (
            itr_atm->CountNonValenceBondsWithProperty
            (
              chemistry::ConfigurationalBondTypeData::e_IsInRing,
              size_t( 1)
            )
          );
          storage::Vector< size_t> order( storage::CreateIndexVector( 0, itr_atm->GetBonds().GetSize()));
          storage::Vector< storage::Pair< size_t, size_t> > keycopy;
          do
          {
            linal::Vector< double> tmp_bond_lengths( order.GetSize() - n_in_ring);
            if( n_in_ring && n_in_ring < order.GetSize())
            {
              storage::Vector< linal::Vector3D> identity( n_in_ring + 1, linal::Vector3D( 0.0));
              util::SiPtrVector< const linal::Vector3D> tmp_ref( n_in_ring + 1);
              for( size_t i( 0); i < n_in_ring; ++i)
              {
                identity( i)( i) = 1.0;
                tmp_ref( i) = util::ToSiPtr( tmp_coords( order( i)));
              }
              tmp_ref( n_in_ring) = util::ToSiPtr( origin);
              math::TransformationMatrix3D tf
              (
                quality::RMSD::SuperimposeCoordinates
                (
                  util::SiPtrVector< const linal::Vector3D>( identity.Begin(), identity.End()),
                  tmp_ref
                )
              );
              linal::Vector3D vo;
              for( size_t i( n_in_ring), nr( order.GetSize()); i < nr; ++i)
              {
                vo = tmp_coords( order( i));
                vo.Transform( tf);
                vo.Normalize();
                coords.GetRow( i - n_in_ring).CopyValues( vo);
                tmp_bond_lengths( i - n_in_ring) = bond_lengths( order( i));
              }
              coords.ShrinkRows( order.GetSize() - n_in_ring);
              if( math::Absolute( coords( 0, 2)) < 0.05 && itr_atm->GetAtomType()->GetNumberBonds() == size_t( 3))
              {
                coords( 0, 2) = 0.0;
                auto row2( coords.GetRow( 0));
                row2.Normalize();
              }
            }
            else
            {
              math::TransformationMatrix3D tf
              (
                coord::LineSegment3D( origin, linal::Vector3D( 1.0, 0.0, 0.0)),
                coord::LineSegment3D( origin, tmp_coords( order( 0)))
              );
              tf.SetTranslation( linal::Vector3D( 0.0, 0.0, 0.0));
              linal::Vector3D vo( tmp_coords( order( 1)));
              vo.Transform( tf);
              vo.Normalize();
              double radian( std::atan2( vo.Z(), vo.Y()));
              math::RotationMatrix3D rot( coord::GetAxes().e_X, radian);
              tf( rot);
              tmp_bond_lengths( 0) = bond_lengths( order( 0));
              for( size_t i( 1), nr( coords.GetNumberRows()); i < nr; ++i)
              {
                vo = tmp_coords( order( i));
                vo.Transform( tf);
                vo.Normalize();
                coords.GetRow( i).CopyValues( vo);
                tmp_bond_lengths( i) = bond_lengths( order( i));
              }
              coords( 1, 2) = 0.0;
              if( math::Absolute( coords( 2, 2)) < 0.05 && itr_atm->GetAtomType()->GetNumberBonds() == size_t( 3))
              {
                coords( 2, 2) = 0.0;
                auto row2( coords.GetRow( 2));
                row2.Normalize();
              }
            }
            auto keycopya( key.Third());
            keycopya.Reorder( order);
            if( keycopya == key.Third())
            {
              InsertMatrixIntoBondAngleMap( coords, tmp_bond_lengths, m_BondAngles[ key]);
            }
          } while( std::next_permutation( order.Begin(), order.End()));
        }
      }
    }

    namespace
    {
      std::string AtomComparisonTypeToString
      (
        const size_t &VALUE,
        const chemistry::ConformationGraphConverter::AtomComparisonType &TYPE
      )
      {
        if( TYPE == chemistry::ConformationGraphConverter::e_AtomType)
        {
          return chemistry::AtomType( VALUE).GetName();
        }
        else if( TYPE == chemistry::ConformationGraphConverter::e_ElementType)
        {
          return chemistry::ElementType( VALUE).GetName();
        }
        return "X";
      }
    }

    //! @brief determine rings and bonds contained in the fragment
    void BuildRotamerLibrary::FindBondAngleFragments() const
    {
      // read in ensemble, don't add any h because they slow down the fragment search and are unnecessary
      // If the user has added h to their molecules, then remove hydrogens
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_StructureDatabaseFlag->GetFirstParameter()->GetValue());
      chemistry::FragmentEnsemble molecules( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( molecules.Begin()), itr_end( molecules.End());
        itr != itr_end;
        ++itr
      )
      {
        // add the small molecules' information to all the histograms
        // Skip molecules with bad 3d coordinates
        GetBondAngleFragments( *itr);
      }

      for( auto itr( m_BondAngles.Begin()), itr_end( m_BondAngles.End()); itr != itr_end; ++itr)
      {
        FilterList( itr->second);
        std::ostringstream o;
        o << chemistry::AtomType( itr->first.Second()).GetName() << ' ';
        o << chemistry::ConformationGraphConverter::GetAtomComparisonType( itr->first.First()) << ' ';
        for( auto itrb( itr->first.Third().Begin()), itrb_end( itr->first.Third().End()); itrb != itrb_end; ++itrb)
        {
          o << itrb->First() << ' ' << AtomComparisonTypeToString( itrb->Second(), itr->first.First()) << ' ';
        }
        o << '\n' << itr->second << '\n';
        m_Output << o.str();
      }
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int BuildRotamerLibrary::Main() const
    {
      // get conformation tolerance value
      m_BinSize = m_BinSizeFlag->GetFirstParameter()->GetNumericalValue< double>();

      if( m_AddSimulatedCounts->GetFlag())
      {
        chemistry::RotamerLibraryFile empty_rotlib;
        m_SampleConformations = chemistry::SampleConformations
        (
          empty_rotlib,
          util::Implementation< chemistry::ConformationComparisonInterface>().GetAlias(),
          0.0,   // tolerance
          10000, // number of conformations
          12000, // number of iterations
          true,  // change chirality
          100000.0, // random dihedral change weight
          false,    // generate 3d
          0.1       // clash tolerance, auto-adjusted if necessary due to clashes
        );
      }

      // open the output file so that we can write rotamer library as we go
      io::File::MustOpenOFStream( m_Output, m_OutputFileFlag->GetFirstParameter()->GetValue());

      sdf::MdlHandler::GetAddAtomMdlLineFlag()->SetFlag();
      if( m_DihedralFragmentFlag->GetFlag())
      {
        FindDihedralFragments();
      }
      if( m_BondAngleFragmentFlag->GetFlag())
      {
        FindBondAngleFragments();
      }
      if( m_DihedralFragmentFlag->GetFlag() || m_BondAngleFragmentFlag->GetFlag())
      {
        io::File::CloseClearFStream( m_Output);
        return 0;
      }

      // read in ensemble, don't add any h because they slow down the fragment search and are unnecessary
      // If the user has added h to their molecules, then remove hydrogens
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_StructureDatabaseFlag->GetFirstParameter()->GetValue());
      chemistry::FragmentEnsemble molecules( input, sdf::e_Saturate);
      io::File::CloseClearFStream( input);

      // read in ensemble, don't add any h because they slow down the scaffold search and are unnecessary
      // If the user has added h to their fragments, then remove hydrogens
      // get a feed to get molecules that need to be fragmented
      chemistry::FragmentFeed feed( "fragments", sdf::e_Maintain);
      chemistry::FragmentEnsemble fragments;

      for( size_t mol_number( 0); feed.NotAtEnd(); ++feed, ++mol_number)
      {
        fragments.PushBack( *feed);
      }

      // For each small molecule, instantiate a graphical representation based on atom type and bond order/aromatic
      m_FragmentSimpleGraphs.Reset();
      m_MoleculeSimpleGraphs.Reset();
      m_EnsembleGraphs.Reset();
      m_EnsembleGraphs.AllocateMemory( molecules.GetSize());

      // create conformation graph maker object
      const chemistry::ConformationGraphConverter graph_maker
      (
        chemistry::ConformationGraphConverter::e_AtomTypeAndDistinguishHydrogens,
        chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
      );

      // create an object that splits molecules to give back only rings.
      const chemistry::FragmentSplitRings split_rings( true);

      // get unique ring id's. This is necessary to find isomorphisms between only whole rings and not parts. This
      // makes sure, for example, that no benzene rings are found in naphthalene.
      m_MapRings = util::ShPtr< chemistry::FragmentGraphMarker>( new chemistry::FragmentGraphMarker( graph_maker, split_rings));

      // update ring mapping class with rings in the whole structure db
      m_FragmentSimpleGraphs = m_MapRings->Insert( fragments);
      m_MoleculeSimpleGraphs = m_MapRings->operator ()( molecules);

      // store atom graphs for molecules
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( molecules.Begin()), itr_end( molecules.End());
        itr != itr_end;
        ++itr
      )
      {
        // create graphs of molecules in structure db
        m_EnsembleGraphs.PushBack
        (
          chemistry::ConformationGraphConverter::CreateGraphWithAtoms( *itr)
        );
      }

      // keep track of which scaffold we are on
      size_t scaffold_number( 0);

      // go through all molecules
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( fragments.Begin()), itr_end( fragments.End());
        itr != itr_end;
        ++itr, ++scaffold_number
      )
      {
        BCL_Message( util::Message::e_Standard, "Molecule number : " + util::Format()( scaffold_number));

        // call the main function that does all the work
        FindConformers( *itr, scaffold_number);
      }

      io::File::CloseClearFStream( m_Output);

      return 0;
    }

    const ApplicationType BuildRotamerLibrary::BuildRotamerLibrary_Instance
    (
      GetAppGroups().AddAppToGroup( new BuildRotamerLibrary(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
