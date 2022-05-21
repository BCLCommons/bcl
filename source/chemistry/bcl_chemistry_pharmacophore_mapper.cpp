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

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_collector_valence.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_molecule_feature_mapper.h"
#include "chemistry/bcl_chemistry_pharmacophore_mapper.h"
#include "descriptor/bcl_descriptor_base.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_linear_least_squares.h"
#include "math/bcl_math_statistics.h"

// external includes
#include <cmath>

namespace bcl
{
  namespace chemistry
  {

    //! @brief constructor
    //! @param SCORE_LABEL the score descriptor to use
    //! @param PROPERTY_LABEL the label of which properties to use
    //! @param ATOM_COMPARISON_TYPE how to compare the scaffold to the molecules (atom info)
    //! @param BOND_TYPE_INFO how to compare the scaffold to molecules (bond info)
    PharmacophoreMapper::PharmacophoreMapper
    (
//      const FragmentEnsemble &MOLECULES,
//      const FragmentComplete &SCAFFOLD,
      const util::ObjectDataLabel &SCORE_LABEL,
      const util::ObjectDataLabel &PROPERTY_LABEL,
      const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO
    ) :
      m_GraphMaker( ATOM_COMPARISON_TYPE, BOND_TYPE_INFO),
      m_ScoreDescriptor(),
      m_Properties(),
      m_ScaffoldData(),
      m_MolData()
    {
      // read the score descriptor
      m_ScoreDescriptor.AssertRead( SCORE_LABEL);

      // molecule-wide descriptor
      m_ScoreDescriptor.SetDimension( 0);

      SetProperties( PROPERTY_LABEL);
    }

    //! @brief set atom properties map; also filters out zero-sized descriptors (Define, ForEach, etc)
    //! @param PROPERTIES descriptor label to use
    void PharmacophoreMapper::SetProperties( const util::ObjectDataLabel &PROPERTIES)
    {

      m_Properties = descriptor::Combine< AtomConformationalInterface, float>();
      m_PropertyStrings = storage::Vector< std::string>();

      // Temporary descriptor
      descriptor::Combine< AtomConformationalInterface, float> read_check;

      // Properties must be molecule-wide; assert they are valid
      read_check.SetDimension( 0);
      read_check.AssertRead( PROPERTIES);

      // Fetch the aliases for each property
      for
      (
        descriptor::Combine< AtomConformationalInterface, float>::const_iterator itr_prop( read_check.Begin()),
          itr_prop_end( read_check.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        // Do not include 0-size descriptors, e.g. Define, ForEach, etc...
        if( ( *itr_prop)->GetSizeOfFeatures())
        {
          BCL_MessageStd( "Mapping " + itr_prop->GetString());
          m_Properties.PushBack( *itr_prop);
          m_PropertyStrings.PushBack( itr_prop->GetString());
        }
      }
    } // SetProperties

    //! @brief get the properties that are going to be mapped
    const descriptor::Combine< AtomConformationalInterface, float> &PharmacophoreMapper::GetProperties() const
    {
      return m_Properties;
    }

    // Get the human-readable of each property that is mapped
    const storage::Vector< std::string> &PharmacophoreMapper::GetPropertyStrings() const
    {
      return m_PropertyStrings;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Scores a set of molecules using a pharmacophore map
    //! @param SCAFFOLD the scaffold to use
    //! @param MOLECULES the molecules to score
    //! @param COEFF_MAP a list of property coefficients for changes at each grow point
    //!                  (calculated by CalculateMaps)
    //! @return a map from molecule indices (keys) to pairs (values) of compared molecule indices (First) and
    //!         values relevant to score calculation.  The first value will always be the changed grow point (scaffold atom)
    //!         and the last value is the activity change going from m1<->m2. Intermediate values are property changes
    storage::Map< size_t, storage::Vector< storage::Pair< size_t, linal::Vector< float> > > > PharmacophoreMapper::ScoreMolecules
    (
      const FragmentComplete &SCAFFOLD,
      const FragmentEnsemble &MOLECULES,
      const storage::Map< size_t, linal::Vector< float> > &COEFF_MAP
    ) const
    {
      return storage::Map< size_t, storage::Vector< storage::Pair< size_t, linal::Vector< float> > > >();
    }

    void PharmacophoreMapper::SetupMolData( const FragmentComplete &SCAFFOLD, const FragmentEnsemble &MOLECULES)
    {
      m_SinglePointChanges.Reset();
      m_MLSSolutions.Reset();

      std::vector< PharmMapMolData> &mol_data( m_MolData);
      mol_data.clear();

      // nominal number of molecules that will be used
      size_t nominal_mols( MOLECULES.GetSize());
      mol_data.reserve( nominal_mols);

      // Set up graph data including scaffold graph
      m_ScaffoldData = PharmMapMolData();
      m_ScaffoldData.m_Molecule = SCAFFOLD;
      m_ScaffoldData.m_MolGraph = m_GraphMaker( SCAFFOLD);

      m_ScaffoldData.m_ScaffoldIsomorphism = storage::Vector< size_t>( m_ScaffoldData.m_MolGraph.GetSize());
      for( size_t v( 0), end_v( m_ScaffoldData.m_ScaffoldIsomorphism.GetSize()); v < end_v; ++v)
      {
        m_ScaffoldData.m_ScaffoldIsomorphism( v) = v;
      }

      FragmentComplete scaff_w_h( SCAFFOLD);
      scaff_w_h.SaturateWithH();
      m_ScaffoldData.m_SummedProperties = m_Properties.SumOverObject( scaff_w_h);
      m_ScaffoldData.m_Score = m_ScoreDescriptor.SumOverObject( scaff_w_h)( 0);

      graph::ConstGraph< size_t, size_t> &scaff_graph( m_ScaffoldData.m_MolGraph);

      graph::SubgraphIsomorphism< size_t, size_t> iso_search;
      iso_search.SetSubgraphExternalOwnership( scaff_graph);

      storage::Set< size_t> changed_scaff_vertices;

      size_t mol_no( 0);
      for
      (
        FragmentEnsemble::const_iterator itr_mol( MOLECULES.Begin()), itr_mol_end( MOLECULES.End());
        itr_mol != itr_mol_end;
        ++itr_mol, ++mol_no
      )
      {
        // for property calculation the molecule must be saturated with hydrogens
        FragmentComplete new_mol_h( *itr_mol);
        new_mol_h.SaturateWithH();

        // Remove hydrogens; don't consider a hydrogen a substitution
        FragmentComplete new_mol( *itr_mol);
        new_mol.RemoveH();

        // molecule graph
        graph::ConstGraph< size_t, size_t> graph( m_GraphMaker( new_mol));

        iso_search.SetGraphExternalOwnership( graph);

        // Find all isomorphisms of the scaffold in the graph
        if( iso_search.FindAllIsomorphisms())
        {
          BCL_MessageVrb( "Molecule number " + util::Format()( mol_no) + " contains the scaffold");
          storage::Vector< storage::Vector< size_t> > isos( iso_search.GetIsomorphisms());

          // check that the scaffold doesn't differ by more than 6 atoms at any point
          {
            util::OwnPtr< const graph::ConstGraph< size_t, size_t> > graph_ptr( &graph, false);
            graph::ConstGraph< size_t, size_t> scaff_complement( graph::Subgraph< size_t, size_t>( graph_ptr, isos( 0)).GetComplement().ToGraph());
            storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( scaff_complement));
            bool keep_mol( true);
            for
            (
              storage::List< storage::Vector< size_t> >::const_iterator itr_comp( components.Begin()), itr_comp_end( components.End());
              itr_comp != itr_comp_end;
              ++itr_comp
            )
            {
              size_t atom_diff_max( 8);
              if( itr_comp->GetSize() > atom_diff_max)
              {
                BCL_MessageVrb( "Molecule contained more than " + util::Format()( atom_diff_max) + " extra atoms at a single point, excluding it from the comparisons");
                keep_mol = false;
                break;
              }
            }
            if( !keep_mol)
            {
              continue;
            }
          }

          size_t extra_points( scaff_graph.GetSize());
          storage::Vector< size_t> picked_iso, picked_substitution;

          // pick the isomorphism that maximizes the overlap with already-seen changed points
          if( changed_scaff_vertices.IsEmpty())
          {
            picked_substitution = GetScaffoldSubstitution( graph, isos( 0));
            extra_points = picked_substitution.GetSize();
            picked_iso = isos( 0);
          }
          else
          {
            for( size_t iso_num( 0), end_num( isos.GetSize()); iso_num < end_num; ++iso_num)
            {
              storage::Vector< size_t> &this_iso( isos( iso_num));
              storage::Vector< size_t> substitution( GetScaffoldSubstitution( graph, this_iso));

              size_t extra_points_this( 0);

              // check how many substituted points don't overlap with seen vertices
              for( size_t v( 0), end_v( substitution.GetSize()); v < end_v; ++v)
              {
                if( !changed_scaff_vertices.Contains( substitution( v)))
                {
                  ++extra_points_this;
                }
              }

              // If we are better than the previous isomorphisms, set the indicators
              if( extra_points_this < extra_points)
              {
                extra_points = extra_points_this;
                picked_iso = this_iso;
                picked_substitution = substitution;
              }
            }
          }

          // Add the substituted points to the changed point set
          for( size_t v( 0), end_v( picked_substitution.GetSize()); v < end_v; ++v)
          {
            changed_scaff_vertices.Insert( picked_substitution( v));
          }

          mol_data.push_back( PharmMapMolData());
          size_t idx( mol_data.size() - 1);

          mol_data[idx].m_Molecule = new_mol;
          mol_data[ idx].m_MolGraph = graph;
          mol_data[ idx].m_ScaffoldIsomorphism = picked_iso;
          mol_data[ idx].m_Substitution = picked_substitution;
          mol_data[ idx].m_SummedProperties = m_Properties.SumOverObject( new_mol_h);
          mol_data[ idx].m_Score = m_ScoreDescriptor.SumOverObject( new_mol_h)( 0);
        }
        else
        {
          BCL_MessageVrb( "Molecule number " + util::Format()( mol_no) + " did not contain the given scaffold, excluding it from consideration");
        }
      }
      storage::Vector< storage::Triplet< size_t, size_t, size_t> > changed_pairs( GetChangedPairs( mol_data));

      // a map of which entries in changed_pairs correspond to which scaffold points
      storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > > &associated_pairs( m_SinglePointChanges);

      // Convert the vector of triplets to a map of pairs
      for( size_t i( 0), end_i( changed_pairs.GetSize()); i < end_i; ++i)
      {
        BCL_MessageVrb
        (
          " Molecules " + util::Format()( changed_pairs( i).First()) + ", " + util::Format()( changed_pairs( i).Second())
          + " changed at " + util::Format()( changed_pairs( i).Third())
        );

        size_t &scaff_pt( changed_pairs( i).Third());
        associated_pairs[ scaff_pt].Append( storage::Pair< size_t, size_t>( changed_pairs( i).First(), changed_pairs( i).Second()));
      }
    }

    void PharmacophoreMapper::CalculateMaps
    (
      const FragmentComplete &SCAFFOLD,
      const FragmentEnsemble &INPUT_MOLECULES,
      const size_t &USE // TODO remove this
    )
    {
      SetupMolData( SCAFFOLD, INPUT_MOLECULES);
      std::vector< PharmMapMolData> &mol_data( m_MolData);
      storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > > &associated_pairs( m_SinglePointChanges);

      for
      (
        storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > >::const_iterator itr_map( associated_pairs.Begin()),
          itr_map_end( associated_pairs.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        util::ShPtr< storage::Pair< linal::Matrix< float>, linal::Vector< float> > > sp_prop_diff_matrix
        (
          CalculatePropertyDifferences( mol_data, itr_map->second)
        );

        if( !sp_prop_diff_matrix.IsDefined())
        {
          continue;
        }

        linal::Matrix< float> &prop_diff_mat( sp_prop_diff_matrix->First());
        linal::Vector< float> &score_diff_v( sp_prop_diff_matrix->Second());

        if( prop_diff_mat.GetNumberRows() <= 2 * m_PropertyStrings.GetSize())
        {
          BCL_MessageStd( "Scaffold point " + util::Format()( itr_map->first) + " is substituted, but only has " + util::Format()( prop_diff_mat.GetNumberRows() / 2) + " data points. cannot perform regression");
          continue;
        }

        BCL_MessageStd
        (
          "Using scaffold atom number " + util::Format()( itr_map->first) + " as a substitution point; linear regression using "
          + util::Format()( prop_diff_mat.GetNumberRows()) + " data points"
        );

        storage::Pair< linal::Vector< float>, float> soln_chi_sq
        (
          math::LinearLeastSquares::SolutionAndChiSquared( prop_diff_mat, score_diff_v)
        );
        //remove( "/tmp/matrix.csv");

        storage::Vector< std::string> coeff_labels( m_PropertyStrings);

        BCL_MessageStd( "  Coefficients for grow point " + util::Format()( itr_map->first));
        for( size_t i( 0); i < coeff_labels.GetSize(); ++i)
        {
          BCL_MessageStd( "    " + coeff_labels( i) + ": " + util::Format()( soln_chi_sq.First()( i)));
        }
        BCL_MessageStd( "  Chi-squared: " + util::Format()( soln_chi_sq.Second()));

        m_MLSSolutions[ itr_map->first] = soln_chi_sq;
      }
    }

    util::ShPtr< storage::Pair< linal::Matrix< float>, linal::Vector< float> > > PharmacophoreMapper::CalculatePropertyDifferences
    (
      const std::vector< PharmacophoreMapper::PharmMapMolData> &DATA,
      const storage::Vector< storage::Pair< size_t, size_t> > &CHANGED_PAIRS
    ) const
    {
      // matrix consisting of 2*(changed pairs) rows (for symmetry), and #cols = number of properties
      util::ShPtr< storage::Pair< linal::Matrix< float>, linal::Vector< float> > > sp_all_prop_diffs
      (
        new storage::Pair< linal::Matrix< float>, linal::Vector< float> >
        (
          linal::Matrix< float>( 2*CHANGED_PAIRS.GetSize(), m_PropertyStrings.GetSize()),
          linal::Vector< float>( 2*CHANGED_PAIRS.GetSize())
        )
      );

      linal::Matrix< float> &all_prop_diffs( sp_all_prop_diffs->First());
      linal::Vector< float> &all_score_diffs( sp_all_prop_diffs->Second());

      for( size_t r( 0), end_r( CHANGED_PAIRS.GetSize()); r < end_r; ++r)
      {
        size_t mol_1( CHANGED_PAIRS( r).First());
        size_t mol_2( CHANGED_PAIRS( r).Second());

        linal::Vector< float> pair_prop_diffs( DATA[ mol_1].m_SummedProperties - DATA[ mol_2].m_SummedProperties);

        // Copy property differences to the matrix
        std::copy( pair_prop_diffs.Begin(), pair_prop_diffs.End(), all_prop_diffs[ ( 2 * r)]);
        all_score_diffs( 2 * r) = DATA[ mol_1].m_Score - DATA[ mol_2].m_Score;

        pair_prop_diffs *= -1;

        // Copy the negative of the results to the matrix, for symmetry
        std::copy( pair_prop_diffs.Begin(), pair_prop_diffs.End(), all_prop_diffs[ ( 2 * r) + 1]);
        all_score_diffs( 2 * r + 1) = DATA[ mol_2].m_Score - DATA[ mol_1].m_Score;
      }
      return sp_all_prop_diffs;
    }

    //! @brief calculates pharmacophore maps given a scaffold and a set of molecules
    //! @param SCAFFOLD the scaffold to use
    //! @param INPUT_MOLECULES the nominal set of molecules to use
    //! @param REPLACE_ZERO_PROPERTIES if true, this will replace zero-difference properties with constant values
    //!                                this is useful to do an MLS regression when some of the properties are nonsensical for
    //!                                certain grow points
    void PharmacophoreMapper::CalculateMaps
    (
      const FragmentComplete &SCAFFOLD,
      const FragmentEnsemble &INPUT_MOLECULES,
      const bool &REPLACE_ZERO_PROPERTIES,
      const float &SCORE_TOLERANCE,
      const bool &BINARY
    )
    {
    }

    //!@brief find common scaffolds for molecules
    //! @param MOLECULES the molecules to look through
    //! @param SAMPLING the percentage of pairwise molecule comparisons that should be done
    //! @return a mapping from scaffolds to which molecules match that scaffold
    storage::Vector< storage::Pair< graph::ConstGraph< size_t, size_t>, storage::Vector< size_t> > >
      PharmacophoreMapper::FindScaffolds
    (
      const FragmentEnsemble &MOLECULES,
      const float &SAMPLING,
      const size_t &MIN_SCAFF_SIZE,
      const bool &PRINT_STATUS
    ) const
    {
      size_t n_mols( MOLECULES.GetSize());
      size_t comparison_stride( std::min( std::max( size_t( float( 1.0) / SAMPLING), size_t( 1)), n_mols - 1));
      size_t n_comparisons( n_mols / comparison_stride);

      BCL_MessageStd( "Comparing every " + util::Format()( comparison_stride) + " molecules");

      // determine which pairs of molecules to compare
      storage::Vector< storage::Pair< size_t, size_t> > pairs;
      pairs.AllocateMemory( n_comparisons);
      for( size_t i( 0); i < n_mols; ++i)
      {
        if( ( i % comparison_stride) != 0)
        {
          continue;
        }
        for( size_t j( i + 1); j < n_mols; ++j)
        {
          if( ( j % comparison_stride) == 0)
          {
            pairs.PushBack( storage::Pair< size_t, size_t>( i, j));
          }
        }
      }

      storage::Vector< graph::ConstGraph< size_t, size_t> > graphs;
      graphs.AllocateMemory( n_mols);

      // make graphs of the molecules
      for
      (
        FragmentEnsemble::const_iterator itr_mol( MOLECULES.Begin()), itr_mol_end( MOLECULES.End());
        itr_mol != itr_mol_end;
        ++itr_mol
      )
      {
        graphs.PushBack( m_GraphMaker( *itr_mol));
      }

      storage::Vector< storage::Pair< graph::ConstGraph< size_t, size_t>, storage::Vector< size_t> > > scaffolds;
      scaffolds.AllocateMemory( n_comparisons);

      // Find common subgraphs between pairs
      graph::CommonSubgraphIsomorphism< size_t, size_t> csi;
      graph::SubgraphIsomorphism< size_t, size_t> iso;
      size_t print_every( pairs.GetSize() * 0.001);
      print_every = std::max( print_every, size_t( 1));
      for( size_t pair_no( 0), end_no( pairs.GetSize()); pair_no < end_no; ++pair_no)
      {
        if( ( pair_no % print_every) == 0)
        {
          util::GetLogger().LogStatus( "Progress: " + util::Format().FFP( 2)( float( pair_no) / float( end_no)) + " done");
        }
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > graph_one( &graphs( pairs( pair_no).First()), false);
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > graph_two( &graphs( pairs( pair_no).Second()), false);

        csi.SetGraphA( graph_one->GetSize() < graph_two->GetSize() ? graph_one : graph_two);
        csi.SetGraphB( graph_one->GetSize() < graph_two->GetSize() ? graph_two : graph_one);

        csi.FindIsomorphism( csi.EstimateUpperBounds(), MIN_SCAFF_SIZE);

        storage::Vector< graph::Subgraph< size_t, size_t> > res( csi.GetSubgraphIsomorphismsOfGraphA());

        if( !res.IsEmpty() && res( 0).GetSize() > 0)
        {
          graph::ConstGraph< size_t, size_t> scaff_graph( res( 0).ToGraph());

          iso.SetGraphExternalOwnership( scaff_graph);
          bool keep( true);
          for( size_t i( 0), end_i( scaffolds.GetSize()); i < end_i; ++i)
          {
            if( scaffolds( i).First().GetSize() != scaff_graph.GetSize())
            {
              continue;
            }
            iso.SetSubgraphExternalOwnership( scaffolds( i).First());
            if( iso.FindIsomorphism())
            {
              keep = false;
              break;
            }
          }

          if( keep)
          {
            scaffolds.PushBack
            (
              storage::Pair< graph::ConstGraph< size_t, size_t>, storage::Vector< size_t> >
              (
                scaff_graph, storage::Vector< size_t>()
              )
            );
          }
        }
      }

      // Determine which molecules match each scaffold
      print_every = n_mols * 0.001;
      for( size_t scaff( 0), end_no( scaffolds.GetSize()); scaff < end_no; ++scaff)
      {
        iso.SetSubgraphExternalOwnership( scaffolds( scaff).First());
        for( size_t mol_no( 0), end_mol( graphs.GetSize()); mol_no < end_mol; ++mol_no)
        {
          iso.SetGraphExternalOwnership( graphs( mol_no));
          if( iso.FindIsomorphism())
          {
            scaffolds( scaff).Second().PushBack( mol_no);
          }
        }
      }
      return scaffolds;

    }

    //! @brief calculates pharmacophore maps by iteratively removing atoms and determining the QSAR score
    void PharmacophoreMapper::MolecularFeatureMap
    (
      const FragmentComplete &SCAFFOLD,
      const FragmentEnsemble &MOLECULES,
      const util::ObjectDataLabel &MODEL,
      const std::string &OUT_PREFIX
    ) const
    {

      descriptor::CheminfoProperty model( MODEL);

      std::string scaff_full_filename( io::File::MakeAbsolutePath( OUT_PREFIX + "_scaffold.sdf"));
      std::string features_full_filename( io::File::MakeAbsolutePath( OUT_PREFIX + "_features.pml"));

      MoleculeFeatureMapper mapper( MODEL);
      std::vector< std::map< size_t, float> > avg_features( mapper.AverageFeatures( SCAFFOLD, MOLECULES));

      if( avg_features.empty())
      {
        return;
      }

      io::OFStream pml;
      io::File::MustOpenOFStream( pml, scaff_full_filename);
      SCAFFOLD.WriteMDL( pml);
      io::File::CloseClearFStream( pml);

      io::File::MustOpenOFStream( pml, features_full_filename);
      pml << "load " << scaff_full_filename << ", molecule" << std::endl;
      pml << "alter molecule and state 1, b=0" << std::endl;
      pml << "alter molecule and state 1, vdw=0" << std::endl;
      pml << "set label_size, -0.4" << std::endl;
      pml << "hide lines " << std::endl;
      pml << "show sticks" << std::endl;
      pml << "set stick_quality, 20" << std::endl;
      pml << "set valence, 1" << std::endl;
//      pml << "spectrum b, blue_white_red, minimum=-1, maximum=1" << std::endl;

      float max( 0);
      float min( 0);

      std::map< size_t, float> &avg_map( avg_features[ 0]);
      for
      (
        std::map< size_t, float>::const_iterator itr_map( avg_map.begin()), itr_map_end( avg_map.end());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        pml << "alter \"molecule\" and state 1 and id " << itr_map->first + 1
            << ", b=" << itr_map->second << std::endl;
        if( itr_map->second < min)
        {
          min = itr_map->second;
        }
        if( itr_map->second > max)
        {
          max = itr_map->second;
        }
      }
      float range( std::max( std::fabs( max), std::fabs( min)));
      pml << "spectrum b, blue_white_red, selection=molecule and state 1, ";
      pml << "minimum=" << ( -1 * range) << ", maximum=" << range << std::endl;

      io::File::CloseClearFStream( pml);
    }

    storage::Vector< math::RunningAverageSD< float> > PharmacophoreMapper::CalculatePropertyDiffStatistics
    (
      const std::vector< PharmacophoreMapper::PharmMapMolData> &DATA,
      const storage::Vector< storage::Pair< size_t, size_t> > &CHANGED_PAIRS
    ) const
    {
      storage::Vector< math::RunningAverageSD< float> > stat_v( m_PropertyStrings.GetSize());
      for( size_t p( 0), end_p( CHANGED_PAIRS.GetSize()); p < end_p; ++p)
      {
        size_t m1( CHANGED_PAIRS( p).First());
        size_t m2( CHANGED_PAIRS( p).Second());

        linal::Vector< float> prop_diffs( DATA[ m1].m_SummedProperties - DATA[ m2].m_SummedProperties);

        for( size_t i( 0), end_i( prop_diffs.GetSize()); i < end_i; ++i)
        {
          stat_v( i) += std::fabs( prop_diffs( i));
          stat_v( i) += -std::fabs( prop_diffs( i));
        }
      }
      return stat_v;
    }

    math::RunningAverageSD< float> PharmacophoreMapper::CalculateScoreDiffStatistics
    (
      const std::vector< PharmacophoreMapper::PharmMapMolData> &DATA,
      const storage::Vector< storage::Pair< size_t, size_t> > &CHANGED_PAIRS
    ) const
    {
      math::RunningAverageSD< float> stat_v;
      for( size_t p( 0), end_p( CHANGED_PAIRS.GetSize()); p < end_p; ++p)
      {
        size_t m1( CHANGED_PAIRS( p).First());
        size_t m2( CHANGED_PAIRS( p).Second());

        float score_diff( DATA[ m1].m_Score - DATA[ m2].m_Score);

        stat_v += std::fabs( score_diff);
        stat_v += -std::fabs( score_diff);
      }
      return stat_v;
    }

    //! @brief write pharmacophore maps as calcualted by CalculateMaps to a PML file
    //! @param SCAFFOLD filename of the SDF containing the scaffold
    //! @param OUTPUT_FILENAME the PML file to write to
    void PharmacophoreMapper::WriteMapsPML( const std::string &SCAFFOLD_FILENAME, const std::string &OUTPUT_FILENAME) const
    {
      if( m_MLSSolutions.IsEmpty())
      {
        // Nothing to do
        BCL_MessageStd( "Least-square regressions have not been performed; there is nothing to write to a PyMol script");
        return;
      }

      storage::Map< size_t, linal::Vector< float> > scaling_factors;

      for
      (
        storage::Map< size_t, storage::Pair< linal::Vector< float>, float> >::const_iterator itr_mls( m_MLSSolutions.Begin()),
          itr_mls_end( m_MLSSolutions.End());
        itr_mls != itr_mls_end;
        ++itr_mls
      )
      {
        size_t changed_point( itr_mls->first);
        const linal::Vector< float> &mls_coeffs( itr_mls->second.First());
        linal::Vector< float> abs_mls_coeffs( mls_coeffs);
        std::transform( abs_mls_coeffs.Begin(), abs_mls_coeffs.End(), abs_mls_coeffs.Begin(), static_cast< float( *)( float)>( &std::fabs));

        const storage::Vector< storage::Pair< size_t, size_t> > &pairs( m_SinglePointChanges.GetValue( changed_point));
        storage::Vector< math::RunningAverageSD< float> > prop_stats( CalculatePropertyDiffStatistics( m_MolData, pairs));

        linal::Vector< float> abs_prop_diffs( prop_stats.GetSize(), float( 0));
        for( size_t i( 0), end_i( abs_prop_diffs.GetSize()); i < end_i; ++i)
        {
          abs_prop_diffs = prop_stats( i).GetStandardDeviation();
        }

        float max_activity_change( abs_mls_coeffs * abs_prop_diffs);
        scaling_factors[ changed_point] = linal::Vector< float>( abs_prop_diffs.GetSize(), float( 0));

        linal::Vector< float> &v( scaling_factors[ changed_point]);

        for( size_t i( 0), end_i( abs_prop_diffs.GetSize()); i < end_i; ++i)
        {
          v( i) = abs_mls_coeffs( i) * abs_prop_diffs( i) / max_activity_change;
        }

        //math::RunningAverageSD< float> score_stats( GetScoreDiffStats( m_MolData, pairs));

      }

      //
      // Write the PML script
      //

      io::OFStream out_pml;
      io::File::MustOpenOFStream( out_pml, OUTPUT_FILENAME);

      // full filename of the scaffold so that this pml can be loaded from anywhere
      std::string full_filename( io::File::MakeAbsolutePath( SCAFFOLD_FILENAME));
      BCL_MessageStd( "Scaffold file: " + full_filename);

      // header part of the pymol script
      out_pml << "#loading mdl/sdf file " << std::endl;
      out_pml << "load " << full_filename << ", molecule " << std::endl;
      out_pml << "alter molecule and state 1, b=0" << std::endl;
      out_pml << "alter molecule and state 1, vdw=0" << std::endl;
      out_pml << "set label_size, -0.4" << std::endl;
      out_pml << "show lines " << std::endl;
      out_pml << "hide sticks" << std::endl;
      out_pml << "set valence, 1" << std::endl;

      // Add PharmMaps for each property
      for( size_t prop_no( 0); prop_no < m_PropertyStrings.GetSize(); ++prop_no)
      {

        std::string prop_name( m_PropertyStrings( prop_no));

        out_pml << "copy " << prop_name << ", molecule" << std::endl;
        //out_pml << "alter \"" << prop_name << "\", vdw=1" << std::endl;
        out_pml << "alter \"" << prop_name << "\", vdw=0" << std::endl;
        out_pml << "color grey, \"" << prop_name << "\"" << std::endl;
        out_pml << "set surface_type, 0" << std::endl;
        out_pml << "set transparency, 0.5" << std::endl;
        //out_pml << "show surface, \"" << prop_name << "\"" << std::endl;
        out_pml << "hide surface, \"" << prop_name << "\"" << std::endl;
        out_pml << "show spheres, \"" << prop_name << "\"" << std::endl;
        out_pml << "hide lines, \"" << prop_name << "\"" << std::endl;

        // Get the coefficient of this property at each grow point
        for
        (
          storage::Map< size_t, storage::Pair< linal::Vector< float>, float> >::const_iterator
            itr_soln( m_MLSSolutions.Begin()), itr_soln_end( m_MLSSolutions.End());
          itr_soln != itr_soln_end;
          ++itr_soln
        )
        {

          const size_t growpt( itr_soln->first);
          size_t atom_no( growpt + 1);

          const storage::Pair< linal::Vector< float>, float> &mls( itr_soln->second);

          float prop_mls_coeff( mls.First()( prop_no));
          linal::Vector< float> &scaling_coeffs( scaling_factors[ growpt]);

          if( prop_mls_coeff != 0)
          {

            //float prop_impact_percent( max_activity_single_prop.GetValue( growpt)( prop_no) / max_activity_all_props.GetValue( growpt));
            float prop_impact_percent( scaling_coeffs( prop_no));
            //float sphere_size( 2 * std::fabs( max_activity_single_prop.GetValue( growpt)( prop_no) / max_activity_all_props.GetValue( growpt)));
            float sphere_size( 2 * scaling_coeffs( prop_no));

            // How impactful changing a property by 1 stddev will be in going from a beneficial to a deleterious effect
            float color_factor( scaling_coeffs( prop_no));

            float red_color( 0.5 + std::min< float>( ( prop_mls_coeff < 0 ? color_factor : 0), 0.5));
            float blue_color( 0.5 + std::min< float>( ( prop_mls_coeff > 0 ? color_factor : 0), 0.5));

            std::string sphere_color( "[" + util::Format()( red_color) + ",0.5," + util::Format()( blue_color) + "]");

            // label
            out_pml << "label \"" << prop_name << "\" and state 1 and id " << atom_no << ", \"" << ( prop_mls_coeff > 0 ? "Increase" : "Decrease") << ", ";
            out_pml << "Impact: " << util::Format().FFP( 2)( std::fabs( prop_impact_percent) * 100) << " %\"" << std::endl;
            out_pml << "alter \"" << prop_name << "\" and state 1 and id " << atom_no << ", vdw=" << sphere_size << std::endl;

            std::string color_name( "colorProp" + util::Format()( prop_no) + "GrowPt" + util::Format()( atom_no));

            out_pml << "set_color " << color_name << ", " << sphere_color << std::endl;
            out_pml << "color " << color_name << ", \"" << prop_name << "\" and id " << atom_no << std::endl;
          }
        }
      }
      io::File::CloseClearFStream( out_pml);
    }

    //! @param DETAILS_PREFIX the prefix to use for files
    void PharmacophoreMapper::WriteDetailsMols( const std::string &DETAILS_PREFIX) const
    {
      std::string out_filename( DETAILS_PREFIX + ".mols.sdf.gz");
      io::OFStream out( out_filename.c_str(), std::ios::out);
      for( size_t i( 0), end_i( m_MolData.size()); i < end_i; ++i)
      {
        m_MolData[ i].m_Molecule.WriteMDL( out);
      }
      io::File::CloseClearFStream( out);
    }

    //! @brief write details of the calculations
    //! @param DETAILS_PREFIX the prefix to use for files
    void PharmacophoreMapper::WriteDetailsCSV( const std::string &DETAILS_PREFIX) const
    {
      std::string out_filename( DETAILS_PREFIX + ".mol_data.csv");
      if( m_MolData.empty())
      {
        return;
      }

      io::OFStream out( out_filename.c_str(), std::ios::out);

      out << "changed_sub,mol_index_1,mol_index_2,";
      for( size_t i( 0); i < m_PropertyStrings.GetSize(); ++i)
      {
        out << "\"" << m_PropertyStrings( i) << "\",";
      }
      out << "pharmmap_score_diff,actual_score_diff\n";

      for
      (
        storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > >::const_iterator
          itr_map( m_SinglePointChanges.Begin()), itr_map_end( m_SinglePointChanges.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        const size_t change_point( itr_map->first);
        if( !m_MLSSolutions.Has( change_point))
        {
          continue;
        }
        const storage::Vector< storage::Pair< size_t, size_t> > &pairs( itr_map->second);

        linal::Vector< float> mls_coeffs( m_MLSSolutions.GetValue( change_point).First());

        for( size_t p( 0), end_p( pairs.GetSize()); p < end_p; ++p)
        {
          size_t m1( pairs( p).First());
          size_t m2( pairs( p).Second());
          linal::Vector< float> prop_diff( m_MolData[ m1].m_SummedProperties - m_MolData[ m2].m_SummedProperties);
          float score_diff( m_MolData[ m1].m_Score - m_MolData[ m2].m_Score);
          float res( prop_diff * mls_coeffs);

          out << change_point << "," << m1 << "," << m2 << ",";
          for( size_t i( 0), end_i( prop_diff.GetSize()); i < end_i; ++i)
          {
            out << prop_diff( i) << ",";
          }
          out << res << "," << score_diff << std::endl;

          out << change_point << "," << m2 << "," << m1 << ",";
          for( size_t i( 0), end_i( prop_diff.GetSize()); i < end_i; ++i)
          {
            out << -prop_diff( i) << ",";
          }
          out << -res << "," << -score_diff << std::endl;
        }
      }
      io::File::CloseClearFStream( out);
    }

    //! @brief collects fragments from a substructure given a substructure
    //! @param MOL_GRAPH graph of the molecule
    //! @param SUBSTRUCT vector of indices representing a substructure to match
    //! @return a map from vertices (keys) to the fragment associated with them
    storage::Map< size_t, storage::Vector< size_t> > PharmacophoreMapper::CollectFragmentsFromSubstructure
    (
      const graph::ConstGraph< size_t, size_t> &MOL_GRAPH,
      const storage::Vector< size_t> &SUBSTRUCTURE
    ) const
    {
      graph::ConstGraph< size_t, size_t> mol_graph( MOL_GRAPH);
      util::OwnPtr< const graph::ConstGraph< size_t, size_t> > graph_ptr( &mol_graph, false);
      graph::Subgraph< size_t, size_t> mol_subgraph( graph_ptr, SUBSTRUCTURE);

      // which atoms should be kept
      storage::Vector< size_t> keep_indices( mol_graph.GetSize(), size_t( 1));

      // which non-scaffold atoms (keys) connect to which scaffold atoms (vectors)
      storage::Map< size_t, storage::Vector< size_t> > to_scaffold;

      // lookup array of which atoms to keep
      for( size_t i( 0); i < SUBSTRUCTURE.GetSize(); ++i)
      {
        keep_indices( SUBSTRUCTURE( i)) = 0;
      }

      // isolate scaffold atoms
      for( size_t i( 0); i < SUBSTRUCTURE.GetSize(); ++i)
      {
        const size_t &mol_atom( SUBSTRUCTURE( i));
        storage::Vector< size_t> neighbors( mol_graph.GetNeighborIndices( mol_atom));
        for( size_t ni( 0); ni < neighbors.GetSize(); ++ni)
        {
          const size_t &neighbor_atom( neighbors( ni));
          mol_graph.RemoveEdge( mol_atom, neighbor_atom);

          // keep track of atom neighbors
          if( keep_indices( neighbor_atom))
          {
            to_scaffold[ neighbor_atom].PushBack( i);
          }
        }
      }

      // which grow points (keys) are associated with which non-scaffold indices (values)
      // note: should be unique
      storage::Map< size_t, storage::Set< size_t> > growpt_fragments;

      storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( mol_graph));

      // Associate indices with each grow point
      for
      (
        storage::Map< size_t, storage::Vector< size_t> >::const_iterator itr_map( to_scaffold.Begin()),
          itr_map_end( to_scaffold.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        const size_t &atom_no( itr_map->first);

        // iterate through each component
        storage::List< storage::Vector< size_t> >::iterator itr_comp( components.Begin()), itr_comp_end( components.End());
        while( itr_comp != itr_comp_end)
        {
          // if true the atom is part of the scaffold.  Remove it and go to the next component
          if( itr_comp->GetSize() == 1 && !( keep_indices( ( *itr_comp)( 0))))
          {
            itr_comp = components.Remove( itr_comp);
            continue;
          }

          // if this atom is in this component, associate all component indices with the grow point
          if( itr_comp->Find( atom_no) < itr_comp->GetSize())
          {
            const storage::Vector< size_t> &associated_growpts( itr_map->second);
            for( size_t i( 0); i < associated_growpts.GetSize(); ++i)
            {
              const size_t &growpt( associated_growpts( i));
              growpt_fragments[ growpt].InsertElements( itr_comp->Begin(), itr_comp->End());
            }
          }
          ++itr_comp;
        }
      }

      //
      // Convert the sets to vectors
      //

      storage::Map< size_t, storage::Vector< size_t> > map;

      for
      (
        storage::Map< size_t, storage::Set< size_t> >::const_iterator itr_map( growpt_fragments.Begin()), itr_map_end( growpt_fragments.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        map[ itr_map->first] = storage::Vector< size_t>( itr_map->second.Begin(), itr_map->second.End());
      }

      return map;
    }

    //! @brief gets pairs of molecules that differ at a single point of a common scaffold
    //! @return a vector of triplets; first=mol 1, second=mol2, third=change point
    storage::Vector< storage::Triplet< size_t, size_t, size_t> > PharmacophoreMapper::GetChangedPairs( const std::vector< PharmMapMolData> &DATA)
    {
      storage::Vector< storage::Triplet< size_t, size_t, size_t> > changed_points;

      size_t n_data( DATA.size());

      // Every molecule matches one other is probably a decent approximation
      changed_points.AllocateMemory( 2 * n_data);

      for( size_t i( 0); i < n_data; ++i)
      {
        for( size_t j( i + 1); j < n_data; ++j)
        {

          // if one of the two molecules is not substituted at all, and the other molecule is only substituted at a single point, then
          // there is also only one point of change
          if
          (
            ( DATA[ i].m_Substitution.IsEmpty() && DATA[ j].m_Substitution.GetSize() == 1)
            || ( DATA[ j].m_Substitution.IsEmpty() && DATA[ i].m_Substitution.GetSize() == 1)
          )
          {
            BCL_MessageVrb( "Molecules " + util::Format()( i) + " and " + util::Format()( j) + " differ trivially");
            changed_points.PushBack
            (
              storage::Triplet< size_t, size_t, size_t>
              (
                i, j,
                DATA[ i].m_Substitution.IsEmpty() ? DATA[ j].m_Substitution( 0) : DATA[ i].m_Substitution( 0)
              )
            );
            continue;
          }

          // Check the common substitutions
          storage::Set< size_t> changed_vertices;
          changed_vertices.InsertElements( DATA[ i].m_Substitution.Begin(), DATA[ i].m_Substitution.End());
          changed_vertices.InsertElements( DATA[ j].m_Substitution.Begin(), DATA[ j].m_Substitution.End());
          size_t present_in_one( 0);
          for( storage::Set< size_t>::const_iterator itr_ch( changed_vertices.Begin()), itr_ch_end( changed_vertices.End()); itr_ch != itr_ch_end; ++itr_ch)
          {
            bool in_i( DATA[ i].m_Substitution.Find( *itr_ch) < DATA[ i].m_Substitution.GetSize());
            bool in_j( DATA[ j].m_Substitution.Find( *itr_ch) < DATA[ j].m_Substitution.GetSize());
            present_in_one += size_t( in_i ^ in_j); // if substitution is only on a single molecule then it is automatically a difference
          }

          // if molecules have different substitution at only a single differing point then they may still qualify as a single-change pair
          // if there are more than 2 groups that are present in only one of the two molecules, though, there is no reason to compare
          // because they will be different automatically
          if( present_in_one >= 2)
          {
            BCL_MessageVrb( "Molecules " + util::Format()( i) + " and " + util::Format()( j) + " differ by at least " + util::Format()( present_in_one) + " points, skipping");
            continue;
          }

          BCL_MessageVrb( "Molecules " + util::Format()( i) + " and " + util::Format()( j) + " may or may not differ at a single point, doing an in-depth search");

          storage::Vector< size_t> changed_subs
          (
            GetSubstitutionDifferences
            (
              DATA[i].m_MolGraph,
              DATA[i].m_ScaffoldIsomorphism,
              DATA[j].m_MolGraph,
              DATA[j].m_ScaffoldIsomorphism
            )
          );

          BCL_MessageVrb( "Molecules changed at " + util::Format()( changed_subs.GetSize()) + " points");
          if( !changed_subs.IsEmpty())
          {
            std::stringstream message;
            message << "  Molecules " + util::Format()( i) + " and " + util::Format()( j) + " changed at points ";
            for( size_t k( 0); k < changed_subs.GetSize(); ++k)
            {
              message << changed_subs( k) << " ";
            }
            BCL_MessageVrb( message.str());
          }

          if( changed_subs.GetSize() == 1)
          {
            changed_points.PushBack( storage::Triplet< size_t, size_t, size_t>( i, j, changed_subs( 0)));
          }
        }
      }
      return changed_points;
    }

    storage::Vector< size_t> PharmacophoreMapper::GetScaffoldSubstitution
    (
      const graph::ConstGraph< size_t, size_t> GRAPH,
      const storage::Vector< size_t> &ISOMORPHISM
    )
    {
      util::OwnPtr< const graph::ConstGraph< size_t, size_t> > op_subgraph( &GRAPH, false);
      graph::Subgraph< size_t, size_t> sub( op_subgraph, ISOMORPHISM);
      storage::List< storage::Pair< size_t, size_t> > adj_edges( sub.GetAdjacentEdgeIndices());
      storage::Set< size_t> subbed_points;
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_adj( adj_edges.Begin()), itr_adj_end( adj_edges.End());
        itr_adj != itr_adj_end;
        ++itr_adj
      )
      {
        size_t iso_index( ISOMORPHISM.Find( itr_adj->First()));
        subbed_points.Insert( iso_index);
      }
      return storage::Vector< size_t>( subbed_points.Begin(), subbed_points.End());
    }

    storage::Vector< size_t> PharmacophoreMapper::GetSubstitutionDifferences
    (
      graph::ConstGraph< size_t, size_t> GRAPH_MOL_1,
      const storage::Vector< size_t> &ISOMORPHISM_MOL_1,
      graph::ConstGraph< size_t, size_t> GRAPH_MOL_2,
      const storage::Vector< size_t> &ISOMORPHISM_MOL_2
    )
    {
      storage::Vector< size_t> changed_points;
      storage::Set< size_t> changed_point_set;
      if( ISOMORPHISM_MOL_1.IsEmpty() || ISOMORPHISM_MOL_1.GetSize() != ISOMORPHISM_MOL_2.GetSize())
      {
        BCL_MessageStd( "    Molecules cannot have the same scaffold, not comparing them");
        return changed_points;
      }

      size_t iso_size( ISOMORPHISM_MOL_1.GetSize());

      graph::CommonSubgraphIsomorphism< size_t, size_t> csi;

      util::OwnPtr< const graph::ConstGraph< size_t, size_t> > op_graph_1( &GRAPH_MOL_1, false);
      util::OwnPtr< const graph::ConstGraph< size_t, size_t> > op_graph_2( &GRAPH_MOL_2, false);

      csi.SetGraphA( op_graph_1);
      csi.SetGraphB( op_graph_2);

      csi.FindIsomorphism( csi.EstimateUpperBounds(), iso_size);

      storage::Vector< graph::Subgraph< size_t, size_t> > csi_mol_1( csi.GetSubgraphIsomorphismsOfGraphA());
      storage::Vector< graph::Subgraph< size_t, size_t> > csi_mol_2( csi.GetSubgraphIsomorphismsOfGraphB());

      // get adjacent indices for the common isomorphisms; these are the changed points between the two molecules
      storage::List< storage::Pair< size_t, size_t> > csi_adj_edges_1( csi_mol_1( 0).GetAdjacentEdgeIndices());
      storage::List< storage::Pair< size_t, size_t> > csi_adj_edges_2( csi_mol_2( 0).GetAdjacentEdgeIndices());

      // If there was no isomorphism, or the ismorphism completely covers both graphs then return nothing
      if( csi_mol_1.IsEmpty())
      {
        BCL_MessageStd( "    Could not find a common isomorphism");
        return changed_points;
      }

      // get scaffold subgraphs of each molecule graph
      graph::Subgraph< size_t, size_t> sg_1( op_graph_1, ISOMORPHISM_MOL_1);
      graph::Subgraph< size_t, size_t> sg_2( op_graph_2, ISOMORPHISM_MOL_2);

      storage::List< storage::Pair< size_t, size_t> > adj_edges_1( sg_1.GetAdjacentEdgeIndices());
      storage::List< storage::Pair< size_t, size_t> > adj_edges_2( sg_2.GetAdjacentEdgeIndices());

      // Remove edges in graph 1 that connect the scaffold to any substituents
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_adj( adj_edges_1.Begin()),
          itr_adj_end( adj_edges_1.End());
        itr_adj != itr_adj_end;
        ++itr_adj
      )
      {
        GRAPH_MOL_1.RemoveEdge( itr_adj->First(), itr_adj->Second());
      }

      // Remove edges in graph 2 that connect the scaffold to substituents
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_adj( adj_edges_2.Begin()),
          itr_adj_end( adj_edges_2.End());
        itr_adj != itr_adj_end;
        ++itr_adj
      )
      {
        GRAPH_MOL_2.RemoveEdge( itr_adj->First(), itr_adj->Second());
      }

      // make an association between scaffold atoms and the fragment vertices attached to them for graph 1
      storage::Map< size_t, storage::Set< size_t> > g_1_fragments;
      storage::Set< size_t> mol_1_scaff_vertices( ISOMORPHISM_MOL_1.Begin(), ISOMORPHISM_MOL_1.End());
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_adj( adj_edges_1.Begin()),
          itr_adj_end( adj_edges_1.End());
        itr_adj != itr_adj_end;
        ++itr_adj
      )
      {
        util::ShPtr< storage::Vector< size_t> > reachable_vertices( graph::Connectivity::GetVerticesReachableFrom( GRAPH_MOL_1, itr_adj->Second()));
        size_t iso_index( ISOMORPHISM_MOL_1.Find( itr_adj->First()));
        g_1_fragments[ iso_index].InsertElements( reachable_vertices->Begin(), reachable_vertices->End());

        /*
        std::cout << "G1 fragment: ";
        for( size_t f( 0), end_f( reachable_vertices->GetSize()); f < end_f; ++f)
        {
          std::cout << ( *reachable_vertices)( f) << " ";
        }
        std::cout << std::endl;
        */
      }

      // make an association between scaffold atoms and the fragment vertices attached to them for graph 1
      storage::Map< size_t, storage::Set< size_t> > g_2_fragments;
      storage::Set< size_t> mol_2_scaff_vertices( ISOMORPHISM_MOL_2.Begin(), ISOMORPHISM_MOL_2.End());
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_adj( adj_edges_2.Begin()),
          itr_adj_end( adj_edges_2.End());
        itr_adj != itr_adj_end;
        ++itr_adj
      )
      {
        util::ShPtr< storage::Vector< size_t> > reachable_vertices( graph::Connectivity::GetVerticesReachableFrom( GRAPH_MOL_2, itr_adj->Second()));
        size_t iso_index( ISOMORPHISM_MOL_2.Find( itr_adj->First()));
        g_2_fragments[ iso_index].InsertElements( reachable_vertices->Begin(), reachable_vertices->End());
      }

      // inspect each change point for the first molecule
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_csi( csi_adj_edges_1.Begin()), itr_csi_end( csi_adj_edges_1.End());
        itr_csi != itr_csi_end;
        ++itr_csi
      )
      {
        // determine which scaffold point each changed vertex is associated with by iterating through the scaff point->fragment vertex map
        for
        (
          storage::Map< size_t, storage::Set< size_t> >::const_iterator itr_v( g_1_fragments.Begin()), itr_v_end( g_1_fragments.End());
          itr_v != itr_v_end;
          ++itr_v
        )
        {
          // check if any fragments associated with this scaffold vertex contains
          // the a change point
          if( itr_v->second.Contains( itr_csi->Second()))
          {
            changed_point_set.Insert( itr_v->first);
          }
        }
      }

      // inspect each change point for the second molecule
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_csi( csi_adj_edges_2.Begin()), itr_csi_end( csi_adj_edges_2.End());
        itr_csi != itr_csi_end;
        ++itr_csi
      )
      {
        // determine which scaffold point each changed vertex is associated with by iterating through the scaff point->fragment vertex map
        for
        (
          storage::Map< size_t, storage::Set< size_t> >::const_iterator itr_v( g_2_fragments.Begin()), itr_v_end( g_2_fragments.End());
          itr_v != itr_v_end;
          ++itr_v
        )
        {
          if( itr_v->second.Contains( itr_csi->Second()))
          {
            changed_point_set.Insert( itr_v->first);
          }
        }
      }

      // return the list of changed scaffold points
      changed_points = storage::Vector< size_t>( changed_point_set.Begin(), changed_point_set.End());
      return changed_points;
    }

    //! @brief gets a list of which grow points from a scaffold were changed
    //! @param MOL_1_FRAGS map from scaffold atom indices (keys) to indices of MOL_1_GRAPH associated with that atom
    //! @param MOL_2_FRAGS same as above, but for molecule 2
    //! @param MOL_1_GRAPH the whole graph from molecule 1
    //! @param MOL_2_GRAPH the whole graph from molecule 2
    //! @param GROW_POINT_INDICES the indices (map keys) which are considered grow points
    //! @return a vector of grow points which changed from mol 1 to mol 2
    storage::Vector< size_t> PharmacophoreMapper::GetChangedGrowPoints
    (
      const storage::Map< size_t, storage::Vector< size_t> > &MOL_1_FRAGS,
      const storage::Map< size_t, storage::Vector< size_t> > &MOL_2_FRAGS,
      const graph::ConstGraph< size_t, size_t> &MOL_1_GRAPH,
      const graph::ConstGraph< size_t, size_t> &MOL_2_GRAPH,
      const storage::Vector< size_t> &GROW_POINT_INDICES
    ) const
    {

      storage::Vector< size_t> changed_growpts;

      graph::SubgraphIsomorphism< size_t, size_t> frag_comparer;

      // Only check designated grow points for changes
      for( size_t i( 0); i < GROW_POINT_INDICES.GetSize(); ++i)
      {
        const size_t &growpt( GROW_POINT_INDICES( i));

        // If one molecule has a fragment but the other doesn't then it counts as a change
        bool has_frag( MOL_1_FRAGS.Has( growpt));
        if( has_frag != MOL_2_FRAGS.Has( growpt))
        {
          changed_growpts.PushBack( growpt);
          BCL_MessageStd( util::Format()( growpt) + ": only one has a fragment");
          continue;
        }

        if( !has_frag)
        {
          // neither mol1 nor mol2 have a fragment at this growpoint
          continue;
        }

        // both molecules have a fragment at this grow point.  compare the two fragments
        // first get the indices of the fragments
        const storage::Vector< size_t> &mol_frag_indices( MOL_1_FRAGS.GetValue( growpt));
        const storage::Vector< size_t> &mol_comp_frag_indices( MOL_2_FRAGS.GetValue( growpt));

        // if fragments have different numbers of atoms then it is definitely a change
        if( mol_frag_indices.GetSize() != mol_comp_frag_indices.GetSize())
        {
          changed_growpts.PushBack( growpt);
          continue;
        }

        // Graphs of the fragments
        graph::ConstGraph< size_t, size_t> mol_frag( MOL_1_GRAPH.GetSubgraph( mol_frag_indices));
        graph::ConstGraph< size_t, size_t> mol_comp_frag( MOL_2_GRAPH.GetSubgraph( mol_comp_frag_indices));

        frag_comparer.SetGraphExternalOwnership( mol_frag);
        frag_comparer.SetSubgraphExternalOwnership( mol_comp_frag);

        // both fragments have the same size, so if there is an isomorphism they are identical
        if( !frag_comparer.FindIsomorphism())
        {
          changed_growpts.PushBack( growpt);
          continue;
        }
      }
      return changed_growpts;
    }

  } // namespace chemistry
} // namespace bcl
