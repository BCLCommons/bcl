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
#include "chemistry/bcl_chemistry_conformation_comparison_property_field_correlation.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_symmetry_rmsd.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_stopwatch.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonPropertyFieldCorrelation::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonPropertyFieldCorrelation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ConformationComparisonPropertyFieldCorrelation::ConformationComparisonPropertyFieldCorrelation() :
      m_Properties(),
      m_PropertyWeights( size_t( 1), double( 1.0)),
      m_MaxAtomDistanceCriterion( 1.0),
      m_OptimizingWeights( false),
      m_MismatchPenalty( 1.0e-2),
      m_HeavyPenaltyFraction( 0.6),
      m_HeavyMismatchPenalty( 2.0),
      m_AnchorWeight( 1.0e-2)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        this->TryRead( util::ObjectDataLabel(), util::GetLogger());
      }
    }

    //! @brief full constructor
    ConformationComparisonPropertyFieldCorrelation::ConformationComparisonPropertyFieldCorrelation
    (
      const descriptor::Combine< AtomConformationalInterface, float> &PROPERTIES,
      const linal::Vector< float> &PROPERTY_WEIGHTS,
      const double MAX_ATOM_DIST,
      const bool OPTI_WEIGHTS,
      const double MISMATCH_PENALTY,
      const double HEAVY_MISMATCH_FRACTION,
      const double HEAVY_MISMATCH_PENALTY,
      const double ANCHOR_WEIGHT
    ) :
      m_Properties( PROPERTIES),
      m_PropertyWeights( PROPERTY_WEIGHTS),
      m_MaxAtomDistanceCriterion( MAX_ATOM_DIST),
      m_OptimizingWeights( OPTI_WEIGHTS),
      m_MismatchPenalty( MISMATCH_PENALTY),
      m_HeavyPenaltyFraction( HEAVY_MISMATCH_FRACTION),
      m_HeavyMismatchPenalty( HEAVY_MISMATCH_PENALTY),
      m_AnchorWeight( ANCHOR_WEIGHT)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        this->TryRead( util::ObjectDataLabel(), util::GetLogger());
      }
    }

    //! virtual copy constructor
    ConformationComparisonPropertyFieldCorrelation *ConformationComparisonPropertyFieldCorrelation::Clone() const
    {
      return new ConformationComparisonPropertyFieldCorrelation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonPropertyFieldCorrelation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonPropertyFieldCorrelation::GetAlias() const
    {
      static const std::string s_name( "PropertyFieldDistance");
      return s_name;
    }

    //! @brief get the max ratio of matched atoms between molecules a and b
    const double ConformationComparisonPropertyFieldCorrelation::GetMaxAtomDistance() const
    {
      return m_MaxAtomDistanceCriterion;
    }

    //! @brief get the linear penalty for unmatched atoms
    const double ConformationComparisonPropertyFieldCorrelation::GetLinearPenalty() const
    {
      return m_MismatchPenalty;
    }

    //! @brief get the fraction below which the heavy mismatch penalty is applied
    const double ConformationComparisonPropertyFieldCorrelation::GetHeavyPenaltyFraction() const
    {
      return m_HeavyPenaltyFraction;
    }

    //! @brief get penalty for having < m_HeavyPenaltyFraction of atoms mutually matched on the maximally matched molecule
    const double ConformationComparisonPropertyFieldCorrelation::GetHeavyPenalty() const
    {
      return m_HeavyMismatchPenalty;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return the RMSD between first and second molecule
    double ConformationComparisonPropertyFieldCorrelation::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      return PropertyDistanceScorer( MOLECULE_A, MOLECULE_B, m_MaxAtomDistanceCriterion).First();
    }

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return RMSDX between first and second molecule and the maximum matched ratio of atom pairs
    storage::Pair< double, double> ConformationComparisonPropertyFieldCorrelation::PropertyDistanceScorer
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B,
      const float &MAX_ATOM_DISTANCE,
      const storage::Vector< size_t> &EXCLUDE_ATOMS_A,
      const storage::Vector< size_t> &EXCLUDE_ATOMS_B
    ) const
    {
      const size_t n_atoms_a( MOLECULE_A.GetNumberAtoms());
      const size_t n_atoms_b( MOLECULE_B.GetNumberAtoms());
      const size_t sum_atoms( n_atoms_a + n_atoms_b);
      const size_t n_properties( m_Properties.GetSizeOfFeatures());

      util::SiPtrVector< const linal::Vector3D> gridpoints;
      gridpoints.AllocateMemory( sum_atoms);
      gridpoints.Append( MOLECULE_A.GetAtomCoordinates());
      gridpoints.Append( MOLECULE_B.GetAtomCoordinates());

      // find the neighbors to atoms in molecule A
      VoxelGridAtom vg_a( MAX_ATOM_DISTANCE), vg_b( MAX_ATOM_DISTANCE);
      vg_a.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( MOLECULE_A.GetAtomsIterator(), MOLECULE_A.GetAtomsIterator().End()));
      vg_b.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( MOLECULE_B.GetAtomsIterator(), MOLECULE_B.GetAtomsIterator().End()));
      auto neighbors_full( vg_a.GetNeighborsIn( vg_b, MAX_ATOM_DISTANCE));

      // get the atom indices for the atoms we will keep
      storage::Vector< size_t> keep_indices_a, keep_indices_b;
      for( size_t a_i( 0), a_sz( MOLECULE_A.GetSize()); a_i < a_sz; ++a_i)
      {
        // do not include the atoms we said we wanted to exclude
        bool exclude( false);
        for( size_t ea_i( 0), ea_sz( EXCLUDE_ATOMS_A.GetSize()); ea_i < ea_sz; ++ea_i)
        {
          if( a_i == EXCLUDE_ATOMS_A( ea_i))
          {
            exclude = true;
            break;
          }
        }
        // if we did not find it in the exclusion indices, add it
        if( !exclude)
        {
          keep_indices_a.PushBack( a_i);
        }
      }
      for( size_t b_i( 0), b_sz( MOLECULE_B.GetSize()); b_i < b_sz; ++b_i)
      {
        // do not include the atoms we said we wanted to exclude
        bool exclude( false);
        for( size_t eb_i( 0), eb_sz( EXCLUDE_ATOMS_B.GetSize()); eb_i < eb_sz; ++eb_i)
        {
          if( b_i == EXCLUDE_ATOMS_B( eb_i))
          {
            exclude = true;
            break;
          }
        }
        // if we did not find it in the exclusion indices, add it
        if( !exclude)
        {
          keep_indices_b.PushBack( b_i);
        }
      }

      // prune neighbor pairs
      auto neighbors( neighbors_full);
      neighbors.Reset();
      for
      (
          auto nbr_itr( neighbors_full.Begin()), nbr_itr_end( neighbors_full.End());
          nbr_itr != nbr_itr_end;
          ++nbr_itr
      )
      {
        // loop over atoms in A that we want
        for( size_t a_i( 0), a_sz( keep_indices_a.GetSize()); a_i < a_sz; ++a_i)
        {
          // loop over atoms in B that we want
          for( size_t b_i( 0), b_sz( keep_indices_b.GetSize()); b_i < b_sz; ++b_i)
          {
            const size_t mol_a_pos( MOLECULE_A.GetAtomIndex( *nbr_itr->First()));
            const size_t mol_b_pos( MOLECULE_B.GetAtomIndex( *nbr_itr->Second()));
            if( mol_a_pos == keep_indices_a( a_i) && mol_b_pos == keep_indices_b( b_i))
            {
              neighbors.InsertElement( *nbr_itr);
            }
          }
        }
      }

      const size_t n_gridpoints( gridpoints.GetSize());
      linal::Matrix< float> properties_gridpoints( n_gridpoints, n_properties, 0.0);
      math::RunningAverage< float> property_distance;

      m_Properties.SetObject( MOLECULE_A);
      iterate::Generic< const AtomConformationalInterface> r_itr_a( MOLECULE_A.GetAtomsIterator());

      storage::Vector< storage::List< storage::Pair< double, size_t> > > mol_a_best_match( n_atoms_a, storage::List< storage::Pair< double, size_t> >());
      storage::Vector< storage::List< storage::Pair< double, size_t> > > mol_b_best_match( n_atoms_b, storage::List< storage::Pair< double, size_t> >());

      // Find neighbors
      for( auto itr_neighbors( neighbors.Begin()), itr_neighbors_end( neighbors.End()); itr_neighbors != itr_neighbors_end; ++itr_neighbors)
      {
        const size_t mol_a_pos( MOLECULE_A.GetAtomIndex( *itr_neighbors->First()));
        const size_t mol_b_pos( MOLECULE_B.GetAtomIndex( *itr_neighbors->Second()));
        const double distance( itr_neighbors->Third());
        mol_a_best_match( mol_a_pos).PushBack( storage::Pair< double, size_t>( distance, mol_b_pos));
        mol_b_best_match( mol_b_pos).PushBack( storage::Pair< double, size_t>( distance, mol_a_pos));
      }

      size_t mol_a_pos = 0;
      math::RunningAverage< double> ave_distance_a, ave_distance_b;
      double normalization( 0);
      size_t counter_a( 0);
      storage::List< double> all_atom_dist_a, all_atom_dist_total_a;

      // compute the value of each row (or property in row) down each column
      for
      (
        descriptor::Iterator< AtomConformationalInterface> itr_a( descriptor::Type( 1, true, descriptor::Type::e_Symmetric), MOLECULE_A);
        itr_a.NotAtEnd();
        ++itr_a, ++mol_a_pos
      )
      {
        auto row( m_Properties( itr_a));
        auto row_gridpoints_aa( properties_gridpoints.GetRow( mol_a_pos));
        for( size_t c( 0); c < n_properties; ++c)
        {
          normalization += math::Sqr( row( c)) * m_PropertyWeightsSqr( c);
        }
        row_gridpoints_aa += row;

        storage::List< storage::Pair< double, size_t> > &list_a( mol_a_best_match( mol_a_pos));
        if( list_a.IsEmpty())
        {
          continue;
        }
        else if( ++list_a.Begin() == list_a.End()) // list of size 1
        {
          auto row_gridpoints_a( properties_gridpoints.GetRow( list_a.FirstElement().Second() + n_atoms_a));
          row_gridpoints_a += row;
          BCL_MessageDbg( "MolA atom: " + util::Format()( mol_a_pos) + " weight: 1.0 to atom: " + util::Format()( list_a.FirstElement().Second()));
        }
        else
        {
          double sum_dist_from_tolerance( 0.0);
          for
          (
            storage::List< storage::Pair< double, size_t> >::const_iterator itr_a_list( list_a.Begin()), itr_a_list_end( list_a.End());
            itr_a_list != itr_a_list_end;
            ++itr_a_list
          )
          {
            sum_dist_from_tolerance += MAX_ATOM_DISTANCE + 0.01 - itr_a_list->First();
          }
          for
          (
            storage::List< storage::Pair< double, size_t> >::const_iterator itr_a_list( list_a.Begin()), itr_a_list_end( list_a.End());
            itr_a_list != itr_a_list_end;
            ++itr_a_list
          )
          {
            auto row_gridpoints_a( properties_gridpoints.GetRow( itr_a_list->Second() + n_atoms_a));
            row_gridpoints_a += row * float( ( MAX_ATOM_DISTANCE + 0.01 - itr_a_list->First()) / sum_dist_from_tolerance);
            BCL_MessageDbg
            (
              "MolA atom: " + util::Format()( mol_a_pos) + " weight: " + util::Format()( float( ( MAX_ATOM_DISTANCE + 0.01 - itr_a_list->First()) / sum_dist_from_tolerance))
              + " to atom: " + util::Format()( itr_a_list->Second())
            );
          }
        }
        if( m_OptimizingWeights)
        {
          for
          (
            auto itr_a_reduced_list( list_a.Begin()), itr_a_reduced_list_end( list_a.End());
            itr_a_reduced_list != itr_a_reduced_list_end; ++itr_a_reduced_list
          )
          {
            all_atom_dist_a.PushBack( itr_a_reduced_list->First());
          }
        }

        //Benchmark data collection
        if( !list_a.IsEmpty() && keep_indices_a.Find( mol_a_pos) < keep_indices_a.GetSize())
        {
          ++counter_a;
        }
      }

      m_Properties.SetObject( MOLECULE_B);
      size_t grid_index( n_atoms_a), mol_b_pos( 0);
      size_t counter_b( 0);
      storage::List< double> all_atom_dist_b, all_atom_dist_total_b;
      for
      (
        descriptor::Iterator< AtomConformationalInterface> itr_b( descriptor::Type( 1, true, descriptor::Type::e_Symmetric), MOLECULE_B);
        itr_b.NotAtEnd();
        ++itr_b, ++grid_index, ++mol_b_pos
      )
      {
        auto row( m_Properties( itr_b));
        auto row_gridpoints_bb( properties_gridpoints.GetRow( grid_index));
        for( size_t c( 0); c < n_properties; ++c)
        {
          normalization += math::Sqr( row( c)) * m_PropertyWeightsSqr( c);
        }
        row_gridpoints_bb -= row;
        storage::List< storage::Pair< double, size_t> > &list_b( mol_b_best_match( mol_b_pos));

        if( list_b.IsEmpty())
        {
          continue;
        }
        else if( ++list_b.Begin() == list_b.End())
        {
          auto row_gridpoints_b( properties_gridpoints.GetRow( list_b.FirstElement().Second()));
          row_gridpoints_b -= row;
          BCL_MessageDbg( "MolB atom: " + util::Format()( mol_b_pos) + " weight: 1.0 to atom: " + util::Format()( list_b.FirstElement().Second()));
        }
        else
        {
          double sum_dist_from_tolerance( 0.0);
          for
          (
            storage::List< storage::Pair< double, size_t> >::const_iterator itr_b_list( list_b.Begin()), itr_b_list_end( list_b.End());
            itr_b_list != itr_b_list_end;
            ++itr_b_list
          )
          {
            sum_dist_from_tolerance += MAX_ATOM_DISTANCE + 0.01 - itr_b_list->First();
          }
          for
          (
            storage::List< storage::Pair< double, size_t> >::const_iterator itr_b_list( list_b.Begin()), itr_b_list_end( list_b.End());
            itr_b_list != itr_b_list_end;
            ++itr_b_list
          )
          {
            auto row_gridpoints_b( properties_gridpoints.GetRow( itr_b_list->Second()));
            row_gridpoints_b -= row * ( float( ( MAX_ATOM_DISTANCE + 0.01 - itr_b_list->First()) / sum_dist_from_tolerance));
            BCL_MessageDbg
            (
              "MolB atom: " + util::Format()( mol_b_pos) + " weight: " + util::Format()( ( float( ( MAX_ATOM_DISTANCE + 0.01 - itr_b_list->First()) / sum_dist_from_tolerance)))
              + " to atom: " + util::Format()( itr_b_list->Second())
            );
          }
        }
        if( m_OptimizingWeights)
        {
          for
          (
            auto itr_b_reduced_list( list_b.Begin()), itr_b_reduced_list_end( list_b.End());
            itr_b_reduced_list != itr_b_reduced_list_end; ++itr_b_reduced_list
          )
          {
            all_atom_dist_b.PushBack( itr_b_reduced_list->First());
          }
        }

        if( !list_b.IsEmpty() && keep_indices_b.Find( mol_b_pos) < keep_indices_b.GetSize())
        {
          ++counter_b;
        }
      }

      //Weight and sum property gridpoints
      float properties_gridpoints_sum( 0);
      double symmetry_rmsd( 0.0);
      util::ObjectDataLabel id_label( "ID");

      linal::Vector< float> unweighted_property_sums( n_properties, 0.0);
//      linal::Vector< float> gridpoint_property_sums( n_gridpoints, 0.0);
      for( size_t r( 0); r < n_gridpoints; ++r)
      {
        for( size_t c( 0); c < n_properties; ++c)
        {
//          gridpoint_property_sums( r) += math::Sqr( properties_gridpoints( r, c)) * m_PropertyWeightsSqr( c);
          properties_gridpoints_sum += math::Sqr( properties_gridpoints( r, c)) * m_PropertyWeightsSqr( c);
          unweighted_property_sums( c) += math::Sqr( properties_gridpoints( r, c));
        }
      }

//      // per gridpoint score
//      for
//      (
//          auto grid_itr( gridpoint_property_sums.Begin()), grid_itr_end( gridpoint_property_sums.End());
//          grid_itr != grid_itr_end;
//          ++grid_itr
//      )
//      {
//        ( *grid_itr) = math::Sqrt( *grid_itr) / std::max( math::Sqrt( normalization), std::numeric_limits< double>::epsilon());
//      }
//      BCL_Debug( gridpoint_property_sums);

      // property normalization
      double rmsdx( math::Sqrt( properties_gridpoints_sum) / std::max( math::Sqrt( normalization), std::numeric_limits< double>::epsilon()));

      // add penalties to property distance score
      double mol_a_matched_ratio( util::GetUndefinedDouble()), mol_b_matched_ratio( util::GetUndefinedDouble());
      double max_matched( util::GetUndefinedDouble());
      if( keep_indices_a.GetSize() && keep_indices_b.GetSize())
      {
        //max ratio of matched atoms between molecules a and b
        mol_a_matched_ratio = static_cast<double>(counter_a) / static_cast<double>(keep_indices_a.GetSize());
        mol_b_matched_ratio = static_cast<double>(counter_b) / static_cast<double>(keep_indices_b.GetSize());
        max_matched = std::max( mol_a_matched_ratio, mol_b_matched_ratio);
        rmsdx += m_HeavyMismatchPenalty *
            ( max_matched < m_HeavyPenaltyFraction ? math::Sqr( ( m_HeavyPenaltyFraction - max_matched) / m_HeavyPenaltyFraction) : 0.0);
        rmsdx += m_MismatchPenalty * ( 1.0 - max_matched);
      }
      double min_ave_distance( std::min( ave_distance_a.GetAverage(), ave_distance_b.GetAverage()));
      rmsdx += m_AnchorWeight * min_ave_distance;

      //collect data in .csv for property weighting
      if( m_OptimizingWeights)
      {
        size_t mol_id_a(MOLECULE_A.GetFromCache(id_label)(0));
        size_t mol_id_b(MOLECULE_B.GetFromCache(id_label)(0));
        //identify the native pose of my current molecule from EnsembleB
        util::SiPtr< const FragmentComplete> native_pose;
        for( FragmentEnsemble::const_iterator itr( m_EnsembleB->Begin()), itr_end( m_EnsembleB->End()); itr != itr_end; ++itr)
        {
          if( itr->GetFromCache( id_label)( 0) == mol_id_a)
          {
            native_pose = &*itr;
            break;
          }
        }
        symmetry_rmsd = m_SymmetryRmsdComparers( mol_id_a)( MOLECULE_A, *native_pose);

        io::OFStream out_csv;
        io::File::MustOpenOFStream( out_csv, m_Path + m_FileName + ".csv", std::ios::app);
        out_csv << mol_id_a << ',' << mol_id_b << ',';
        out_csv << n_atoms_a << ',' << n_atoms_b << ','
            << counter_a << ',' << counter_b << ','
            << mol_a_matched_ratio << ',' << mol_b_matched_ratio << ','
            << rmsdx << "," << min_ave_distance << ","
            << ave_distance_a.GetAverage() << ',' << ave_distance_b.GetAverage() << ','
            << std::min( ave_distance_a.GetAverage(), ave_distance_b.GetAverage()) << ',';
        for( size_t output_index( 0); output_index < unweighted_property_sums.GetSize(); ++output_index)
        {
          out_csv << math::Sqrt( unweighted_property_sums( output_index) / double( sum_atoms)) << ',';
        }
        out_csv << symmetry_rmsd << '\n';
        io::File::CloseClearFStream( out_csv);

        io::OFStream out;
        io::File::MustOpenOFStream( out, m_Path + m_FileName + "_all_atom_dist_a.csv", std::ios::app);
        out << all_atom_dist_a;
        io::File::CloseClearFStream( out);

        io::File::MustOpenOFStream( out, m_Path + m_FileName + "_total_all_atom_dist_a.csv", std::ios::app);
        out << all_atom_dist_total_a;
        io::File::CloseClearFStream( out);

        io::File::MustOpenOFStream( out, m_Path + m_FileName + "_all_atom_dist_b.csv", std::ios::app);
        out << all_atom_dist_b;
        io::File::CloseClearFStream( out);

        io::File::MustOpenOFStream( out, m_Path + m_FileName + "_total_all_atom_dist_b.csv", std::ios::app);
        out << all_atom_dist_total_b;
        io::File::CloseClearFStream( out);

        io::File::MustOpenOFStream( out, m_Path + m_FileName + "_total_all_atom_dist_min.csv", std::ios::app);
        out << ( ave_distance_a.GetAverage() < ave_distance_b.GetAverage() ? all_atom_dist_total_a : all_atom_dist_total_b);
        io::File::CloseClearFStream( out);
      }
      //done
      storage::Pair< double, double> output( rmsdx, max_matched);
      return output;
    }

    storage::Pair< storage::Vector< size_t>, storage::Vector< size_t> > ConformationComparisonPropertyFieldCorrelation::GetAlignedAtoms
    (
      const FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const storage::Vector< size_t> &KEEP_INDICES_A,
      const storage::Vector< size_t> &KEEP_INDICES_B
    )
    {
      storage::Vector< size_t> atom_indices_aligned_a, atom_indices_aligned_b;
      const size_t n_atoms_a( MOLECULE_A.GetNumberAtoms());
      const size_t n_atoms_b( MOLECULE_B.GetNumberAtoms());

      iterate::Generic< const AtomConformationalInterface> r_itr_a( MOLECULE_A.GetAtomsIterator());
      storage::Vector< storage::Pair< double, size_t> > mol_a_best_match( n_atoms_a, storage::Pair< double, size_t>( util::GetUndefined< double>(), n_atoms_b));
      storage::Vector< storage::Pair< double, size_t> > mol_b_best_match( n_atoms_b, storage::Pair< double, size_t>( util::GetUndefined< double>(), n_atoms_a));

      // for each atom in MOLECULE_A find the nearest atom in MOLECULE_B
      for( size_t mol_a_pos( 0); mol_a_pos < n_atoms_a; ++r_itr_a, ++mol_a_pos)
      {
        // make sure this is an allowed molecule a index
        bool keep( false);
        for( size_t keep_a_i( 0), keep_a_sz( KEEP_INDICES_A.GetSize()); keep_a_i < keep_a_sz; ++keep_a_i)
        {
          if( mol_a_pos == KEEP_INDICES_A( keep_a_i))
          {
            keep = true;
            break;
          }
        }
        if( !keep)
        {
          continue;
        }

        size_t mol_b_pos( 0);
        const linal::Vector3D &a_pos( r_itr_a->GetPosition());
        for
        (
          iterate::Generic< const AtomConformationalInterface> r_itr_b( MOLECULE_B.GetAtomsIterator());
          mol_b_pos < n_atoms_b;
          ++r_itr_b, ++mol_b_pos
        )
        {
          // make sure this is an allowed molecule b index
          bool keep( false);
          for( size_t keep_b_i( 0), keep_b_sz( KEEP_INDICES_B.GetSize()); keep_b_i < keep_b_sz; ++keep_b_i)
          {
            if( mol_b_pos == KEEP_INDICES_B( keep_b_i))
            {
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }

          const float distance( linal::SquareDistance( r_itr_b->GetPosition(), a_pos));
          if( mol_a_best_match( mol_a_pos).Second() == n_atoms_b || distance < mol_a_best_match( mol_a_pos).First())
          {
            mol_a_best_match( mol_a_pos).Second() = mol_b_pos;
            mol_a_best_match( mol_a_pos).First() = distance;
          }
          if( mol_b_best_match( mol_b_pos).Second() == n_atoms_a || distance < mol_b_best_match( mol_b_pos).First())
          {
            mol_b_best_match( mol_b_pos).Second() = mol_a_pos;
            mol_b_best_match( mol_b_pos).First() = distance;
          }
        }
      }
      atom_indices_aligned_a.Reset();
      atom_indices_aligned_b.Reset();

      // save the indices of the mutually closest atom pairs from MOLECULE_A and MOLECULE_B
      for( size_t mol_a_pos( 0); mol_a_pos < n_atoms_a; ++mol_a_pos)
      {
        if( mol_a_best_match( mol_a_pos).Second() < n_atoms_b && mol_b_best_match( mol_a_best_match( mol_a_pos).Second()).Second() == mol_a_pos)
        {
          atom_indices_aligned_a.PushBack( mol_a_pos);
          atom_indices_aligned_b.PushBack( mol_a_best_match( mol_a_pos).Second());
        }
      }

      // return the final sets of mutually matched atoms
      return storage::Pair< storage::Vector< size_t>, storage::Vector< size_t>>( std::make_pair( atom_indices_aligned_a, atom_indices_aligned_b));
    }

    //! @brief prepare the class for comparing a conformation
    //! @param MOLECULE the molecule to prepare to compare
    void ConformationComparisonPropertyFieldCorrelation::Prepare( const ConformationInterface &MOLECULE) const
    {
      // check that the properties are initialized
      if( !m_Properties.GetSizeOfFeatures())
      {
        // alignment properties and weights
        util::ObjectDataLabel property_labels
        (
          "("
          "Combine(Define(ChargeNegative=Multiply(Less(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
          "Define(ChargePositive=Multiply(Greater(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
          "Define(HbondAcceptorsStrict=Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))),"
          "Define(ENegOffset=Subtract(lhs=Atom_ElectroNegativity,rhs=Constant(2.5))),"
          "Define(IsENeg=Greater(lhs=ENegOffset,rhs=Constant(0))),"
          "Define(PolarTernary=Subtract(lhs=Multiply(Add(HbondAcceptorsStrict, Atom_HbondDonors), Constant(2)), rhs=Constant(1))),"
          "Define(Atom_Hydrophobic=Not(DescriptorSum(Abs(Partial(indices(1,2,5,8),2DASign(property=PolarTernary, steps=3)))))),"
          "Atom_SigmaCharge,"
          "HbondAcceptorsStrict,"
          "Atom_HbondDonors,"
          "Atom_Polarizability,"
          "Atom_AromaticityAxes,"
          "Multiply(IsENeg,ENegOffset),"
          "Atom_Hydrophobic,"
          "Atom_VDWVolume)"
          ")"
        );
        m_Properties = descriptor::Combine< AtomConformationalInterface, float>( PrepareProperties( property_labels));
      }

      // calculate all properties on the molecule; this avoids potential clashes if multiple threads calculate properties
      // for the same molecule and try to update the cache simultaneously
      m_Properties.CollectValuesOnEachElementOfObject( MOLECULE);

      for
      (
          storage::Vector< descriptor::CheminfoProperty>::iterator
          itr( m_Properties.Begin()), itr_end( m_Properties.End());
          itr != itr_end;
          ++itr
      )
      {
        MOLECULE.Uncache( ( *itr)->GetCacheLabel());
      }
    }

    void ConformationComparisonPropertyFieldCorrelation::PrepareEnsemble( const FragmentEnsemble &ENSEMBLE) const
    {
      // check that the properties are initialized
      if( !m_Properties.GetSizeOfFeatures())
      {
        // alignment properties and weights
        util::ObjectDataLabel property_labels
        (
          "("
          "Combine(Define(ChargeNegative=Multiply(Less(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
          "Define(ChargePositive=Multiply(Greater(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
          "Define(HbondAcceptorsStrict=Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))),"
          "Define(ENegOffset=Subtract(lhs=Atom_ElectroNegativity,rhs=Constant(2.5))),"
          "Define(IsENeg=Greater(lhs=ENegOffset,rhs=Constant(0))),"
          "Define(PolarTernary=Subtract(lhs=Multiply(Add(HbondAcceptorsStrict, Atom_HbondDonors), Constant(2)), rhs=Constant(1))),"
          "Define(Atom_Hydrophobic=Not(DescriptorSum(Abs(Partial(indices(1,2,5,8),2DASign(property=PolarTernary, steps=3)))))),"
          "Atom_SigmaCharge,"
          "HbondAcceptorsStrict,"
          "Atom_HbondDonors,"
          "Atom_Polarizability,"
          "Atom_AromaticityAxes,"
          "Multiply(IsENeg,ENegOffset),"
          "Atom_Hydrophobic,"
          "Atom_VDWVolume)"
          ")"
        );
        m_Properties = descriptor::Combine< AtomConformationalInterface, float>( PrepareProperties( property_labels));
      }

      //molecules from ensemble B
      //used for weight optimization when ensemble B contains native poses
      m_EnsembleB = util::ToSiPtr( ENSEMBLE);
      m_SymmetryRmsdComparers.Resize( ENSEMBLE.GetSize(), ConformationComparisonBySymmetryRmsd( false));
      if( m_OptimizingWeights)
      {
        linal::Vector< float> ensemble_index( 1, 0.0);
        for( FragmentEnsemble::const_iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End()); itr != itr_end; ++itr, ensemble_index( 0) += 1.0)
        {
          itr->Cache( util::ObjectDataLabel( "ID"), ensemble_index);
        }
      }
      ConformationComparisonInterface::PrepareEnsemble( ENSEMBLE);
    }

    //! @brief set atom properties and remove size zero properties
    //! @param PROPERTIES descriptor label to use
    descriptor::Combine< AtomConformationalInterface, float> ConformationComparisonPropertyFieldCorrelation::PrepareProperties( const util::ObjectDataLabel &PROPERTIES) const
    {
      descriptor::Combine< AtomConformationalInterface, float> properties, read_check;

      // Properties must be molecule-wide; assert they are valid
      read_check.SetDimension( 0);
      read_check.AssertRead( PROPERTIES);

      // Fetch the aliases for each property
      for
      (
        descriptor::Combine< AtomConformationalInterface, float>::const_iterator
        itr_prop( read_check.Begin()), itr_prop_end( read_check.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        // Do not include 0-size descriptors, e.g. Define, ForEach, etc...
        if( ( *itr_prop)->GetSizeOfFeatures())
        {
          BCL_MessageVrb( "Preparing properties: " + itr_prop->GetString());
          properties.PushBack( *itr_prop);
        }
      }

      return properties;
    } // PrepareProperties

    void ConformationComparisonPropertyFieldCorrelation::SetProperties( const descriptor::Combine< AtomConformationalInterface, float> &PROPERTIES)
    {
      m_Properties = PROPERTIES;
    }

    void ConformationComparisonPropertyFieldCorrelation::SetPropertyWeights( const linal::Vector< float> &PROPERTY_WEIGHTS)
    {
      m_PropertyWeights = PROPERTY_WEIGHTS;
      m_PropertyWeightsSqr = m_PropertyWeights;
      linal::ElementwiseMultiply( m_PropertyWeightsSqr, m_PropertyWeightsSqr);
    }

    //! @brief get the max ratio of matched atoms between molecules a and b
    void ConformationComparisonPropertyFieldCorrelation::SetMaxAtomDistance( const double MAX_ATOM_DISTANCE)
    {
      m_MaxAtomDistanceCriterion = MAX_ATOM_DISTANCE;
    }

    //! @brief get the linear penalty for unmatched atoms
    void ConformationComparisonPropertyFieldCorrelation::SetLinearPenalty( const double MISMATCH_PENALTY)
    {
      m_MismatchPenalty = MISMATCH_PENALTY;
    }

    //! @brief get the fraction below which the heavy mismatch penalty is applied
    void ConformationComparisonPropertyFieldCorrelation::SetHeavyPenaltyFraction( const double HEAVY_PENALTY_FRACTION)
    {
      m_HeavyPenaltyFraction = HEAVY_PENALTY_FRACTION;
    }

    //! @brief get penalty for having < m_HeavyPenaltyFraction of atoms mutually matched on the maximally matched molecule
    void ConformationComparisonPropertyFieldCorrelation::SetHeavyPenalty( const double HEAVY_MISMATCH_PENALTY)
    {
      m_HeavyMismatchPenalty = HEAVY_MISMATCH_PENALTY;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonPropertyFieldCorrelation::GetSerializer() const
    {

      io::Serializer parameters;

      parameters.SetClassDescription
      (
        "Computes distance or correlation between molecular fields"
      );

      parameters.AddInitializer
      (
        "properties",
        "atom properties to consider, use multiply(Constant(X),property y) for weighting",
        io::Serialization::GetAgent( &m_Properties),
        util::ObjectDataLabel
        (
          "Combine(Define(ChargeNegative=Multiply(Less(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
          "Define(ChargePositive=Multiply(Greater(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
          "Define(HbondAcceptorsStrict=Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))),"
          "Define(ENegOffset=Subtract(lhs=Atom_ElectroNegativity,rhs=Constant(2.5))),"
          "Define(IsENeg=Greater(lhs=ENegOffset,rhs=Constant(0))),"
          "Define(PolarTernary=Subtract(lhs=Multiply(Add(HbondAcceptorsStrict, Atom_HbondDonors), Constant(2)), rhs=Constant(1))),"
          "Define(Atom_Hydrophobic=Not(DescriptorSum(Abs(Partial(indices(1,2,5,8),2DASign(property=PolarTernary, steps=3)))))),"
          "Atom_SigmaCharge,"
          "HbondAcceptorsStrict,"
          "Atom_HbondDonors,"
          "Atom_Polarizability,"
          "Atom_AromaticityAxes,"
          "Multiply(IsENeg,ENegOffset),"
          "Atom_Hydrophobic,"
          "Atom_VDWVolume)"
        )
      );

      parameters.AddInitializer
      (
        "property_weights",
        "Weighting to give the properties. Typically these are the inverse standard deviations of the properties over a representative set of molecules, "
        "which can be computed w/ molecule:Properties -statistics",
        io::Serialization::GetAgent( &m_PropertyWeights),
        util::ObjectDataLabel( "(5, 3.5, 7, 1.71, 2.41, 2.41, 2.41, 3.0, 2.0, 0.714)")
      );

      parameters.AddInitializer
      (
        "max_atom_distance_criterion",
        "independent scoring: max atom distance when directly calling this class for scoring",
        io::Serialization::GetAgent( &m_MaxAtomDistanceCriterion),
        "1.0"
      );

      parameters.AddInitializer
      (
        "optimizing_weights",
        "string - perform molecule:Compare w/ SymmetryRMSD to native binding poses in EnsembleB",
        io::Serialization::GetAgent( &m_OptimizingWeights),
        "false"
      );

      parameters.AddInitializer
      (
        "anchor_weight",
        "deprecated; "
        "enter in scientific notation",
        io::Serialization::GetAgent( &m_AnchorWeight),
        "1.0e-3"
      );

      parameters.AddInitializer
      (
        "linear_mismatch_penalty",
        "this penalty is applied as 'penalty = linear_mismatch_penalty * (1.0 - max_fraction_matched)'; "
        "enter in scientific notation",
        io::Serialization::GetAgent( &m_MismatchPenalty),
        "1.0e-2"
      );

      parameters.AddInitializer
      (
        "heavy_mismatch_penalty_fraction",
        "molecule alignments where at least one molecule has fewer than this fraction of its atoms "
        "mutually matched to atoms in the other molecule suffer the 'heavy_mismatch_penalty'",
        io::Serialization::GetAgent( &m_HeavyPenaltyFraction),
        "0.6"
      );

      parameters.AddInitializer
      (
        "heavy_mismatch_penalty",
        "this penalty is applied as "
        "'heavy_mismatch_penalty * ( Sqr(( 0.6 - max_matched) / 0.6))' "
        "when the fraction of matched atoms in the maximally matched molecule is less than heavy_mismatch_penalty_fraction; "
        "enter in scientific notation",
        io::Serialization::GetAgent( &m_HeavyMismatchPenalty),
        "2"
      );

      parameters.AddInitializer
      (
        "weight_file_path",
        "path for outputting .csv files for weight optimization",
        io::Serialization::GetAgent( &m_Path),
        "/tmp/"
      );

      parameters.AddInitializer
      (
        "data_file_name",
        "file tag for unweighted property data in ,csv",
        io::Serialization::GetAgent( &m_FileName),
        "unweighted_files"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ConformationComparisonPropertyFieldCorrelation::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      if( m_PropertyWeights.GetSize() != m_Properties.GetSizeOfFeatures())
      {
        ERR_STREAM << "Wrong number of property weights; had " << m_PropertyWeights.GetSize() << " but had " << m_Properties.GetSizeOfFeatures() << " descriptors";
        return false;
      }
      if( m_OptimizingWeights)
      {
        io::OFStream out_csv;
        io::File::MustOpenOFStream( out_csv, m_Path + m_FileName + ".csv");
        io::File::CloseClearFStream( out_csv);
      }
      m_PropertyWeightsSqr = m_PropertyWeights;
      // compute square property distances. This is b/c the user-given weights are typically inverse standard deviations of the given property,
      // and so to compute the Mahalanobis distance, we we need inverse variances
      linal::ElementwiseMultiply( m_PropertyWeightsSqr, m_PropertyWeightsSqr);
      return true;
    }
  } // namespace chemistry
} // namespace bcl
