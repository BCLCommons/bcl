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
#include "chemistry/bcl_chemistry_ligand_pocket_fit_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_multi_align.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_linear_function.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_pair.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
//    const util::SiPtr< const util::ObjectInterface> LigandPocketFitScore::s_Instance
//    (
//      util::Enumerated< ConformationComparisonInterface>::AddInstance( new LigandPocketFitScore())
//    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LigandPocketFitScore::LigandPocketFitScore
    (
    ) :
      m_CentroidWeight( 0.20),
      m_StartPosition( util::GetUndefinedDouble()),
      m_Properties
      (
        size_t( 1),
        descriptor::CheminfoProperty( descriptor::Constants< AtomConformationalInterface, float>( 1.0))
      ),
      m_NumberOfAtomsToAlign
      (
        size_t( 100)
      ),
      m_PropertyWeights
      (
        size_t( 1),
        double( 1.0)
      ),
      m_PropertyCorrelationLengths
      (
        size_t( 1),
        double( 0.5)
      )
    {
    }

    //! @brief pocket constructor
    LigandPocketFitScore::LigandPocketFitScore
    (
      const FragmentComplete &POCKET
    ) :
      m_BFactor( POCKET.GetSize(), 1.0),
      m_CentroidWeight( 0.20),
      m_Pocket( POCKET),
      m_StartPosition( POCKET.GetCenter()),
      m_Properties
      (
        size_t( 1),
        descriptor::CheminfoProperty( descriptor::Constants< AtomConformationalInterface, float>( 1.0))
      ),
      m_NumberOfAtomsToAlign
      (
        size_t( 100)
      ),
      m_PropertyWeights
      (
        size_t( 1),
        double( 1.0)
      ),
      m_PropertyCorrelationLengths
      (
        size_t( 1),
        double( 0.5)
      )
    {
    }

    //! @brief pocket constructor
    LigandPocketFitScore::LigandPocketFitScore
    (
      const FragmentComplete &POCKET,
      const linal::Vector3D &START_POSITION
    ) :
      m_BFactor( POCKET.GetSize(), 1.0),
      m_CentroidWeight( 0.20),
      m_Pocket( POCKET),
      m_StartPosition( START_POSITION),
      m_Properties
      (
        size_t( 1),
        descriptor::CheminfoProperty( descriptor::Constants< AtomConformationalInterface, float>( 1.0))
      ),
      m_NumberOfAtomsToAlign
      (
        size_t( 100)
      ),
      m_PropertyWeights
      (
        size_t( 1),
        double( 1.0)
      ),
      m_PropertyCorrelationLengths
      (
        size_t( 1),
        double( 0.5)
      )
    {
    }

    //! virtual copy constructor
    LigandPocketFitScore *LigandPocketFitScore::Clone() const
    {
      return new LigandPocketFitScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LigandPocketFitScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LigandPocketFitScore::GetAlias() const
    {
      static std::string s_name( "LigPocketScore");
      return s_name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the property RMSD
    //! @param MOLECULE - molecule to be fit into pocket
    //! @return the molecule in its new pose relative to the static pocket
    double LigandPocketFitScore::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      return this->operator ()( MOLECULE, m_Pocket);
    }

    //! @brief align two small molecule objects and find the property RMSD
    //! @param MOLECULE - molecule to be fit into pocket
    //! @param POCKET - pocket into which MOLECULE is being geometrically fit
    //! @return the molecule in its new pose relative to the static pocket
    double LigandPocketFitScore::operator()
    (
      const FragmentComplete &MOLECULE,
      const FragmentComplete &POCKET
    ) const
    {
      double score( math::GetHighestBoundedValue< double>());
      if( !m_BFactor.IsEmpty())
      {
        BCL_Assert
        (
          m_BFactor.GetSize() == POCKET.GetHeavyAtomCoordinates().GetSize(),
          "The number of values in the B-factor list must be equal to the number of heavy atoms in the pocket file."
        );
//        score = LigandPocketFitScore::CollisionScore(MOLECULE,POCKET,m_BFactor,false);
        score = LigandPocketFitScore::PropertyCorrelationScore( MOLECULE, POCKET);
      }
      else
      {
        for( size_t index( 0); index < POCKET.GetHeavyAtomCoordinates().GetSize(); ++index)
        {
          m_BFactor.PushBack( 1.0);
        }
//        score = LigandPocketFitScore::CollisionScore(MOLECULE,POCKET,m_BFactor,false);
        score = LigandPocketFitScore::PropertyCorrelationScore( MOLECULE, POCKET);
      }
      BCL_MessageStd( "LigandPocketFitScore_operator: " + util::Format()( score));
      return score;
    }

    double LigandPocketFitScore::CollisionScore
    (
      const FragmentComplete &MOL,
      const FragmentComplete &POCK,
      const storage::Vector< double> &BFACTOR,
      const bool &STARTING_POSE_BIAS
    ) const
    {
      //Build voxel grids for the molecule and the pocket
      VoxelGridAtom voxel_grid_pocket( 4.0), voxel_grid_mol( 4.0);
      voxel_grid_mol.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( MOL.GetAtomsIterator().Begin(), MOL.GetAtomsIterator().End()));
      voxel_grid_pocket.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( POCK.GetAtomsIterator().Begin(), POCK.GetAtomsIterator().End()));
      auto neighbors( voxel_grid_mol.GetNeighborsIn( voxel_grid_pocket, 4.0));

      storage::Vector< linal::Vector< double>> radii( neighbors.GetSize(), linal::Vector< double>( size_t( 4)));
      size_t neighbors_index( 0);
      double score( 0.0);

      //Compute distances between geometric centroid of the molecule
      linal::Vector3D mol_centroid( MOL.GetCenter());

      //recenter molecule in binding pocket center
      BCL_Assert( m_StartPosition.IsDefined(), "No centroid coordinates specified");
      double centroid_distance( linal::Distance( mol_centroid, m_StartPosition));
      score += centroid_distance * m_CentroidWeight;

      for( auto itr_neighbors( neighbors.Begin()), itr_neighbors_end( neighbors.End()); itr_neighbors != itr_neighbors_end; ++itr_neighbors, ++neighbors_index)
      {
        //Get distance between atoms
        const double distance( itr_neighbors->Third());

        //Get b-factor
        double bfactor( BFACTOR( POCK.GetAtomIndex( *itr_neighbors->Second())));

        //Compute covalent radius of molecule atom
        const double mol_atom_type_covalent_radius( BondLengths::GetAverageCovalentRadius( *itr_neighbors->First()));
        const double mol_covalent_radius
        (
          util::IsDefined( mol_atom_type_covalent_radius)
          ? mol_atom_type_covalent_radius
            : itr_neighbors->First()->GetElementType()->GetProperty( ElementTypeData::e_CovalentRadius)
        );
        radii( neighbors_index)( 0) = mol_covalent_radius;

        //Compute csd_vdw radius of molecule atom
        const double mol_csd_vdw_radius( itr_neighbors->First()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
        radii( neighbors_index)( 1) = mol_csd_vdw_radius;

        //Compute covalent radius of pocket atom
        const double pocket_atom_type_covalent_radius( BondLengths::GetAverageCovalentRadius( *itr_neighbors->Second()));
        const double pocket_covalent_radius
        (
          util::IsDefined( pocket_atom_type_covalent_radius)
          ? pocket_atom_type_covalent_radius
            : itr_neighbors->Second()->GetElementType()->GetProperty( ElementTypeData::e_CovalentRadius)
        );
        radii( neighbors_index)( 2) = pocket_atom_type_covalent_radius;

        //Compute csd_vdw radius of pocket atom
        const double pocket_csd_vdw_radius( itr_neighbors->Second()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
        radii( neighbors_index)( 3) = pocket_csd_vdw_radius;

        //Choose the minimum between the sum of the covalent radii of the molecule and pocket atoms or sum of the vdw radii
        const double min_radius
        (
          std::min( mol_covalent_radius + pocket_covalent_radius, mol_csd_vdw_radius + pocket_csd_vdw_radius)
          ? mol_covalent_radius + pocket_covalent_radius
            : mol_csd_vdw_radius + pocket_csd_vdw_radius
        );

        //If distance is less than min_radius then there is a collision
//        double clash_penalty
//        (
//          std::min( distance, min_radius)
//          ? 1.0
//          : 0.0
//        );

        if( mol_csd_vdw_radius + pocket_csd_vdw_radius > distance)
        {
          //Compute BFactor-weighted distance of molecule atom from neighbor pocket atom
          score += ( mol_csd_vdw_radius + pocket_csd_vdw_radius - distance) / bfactor;
        }
      }

      double normalized_score(score/neighbors.GetSize());
      BCL_MessageStd( "LigandPocketFitScore: " + util::Format()( normalized_score));
      return normalized_score;
    }

    double LigandPocketFitScore::PropertyCorrelationScore
    (
      const FragmentComplete &MOL,
      const FragmentComplete &POCK
    ) const
    {
      if( MOL.GetNumberAtoms() == size_t( 0) || POCK.GetNumberAtoms() == size_t( 0))
      {
        return math::GetLowestBoundedValue< double>();
      }

      //Initialize objects
      const size_t n_atoms_a( MOL.GetNumberAtoms());
      const size_t n_atoms_b( POCK.GetNumberAtoms());
      const size_t sum_atoms( n_atoms_a + n_atoms_b);
      util::SiPtrVector< const linal::Vector3D> gridpoints;
      gridpoints.AllocateMemory( sum_atoms);
      gridpoints.Append( MOL.GetAtomCoordinates());
      gridpoints.Append( POCK.GetAtomCoordinates());
      linal::Vector< float> properties_a( n_atoms_a), properties_b( n_atoms_b);
      linal::Vector< float> properties_gridpoints_a( sum_atoms), properties_gridpoints_b( sum_atoms);
      linal::Matrix< float> property_matrix_a( sum_atoms, m_Properties.GetSize()), property_matrix_b( sum_atoms, m_Properties.GetSize());
      math::RunningAverage< float> average_correlation;

      // set properties on molecule and pocket objects
      storage::Vector< descriptor::CheminfoProperty> mol_properties( 0);
      for( size_t i( 0), n_properties( m_Properties.GetSize()); i < n_properties; ++i)
      {
        auto &descrptr( *m_Properties( i));
        descrptr.SetObject( MOL);
        for
        (
          descriptor::Iterator< AtomConformationalInterface> itr_a( descriptor::Type( 1, true, descriptor::Type::e_Symmetric), MOL);
          itr_a.NotAtEnd();
          ++itr_a
        )
        {
          properties_a( itr_a.GetPosition()) = descrptr( itr_a)( 0);
        }
        descrptr.SetObject( POCK);
        for
        (
          descriptor::Iterator< AtomConformationalInterface> itr_b( descriptor::Type( 1, true, descriptor::Type::e_Symmetric), POCK);
          itr_b.NotAtEnd();
          ++itr_b
        )
        {
          properties_b( itr_b.GetPosition()) = descrptr( itr_b)( 0);
        }

        math::RunningAverageSD< float> ave_sd_a, ave_sd_b;

        properties_gridpoints_a = 0.0;
        properties_gridpoints_b = 0.0;
        size_t gridpoint_id( 0);
        const double correlation_length( m_PropertyCorrelationLengths( i));
        const double neg_inv_sqr_correlation_length( -1.0 / ( 2 * math::Sqr( correlation_length)));
        for
        (
          util::SiPtrVector< const linal::Vector3D>::const_iterator
            itr_grid( gridpoints.Begin()), itr_grid_end( gridpoints.End());
          itr_grid != itr_grid_end;
          ++itr_grid, ++gridpoint_id
        )
        {
          // calculate value of point for molecule A
          double field_value_a( 0.0);
          const float *itr_prop_a( properties_a.Begin()), *itr_prop_b( properties_b.Begin());
          for
          (
            iterate::Generic< const AtomConformationalInterface> itr_a( MOL.GetAtomsIterator());
            itr_a.NotAtEnd();
            ++itr_a, ++itr_prop_a
          )
          {
            if( float gaussian_arg = linal::SquareDistance( **itr_grid, itr_a->GetPosition()) * ( neg_inv_sqr_correlation_length) > -7.0)
            {
              field_value_a += *itr_prop_a * std::exp( gaussian_arg);
            }

          }
          properties_gridpoints_a( gridpoint_id) = field_value_a;
          ave_sd_a += field_value_a;

          // calculate value of point for molecule B
          double field_value_b( 0.0);
          for
          (
            iterate::Generic< const AtomConformationalInterface> itr_b( POCK.GetAtomsIterator());
            itr_b.NotAtEnd();
            ++itr_b, ++itr_prop_b //, ++itr_row
          )
          {
            if( float gaussian_arg = linal::SquareDistance( **itr_grid, itr_b->GetPosition()) * ( neg_inv_sqr_correlation_length) > -7.0)
            {
              field_value_b += *itr_prop_b * std::exp( gaussian_arg);
            }

          }
          properties_gridpoints_b( gridpoint_id) = field_value_b;
          ave_sd_b += field_value_b;
        }

        //compute matrix of normalized property values for each molecule
        for( size_t grid_index( 0); grid_index < sum_atoms; ++grid_index)
        {
          if( ave_sd_a.GetStandardDeviation())
          {
            property_matrix_a( grid_index, i) = properties_gridpoints_a( grid_index) / ave_sd_a.GetStandardDeviation();
          }
          if( ave_sd_b.GetStandardDeviation())
          {
            property_matrix_b( grid_index, i) = properties_gridpoints_b( grid_index) / ave_sd_b.GetStandardDeviation();
          }
        }
      }

      //identify grid indices that have the biggest property differences between molecule and pocket
      storage::Vector< storage::Pair< size_t, float> > atom_inclusion_index( sum_atoms);
      for (size_t r(0); r < sum_atoms; ++r)
      {
        atom_inclusion_index(r).First() = r;
        for( size_t c( 0); c < m_Properties.GetSize(); ++c)
        {
          atom_inclusion_index( r).Second() = m_PropertyWeights( c) * ( property_matrix_a( r, c) - property_matrix_b( r, c));
        }
        atom_inclusion_index( r).Second() = math::Absolute( atom_inclusion_index( r).Second());
      }

      //do not include values with largest differences in the subsequent correlation
      atom_inclusion_index.Sort( storage::PairBinaryPredicateSecond< size_t, float>( **math::Comparisons< float>::GetEnums().e_Less));
      size_t min_atoms( std::min( m_NumberOfAtomsToAlign, sum_atoms));

      //compute individual property correlations on the reduced sets of atoms
      linal::Vector< float> reduced_weighted_sum_prop_a( min_atoms), reduced_weighted_sum_prop_b( min_atoms);
      for( size_t c( 0), n_properties( m_Properties.GetSize()); c < n_properties; ++c)
      {
        math::RunningAverageSD< float> ave_sd_a, ave_sd_b;
        for (size_t i(0); i < min_atoms; ++i)
        {
          reduced_weighted_sum_prop_a(i) = property_matrix_a(atom_inclusion_index(i).First(), c);
          ave_sd_a += property_matrix_a(atom_inclusion_index(i).First(), c);
          reduced_weighted_sum_prop_b(i) = property_matrix_b(atom_inclusion_index(i).First(), c);
          ave_sd_b += property_matrix_b(atom_inclusion_index(i).First(), c);
        }

        if( ave_sd_a.GetStandardDeviation() <= 1.0e-8 || ave_sd_b.GetStandardDeviation() <= 1.0e-8)
        {
          if
          (
             !math::EqualWithinAbsoluteTolerance( ave_sd_a.GetAverage(), ave_sd_b.GetAverage(), 1.0e-6)
            || ave_sd_a.GetStandardDeviation() > 1.0e-8
            || ave_sd_b.GetStandardDeviation() > 1.0e-8
          )
          {
            average_correlation.AddWeightedObservation( 0.0, m_PropertyWeights( c));
          }
          continue;
        }
          reduced_weighted_sum_prop_a -= ave_sd_a.GetAverage();
          reduced_weighted_sum_prop_b -= ave_sd_b.GetAverage();

          const double covar
          (
            math::Absolute( linal::ScalarProduct( reduced_weighted_sum_prop_a, reduced_weighted_sum_prop_b)) / double( min_atoms)
          );
          const double correl = covar / ave_sd_a.GetStandardDeviation() / ave_sd_b.GetStandardDeviation();
          average_correlation.AddWeightedObservation( correl, m_PropertyWeights( c));
      }
      return average_correlation.GetAverage();
    }

    //! @brief prepare the class for comparing a conformation
    //! @param MOLECULE the molecule to prepare to compare
    void LigandPocketFitScore::Prepare( const ConformationInterface &MOLECULE) const
    {
      // calculate all properties on the molecule; this avoids potential clashes if multiple threads calculate properties
      // for the same molecule and try to update the cache simultaneously
      for
      (
        storage::Vector< descriptor::CheminfoProperty>::iterator
          itr( m_Properties.Begin()), itr_end( m_Properties.End());
        itr != itr_end;
        ++itr
      )
      {
        ( *itr)->CollectValuesOnEachElementOfObject( MOLECULE);
      }
    }

     //! @brief return parameters for member data that are set up from the labels
     //! @return parameters for member data that are set up from the labels
    io::Serializer LigandPocketFitScore::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Orients a small molecule in a pocket cavity by minimizing geometric overlap of matched atoms");
      member_data.AddOptionalInitializer
      (
        "b factors",
        "modify collision tolerance at each residue based on per-residue b-factor or RMSF",
        io::Serialization::GetAgent( &m_BFactor)
      );
      member_data.AddOptionalInitializer
      (
        "start position",
        "start ligand centered at this position instead of pocket center",
        io::Serialization::GetAgent( &m_StartPosition)
      );
      member_data.AddInitializer
      (
        "aligned atoms",
        "number of atoms scored in correlation",
        io::Serialization::GetAgent( &m_NumberOfAtomsToAlign),
        "100"
      );
      member_data.AddInitializer
      (
        "properties",
        "atom properties to consider, use multiply(Constant(X),property y) for weighting",
        io::Serialization::GetAgent( &m_Properties),
        util::ObjectDataLabel
        (
          "Combine("
          "Atom_SigmaCharge,"
          "Atom_HbondAcceptors,"
          "Atom_HbondDonors,"
          "Atom_Polarizability,"
          "Atom_VDWVolume)"
        )
      );
      member_data.AddInitializer
      (
        "property weights",
        "Weighting to give the properties",
        io::Serialization::GetAgent( &m_PropertyWeights),
        util::ObjectDataLabel( "(5, 3.5, 7, 1.71, 0.714)")
      );
      member_data.AddInitializer
      (
        "property lengths",
        "decay lengths for each property. Longer range interactions, such as H-bonding, should use longer decay lengths",
        io::Serialization::GetAgent( &m_PropertyCorrelationLengths),
        "0.5"
      );
      return member_data;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool LigandPocketFitScore::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

  } // namespace chemistry
} // namespace bcl

