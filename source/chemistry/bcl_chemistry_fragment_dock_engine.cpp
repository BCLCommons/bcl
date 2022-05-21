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
#include "chemistry/bcl_chemistry_fragment_dock_engine.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_pair.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentDockEngine::s_Instance
    (
      GetObjectInstances().AddInstance( new FragmentDockEngine())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentDockEngine::FragmentDockEngine() :
        m_SampleConfs(),
        m_Scorer()
    {
    }

    //! @brief movie constructor
    FragmentDockEngine::FragmentDockEngine
    (
      const std::string &MOVIE_FILENAME
    ) :
      m_SampleConfs(),
      m_Scorer(),
      m_Movie( MOVIE_FILENAME)
    {
    }

    //! @brief sample confs constructor
    FragmentDockEngine::FragmentDockEngine
    (
      const SampleConformations &SAMPLER
    ) :
        m_SampleConfs( SAMPLER),
        m_Scorer()
    {
    }

    //! @brief sample confs and movie constructor
    FragmentDockEngine::FragmentDockEngine
    (
      const SampleConformations &SAMPLER,
      const std::string &MOVIE_FILENAME
    ) :
        m_SampleConfs( SAMPLER),
        m_Scorer(),
        m_Movie( MOVIE_FILENAME)
    {
    }

    //! @brief constructor
    FragmentDockEngine::FragmentDockEngine
    (
      const SampleConformations &SAMPLER,
      const descriptor::CheminfoProperty &SCORER
    ) :
        m_SampleConfs( SAMPLER),
        m_Scorer( SCORER)
    {
    }

    //! virtual copy constructor
    FragmentDockEngine *FragmentDockEngine::Clone() const
    {
      return new FragmentDockEngine( *this);
    }

  /////////////////
  // data access // 
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentDockEngine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &FragmentDockEngine::GetAlias() const
    {
      static std::string s_name( "FragmentDockEngine");
      return s_name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the property RMSD
    //! @param MOLECULE - molecule to be fit into pocket
    //! @param POCKET - pocket into which MOLECULE is being geometrically fit
    //! @return the molecule in its new pose relative to the static pocket
    math::MutateResult< FragmentComplete> FragmentDockEngine::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      // Do something useful

      // After having done something useful, return the mutate
      util::ShPtr< FragmentComplete> mutant( new FragmentComplete( MOLECULE));
      math::MutateResult< FragmentComplete> mutate_result( mutant, *this);

      // Shit out each frame as a movie
      if( !m_Movie.empty())
      {
        io::OFStream out;
        io::File::MustOpenOFStream( out, m_Movie + ".sdf", std::ios::app);
        mutant->WriteMDL( out);
      }
      return mutate_result;
    }

    //! @brief apply rotation/translation transformation matrix to molecule
    //! @param MOLECULE - molecule to be fit into pocket
    //! @return the molecule in its new pose
    FragmentComplete FragmentDockEngine::TransformCoordinates
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      FragmentComplete mol( MOLECULE);
      double prob( random::GetGlobalRandom().Random< double>( 1));

      // rotate molecule randomly with magnitude between 0 and 5 degrees
      if( prob < 0.25)
      {
        math::RotationMatrix3D rot;
        rot.SetRand( 0.08);
        linal::Vector3D mol_centered( mol.GetCenter());
        mol.Translate( -mol_centered);
        mol.Rotate( rot);
        mol.Translate( mol_centered);
      }

      // translate molecule a random distance between 0 and 1 angstroms
      else if( prob < 0.50)
      {
        const double mag( random::GetGlobalRandom().Random( 1.0));
        linal::Vector3D trans( mag, 0.0, 0.0);
        mol.Translate( trans.Rotate( math::RotationMatrix3D().SetRand()));
      }

      // rotate molecule randomly with magnitude between 0 and 180 degrees
      else if( prob < 0.75)
      {
        math::RotationMatrix3D rot;
        rot.SetRand();
        linal::Vector3D mol_centered( mol.GetCenter());
        mol.Translate( -mol_centered);
        mol.Rotate( rot);
        mol.Translate( mol_centered);
      }
      else
      {
        const storage::Vector< float> flip_angles( storage::Vector< float>::Create( 90.0, 180.0));
        linal::Vector3D axis;
        axis.SetRandomTranslation( 1.0);
        math::RotationMatrix3D flip( axis, flip_angles( size_t( ( prob - 0.75) / 0.125)));
        linal::Vector3D mol_centered( mol.GetCenter());
        mol.Translate( -mol_centered);
        mol.Rotate( flip);
        mol.Translate( mol_centered);
      }

      // return the molecule with transformed coordinates
      return mol;
    }

    //! @brief locally sample conformations of input molecule
    //! @param MOLECULE - molecule for which to generate conformers
    //! @return the local conformational ensemble
    FragmentEnsemble FragmentDockEngine::SampleLigandConfsLocal
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      FragmentComplete mol( MOLECULE);
      mol.RemoveH();

      // initialize ensembles
      FragmentEnsemble local_ensemble;
      util::SiPtrVector< const linal::Vector3D> coords( mol.GetAtomCoordinates());

      // get local ensemble
      m_SampleConfs.SetSamplingPreferences( false, false, true, false);
      local_ensemble = m_SampleConfs( mol).First();

      // Realign local ensemble
      for( auto itr( local_ensemble.Begin()), itr_end( local_ensemble.End()); itr != itr_end; ++itr)
      {
        auto rotamer_coords( itr->GetAtomCoordinates());
        auto transform( quality::RMSD::SuperimposeCoordinates( coords, rotamer_coords));
        itr->Transform( transform);
      }
      return local_ensemble;
    } // end SampleLigandConfsLocal

    //! @brief broadly sample conformations of input molecule
    //! @param MOLECULE - molecule for which to generate conformers
    //! @return the global conformational ensemble
    FragmentEnsemble FragmentDockEngine::SampleLigandConfsGlobal
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      FragmentComplete mol( MOLECULE);
      mol.RemoveH();

      // initialize ensembles
      FragmentEnsemble global_ensemble;
      util::SiPtrVector< const linal::Vector3D> coords( mol.GetAtomCoordinates());

      // Get global ensemble
      m_SampleConfs.SetSamplingPreferences( true, true, true, false);
      global_ensemble = m_SampleConfs( mol).First();

      // Realign global ensemble
      for( auto itr( global_ensemble.Begin()), itr_end( global_ensemble.End()); itr != itr_end; ++itr)
      {
        auto rotamer_coords( itr->GetAtomCoordinates());
        auto transform( quality::RMSD::SuperimposeCoordinates( coords, rotamer_coords));
        itr->Transform( transform);
      }
      return global_ensemble;
    } // end SampleLigandConfsGlobal

    //! @brief detect protein-ligand collisions
    //! @param MOLECULE - molecule to be fit into pocket
    //! @param POCKET - pocket into which MOLECULE is being fit
    //! @param BFACTOR - per-atom b factors
    //! @return the score reflecting MOLECULE clashes with POCKET
    double FragmentDockEngine::CollisionScore
    (
      const FragmentComplete &MOLECULE,
      const FragmentEnsemble &POCKET,
      const linal::Vector< double> &BFACTOR
    ) const
    {
      // only concerned with heavy atoms
      FragmentComplete mol( MOLECULE), pocket( *( POCKET.GetMolecules().Begin()));
      mol.RemoveH();
      pocket.RemoveH();

      // get normalized b factors
      linal::Vector< double> bfactors_normed( pocket.GetSize());
      if( BFACTOR.IsDefined() && BFACTOR.GetSize() == POCKET.GetMolecules().Begin()->GetSize())
      {
        math::RunningAverageSD< double> b_ave_sd;
        math::RunningMinMax< double> b_min_max;
        for
        (
            auto b_itr( BFACTOR.Begin()), b_itr_end( BFACTOR.End());
            b_itr != b_itr_end;
            ++b_itr
        )
        {
          b_ave_sd += *b_itr;
          b_min_max += *b_itr;
        }
        size_t b_index( 0);
        for
        (
            auto b_itr( BFACTOR.Begin()), b_itr_end( BFACTOR.End());
            b_itr != b_itr_end;
            ++b_itr, ++b_index
        )
        {
          double ave( b_ave_sd.GetAverage()), min( b_min_max.GetMin());
          bfactors_normed( b_index) = ( *b_itr - ave - ( min - ave)) / b_ave_sd.GetStandardDeviation();
        }
      }
      else
      {
        bfactors_normed = 1.0;
      }

      // Recenter molecule centroid to starting point
      if( m_StartPosition.IsDefined())
      {
        mol.Translate( -linal::Vector3D( mol.GetCenter()));
        mol.Translate( m_StartPosition);
      }

      // Build voxel grids for the molecule and the pocket
      VoxelGridAtom voxel_grid_pocket( 4.0), voxel_grid_mol( 4.0);
      voxel_grid_mol.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( mol.GetAtomsIterator().Begin(), mol.GetAtomsIterator().End()));
      voxel_grid_pocket.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( pocket.GetAtomsIterator().Begin(), pocket.GetAtomsIterator().End()));
      auto neighbors( voxel_grid_mol.GetNeighborsIn( voxel_grid_pocket, 4.0));

      // Prepare to collect radii from each molecule atom and receptor neighbor atom
      storage::Vector< linal::Vector< double>> radii( neighbors.GetSize(), linal::Vector< double>( size_t( 4)));
      size_t neighbors_index( 0);
      double score( 0.0);

      // Look for clashes
      for
      (
          auto itr_neighbors( neighbors.Begin()), itr_neighbors_end( neighbors.End());
          itr_neighbors != itr_neighbors_end; ++itr_neighbors, ++neighbors_index
      )
      {
        //Get distance between atoms
        const double distance( itr_neighbors->Third());

        //Get b-factor
        double b_weight( bfactors_normed( pocket.GetAtomIndex( *itr_neighbors->Second())));

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
        const double mol_csd_vdw_radius( itr_neighbors->First()->GetAtomType()->GetElementType()->GetProperty( ElementTypeData::e_DaltonVdwRadius));
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
        const double pocket_csd_vdw_radius( itr_neighbors->Second()->GetAtomType()->GetElementType()->GetProperty( ElementTypeData::e_DaltonVdwRadius));
        radii( neighbors_index)( 3) = pocket_csd_vdw_radius;

        //Choose the maximum between the sum of the covalent radii of the molecule and pocket atoms or sum of the vdw radii
        const double max_radius
        (
          std::max( mol_covalent_radius + pocket_covalent_radius, mol_csd_vdw_radius + pocket_csd_vdw_radius)
          ? mol_covalent_radius + pocket_covalent_radius
            : mol_csd_vdw_radius + pocket_csd_vdw_radius
        );

        //If distance is less than max_radius then there is a collision
        //       double clash_penalty
        //       (
        //         std::min( distance, max_radius)
        //         ? 1.0
        //         : 0.0
        //       );

        if( distance < max_radius)
        {
          //Compute BFactor-weighted distance of molecule atom from neighbor pocket atom
          score += ( max_radius - distance) * b_weight; // low b_weight means high bfactor
        }
      }

      // normalize by the number of neighbors
      double normalized_score(score/neighbors.GetSize());
      BCL_MessageStd( "FragmentDockEngine Collision Score: " + util::Format()( normalized_score));
      return normalized_score;
    }

    //! @brief score protein-ligand interactions
    //! @param MOLECULE - molecule to be fit into pocket
    //! @param POCKET - pocket into which MOLECULE is being fit
    //! @param MODEL - the descriptor / model with which to compute the score
    //! @return the score reflecting MOLECULE's interaction with POCKET
    double FragmentDockEngine::InteractionScore
    (
      const FragmentComplete &MOLECULE,
      const FragmentComplete &POCKET,
      const descriptor::CheminfoProperty &MODEL
    ) const
    {
      // enforce saturation of the molecule
      FragmentComplete saturated_mol( MOLECULE);
      saturated_mol.SaturateWithH();

      // score with model
      linal::Vector< double> res( m_Scorer->SumOverObject( saturated_mol));
      BCL_MessageStd( "FragmentDockEngine Interaction Score: " + util::Format()( res.Sum() / res.GetSize()));
      return res.Sum() / res.GetSize();
    }

     //! @brief return parameters for member data that are set up from the labels
     //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentDockEngine::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Orients a small molecule in a pocket cavity by minimizing geometric overlap of matched atoms");
      parameters.AddInitializer
      (
        "sampler",
        "sample configurational space across dihedral bins",
        io::Serialization::GetAgent( &m_SampleConfs),
        util::ObjectDataLabel
        ( "("
          "conformation_comparer=bcl::chemistry::ConformationComparisonInterface,tolerance=0.0,"
          "generate_3D=0,relative_random_dihedral_change_weight=0.0,cluster=false,"
          "max_iterations=250, max_conformations=100"
          ")"
        )
      );
      parameters.AddInitializer
      (
        "score_function",
        "cheminfo property to be used as score function",
        io::Serialization::GetAgent( &m_Scorer),
        util::ObjectDataLabel
        ( "("
            "Limit(PredictionMean(storage="
            "File("
            "directory=/dors/meilerlab/home/brownbp1/hybrid_qsar/qsar_dG_pred/pdb_bind/CASF-2016/pdbbind_v2016/"
            "refined_general_combined_minus_casf16/groups/3DAPairRealSpaceAsym_050.final.080.fast.interp_gauss_050/"
            "models/combined.rand.3DAPairRealSpaceAsym_050.final.080.fast.interp_gauss_050.2x512-32_005_025_005/"
            ",prefix=model)"
            "),max=10.0,min=3.0)"
          ")"
        )
      );
      parameters.AddOptionalInitializer
      (
        "b_factors",
        "modify collision tolerance at each residue based on per-residue b-factor or RMSF",
        io::Serialization::GetAgent( &m_BFactors)
      );
      parameters.AddOptionalInitializer
      (
        "start_position",
        "start ligand centered at this position instead of at initial coordinates",
        io::Serialization::GetAgent( &m_StartPosition)
      );
      parameters.AddOptionalInitializer
      (
        "output_movie",
        "output each step of the attempted orientation of the ligand into the pocket - WARNING - large file size easily accumulates",
        io::Serialization::GetAgent( &m_Movie)
      );
      return parameters;
    }
    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool FragmentDockEngine::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return math::MutateInterface< FragmentComplete>::ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }

  } // namespace chemistry
} // namespace bcl

