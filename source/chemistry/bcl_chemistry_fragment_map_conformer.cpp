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
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"

// include bcl headers
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_rmsd.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_fragment_split_rings_with_unsaturated_substituents.h"
#include "chemistry/bcl_chemistry_fragment_stochastic_pose_optimizer.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_molecule_druglike.h"
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "mm/bcl_mm_rdkit_energy_minimize_mmff94.h"
#include "pdb/bcl_pdb_factory.h"
#include "quality/bcl_quality_rmsd.h"
#include "random/bcl_random_distribution_interface.h"
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

#undef AddAtom
#undef RemoveAtom
#undef ATOMS

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // initialize static
    sched::Mutex &FragmentMapConformer::GetMutex()
    {
      static sched::Mutex s_Mutex;
      return s_Mutex;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMapConformer::FragmentMapConformer() :
      m_DrugLikenessType( "IsConstitutionDruglike"),
      m_ResolveClashes( false),
      m_VDWClashCutoff( 5.0),
      m_Corina( false),
      m_MoveableIndices( storage::Vector< size_t>()),
      m_ChooseBestAlignedConf( false),
      m_FixGeometry( true),
      m_AdjacentNeighbors( size_t( 1)),
      m_MapSubgraphRingAtoms( true),
      m_SkipConfGen( false),
      m_GeoOpt( false)
    {
    }

    //! @brief druglikeness constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMapConformer::FragmentMapConformer
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const bool CORINA_CONFS,
      const storage::Vector< size_t> &MOVEABLE_INDICES
    ) :
      m_DrugLikenessType( "IsConstitutionDruglike"),
      m_ResolveClashes( false),
      m_VDWClashCutoff( 5.0),
      m_Corina( CORINA_CONFS),
      m_MoveableIndices( MOVEABLE_INDICES),
      m_ChooseBestAlignedConf( false),
      m_FixGeometry( true),
      m_AdjacentNeighbors( size_t( 1)),
      m_MapSubgraphRingAtoms( true),
      m_SkipConfGen( false),
      m_GeoOpt( false)
    {
    }

    //! @brief pose resolver constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param MDL property label containing path to protein binding pocket PDB file
    //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
    //! @param PROPERTY_SCORER enables pose-dependent optimization of score
    //! @param RESOLVE_CLASHES true will seek to resolve clashes with pocket by changing the ligand conformer
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    FragmentMapConformer::FragmentMapConformer
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const std::string &MDL,
      const std::string &BINDING_POCKET_FILENAME,
      const descriptor::CheminfoProperty &PROPERTY_SCORER,
      const bool RESOLVE_CLASHES,
      const storage::Vector< float> &BFACTORS,
      const bool CORINA_CONFS,
      const storage::Vector< size_t> &MOVEABLE_INDICES,
      const bool CHOOSE_BEST_ALIGNED_CONF,
      const bool FIX_GEOMETRY,
      const size_t ADJACENT_NBRS,
      const bool MAP_SUBGRAPH_RINGS,
      const bool SKIP_CONFGEN,
      const bool GEO_OPT
    ) :
      m_DrugLikenessType( DRUG_LIKENESS_TYPE),
      m_MDL( MDL),
      m_BindingPocketFilename( util::Strip( BINDING_POCKET_FILENAME, " \t\n\r")),
      m_PropertyScorer( PROPERTY_SCORER),
      m_ResolveClashes( RESOLVE_CLASHES),
      m_BFactors( BFACTORS),
      m_VDWClashCutoff( 5.0),
      m_Corina( CORINA_CONFS),
      m_MoveableIndices( MOVEABLE_INDICES),
      m_ChooseBestAlignedConf( CHOOSE_BEST_ALIGNED_CONF),
      m_FixGeometry( FIX_GEOMETRY),
      m_AdjacentNeighbors( ADJACENT_NBRS),
      m_MapSubgraphRingAtoms( MAP_SUBGRAPH_RINGS),
      m_SkipConfGen( SKIP_CONFGEN),
      m_GeoOpt( GEO_OPT)
    {
    }

    //! @brief clone constructor
    FragmentMapConformer *FragmentMapConformer::Clone() const
    {
      return new FragmentMapConformer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMapConformer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the MDL property associated with this object
    //! @return the MDL string
    const std::string &FragmentMapConformer::GetMDL() const
    {
      return m_MDL;
    }

    //! @brief return the pocket filename associated with this object
    //! @return the PDB filename
    const std::string &FragmentMapConformer::GetPocketFilename() const
    {
      return m_BindingPocketFilename;
    }

    //! brief return the bfactors associated with this object
    //! return the bfactors for each atom in the pocket
    const storage::Vector< float> &FragmentMapConformer::GetBFactors() const
    {
      return m_BFactors;
    }

    //! @brief return whether to skip conformer generation
    bool FragmentMapConformer::GetSkipConfGen() const
    {
      return m_SkipConfGen;
    }

    //! @brief return whether to run MMFF94s minimization
    bool FragmentMapConformer::GetGeoOpt() const
    {
      return m_GeoOpt;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set skip conformer generation
    void FragmentMapConformer::SetSkipConfGen( const bool SKIP)
    {
      m_SkipConfGen = SKIP;
    }

    //! @brief set geometry optimization
    void FragmentMapConformer::SetGeoOpt( const bool OPT)
    {
      m_GeoOpt = OPT;
    }

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @param REFERENCE_MOL the scaffold molecule for the substructure-based alignment of the new 3D conformer
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @return util::ShPtr< FragmentComplete> of a new 3D molecule following the mutation
    util::ShPtr< FragmentComplete> FragmentMapConformer::Clean
    (
      const AtomVector< AtomComplete> &ATOM_VEC,
      const FragmentComplete &REFERENCE_MOL,
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &SKIP_NEUT
    ) const
    {
      // clean atoms
      AtomVector< AtomComplete> new_mol_atoms_noh
      (
        CleanAtoms( ATOM_VEC, DRUG_LIKENESS_TYPE, SKIP_NEUT)
      );

      // exit if we failed the atom cleaning
      if( !new_mol_atoms_noh.GetSize())
      {
        return util::ShPtr< FragmentComplete>();
      }

      // isomorphism search is too slow if we do not remove hydrogen atoms
      FragmentComplete new_mol( new_mol_atoms_noh, "");
      HydrogensHandler::Remove( new_mol_atoms_noh);
      FragmentComplete new_mol_noh( new_mol_atoms_noh, "");

      // remove hydrogen atoms from reference, too
      AtomVector< AtomComplete> ref_mol_atoms_noh( REFERENCE_MOL.GetAtomVector());
      HydrogensHandler::Remove( ref_mol_atoms_noh);
      FragmentComplete ref_mol_noh( ref_mol_atoms_noh, "");

      // make 3d conformer
      if( !m_MoveableIndices.GetSize())
      {
        // obtain moveable indices from substructure comparison to reference
        BCL_MessageStd( "Getting atom indices for conformer sampling...");
        storage::Pair< storage::Set< size_t>, storage::Set< size_t>> moveable_indices( MapAtoms( ref_mol_noh, new_mol_noh));

        //convert unique atoms in set to vector
        m_MoveableIndices = storage::Vector< size_t>( moveable_indices.First().Begin(), moveable_indices.First().End());
        m_RestrainedIndices = storage::Vector< size_t>( moveable_indices.Second().Begin(), moveable_indices.Second().End());
      }
      // We should have moveable indices now either from MapAtoms or from input;
      // If we do not, then our molecule must not have been mutated
      if( m_MoveableIndices.GetSize())
      {
        // TODO this only works reliably if hydrogen atoms are all indexed higher than heavy atoms;
        // bake that into the protocol by force removing hydrogen atoms from molecules in CleanAtoms, maybe?
        // hmm... think about this.
        new_mol.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", m_MoveableIndices);
        BCL_MessageStd( "SampleByParts with the following atom indices: " + util::Format()( new_mol.GetMDLProperty( "SampleByParts")));
      }
      // no moveable indices, just return existing conformer, after fixing any issues with bond lengths or H
      if( m_MoveableIndices.IsEmpty())
      {
        BCL_MessageStd( "FragmentMapConformer::Clean no mobile atom indices");
        new_mol.StandardizeBondLengths();
        new_mol.UpdateH();
      }

      // if we have moveable indices, then consider optimizing the conformer
      // TODO refactor - require that FragmentMapConformer functions that need to generate conformers
      // are given a SampleConformations object as an argument. Then for the drug design framework make
      // the rotamer library and sample confs objects static member data in the FragmentMutateInterface
      // base class. Then options can be controlled at the derived class level through setters.
      // TODO 2: double check that the above-described refactor will actually work if this is run in parallel,
      // or if I should just stick with this hard-coded confgen stuff
      static RotamerLibraryFile rotamer_library;
      util::ShPtr< FragmentComplete> gen_mol_3d_sp;

      // conformer-dependent score-based refinement
      if( m_PropertyScorer.IsDefined() && !m_ChooseBestAlignedConf && !m_SkipConfGen)
      {
        // make conformers
        static FragmentMakeConformers make_confs
        (
          rotamer_library,
          "",
          0.0,
          true,
          2000,
          0.50,
          false,
          false
        );
        util::ShPtr< FragmentEnsemble> conf_ens_opti;
        util::ShPtr< FragmentEnsemble> conf_ens( make_confs.MakeConformers( new_mol));

        // make sure we have conformers, then align them to common substructure
        if( conf_ens->GetSize())
        {
          for
          (
              auto conf_itr( conf_ens->Begin()), conf_itr_end( conf_ens->End());
              conf_itr != conf_itr_end;
              ++conf_itr
          )
          {
            // set MDL property in our molecule
            if( m_MDL.size() && m_BindingPocketFilename.size())
            {
              conf_itr->GetStoredPropertiesNonConst().SetMDLProperty( m_MDL, m_BindingPocketFilename);
            }

            // perform substructure-based alignment to reference
            if( REFERENCE_MOL.GetSize())
            {
              m_Aligner.AlignToScaffold
              (
                *conf_itr,
                REFERENCE_MOL
              );
            }
          }
        }
        // return if there are no 3D confs
        else
        {
          return util::ShPtr< FragmentComplete>();
        }

        // make a pose optimizer object
        // TODO: never initialize static variables with member data like this
        // refactor to make the pose optimization object part of member data, build in class construction
        static FragmentStochasticPoseOptimizer pose_optimizer
        (
          m_PropertyScorer,
          m_BFactors,
          m_MDL,
          m_BindingPocketFilename,
          size_t( 100),
          size_t( 1),
          size_t( 100),
          float( 1.0),
          opti::e_LargerIsBetter,
          float( 5.0)
        );
        gen_mol_3d_sp =
            pose_optimizer.RunMCMDocker( *conf_ens, true, false);

        if( gen_mol_3d_sp.IsDefined() && gen_mol_3d_sp->GetSize() && !gen_mol_3d_sp->HasBadGeometry())
        {
          gen_mol_3d_sp->GetStoredPropertiesNonConst().SetMDLProperty( m_MDL, m_BindingPocketFilename);
          BCL_MessageStd( "Returning optimized pose!");
          return gen_mol_3d_sp;
        }
        else
        {
          return util::ShPtr< FragmentComplete>();
        }
      }
      // no optimization, just get a conformer
      else if( !m_SkipConfGen)
      {
        if( m_ChooseBestAlignedConf)
        {
          util::ShPtr< FragmentEnsemble> conf_ens;
          static FragmentMakeConformers make_confs
          (
            rotamer_library,
            "",
            0.25,  // was 0.25
            true,
            2000,   // was 2000
            0.1,
            false,  // was true
            false,
            false
          );
          if( m_MoveableIndices.IsEmpty())
          {
            gen_mol_3d_sp = util::CloneToShPtr( new_mol); // use input conformer
          }
          else
          {
            conf_ens = make_confs.MakeConformers( new_mol); // make a single new conformer
          }

          // bail if no conformers could be generated and if input conformer is undefined
          if( !conf_ens.IsDefined() && !gen_mol_3d_sp.IsDefined())
          {
            BCL_MessageStd( "Could not generate a valid conformer ensemble, returning null");
            return util::ShPtr< FragmentComplete>();
          }

          if( conf_ens->GetSize() && REFERENCE_MOL.GetSize())
          {
            storage::Vector< storage::Pair< bool, float> > aets
            (
              m_Aligner.AlignEnsembleToScaffold
              (
                *conf_ens,
                REFERENCE_MOL
              )
            );

            // sort by alignment score
            conf_ens->Sort( "PropertyFieldDistance");

            // output best by MolAlign score
            size_t best_conf( ( conf_ens->GetMolecules().FirstElement().GetStoredProperties().GetMDLPropertyAsVector( "ConformerIndex"))( 0));
            if( aets( best_conf).First())
            {
              gen_mol_3d_sp = util::CloneToShPtr( conf_ens->GetMolecules().FirstElement());
            }
            else
            {
              gen_mol_3d_sp = util::CloneToShPtr( new_mol);
            }
          }
          // score the ensemble without re-alignment and choose best by score
          else if( conf_ens->GetSize() && m_PropertyScorer.IsDefined())
          {
            for( auto conf_itr( conf_ens->Begin()), conf_itr_end( conf_ens->End()); conf_itr != conf_itr_end; ++conf_itr)
            {
              linal::Vector<float> alignment_score( m_PropertyScorer->SumOverObject( *conf_itr));
              conf_itr->GetStoredPropertiesNonConst().SetMDLProperty( m_PropertyScorer->GetAlias(), alignment_score);
            }
            // output best by MolAlign score
            conf_ens->Sort( m_PropertyScorer->GetAlias());
            gen_mol_3d_sp = util::CloneToShPtr( conf_ens->GetMolecules().FirstElement());
          }
        }
        // just make a single conformer and choose by conf score
        else
        {
          BCL_MessageStd( "Generating single conformer");
          static FragmentMakeConformers make_one_conf
          (
            rotamer_library,
            1000,
            1.0,
            true,
            m_Corina
          );

          // I hate everything about this class, but I wrote it and continue to grow its filth, so it's my fault - BPB
          // Until this gets refactored, it'll be used as a teaching tool. This is what happens when a hack grows too large.
          GetMutex().Lock();
          m_MoveableIndices.IsEmpty() ?
              gen_mol_3d_sp = util::CloneToShPtr( new_mol) : // use input conformer
              gen_mol_3d_sp = make_one_conf.MakeConformer( new_mol); // make a single new conformer
          GetMutex().Unlock();
        }
      }
      else
      {
        gen_mol_3d_sp = util::CloneToShPtr( new_mol);
        BCL_MessageStd( "FragmentMapConformer::Clean skipping conformer generation");
      }

      // Stick MMFF94s geometry optimization HERE (because why the fuck not? it all needs refactoring anyway D: )
      if( m_GeoOpt)
      {
        // DEBUG
        const FragmentComplete &mol( *gen_mol_3d_sp);
        io::OFStream debug_out;
        io::File::MustOpenOFStream( debug_out, "pre_opt.sdf");
        mol.WriteMDL( debug_out);
        io::File::CloseClearFStream( debug_out);

        // setup our restraint terms in the minimizer
        static mm::RdkitEnergyMinimizeMmff94 minimizer;
        minimizer.SetPositionalRestraintAtoms( m_RestrainedIndices);
        storage::Vector< double> restraint_forces( m_RestrainedIndices.GetSize(), 1000.0);
        minimizer.SetRestraintForce( restraint_forces);
        storage::Vector< double> max_displacement( m_RestrainedIndices.GetSize(), 0.0);
        minimizer.SetMaxUnrestrainedDisplacement( max_displacement);

        // run minimization
        storage::Triplet< FragmentComplete, int, double> opt_result( minimizer.OptimizeGeometry( mol));
        io::File::MustOpenOFStream( debug_out, "post_opt.sdf");
        opt_result.First().WriteMDL( debug_out);
        io::File::CloseClearFStream( debug_out);
        if( opt_result.Second() == 0)
        {
          BCL_MessageStd( "Minimization was successful! Conformer potential energy = " + util::Format()( opt_result.Third()) + " kcal/mol");
          gen_mol_3d_sp = util::CloneToShPtr( opt_result.First());
        }
        else if( opt_result.Second() == 1)
        {
          BCL_MessageStd( "Minimization proceeded but did not converge! Conformer potential energy = " + util::Format()( opt_result.Third()) + " kcal/mol");
          gen_mol_3d_sp = util::CloneToShPtr( opt_result.First());
        }
        else if( opt_result.Second() == -1)
        {
          BCL_MessageStd( "Minimization failed due to missing parameters!");
        }
      }

      // check if defined with atoms and reasonable geometry
      if( gen_mol_3d_sp.IsDefined() && gen_mol_3d_sp->GetSize() && !gen_mol_3d_sp->HasBadGeometry())
      {
        bool good_conf( true);
        // filter molecules with strained 3D conformers
//        GetMutex().Lock();
//        (
//            descriptor::GetCheminfoProperties().calc_MoleculeVdwScore->SumOverObject( *gen_mol_3d_sp)( 0) < m_VDWClashCutoff ?
//                true :
//                false
//        ); // TODO this is not doing anything
//        GetMutex().Unlock();
        if( !good_conf)
        {
          return util::ShPtr< FragmentComplete>();
        }
        // set MDL property in our molecule
        if( m_MDL.size() && m_BindingPocketFilename.size())
        {
          gen_mol_3d_sp->GetStoredPropertiesNonConst().SetMDLProperty( m_MDL, m_BindingPocketFilename);
        }

        // perform substructure-based alignment
        if( REFERENCE_MOL.GetSize())
        {
        bool success_aln
        (
          m_Aligner.AlignToScaffold
           (
            *gen_mol_3d_sp,
            REFERENCE_MOL,
            m_RestrainedIndices
           )
          );
          BCL_MessageStd( "Final realignment of scaffold successful? " + util::Format()( success_aln ? "Yes" : "No"));
        }

        // return final molecule
        return gen_mol_3d_sp;
      }
      else if( !gen_mol_3d_sp.IsDefined())
      {
        BCL_MessageStd("Molecule cleaning failed to generate a valid 3D conformer!")
        BCL_MessageStd("Defined: " + util::Format()( gen_mol_3d_sp.IsDefined() ? "true" : "false"));
      }
      else
      {
        BCL_MessageStd("Molecule cleaning failed to generate a valid 3D conformer!")
        BCL_MessageStd("Defined: " + util::Format()( gen_mol_3d_sp.IsDefined() ? "true" : "false"));
        BCL_MessageStd("Has good geometry: " + util::Format()( gen_mol_3d_sp->HasBadGeometry() ? "false" : "true" ));
        BCL_MessageStd("Final molecule size: " + util::Format()( gen_mol_3d_sp->GetSize()));
        BCL_MessageStd("Returning null...")
      }
      return util::ShPtr< FragmentComplete>();
    }

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @return AtomVector< AtomComplete> of a new set of clean atoms following the mutation
    AtomVector< AtomComplete> FragmentMapConformer::CleanAtoms
    (
      const AtomVector< AtomComplete> &ATOM_VEC,
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &SKIP_NEUT,
      const bool &SKIP_SATURATE_H
    )
    {
      // make sure we have atoms
      if( !ATOM_VEC.GetSize())
      {
        BCL_MessageStd( "Failed clean!");
        return AtomVector< AtomComplete>();
      }

      // work from here
      AtomVector< AtomComplete> atom_vec( ATOM_VEC);

      // clean at atom_vector level
      AtomsCompleteStandardizer standardizer( atom_vec, "", true);
      if( !SKIP_NEUT)
      {
        standardizer.TryNeutralize( atom_vec, sdf::e_BondOrderAndpHAromaticityLossOk);
      }
      standardizer.SetConjugationOfBondTypes( atom_vec);

      // add isometry info, which may be different for a subgraph
      BondIsometryHandler::AddIsometryInformation( atom_vec, true);

      // add stereocenter information
      StereocentersHandler::AddChiralityFromConformation( atom_vec);

      // split out small pieces that may have broken off during atom swap / bond break / etc.
      FragmentComplete new_mol( atom_vec, "");
      // if( new_mol.HasNonGasteigerAtomTypes()) gives ICE error on windows compiler, so just inlining it here
      for( auto itr( new_mol.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        if( !itr->GetAtomType()->IsGasteigerAtomType())
        {
          return AtomVector< AtomComplete>();
        }
      }
      FragmentSplitLargestComponent splitter;
      FragmentEnsemble largest_component( splitter( new_mol));
      new_mol = largest_component.GetMolecules().FirstElement();
      if( !SKIP_SATURATE_H)
      {
        new_mol.SaturateWithH();
      }

      // check drug-likeness
      if( DRUG_LIKENESS_TYPE.size() && DRUG_LIKENESS_TYPE != "None")
      {
        GetMutex().Lock();
        static descriptor::CheminfoProperty drug_likeness_filter( DRUG_LIKENESS_TYPE);
        bool druglike( drug_likeness_filter->SumOverObject( new_mol)( 0));
        GetMutex().Unlock();
        if( druglike)
        {
          BCL_MessageStd( "Pass drug-likeness filter!");
        }
        else
        {
          BCL_MessageStd( "Fail drug-likeness filter!");
          return AtomVector< AtomComplete>();
        }
      }
      else
      {
        BCL_MessageVrb( "Skipping drug-likeness filter!");
      }
      return new_mol.GetAtomVector();
    }

    //! @brief preserve conformational information from starting molecule in new molecule
    //! @param STARTING_MOL the starting molecule
    //! @param NEW_MOL the resulting molecule post-design
    //! @param ATOM_COMPARISON the atom comaprison type to use for subgraph isomorphism
    //! @param BOND_COMPARISON the bond comaprison type to use for subgraph isomorphism
    //! @param COMPLEMENT if true, the mapped atom indices returned are the subgraph complement,
    //! if false then the returned indices are of the common subgraph
    //! @return the NEW_MOL indices mapped to STARTING_MOL (first) and the unmapped indices (second)
    storage::Pair< storage::Set< size_t>, storage::Set< size_t>> FragmentMapConformer::MapAtoms
        (
          const FragmentComplete &STARTING_MOL,
          const FragmentComplete &NEW_MOL,
          const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON,
          const ConfigurationalBondTypeData::Data &BOND_COMPARISON,
          const bool &COMPLEMENT
        ) const
        {
          // prepare dehydrogenated copy
          ConformationGraphConverter graph_maker
          (
            ATOM_COMPARISON,
            BOND_COMPARISON,
            false
          );

      // set up the isomorphisms with pointers to the graphs
      graph::CommonSubgraphIsomorphism< size_t, size_t> isomorphism( graph::CommonSubgraphIsomorphismBase::SolutionType::e_GreedyUnconnected);
      isomorphism.SetGraphA( graph_maker( STARTING_MOL));
      isomorphism.SetGraphB( graph_maker( NEW_MOL));

      // get the isomorphism - allocate array based on maximum possible atom matches without considering connectivity
      isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds());

      // get graph of the new mol
      graph::Subgraph< size_t, size_t> subgraph_b
      (
        util::OwnPtr< const graph::ConstGraph< size_t, size_t> >( &isomorphism.GetGraphB(), false),
        isomorphism.GetIsomorphism().GetMappedValues()
      );

      // get inverted subgraph of the new mol
      graph::Subgraph< size_t, size_t> complement_subgraph_b;
      if( COMPLEMENT)
      {
        complement_subgraph_b = subgraph_b.GetComplement();
      }
      else
      {
        complement_subgraph_b = subgraph_b;
      }

      // these will be our moveable indices
      storage::Set< size_t> subgraph_b_vertices, remaining_vertices;
      for( size_t a_i( 0), a_sz( NEW_MOL.GetAtomVector().GetSize()); a_i < a_sz; ++a_i)
      {
        remaining_vertices.Insert( a_i);
      }

      // allow conformational sampling of the atoms in the new mol that are not in the original mol
      subgraph_b_vertices.InsertElements( complement_subgraph_b.GetVertexIndices().Begin(), complement_subgraph_b.GetVertexIndices().End());
      remaining_vertices.EraseKeys( complement_subgraph_b.GetVertexIndices().Begin(), complement_subgraph_b.GetVertexIndices().End());

      // get ring components
      if( m_MapSubgraphRingAtoms)
      {
        ConformationGraphConverter::t_AtomGraph molecule_graph( ConformationGraphConverter::CreateGraphWithAtoms( NEW_MOL));
        FragmentSplitRings splitter( true, size_t( 3));
        storage::List< storage::Vector< size_t> > ring_components( splitter.GetComponentVertices( NEW_MOL, molecule_graph));
        storage::Vector< storage::Vector< size_t>> ring_components_vec( ring_components.Begin(), ring_components.End());

        // decide which rings matter, as defined by rings that contain perturbed atoms
        storage::Set< size_t> important_ring_atoms( MapSubgraphRingAtoms( complement_subgraph_b, ring_components));

        // add the atoms from the important rings to the moveable indices;
        // allows re-sampling perturbed rings, very important
        subgraph_b_vertices.InsertElements( important_ring_atoms.Begin(), important_ring_atoms.End());
        remaining_vertices.EraseKeys( important_ring_atoms.Begin(), important_ring_atoms.End());
      }

      // allow conformational sampling of atoms involved in bad geometry
      if( m_FixGeometry)
      {
        storage::Vector< size_t> bad_geo_atoms( NEW_MOL.GetAtomsWithBadGeometry());
        subgraph_b_vertices.InsertElements( bad_geo_atoms.Begin(), bad_geo_atoms.End());
        remaining_vertices.EraseKeys( bad_geo_atoms.Begin(), bad_geo_atoms.End());
      }

      // include atoms adjacent to the new mol edges; MUST do this last since we modify the complement subgraph
      if( m_AdjacentNeighbors)
      {
        // add all the adjacent edges to our unique subgraph vertices
        storage::Set< size_t> all_adjacent_indices( MapSubgraphAdjacentAtoms( complement_subgraph_b, m_AdjacentNeighbors));
        subgraph_b_vertices.InsertElements( all_adjacent_indices.Begin(), all_adjacent_indices.End());
        remaining_vertices.EraseKeys( all_adjacent_indices.Begin(), all_adjacent_indices.End());
      }

      // end
      return storage::Pair< storage::Set< size_t>, storage::Set< size_t>>( subgraph_b_vertices, remaining_vertices);
    }

    //! @brief preserve conformational information from starting molecule in new molecule
    //! @param STARTING_MOL the starting molecule
    //! @param NEW_MOL the resulting molecule post-design
    //! @param ATOM_COMPARISON the atom comaprison type to use for subgraph isomorphism
    //! @param BOND_COMPARISON the bond comaprison type to use for subgraph isomorphism
    //! @param COMPLEMENT if true, the mapped atom indices returned are the subgraph complement,
    //! if false then the returned indices are of the common subgraph
    //! @return the NEW_MOL indices mapped to STARTING_MOL
    storage::Set< size_t> FragmentMapConformer::MapSubgraphRingAtoms
    (
      const FragmentComplete &STARTING_MOL,
      const FragmentComplete &NEW_MOL,
      const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON,
      const ConfigurationalBondTypeData::Data &BOND_COMPARISON,
      const bool &COMPLEMENT
    ) const
    {
      // prepare dehydrogenated copy
      ConformationGraphConverter graph_maker
      (
        ATOM_COMPARISON,
        BOND_COMPARISON,
        false
      );

      // set up the isomorphisms with pointers to the graphs
      graph::CommonSubgraphIsomorphism< size_t, size_t> isomorphism( graph::CommonSubgraphIsomorphismBase::SolutionType::e_GreedyUnconnected);
      isomorphism.SetGraphA( graph_maker( STARTING_MOL));
      isomorphism.SetGraphB( graph_maker( NEW_MOL));

      // get the isomorphism - allocate array based on maximum possible atom matches without considering connectivity
      isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds());

      // get graph of the new mol
      graph::Subgraph< size_t, size_t> subgraph_b
      (
        util::OwnPtr< const graph::ConstGraph< size_t, size_t> >( &isomorphism.GetGraphB(), false),
        isomorphism.GetIsomorphism().GetMappedValues()
      );

      // get inverted subgraph of the new mol
      graph::Subgraph< size_t, size_t> complement_subgraph_b;
      if( COMPLEMENT)
      {
        complement_subgraph_b = subgraph_b.GetComplement();
      }
      else
      {
        complement_subgraph_b = subgraph_b;
      }

      // get ring components
      ConformationGraphConverter::t_AtomGraph molecule_graph( ConformationGraphConverter::CreateGraphWithAtoms( NEW_MOL));
      FragmentSplitRings splitter( true, size_t( 3));
      storage::List< storage::Vector< size_t> > ring_components( splitter.GetComponentVertices( NEW_MOL, molecule_graph));
      return MapSubgraphRingAtoms( complement_subgraph_b, ring_components);
    }

    //! @brief preserve conformational information from starting molecule in new molecule
    //! @param STARTING_MOL the starting molecule
    //! @param NEW_MOL the resulting molecule post-design
    //! @param ATOM_COMPARISON the atom comaprison type to use for subgraph isomorphism
    //! @param BOND_COMPARISON the bond comaprison type to use for subgraph isomorphism
    //! @param COMPLEMENT if true, the mapped atom indices returned are the subgraph complement,
    //! if false then the returned indices are of the common subgraph
    //! @return the NEW_MOL indices mapped to STARTING_MOL
    storage::Set< size_t> FragmentMapConformer::MapSubgraphRingAtoms
    (
      const graph::Subgraph< size_t, size_t> &SUBGRAPH,
      const storage::List< storage::Vector< size_t> > &RING_COMPONENTS
    ) const
    {
      // decide which rings matter, as defined by rings that contain perturbed atoms
      storage::Vector< size_t> initial_indices( SUBGRAPH.GetVertexIndices().Begin(), SUBGRAPH.GetVertexIndices().End());
      storage::Set< size_t> rings_that_matter;
      size_t ring_index( 0);
      for
      (
          auto ring_comp_itr( RING_COMPONENTS.Begin()), ring_comp_itr_end( RING_COMPONENTS.End());
          ring_comp_itr != ring_comp_itr_end;
          ++ring_comp_itr, ++ring_index
      )
      {
        for
        (
            auto atom_itr( ring_comp_itr->Begin()), atom_itr_end( ring_comp_itr->End());
            atom_itr != atom_itr_end;
            ++atom_itr
        )
        {
          // check against all of my perturbed atoms
          for( size_t i( 0); i < initial_indices.GetSize(); ++i)
          {
            if( *atom_itr == initial_indices( i))
            {
              rings_that_matter.InsertElement( ring_index);
              break;
            }
          }
        }
      }

      // get atoms of rings that matter
      storage::Set< size_t> important_ring_atoms;
      storage::Vector< storage::Vector< size_t> > ring_components_vec
      (
        RING_COMPONENTS.Begin(),
        RING_COMPONENTS.End()
      );
      for
      (
        auto important_rings_itr( rings_that_matter.Begin()),
          import_rings_itr_end( rings_that_matter.End());
        important_rings_itr != import_rings_itr_end;
        ++important_rings_itr
      )
      {
        for( size_t atom( 0); atom < ring_components_vec( *important_rings_itr).GetSize(); ++atom)
        {
          important_ring_atoms.InsertElement( ring_components_vec( *important_rings_itr)( atom));
        }
      }

      return important_ring_atoms;
    }

    //! @brief preserve conformational information from starting molecule in new molecule
    //! @param SUBGRAPH the subgraph containing rings
    //! @param N_ADJACENT_NBRS the number of adjacent neighbors to include
    //! @return the indices of adjacent atoms in the subgraph
    storage::Set< size_t> FragmentMapConformer::MapSubgraphAdjacentAtoms
    (
      const graph::Subgraph< size_t, size_t> &SUBGRAPH,
      const size_t N_ADJACENT_NBRS
    ) const
    {
      if( !N_ADJACENT_NBRS)
      {
        return storage::Set< size_t>();
      }

      graph::Subgraph< size_t, size_t> subgraph( SUBGRAPH);
      storage::Set< size_t> all_adjacent_indices;
      for( size_t i( 0); i < N_ADJACENT_NBRS; ++i)
      {
        storage::List< storage::Pair< size_t, size_t> > b_adjacent_indices( subgraph.GetAdjacentEdgeIndices());
        for
        (
            auto itr( b_adjacent_indices.Begin()), itr_end( b_adjacent_indices.End());
            itr != itr_end;
            ++itr
        )
        {
          all_adjacent_indices.Insert( itr->First());
          all_adjacent_indices.Insert( itr->Second());
        }
        subgraph.SetVertexIndices( storage::Vector< size_t>( all_adjacent_indices.Begin(), all_adjacent_indices.End()));
      }
      return all_adjacent_indices;
    }

    //! @brief clean and generate a 3D structure of a molecule without perturbing atom indices
    //! @param MOL the molecule of interest
    //! @param ATOM_COMPARISON the atom comaprison type to use for subgraph isomorphism
    //! @param BOND_COMPARISON the bond comaprison type to use for subgraph isomorphism
    //! @return a 3D conformer of the original molecule with preserved atom indices
    FragmentComplete FragmentMapConformer::Clean3DCoords
    (
      const FragmentComplete &MOL,
      const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON,
      const ConfigurationalBondTypeData::Data &BOND_COMPARISON
    ) const
    {
      // generate 3D conformer
      FragmentComplete mol( MOL);
      util::ShPtr< FragmentComplete> clean_mol( Clean( MOL.GetAtomVector(), MOL, "None", true));
      if( !util::IsDefined( clean_mol))
      {
        return FragmentComplete();
      }

      // prepare dehydrogenated copy
      ConformationGraphConverter graph_maker
      (
        ATOM_COMPARISON,
        BOND_COMPARISON,
        false
      );

      // set up the isomorphisms with pointers to the graphs
      graph::CommonSubgraphIsomorphism< size_t, size_t> isomorphism( graph::CommonSubgraphIsomorphismBase::SolutionType::e_GreedyUnconnected);
      isomorphism.SetGraphA( graph_maker( mol));
      isomorphism.SetGraphB( graph_maker( *clean_mol));

      // get the isomorphism - allocate array based on maximum possible atom matches without considering connectivity
      isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds());
      const storage::Map< size_t, size_t> &iso_map( isomorphism.GetIsomorphism());

      // set the positions of the atoms in target molecule
      size_t atom_index( 0);
      AtomVector< AtomComplete> atoms( mol.GetAtomVector());
      for
      (
          auto atom_itr( atoms.Begin()), atom_itr_end( atoms.End());
          atom_itr != atom_itr_end;
          ++atom_itr, ++atom_index
      )
      {
        const linal::Vector3D position( clean_mol->GetAtomVector()( iso_map.Find( atom_index)->second).GetPosition());
        atom_itr->SetPosition( position);
      }
      return FragmentComplete( atoms, mol.GetName());
    }

    //! @brief reset the moveable indices member data to an empty vector
    void FragmentMapConformer::ResetMoveableIndices() const
    {
      m_MoveableIndices = storage::Vector< size_t>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentMapConformer::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FragmentMapConformer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
