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
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_interface_retrieve_from_file.h"
#include "model/bcl_model_kappa_nearest_neighbor.h"
#include "model/bcl_model_retrieve_interface.h"
#include "util/bcl_util_assert.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonPsiFlexField::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonPsiFlexField())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    ConformationComparisonPsiFlexField::ConformationComparisonPsiFlexField() :
      ConformationComparisonPsiField() ,
      m_PairNumber( 100),
      m_NumberBestPairsFromEachTraj( 3),
      m_FilterIterations( 600),
      m_FilterLimit( 240),
      m_RefinementIterations( 200),
      m_RefinementLimit( 80),
      m_NumberFlexibleTrajectories( 5),
      m_FractionFilterInitial( 0.25),
      m_FractionFilterIterations( 0.50),
      m_RigidMoleculeB( true),
      m_SampleConformations()
    {
    }

    //! @brief full constructor
    ConformationComparisonPsiFlexField::ConformationComparisonPsiFlexField
    (
      const ConformationComparisonPsiField &PSIFIELD,
      const size_t &N_PAIRS,
      const size_t &N_BEST_PAIRS_EACH_TRAJ,
      const size_t &FILTER_ITERATIONS,
      const size_t &FILTER_LIMIT,
      const size_t &REFINEMENT_ITERATIONS,
      const size_t &REFINMENT_LIMIT,
      const size_t &N_FLEXIBLE_TRAJ,
      const float &FRACTION_FILTER_INITIAL,
      const float &FRACTION_FILTER_ITERATIONS,
      const bool &RIGID_MOL_B,
      const SampleConformations &SAMPLE_CONFS
    ) :
      ConformationComparisonPsiField( PSIFIELD),
      m_PairNumber( N_PAIRS),
      m_NumberBestPairsFromEachTraj( N_BEST_PAIRS_EACH_TRAJ),
      m_FilterIterations( FILTER_ITERATIONS),
      m_FilterLimit( FILTER_LIMIT),
      m_RefinementIterations( REFINEMENT_ITERATIONS),
      m_RefinementLimit( REFINMENT_LIMIT),
      m_NumberFlexibleTrajectories( N_FLEXIBLE_TRAJ),
      m_FractionFilterInitial( FRACTION_FILTER_INITIAL),
      m_FractionFilterIterations( FRACTION_FILTER_ITERATIONS),
      m_RigidMoleculeB( RIGID_MOL_B),
      m_SampleConformations( SAMPLE_CONFS)
    {
    }

    //! virtual copy constructor
    ConformationComparisonPsiFlexField *ConformationComparisonPsiFlexField::Clone() const
    {
      return new ConformationComparisonPsiFlexField( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConformationComparisonPsiFlexField::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ConformationComparisonPsiFlexField::GetAlias() const
    {
      static std::string s_name( "PsiFlexField");
      return s_name;
    }

    //! @brief return the sample conformations object used for flexible alignment
    //! @return the sample conformers object
    const SampleConformations ConformationComparisonPsiFlexField::GetSampleConfs() const
    {
      return m_SampleConformations;
    }

    //! @brief return the number of conformer pairs
    //! @return the number of conformer pairs
    const size_t ConformationComparisonPsiFlexField::GetNumberConformerPairs() const
    {
      return m_PairNumber;
    }

    //! @brief return the number of iterations during iterative filter steps
    //! @return the number of filter iterations
    const size_t ConformationComparisonPsiFlexField::GetNumberFilterIterations() const
    {
      return m_FilterIterations;
    }

    //! @brief return the number of max unimproved iterations allowed in the MCM during filter steps
    //! @return the number of max unimproved filter iterations
    const size_t ConformationComparisonPsiFlexField::GetNumberFilterMaxUnimproved() const
    {
      return m_FilterLimit;
    }

    //! @brief return the number of refinement iterations post-filtering
    //! @return the number of refinement iterations
    const size_t ConformationComparisonPsiFlexField::GetNumberRefinementIterations() const
    {
      return m_RefinementIterations;
    }

    //! @brief return the number of max unimproved iterations allowed during the refinement MCM
    //! @return the number of max unimproved refinement iterations
    const size_t ConformationComparisonPsiFlexField::GetNumberRefinementMaxUnimproved() const
    {
      return m_RefinementLimit;
    }

    //! @brief return the fraction of conformer pairs filtered out after first alignment round
    //! @return the fraction filtered initially
    const float ConformationComparisonPsiFlexField::GetFractionFilterInitial() const
    {
      return m_FractionFilterInitial;
    }

    //! @brief return the fraction of conformer pairs filtered out at each successive alignment iteration
    //! @return the fraction filtered iteratively
    const float ConformationComparisonPsiFlexField::GetFractionFilterIterative() const
    {
      return m_FractionFilterIterations;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief make conformers of input molecules, randomly pair conformers
    //! @brief align and find best property rmsd between input molecules
    //! @param MOLECULE_A - first molecule from which conformational library will be generated
    //! @param MOLECULE_B - second molecule from which conformational library will be generated

    double ConformationComparisonPsiFlexField::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      storage::Vector< size_t> keep_indices_a, keep_indices_b;
      GetNonMaskedAtoms( MOLECULE_A, MOLECULE_B, m_ExclusionIndicesA, m_ExclusionIndicesB, keep_indices_a, keep_indices_b);
      m_KeepIndicesA = keep_indices_a;
      m_KeepIndicesB = keep_indices_b;

      // make conformers of input molecules
      FragmentEnsemble mol_ensA(MakeConformers(MOLECULE_A));
      FragmentEnsemble mol_ensB
      (
        !m_RigidMoleculeB
        ? MakeConformers( MOLECULE_B)
        : FragmentEnsemble( storage::List< FragmentComplete>( size_t( 1), MOLECULE_B))
      );

      // perform alignment
      return FieldOptimizeOrientationFlex
          (
            mol_ensA,
            mol_ensB,
            m_PairNumber,
            m_IterationNumber,
            m_LimitNumber,
            m_FilterIterations,
            m_FilterLimit,
            m_RefinementIterations,
            m_RefinementLimit,
            m_FractionFilterInitial,
            m_FractionFilterIterations,
            true
          )( 0).Third();
    }

    storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > ConformationComparisonPsiFlexField::FieldOptimizeOrientationFlex
    (
      const FragmentEnsemble &CONFORMERS_MOL_A,
      const FragmentEnsemble &CONFORMERS_MOL_B,
      const size_t &PAIRS,
      const size_t &ITERATIONS,
      const size_t &MAX_UNIMPROVED,
      const size_t &FILTER_ITERATIONS,
      const size_t &FILTER_LIMIT,
      const size_t &REFINEMENT_ITERATIONS,
      const size_t &REFINEMENT_LIMIT,
      const float &FRACTION_FILTER_INITIAL,
      const float &FRACTION_FILTER_ITERATIONS,
      const bool &CENTER_ACONFS_ON_BCONFS,
      const FragmentEnsemble &POCKETS
    ) const

    {
      FragmentEnsemble mol_ensA( CONFORMERS_MOL_A);
      size_t mol_ensA_size( mol_ensA.GetSize());
      size_t mol_ensB_size( CONFORMERS_MOL_B.GetSize());
      // total possible conformer pairings
      size_t possible_pairs = mol_ensA_size * mol_ensB_size;

      // use a KNN to identify property-matched conformer pairs between molecules A and B
//      BCL_Debug( m_DatasetBuilder.GetFeatureCode());
      descriptor::Dataset dataset_a( m_DatasetBuilder( iterate::Generic< const descriptor::SequenceInterface< AtomConformationalInterface> >( mol_ensA.Begin(), mol_ensA.End())));
      descriptor::Dataset dataset_b( m_DatasetBuilder( iterate::Generic< const descriptor::SequenceInterface< AtomConformationalInterface> >( CONFORMERS_MOL_B.Begin(), CONFORMERS_MOL_B.End())));
      auto closest_by_properties( model::KappaNearestNeighbor::FindWithoutRescaling( dataset_a.GetFeatures(), PAIRS, dataset_b.GetFeatures()));

      // select conformer pairings
      storage::Vector< storage::Pair< size_t, size_t> > picked_pairs;
      const size_t find_pairs( std::min( m_PairNumber, possible_pairs));
      for( auto itr( closest_by_properties.Begin()), itr_end( closest_by_properties.End()); itr != itr_end; ++itr)
      {
        picked_pairs.PushBack( storage::Pair< size_t, size_t>( itr->Second(), itr->Third()));
      }

      // storage for conformer pairings
      storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > cumulative_pair_scores;

      // perform m_NumberFlexibleTrajectories independent trajectories
      // These independent trajectories all go through 3 tiers of MCM optimization
      // At the end, the best m_NumberBestPairsFromEachTraj pairs from each independent trajectory are compared by
      // re-scoring using common max atom distances.
      storage::Vector< double> distances;
      float dmax_increment( ( m_AtomDistanceUpperLimit - m_AtomDistanceLowerLimit) / ( m_NumberFlexibleTrajectories - 1));
      for( size_t traj( 0); traj < m_NumberFlexibleTrajectories; ++traj)
      {
        // MC-metropolis optimization tier 1 - all pairs
        storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > pair_scores;
        float max_atom_distance( m_AtomDistanceLowerLimit + ( dmax_increment * traj));
        distances.PushBack( max_atom_distance);
        for( size_t pair_index = 0; pair_index < picked_pairs.GetSize(); ++pair_index)
        {
          // center each target molecule (molecule A) conformer on the scaffold molecule (molecule B) conformer
          size_t posA( picked_pairs( pair_index).First());
          size_t posB( picked_pairs( pair_index).Second());

          FragmentEnsemble::const_iterator mol_b_itr( CONFORMERS_MOL_B.Begin());
          std::advance( mol_b_itr, posB);

          FragmentEnsemble::iterator mol_a_itr( mol_ensA.Begin());
          std::advance( mol_a_itr, posA);

          // Should always be true for general flexible alignment, but is specified false in MultiAlign
          if( CENTER_ACONFS_ON_BCONFS)
          {
            mol_a_itr->Translate( mol_b_itr->GetCenter() - mol_a_itr->GetCenter());
          }

          // MCM optimization #1
          auto opti_pair
          (
            this->FieldOptimizeOrientation
            (
              *mol_a_itr,
              *mol_b_itr,
              ITERATIONS,
              MAX_UNIMPROVED,
              false,
              max_atom_distance,
              mol_ensA//,
//              POCKETS
            )
          );
          pair_scores.PushBack( storage::Triplet< FragmentComplete, FragmentComplete, double>( opti_pair.First(), *mol_b_itr, opti_pair.Second()));
        } // end MC-metropolis optimization tier 1

        // sort pair_scores to put best scoring pairs first
        TripletComparison triplet_compare;
        pair_scores.Sort< TripletComparison>( triplet_compare);

        // keep only best m_FractionFilterInitial x 100 % of conformer pairs; delete the rest
        pair_scores.Resize( std::min( pair_scores.GetSize(), std::max< size_t>( size_t( m_NumberBestPairsFromEachTraj), size_t( find_pairs * FRACTION_FILTER_INITIAL))));

        // MCM optimization tier 2
        // repeat optimization saving best m_FractionFilterIterations x 100 % conformer pairs each time until
        // only the best conformer pairing(s) for each pair of molecules is left
        while( pair_scores.GetSize() > m_NumberBestPairsFromEachTraj)
        {
          size_t itr_index( 0);
          for
          (
              storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> >::iterator
              itr = pair_scores.Begin(), itr_end = pair_scores.End();
              itr != itr_end;
              ++itr, ++itr_index
          )
          {
            auto opti_pair( this->FieldOptimizeOrientation
              (
                itr->First(),
                itr->Second(),
                FILTER_ITERATIONS,
                FILTER_LIMIT,
                false,
                max_atom_distance,
                mol_ensA//,
//                POCKETS
              ));
            itr->First() = opti_pair.First();
            itr->Third() = opti_pair.Second();
          }
          pair_scores.Sort< TripletComparison>( triplet_compare);
          pair_scores.Resize( std::min( pair_scores.GetSize(), std::max( size_t( pair_scores.GetSize() * FRACTION_FILTER_ITERATIONS), m_NumberBestPairsFromEachTraj)));
        } // end MC-metropolis optimization tier 2

        // MCM Optimization tier 3 - final optimization
        for( size_t best_pairs_index( 0), mx( std::min( m_NumberBestPairsFromEachTraj, pair_scores.GetSize())); best_pairs_index < mx; ++best_pairs_index)
        {
          auto opti_pair
          (
            this->FieldOptimizeOrientation
            (
              pair_scores( best_pairs_index).First(),
              pair_scores( best_pairs_index).Second(),
              REFINEMENT_ITERATIONS,
              REFINEMENT_LIMIT,
              false,
              max_atom_distance,
              mol_ensA//,
//              POCKETS
            )
          );
          pair_scores( best_pairs_index).First() = opti_pair.First();
          pair_scores( best_pairs_index).Third() = opti_pair.Second();
        } // end MC-metropolis optimization tier 3

        // Performs substructure-based alignment with best pair
        if( m_AlignToScaffold)
        {
          // cache properties
          FragmentComplete mol_a( pair_scores( 0).First());
          mol_a.ShareCache( pair_scores( 0).First());
          FragmentComplete mol_b( pair_scores( 0).Second());
          mol_b.ShareCache( pair_scores( 0).Second());

          storage::Vector< storage::Pair< ConformationGraphConverter::AtomComparisonTypeEnum, ConfigurationalBondTypeData::DataEnum>> ats_options;
          for( size_t e_atom( 0); e_atom < ConformationGraphConverter::s_NumberAtomComparisonTypes; ++e_atom)
          {
            for( size_t e_bond( 1); e_bond < ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness; ++e_bond)
            {
              ArbitraryScaffoldAlignment
              (
                mol_a,
                mol_b,
                static_cast< ConformationGraphConverter::AtomComparisonType>( e_atom),
                static_cast< ConfigurationalBondTypeData::Data>( e_bond)
              );

              pair_scores.PushBack();
              pair_scores.LastElement().First() = mol_a;
              pair_scores.LastElement().Second() = mol_b;
              pair_scores.LastElement().Third() = this->PropertyDistanceScorer( mol_a, mol_b, max_atom_distance).First();
            }
          }
        } // end optional substructure-based alignment

        // at the end of each independent trajectory save the best pair(s)
        cumulative_pair_scores.Append( pair_scores);
      }

      // Average scores over all max atom distances; naive but worked originally
      size_t tracker_index( 0);
      storage::Vector< storage::Pair< double, size_t> > sorter( cumulative_pair_scores.GetSize());
      for( auto itr( cumulative_pair_scores.Begin()), itr_end( cumulative_pair_scores.End()); itr < itr_end; ++itr, ++tracker_index)
      {
        math::RunningAverage< double> sum( 0.0);
        for( size_t i( 0); i < distances.GetSize(); ++i)
        {
          sum += PropertyDistanceScorer( itr->First(), itr->Second(), distances( i)).First();
        }
        sorter( tracker_index).First() = sum.GetAverage();
        itr->First().StoreProperty( "PropertyFieldDistance", util::Format()( sum.GetAverage()));
        itr->Second().StoreProperty( "PropertyFieldDistance", util::Format()( sum.GetAverage()));
        sorter( tracker_index).Second() = tracker_index;
      }
      // Sort molecules by best
      sorter.Sort( std::less< storage::Pair< double, size_t>>());
      storage::Vector< size_t> best_order( sorter.GetSize());
      for( size_t i( 0); i < sorter.GetSize(); ++i)
      {
        best_order( i) = sorter( i).Second();
      }
      cumulative_pair_scores.Reorder( best_order);

      // prune similar molecules
      ConformationComparisonBySymmetryRmsd compare_symrmsd( false, 100); // consider real space difference
      storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > pruned_pair_scores;

      // add first molecule to new molecules collection
      pruned_pair_scores.PushBack( cumulative_pair_scores( 0));

      // check other alignment poses; if they deviate by > 0.50 angstroms then they are different and save
      double tolerance( 1.0);

      // loop over all poses
      for( size_t i( 0); i < cumulative_pair_scores.GetSize(); ++i)
      {
        // loop over the saved poses
        bool same( false);
        for( size_t j( 0); j < pruned_pair_scores.GetSize(); ++j)
        {
          if( compare_symrmsd( cumulative_pair_scores( i).First(), pruned_pair_scores( j).First()) <= tolerance)
          {
            // outer loop contains something same as already in new collection
            same = true;
            break;
          }
        }
        if( !same)
        {
          pruned_pair_scores.PushBack( cumulative_pair_scores( i));
        }
      }

      // output the best alignment pairs and their rmsdx values
      storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > final_pair_scores;

      // this is fundamentally flawed because it does not inherently preserve order if multiple alignments are being run simultaneously
      // this is just a hack that allows molecules to be output via this function directly so that we can get molecules via molecule:Compare
      // molecule:Compare is really only supposed to return the double from operator()
      // TODO: deprecate this with cheminfo:MoleculeFit when possible
      static sched::Mutex s_mutex;
      s_mutex.Lock();
      for( size_t best_pairs_index( 0); best_pairs_index < m_NumberOutputSDF; ++best_pairs_index)
      {
        if( best_pairs_index < pruned_pair_scores.GetSize())
        {
          // save best for return
          final_pair_scores.PushBack( pruned_pair_scores( best_pairs_index));

          // WriteMDL
          if( !m_OutputA.empty())
          {
            io::OFStream out;
            io::File::MustOpenOFStream( out, m_OutputA + "_" + std::to_string( best_pairs_index) + ".sdf", std::ios::app);
            pruned_pair_scores( best_pairs_index).First().StoreProperty( "PropertyFieldDistance", util::Format()( sorter( best_pairs_index).First()));
            pruned_pair_scores( best_pairs_index).First().WriteMDL( out);
            io::File::CloseClearFStream( out);
          }
          if( !m_OutputB.empty())
          {
            io::OFStream out;
            io::File::MustOpenOFStream( out, m_OutputB + "_" + std::to_string( best_pairs_index) + ".sdf", std::ios::app);
            pruned_pair_scores( best_pairs_index).Second().StoreProperty( "PropertyFieldDistance", util::Format()( sorter( best_pairs_index).First()));
            pruned_pair_scores( best_pairs_index).Second().WriteMDL( out);
            io::File::CloseClearFStream( out);
          }
        }
        else
        {
          break;
        }
      }
      s_mutex.Unlock();

      // return final structures/score
      return final_pair_scores;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    // Create conformers from input molecules
    FragmentEnsemble ConformationComparisonPsiFlexField::MakeConformers ( const FragmentComplete &MOL) const
    {
      return m_SampleConformations( MOL).First();
    }

    namespace
    {
      util::Implementation< descriptor::Base< AtomConformationalInterface, float> > &GetDescriptorsPlease()
      {
        static util::Implementation< descriptor::Base< AtomConformationalInterface, float> > s_implementation
        (
          "Combine(                                                                                                                         "
          "  Define(BondGirth=DescriptorSum(2DAMax(steps=96,property=Atom_Identity,substitution_value=nan))),                               "
          "  Define(IsHTernary=Add(Constant(-1),Multiply(IsH,Constant(2)))),                                                                "
          "  Define(Atom_IsInAromaticRing=GreaterEqual(lhs=BondTypeCount(property=IsAromatic,value=1),rhs=2)),                              "
          "  Define(Atom_IsInAromaticRingTernary=Add(Constant(-1),Multiply(Atom_IsInAromaticRing,Constant(2)))),                            "
          "  Define(Atom_InAromaticRingIntersection=GreaterEqual(lhs=BondTypeCount(property=IsAromatic,value=1),rhs=3)),                    "
          "  Define(Atom_InRingIntersection=GreaterEqual(lhs=BondTypeCount(property=IsInRing,value=1),rhs=3)),                              "
          "  Define(ChargeNegative=Multiply(Less(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),                                  "
          "  Define(ChargePositive=Multiply(Greater(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),                               "
          "  Define(HbondAcceptorsStrict=Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))), "
          "  Define(ENegOffset=Subtract(lhs=Atom_ElectroNegativity,rhs=Constant(2.5))),                                                     "
          "  Define(IsENeg=Greater(lhs=ENegOffset,rhs=Constant(0))),                                                                        "
          "  Define(PolarTernary=Subtract(lhs=Multiply(Add(HbondAcceptorsStrict, Atom_HbondDonors), Constant(2)), rhs=Constant(1))),        "
          "  Define(Atom_HBondDonorTernary=Subtract(lhs=Multiply(Atom_HbondDonors, Constant(2)), rhs=HbondAcceptorsStrict)),                "
          "  Define(Atom_Hydrophobic=Not(DescriptorSum(Abs(Partial(indices(1,2,5,8),2DASign(property=PolarTernary, steps=3)))))),           "
          "  Define(Atom_HydrophobicTernary=Subtract(lhs=Multiply(Atom_Hydrophobic, Constant(2)), rhs=Constant(1))),                        "
          "  Weight,                                                                                                                        "
          "  HbondDonor,                                                                                                                    "
          "  HbondAcceptor,                                                                                                                 "
          "  LogP,                                                                                                                          "
          "  TotalCharge,                                                                                                                   "
          "  NRotBond,                                                                                                                      "
          "  NAromaticRings,                                                                                                                "
          "  NRings,                                                                                                                        "
          "  TopologicalPolarSurfaceArea,                                                                                                   "
          "  Girth,                                                                                                                         "
          "  BondGirth,                                                                                                                     "
          "  MaxRingSize,                                                                                                                   "
          "  Limit(MinRingSize,max=8,min=0),                                                                                                "
          "  MoleculeSum(Atom_InAromaticRingIntersection),                                                                                  "
          "  MoleculeSum(Atom_InRingIntersection),                                                                                          "
          "  MoleculeStandardDeviation(Atom_Vcharge),                                                                                       "
          "  MoleculeStandardDeviation(Atom_SigmaCharge),                                                                                   "
          "  MoleculeMax(Atom_Vcharge),                                                                                                     "
          "  MoleculeMax(Atom_SigmaCharge),                                                                                                 "
          "  MoleculeMin(Atom_Vcharge),                                                                                                     "
          "  MoleculeMin(Atom_SigmaCharge),                                                                                                 "
          "  MoleculeSum(Abs(Atom_Vcharge)),                                                                                                "
          "  MoleculeSum(Abs(Atom_SigmaCharge)),                                                                                            "
          "  ConfScore,                                                                                                                     "
          "  ConfScorePBest,                                                                                                                "
          "  Template(                                                                                                                      "
          "    signature=2DASign11(X),                                                                                                      "
          "    Partial(                                                                                                                     "
          "      2DASign(property=X,steps=11),                                                                                              "
          "      indices(0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)                            "
          "    )                                                                                                                            "
          "  ),                                                                                                                             "
          "  Template(                                                                                                                      "
          "    signature=3DASign24(X),                                                                                                      "
          "    Partial(                                                                                                                     "
          "      3daSmoothSign(property=X,step size=0.25,temperature=100,steps=24,gaussian=False,interpolate=True),                         "
          "      indices(                                                                                                                   "
          "        12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,                                              "
          "        37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,                                              "
          "        62,63,64,65,66,67,68,69,70,71                                                                                            "
          "      )                                                                                                                          "
          "    )                                                                                                                            "
          "  ),                                                                                                                             "
          "  2DASign11(Atom_SigmaCharge),                                                                                                   "
          "  2DASign11(Atom_Vcharge),                                                                                                       "
          "  2DASign11(IsHTernary),                                                                                                         "
          "  2DASign11(Atom_IsInAromaticRingTernary),                                                                                       "
          "  2DASign11(Atom_HydrophobicTernary),                                                                                            "
          "  2DASign11(Atom_HBondDonorTernary),                                                                                             "
          "  3DASign24(Atom_SigmaCharge),                                                                                                   "
          "  3DASign24(Atom_Vcharge),                                                                                                       "
          "  3DASign24(IsHTernary),                                                                                                         "
          "  3DASign24(Atom_IsInAromaticRingTernary),                                                                                       "
          "  3DASign24(Atom_HydrophobicTernary),                                                                                            "
          "  3DASign24(Atom_HBondDonorTernary)                                                                                              "
          ")                                                                                                                                "
        );
        return s_implementation;
      }

      const model::RetrieveInterface::t_Container &GetModelsPlease()
      {
        static const model::InterfaceRetrieveFromFile storage( model::Model::AddModelPath( "conf_pair_prediction"), "model");
        static const model::RetrieveInterface::t_Container anns( storage.RetrieveEnsemble());
        return anns;
      }
    }

    // TODO: use an ANN to pick the best conformer pairs for a given alignment
    // Pick best conformer pair
    linal::Vector< double> ConformationComparisonPsiFlexField::PickConformerPairs( const FragmentComplete &MOL, const FragmentEnsemble &ENS)
    {
      return linal::Vector< double>();
    }

    //! @brief set the sample conformations object
    void ConformationComparisonPsiFlexField::SetSampleConformations( const SampleConformations &SAMPLE_CONFS)
    {
      m_SampleConformations = SAMPLE_CONFS;
    }

    //! @brief set the number of conformer pairs
    void ConformationComparisonPsiFlexField::SetNumberConformerPairs( const size_t N_CONF_PAIRS)
    {
      m_PairNumber = N_CONF_PAIRS;
    }

    //! @brief set the number of iterations during iterative filter steps
    void ConformationComparisonPsiFlexField::SetNumberFilterIterations( const size_t N_FILTER_ITERATIONS)
    {
      m_FilterIterations = N_FILTER_ITERATIONS;
    }

    //! @brief set the number of max unimproved iterations allowed in the MCM during filter steps
    void ConformationComparisonPsiFlexField::SetNumberFilterMaxUnimproved( const size_t N_FILTER_MAX_UNIMPROVED)
    {
      m_FilterLimit = N_FILTER_MAX_UNIMPROVED;
    }

    //! @brief set the number of refinement iterations post-filtering
    void ConformationComparisonPsiFlexField::SetNumberRefinementIterations( const size_t N_REFINE_ITERATIONS)
    {
      m_RefinementIterations = N_REFINE_ITERATIONS;
    }

    //! @brief set the number of max unimproved iterations allowed during the refinement MCM
    void ConformationComparisonPsiFlexField::SetNumberRefinementMaxUnimproved( const size_t N_REFINE_MAX_UNIMPROVED)
    {
      m_RefinementLimit = N_REFINE_MAX_UNIMPROVED;
    }

    //! @brief set the fraction of conformer pairs filtered out after first alignment round
    void ConformationComparisonPsiFlexField::SetFractionFilterInitial( const float FILTER_FRACTION_INITIAL)
    {
      m_FractionFilterInitial = FILTER_FRACTION_INITIAL;
    }

    //! @brief set the fraction of conformer pairs filtered out at each successive alignment iteration
    void ConformationComparisonPsiFlexField::SetFractionFilterIterative( const float FILTER_FRACTION_ITERATIVE)
    {
      m_FractionFilterIterations = FILTER_FRACTION_ITERATIVE;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonPsiFlexField::GetSerializer() const
    {
      io::Serializer parameters( ConformationComparisonPsiField::GetSerializer());
      parameters.SetClassDescription( "Aligns and computes property RMSD for molecule conformer libraries");

      parameters.AddInitializer
      (
        "rigid_mol_b",
        "do not perform conformational sampling for molecule b",
        io::Serialization::GetAgent( &m_RigidMoleculeB),
        "true"
      );

      parameters.AddInitializer
      (
        "conformer_pairs",
        "number of pairs of conformers sampled",
        io::Serialization::GetAgent( &m_PairNumber),
        "100"
      );

      parameters.AddInitializer
      (
        "top_pairs",
        "number of best pairs within each trajectory included during inter-trajectory scoring",
        io::Serialization::GetAgent( &m_NumberBestPairsFromEachTraj),
        "5"
      );

      parameters.AddInitializer
      (
        "filter_iterations",
        "number of mc iterations after initial mc",
        io::Serialization::GetAgent( &m_FilterIterations),
        "600"
      );

      parameters.AddInitializer
      (
        "filter_limit",
        "max number of unimproved mc iterations after initial mc",
        io::Serialization::GetAgent( &m_FilterLimit),
        "240"
      );

      parameters.AddInitializer
      (
        "refinement_iterations",
        "number of mc iterations on the final molecule(s)",
        io::Serialization::GetAgent( &m_RefinementIterations),
        "200"
      );

      parameters.AddInitializer
      (
        "refinement_limit",
        "max number of unimproved mc iterations during final refinement",
        io::Serialization::GetAgent( &m_RefinementLimit),
        "80"
      );

      parameters.AddInitializer
      (
        "number_flexible_trajectories",
        "number of independent flexible trajectories to run",
        io::Serialization::GetAgent( &m_NumberFlexibleTrajectories),
        "5"
      );

      parameters.AddInitializer
      (
        "fraction_filtered_initially",
        "set the fraction of best conformer pairs to be saved after first pass",
        io::Serialization::GetAgent( &m_FractionFilterInitial),
        "0.25"
      );

      parameters.AddInitializer
      (
        "fraction_filtered_iteratively",
        "fraction of best pairs saved with each subsequent iteration until best pair",
        io::Serialization::GetAgent( &m_FractionFilterIterations),
        "0.50"
      );

      parameters.AddInitializer
      (
        "sample_conformers",
        "makes conformers of molecules",
        io::Serialization::GetAgent( &m_SampleConformations),
        util::ObjectDataLabel
        (
          "(conformation_comparer=SymmetryRMSD,tolerance=0.25,"
          "generate_3D=0, cluster=true,"
          "max_iterations=2000, max_conformations=200,"
          "change_chirality=0)"
        )
      );

      parameters.AddInitializer
      (
        "conf_properties",
        "atom properties to consider, use multiply(Constant(X),property y) for weighting",
        io::Serialization::GetAgent( &m_ConfPairProperties),
        "(1)"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ConformationComparisonPsiFlexField::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      for( size_t best_pairs_index( 0); 1; ++best_pairs_index)
      {
        const std::string filename_A( m_OutputA + "_" + std::to_string( best_pairs_index) + ".sdf");
        io::DirectoryEntry entry_A( filename_A);
        if( entry_A.DoesExist())
        {
          entry_A.Remove();
        }
        else
        {
          break;
        }
        const std::string filename_B( m_OutputB + "_" + std::to_string( best_pairs_index) + ".sdf");
        io::DirectoryEntry entry_B( filename_B);
        if( entry_B.DoesExist())
        {
          entry_B.Remove();
        }
        else
        {
          break;
        }
      }

//      BCL_Debug(m_ConfPairProperties.GetCacheLabel());
      m_DatasetBuilder.SetFeatureCode( m_ConfPairProperties.GetCacheLabel());
      return ConformationComparisonPsiField::ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }
  } // namespace chemistry
} // namespace bcl
