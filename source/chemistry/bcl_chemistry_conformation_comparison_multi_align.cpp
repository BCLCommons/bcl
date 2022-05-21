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
#include "chemistry/bcl_chemistry_conformation_comparison_multi_align.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_volume.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "math/bcl_math_running_min_max.h"
#include "model/bcl_model_data_set_reduced_to_cluster_centers.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonMultiAlign::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonMultiAlign())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    ConformationComparisonMultiAlign::ConformationComparisonMultiAlign() :
      ConformationComparisonPsiFlexField(),
      m_MovieLigFitFilename( ""),
      m_StartPos( util::GetUndefinedDouble())
    {
    }

    //! virtual copy constructor
    ConformationComparisonMultiAlign *ConformationComparisonMultiAlign::Clone() const
    {
      return new ConformationComparisonMultiAlign( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConformationComparisonMultiAlign::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ConformationComparisonMultiAlign::GetAlias() const
    {
      static std::string s_name( "MultiAlign");
      return s_name;
    }

  /////////////////
  //  operations //
  /////////////////

    double ConformationComparisonMultiAlign::operator ()
        (
          const FragmentEnsemble &MOLECULE_LIST,
          const FragmentEnsemble &POCKET_LIST
        ) const
    {

      //Check if list has only one molecule or is empty.
      if( MOLECULE_LIST.GetSize() <= 1)
      {
        BCL_MessageStd("0 or 1 molecule in vector");
        BCL_MessageStd("It is difficult to align less than 2 molecules.");
//        return 0;
      }

      //Small molecule setup first

      //copy input molecules to modifiable variable
      storage::Vector< FragmentComplete> mol_vec( MOLECULE_LIST.Begin(), MOLECULE_LIST.End());

      //Make pointer to molecules to allow easy reordering without copying
      util::SiPtrVector< FragmentComplete> mol_vec_ptr( mol_vec.Begin(), mol_vec.End());
      size_t mol_vec_size( mol_vec_ptr.GetSize());

      //Generate conformers of each molecule
      storage::Vector< FragmentEnsemble> conformers( mol_vec_size);
      for( size_t index_f( 0); index_f < mol_vec_size; ++index_f)
      {
        conformers( index_f) = MakeConformers( *mol_vec_ptr( index_f));
      }

      //Make pointer to conformers to allow easy reordering without copying
      util::SiPtrVector< FragmentEnsemble> conformers_ptr( conformers.Begin(), conformers.End());

      //Calculate number of rotatable bonds and heavy atoms of each molecule
      size_t n_rotatable_bonds( 0); // n_atoms( 0);
      storage::Vector< storage::Pair< size_t, size_t>> mol_nrb_tracker( mol_vec_size);
      for( size_t mol_index( 0); mol_index < mol_vec_size; ++mol_index)
      {
        n_rotatable_bonds = descriptor::GetCheminfoProperties().calc_NRotBond->SumOverObject( *mol_vec_ptr( mol_index))( 0);
        mol_nrb_tracker( mol_index).First() = n_rotatable_bonds;
        mol_nrb_tracker( mol_index).Second() = mol_index;
      }

      //Reorder the molecules from least to most number of rotatable bonds
      mol_nrb_tracker.Sort( std::less< storage::Pair< size_t, size_t> >());

      storage::Vector< size_t> nrb_order( mol_vec_size);
      for( size_t i( 0); i < mol_vec_size; ++i)
      {
        nrb_order( i) = mol_nrb_tracker( i).Second();
      }
      mol_vec_ptr.Reorder( nrb_order);
      conformers_ptr.Reorder( nrb_order);

      return 0;

      //ANN to choose best conformer pair

      //Fully flexible alignment of conformers
      storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double>> best_pairwise_alignments( mol_vec_size - 1);
      size_t alignment_index( 0);
      for( size_t molecule_index( 1); molecule_index < mol_vec_size; ++molecule_index, ++alignment_index)
      {
        storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > sequential_alignment
        (
          FieldOptimizeOrientationFlex
          (
            *conformers_ptr(molecule_index),
            *conformers_ptr(0),
            m_PairNumber,
            m_IterationNumber,
            m_LimitNumber,
            m_FilterIterations,
            m_FilterLimit,
            m_RefinementIterations,
            m_RefinementLimit,
            m_FractionFilterInitial,
            m_FractionFilterIterations,
            true,
            POCKET_LIST
          )
        );
        best_pairwise_alignments( alignment_index) = sequential_alignment( 0);
      }

      //Get the best pose of each molecule
      FragmentEnsemble best_poses;
      best_poses.PushBack( best_pairwise_alignments( 0).Second());
      math::RunningAverage< double> multi_align_score( 0.0);
      for( size_t i( 0); i < best_pairwise_alignments.GetSize(); ++i)
      {
        best_poses.PushBack( best_pairwise_alignments( i).First());
        multi_align_score += best_pairwise_alignments( i).Third();
      }

      // output the best overall multiple alignment
      static sched::Mutex s_mutex;
      s_mutex.Lock();
      io::OFStream out;
      io::File::MustOpenOFStream( out, "MultiAlign.sdf");
      best_poses.WriteMDL( out);
      io::File::CloseClearFStream( out);
      s_mutex.Unlock();

      return multi_align_score.GetAverage();
    } // MultiAlign operator

  ////////////////////////////
  // Helper Classes         //
  ////////////////////////////

    class CompareTriplets
    {
    public:
      bool operator()( const storage::Triplet< size_t, size_t, float> &FIRST, const storage::Triplet< size_t, size_t, float> &SECOND) const
      {
        return FIRST.Third() > SECOND.Third();
      }
    };

  ////////////////////////////
  // Helper Functions       //
  ////////////////////////////

    // Create 3D bounding box
    // for ben: type<blah> foo::bar(const ArgumentType &INPUT) const
    // the const near the argument input says that I will not modify the argument object
    // the const at the end of the declaration means I promise not to modfify any member variables in the class bar belongs to (foo)
    linal::Vector< double> ConformationComparisonMultiAlign::GetAxisAlignedBoundingBox( FragmentComplete &TARGET) const
    {
      math::RunningMinMax< double> x, y, z;
      linal::Vector< double> bounding_box( 6);
      AtomVector< AtomComplete> atoms( TARGET.GetAtomVector());
      for( AtomVector< AtomComplete>::iterator itr( atoms.Begin()), itr_end( atoms.End()); itr != itr_end; ++itr)
      {
        x += itr->GetPosition().X();
        y += itr->GetPosition().Y();
        z += itr->GetPosition().Z();
      }

      // These points will define the edges of the bounding box
      bounding_box(0) = x.GetMin();
      bounding_box(1) = x.GetMax();
      bounding_box(2) = y.GetMin();
      bounding_box(3) = y.GetMax();
      bounding_box(4) = z.GetMin();
      bounding_box(5) = z.GetMax();
      return bounding_box;
    }

    FragmentComplete ConformationComparisonMultiAlign::BondAlignInfEns
    (
      FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const storage::Vector< storage::Triplet< size_t, size_t, float> > &INDICES
        ) const
        {
          storage::Vector< size_t> atom_indices_aligned_a, atom_indices_aligned_b;
      for( size_t index( 0); index < INDICES.GetSize(); ++index)
      {
        atom_indices_aligned_a.PushBack( INDICES( index).First());
        atom_indices_aligned_b.PushBack( INDICES( index).Second());
      }

      util::SiPtrVector< const linal::Vector3D> matched_atom_pos_a, matched_atom_pos_b;
      for( size_t aligned_pos( 0), n_aligned( atom_indices_aligned_a.GetSize()); aligned_pos < n_aligned; ++aligned_pos)
      {
        matched_atom_pos_a.PushBack( MOLECULE_A.GetAtomVector()( atom_indices_aligned_a( aligned_pos)).GetPosition());
        matched_atom_pos_b.PushBack( MOLECULE_B.GetAtomVector()( atom_indices_aligned_b( aligned_pos)).GetPosition());
      }

      quality::RMSD calculator;
      math::TransformationMatrix3D transformed_coords( calculator.CalculateSuperimposition( matched_atom_pos_a, matched_atom_pos_b));
      MOLECULE_A.Transform( transformed_coords);
      return MOLECULE_A;
    }

    storage::Pair< storage::Vector< storage::Triplet< size_t, size_t, float> >, storage::Pair< FragmentComplete, FragmentComplete> >
    ConformationComparisonMultiAlign::AlignMolecularEnsemble( const FragmentEnsemble &ENSEMBLE_A, const FragmentEnsemble &ENSEMBLE_B) const
    {
      // Make pseudo molecules from each ensemble
      storage::Vector< sdf::BondInfo> empty_bonds( 0);

      // molecule a
      AtomVector< AtomComplete> vector_a( ENSEMBLE_A.GetMolecules().FirstElement().GetAtomVector());
      for( FragmentEnsemble::const_iterator ens_itr_a( ENSEMBLE_A.Begin()), ens_itr_a_end( ENSEMBLE_A.End()); ens_itr_a != ens_itr_a_end; ++ens_itr_a)
      {
        if( ens_itr_a == ENSEMBLE_A.Begin())
        {
          continue;
        }
        vector_a.AddAtomsWithConnectivity( ens_itr_a->GetAtomVector(), empty_bonds);
      }
      FragmentComplete pseudo_mol_a( vector_a, "");

      // molecule b
      AtomVector< AtomComplete> vector_b( ENSEMBLE_B.GetMolecules().FirstElement().GetAtomVector());
      for( FragmentEnsemble::const_iterator ens_itr_b( ENSEMBLE_B.Begin()), ens_itr_b_end( ENSEMBLE_B.End()); ens_itr_b != ens_itr_b_end; ++ens_itr_b)
      {
        if( ens_itr_b == ENSEMBLE_B.Begin())
        {
          continue;
        }
        vector_b.AddAtomsWithConnectivity( ens_itr_b->GetAtomVector(), empty_bonds);
      }
      FragmentComplete pseudo_mol_b( vector_b, "");

      storage::Pair< FragmentComplete, FragmentComplete> pseudo_molecules;
      pseudo_molecules.First() = pseudo_mol_a;
      pseudo_molecules.Second() = pseudo_mol_b;

      // Compute atom volumes for each molecule
      size_t n_atoms_a( pseudo_mol_a.GetSize()), n_atoms_b( pseudo_mol_b.GetSize());
      linal::Vector< float> properties_a( n_atoms_a), properties_b( n_atoms_b);

      descriptor::AtomVolume atom_volume;
      atom_volume.SetObject( pseudo_mol_a);
      for
      (
          descriptor::Iterator< AtomConformationalInterface> itr_a( descriptor::Type( 1, true, descriptor::Type::e_Symmetric), pseudo_mol_a);
          itr_a.NotAtEnd();
          ++itr_a
      )
      {
        properties_a( itr_a.GetPosition()) = atom_volume( itr_a)( 0);
      }
      atom_volume.SetObject( pseudo_mol_b);
      for
      (
          descriptor::Iterator< AtomConformationalInterface> itr_b( descriptor::Type( 1, true, descriptor::Type::e_Symmetric), pseudo_mol_b);
          itr_b.NotAtEnd();
          ++itr_b
      )
      {
        properties_b( itr_b.GetPosition()) = atom_volume( itr_b)( 0);
      }

      // Compute maximum value for all row elements
      storage::Vector< storage::Triplet< size_t, size_t, float>> max_volume_atom_pairs;
      linal::Matrix< float> volume_sum_matrix( n_atoms_a, n_atoms_b);
      for( size_t a( 0); a < n_atoms_a; ++a)
      {
        storage::Triplet< size_t, size_t, float> max_value( util::GetUndefinedSize_t(), util::GetUndefinedSize_t(), math::GetLowestBoundedValue< float>());
        for( size_t b( 0); b < n_atoms_b; ++b)
        {
          volume_sum_matrix( a, b) = properties_a( a) + properties_b( b);
          if( volume_sum_matrix( a, b) > max_value.Third())
          {
            max_value.First() = a;
            max_value.Second() = b;
            max_value.Third() = volume_sum_matrix(a, b);
          }
        }
        max_volume_atom_pairs.PushBack( max_value);
      }
      // Sort by max sum
      CompareTriplets triplet_compare;
      max_volume_atom_pairs.Sort< CompareTriplets>( triplet_compare);

      // Identify max values for each row, remove columns redundant with max pairs
      storage::Vector< size_t> col_seen_before;
      storage::Vector< storage::Triplet< size_t, size_t, float>> best_pairs;
      for
      (
        storage::Vector< storage::Triplet< size_t, size_t, float> >::iterator
        atom_pair_itr( max_volume_atom_pairs.Begin()), atom_pair_itr_end( max_volume_atom_pairs.End());
        atom_pair_itr != atom_pair_itr_end;
        ++atom_pair_itr
      )
      {
        // check to see if the max value from a column has already been used or not
        if( !col_seen_before.IsEmpty())
        {
          for( size_t i( 0); i < col_seen_before.GetSize(); ++i)
          {
            if( atom_pair_itr->Second() == col_seen_before( i))
            {
              // max value already associated with this column
              // move on to the next triplet
              break;
            }
            else
            {
              // check next column value in col_seen_before
              continue;
            }
          }

          //Column value has not been see yet; save
          col_seen_before.PushBack( atom_pair_itr->Second());
          best_pairs.PushBack( *atom_pair_itr);
        }
        else
        {
          col_seen_before.PushBack( atom_pair_itr->Second());
          best_pairs.PushBack( *atom_pair_itr);

          // move on to the next triplet
          continue;
        }
      }

      // Now we need to repeat this step to get the next-most-maximum value for the previous redundant columns values
      // For the time being, though, let's just use what we got for bondalign3 rigid alignment
      storage::Pair< storage::Vector< storage::Triplet< size_t, size_t, float>>, storage::Pair< FragmentComplete, FragmentComplete>> holy_cow_output;
      holy_cow_output.Second() = pseudo_molecules;
      holy_cow_output.First() = best_pairs;

      return holy_cow_output;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels

    io::Serializer ConformationComparisonMultiAlign::GetSerializer() const
    {
      io::Serializer member_data( ConformationComparisonPsiFlexField::GetSerializer());
      member_data.SetClassDescription( "Align multiple molecules simultaneously");
      member_data.AddInitializer
       (
         "file tag",
         "tag in the output filename",
         io::Serialization::GetAgent( &m_FileTagMulti),
         "MultiAlign"
       );
      member_data.AddOptionalInitializer
      (
        "start position",
        "start ligand centered at this position instead of pocket center",
        io::Serialization::GetAgent( &m_StartPos)
      );
      member_data.AddOptionalInitializer
      (
        "ligfit movie filename",
        "supply a filename to activate output of each frame of the ligand-pocket geometric minimization",
        io::Serialization::GetAgent( &m_MovieLigFitFilename)
      );
      return member_data;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ConformationComparisonMultiAlign::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return ConformationComparisonPsiFlexField::ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }

  } // namespace chemistry
} // namespace bcl
