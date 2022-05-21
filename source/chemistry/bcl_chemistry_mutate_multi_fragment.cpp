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
#include "math/bcl_math_statistics.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_mutate_multi_fragment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_dihedral_angles.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_mutate_chirality.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "chemistry/bcl_chemistry_valence_handler.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_histogram.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateMultiFragment::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateMultiFragment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateMultiFragment::MutateMultiFragment()
    {
    }

    //! @brief constructor taking the member variable parameters
    //! @param FRAGMENT fragment which the mutate mutates
    //! @param GRAPH constitution graph of the molecule of interest
    //! @param MOLECULE molecule of interest
    MutateMultiFragment::MutateMultiFragment
    (
      const util::ShPtrVector< MutateDihedralsInterface> &FRAGMENT_MUTATES,
      const util::ShPtrVector< MutateDihedralsInterface> &DIHEDRAL_BOND_MUTATES_RANDOM,
      const double &RANDOM_DIHEDRAL_MUTATE_FREQUENCY,
      const size_t &N_DIHEDRALS
    ) :
      m_DihedralMutators( N_DIHEDRALS),
      m_ProbabilityDistributions( N_DIHEDRALS),
      m_ProbabilityDistributionSums( N_DIHEDRALS, double( 0.0))
    {
      storage::Vector< double> n_single( N_DIHEDRALS, double( 0));
      storage::Vector< double> sum_multi( N_DIHEDRALS, double( 0.0));
      for( auto itr( FRAGMENT_MUTATES.Begin()), itr_end( FRAGMENT_MUTATES.End()); itr != itr_end; ++itr)
      {
        if( ( *itr)->GetMoleculeDihedralBondIndices().IsEmpty())
        {
          continue;
        }
        double fract( 1.0 / double( ( *itr)->GetMoleculeDihedralBondIndices().GetSize()));
        for
        (
          auto itr_bnd( ( *itr)->GetMoleculeDihedralBondIndices().Begin()),
               itr_bnd_end( ( *itr)->GetMoleculeDihedralBondIndices().End());
          itr_bnd != itr_bnd_end;
          ++itr_bnd
        )
        {
          ( ( *itr)->GetMoleculeDihedralBondIndices().GetSize() == 1 ? n_single : sum_multi)( *itr_bnd) += fract;
          m_DihedralMutators( *itr_bnd).PushBack( *itr);
          m_ProbabilityDistributions( *itr_bnd).PushBack( fract);
          m_ProbabilityDistributionSums( *itr_bnd) += fract;
        }
      }
      if( RANDOM_DIHEDRAL_MUTATE_FREQUENCY)
      {
        for( auto itr( DIHEDRAL_BOND_MUTATES_RANDOM.Begin()), itr_end( DIHEDRAL_BOND_MUTATES_RANDOM.End()); itr != itr_end; ++itr)
        {
          if( ( *itr)->GetMoleculeDihedralBondIndices().IsEmpty())
          {
            continue;
          }
          double fract( 1.0 / double( ( *itr)->GetMoleculeDihedralBondIndices().GetSize()));
          for
          (
            auto itr_bnd( ( *itr)->GetMoleculeDihedralBondIndices().Begin()),
                 itr_bnd_end( ( *itr)->GetMoleculeDihedralBondIndices().End());
            itr_bnd != itr_bnd_end;
            ++itr_bnd
          )
          {
            m_DihedralMutators( *itr_bnd).PushBack( *itr);
            if( !m_ProbabilityDistributionSums( *itr_bnd))
            {
              m_ProbabilityDistributions( *itr_bnd).PushBack( fract);
              m_ProbabilityDistributionSums( *itr_bnd) += fract;
            }
            else
            {
              double val
              (
                m_ProbabilityDistributionSums( *itr_bnd) * RANDOM_DIHEDRAL_MUTATE_FREQUENCY
                / ( 1.0 - RANDOM_DIHEDRAL_MUTATE_FREQUENCY)
              );
              m_ProbabilityDistributions( *itr_bnd).PushBack( val);
              m_ProbabilityDistributionSums( *itr_bnd) += val;
            }
          }
        }
      }
      for( size_t dihedral_id( 0); dihedral_id < N_DIHEDRALS; ++dihedral_id)
      {
        if( !m_ProbabilityDistributions( dihedral_id).IsEmpty())
        {
          if( n_single( dihedral_id))
          {
            m_ChainBondIndices.PushBack( dihedral_id);
          }
          else
          {
            m_RingBondIndices.PushBack( dihedral_id);
          }
        }
      }
    }

    //! @brief Clone function
    //! @return pointer to new MutateMultiFragment
    MutateMultiFragment *MutateMultiFragment::Clone() const
    {
      return new MutateMultiFragment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateMultiFragment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateMultiFragment::GetScheme() const
    {
      static const std::string s_name( "MutateMolecule");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a conformation and returns a mutated conformation
    //! @param MOLECULE conformation of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< FragmentComplete> MutateMultiFragment::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      size_t n_dihedrals( m_DihedralMutators.GetSize());
      storage::Vector< size_t> index_rings( m_RingBondIndices);
      storage::Vector< size_t> index_chains( m_ChainBondIndices);
      index_rings.Shuffle();
      index_chains.Shuffle();
      index_rings.Append( index_chains);
      auto index_v( index_rings);
      storage::Vector< size_t> need_to_change( n_dihedrals, size_t( 1));
      util::ShPtr< FragmentComplete> new_molecule( MOLECULE.Clone());
      for
      (
        auto itr_index( index_v.Begin()), itr_index_end( index_v.End());
        itr_index != itr_index_end;
        ++itr_index
      )
      {
        if( !need_to_change( *itr_index) || !m_ProbabilityDistributionSums( *itr_index))
        {
          continue;
        }

        // initialize the position holder
        auto itr_mutate( m_DihedralMutators( *itr_index).Begin());
        bool all_need_change( false);
        // try to find a dihedral mutate up to 20 times which doesn't have all dihedral bonds already covered.
        for( size_t try_iteration( 0), max_iterations( 20); try_iteration < max_iterations; ++try_iteration)
        {
          itr_mutate = m_DihedralMutators( *itr_index).Begin();
          // initiate a random number from 0 to 1.0
          double temp( random::GetGlobalRandom().Random< double>( m_ProbabilityDistributionSums( *itr_index)));

          for
          (
            // initialize the pointers inside the distribution
            auto ptr( m_ProbabilityDistributions( *itr_index).Begin()),
                 ptr_end( m_ProbabilityDistributions( *itr_index).End());
            ptr != ptr_end;
            ++ptr, ++itr_mutate
          )
          {
            // subtract the random number generated from the value pointed to in m_Distribution1D
            temp -= *ptr;

            // if the difference is negative or zero, then return the current position in the vector
            // because the larger the probability the greater the chance of causing this difference to be negative
            if( temp <= 0.0)
            {
              break;
            }
          }
          all_need_change = true;
          for
          (
            auto itr_bnd( ( *itr_mutate)->GetMoleculeDihedralBondIndices().Begin()),
                 itr_bnd_end( ( *itr_mutate)->GetMoleculeDihedralBondIndices().End());
            itr_bnd != itr_bnd_end && all_need_change;
            ++itr_bnd
          )
          {
            all_need_change = need_to_change( *itr_bnd);
          }

          if( all_need_change || try_iteration > size_t( 10))
          {
            auto result( ( *itr_mutate)->operator ()( *new_molecule));
            if( result.GetArgument().IsDefined() && result.GetArgument()->GetSize() != size_t( 0))
            {
              new_molecule = result.GetArgument();
              for
              (
                auto itr_bnd( ( *itr_mutate)->GetMoleculeDihedralBondIndices().Begin()),
                     itr_bnd_end( ( *itr_mutate)->GetMoleculeDihedralBondIndices().End());
                itr_bnd != itr_bnd_end;
                ++itr_bnd
              )
              {
                need_to_change( *itr_bnd) = size_t( 0);
              }
              break;
            }
          }
        }
      }
      return math::MutateResult< FragmentComplete>( new_molecule, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateMultiFragment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DihedralMutators, ISTREAM);
      io::Serialize::Read( m_ProbabilityDistributions, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateMultiFragment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DihedralMutators, OSTREAM, INDENT);
      io::Serialize::Write( m_ProbabilityDistributions, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
