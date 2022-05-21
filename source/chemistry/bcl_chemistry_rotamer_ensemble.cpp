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
#include "chemistry/bcl_chemistry_rotamer_ensemble.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "quality/bcl_quality_rmsd_preprocessor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief constructor
    //! @param MOLECULE Molecule that needs to be added to the cluster
    //! @param GRAPH method to be used for finding cluster center
    //! @param BIN_SIZE dihedral bin size used for conformation determination
    RotamerEnsemble::RotamerEnsemble
    (
      const FragmentConformationShared &MOLECULE,
      const graph::ConstGraph< size_t, size_t> &GRAPH,
      const double WEIGHT,
      double BIN_SIZE,
      const bool &SIMULATED
    ) :
      m_Graph( GRAPH),
      m_BinSize( BIN_SIZE),
      m_Instances( size_t( 0)),
      m_SimulatedInstances( 0.0)
    {
      operator ()( MOLECULE, WEIGHT, SIMULATED);
    }

    //! @brief Clone function
    //! @return pointer to new RotamerEnsemble
    RotamerEnsemble *RotamerEnsemble::Clone() const
    {
      return new RotamerEnsemble( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RotamerEnsemble::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns average dihedral angles for conformations that belong to this cluster
    //! @return average dihedral angles
    const linal::Vector< double> &RotamerEnsemble::GetAverage() const
    {
      return m_AverageBinAngles.GetAverage();
    }

    //! @brief returns the size of this cluster
    //! @return the size of the cluster
    const double &RotamerEnsemble::GetWeights() const
    {
      return m_AverageBinAngles.GetWeight();
    }

    //! @brief returns the number of times fragment has been seen in the structure database
    //! @return the number of times fragment has been seen in the structure database
    const size_t &RotamerEnsemble::GetNumberInstances() const
    {
      return m_Instances;
    }

    //! @brief return the conformation graph colored by dihedral bins for this cluster
    //! @return the conformatin graph for this cluster
    const graph::ConstGraph< size_t, size_t> &RotamerEnsemble::GetGraph() const
    {
      return m_Graph;
    }

    //! @brief return the conformation which is the center of this cluster
    //! @return the center of this cluster
    const FragmentConformationShared &RotamerEnsemble::GetClusterCenter() const
    {
      return m_Center;
    }

    //! @brief return members of the rotamer ensemble
    //! @return all members of rotamer ensemble
    const storage::List< FragmentConformationShared> &RotamerEnsemble::GetRotamers() const
    {
      return m_Conformers;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief add a conformation to the cluster info
    //! @param FRAGMENT the conformation to add
    void RotamerEnsemble::operator()( const FragmentConformationShared &FRAGMENT, const double WEIGHT, const bool &SIMULATED)
    {
      if( !SIMULATED)
      {
        ++m_Instances;
        if( m_Conformers.IsEmpty() || random::GetGlobalRandom().Double() < WEIGHT)
        {
          // avoid retaining excessive numbers of redundant rotamers by only stochastically adding them based on weight.
          m_Conformers.PushBack( FRAGMENT);
        }
        m_AverageBinAngles.AddWeightedObservation
        (
          ConformationComparisonByDihedralBins( m_BinSize).DetermineWrappedAroundAngles
          (
            linal::Vector< double>( PriorityDihedralAngles()( FRAGMENT).First()), false
          ),
          WEIGHT
        );
      }
      if( SIMULATED)
      {
        m_SimulatedInstances += WEIGHT;
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the center of this cluster from the cluster info
    //! @param MAX_COUNTS maximum number of examples over which to calculate cluster center
    //! @return void
    void RotamerEnsemble::CalculateCenter( size_t MAX_COUNTS)
    {
      // create a object to store RMSD
      {
        size_t n_conformers( m_Conformers.GetSize());
        if( n_conformers > MAX_COUNTS)
        {
          std::vector< storage::List< FragmentConformationShared>::iterator> ens_list;
          ens_list.reserve( n_conformers);
          for( auto itr( m_Conformers.Begin()), itr_end( m_Conformers.End()); itr != itr_end; ++itr)
          {
            ens_list.push_back( itr);
          }
          random::RandomShuffle( ens_list.begin(), ens_list.end());
          storage::List< FragmentConformationShared> new_frag_ensemble;
          for( auto itr( ens_list.begin()), itr_end( ens_list.begin() + MAX_COUNTS); itr != itr_end; ++itr)
          {
            new_frag_ensemble.InternalData().splice( new_frag_ensemble.End(), m_Conformers.InternalData(), *itr);
          }
          std::swap( new_frag_ensemble, m_Conformers);
        }
      }

      size_t row_number( 0);

      // create a vector to store parent molecule indices
      storage::Vector< double> molecule_indices;
      linal::Vector< double> rmsd_sum( m_Conformers.GetSize(), double( 0.0));

      // iterate through the conformations that belong to this cluster
      size_t outer_count( 0);
      storage::List< quality::RMSDPreprocessor> rmsd_helpers;
      for
      (
        storage::List< FragmentConformationShared>::const_iterator
          itr_outer( m_Conformers.Begin()), itr_outer_end( m_Conformers.End());
        itr_outer != itr_outer_end;
        ++itr_outer
      )
      {
        rmsd_helpers.PushBack( quality::RMSDPreprocessor( itr_outer->GetAtomCoordinates(), true));
      }
      auto itr_rmsd_helpers( rmsd_helpers.Begin());
      for
      (
        storage::List< FragmentConformationShared>::const_iterator
          itr_outer( m_Conformers.Begin()), itr_outer_end( m_Conformers.End());
        itr_outer != itr_outer_end;
        ++itr_outer, ++row_number, ++outer_count, ++itr_rmsd_helpers
      )
      {
        // if there are large counts then dont bother averaging over all of them as bin width is narrow, so just take a sample
        if( outer_count == MAX_COUNTS)
        {
          break;
        }

        // start another iterations through the conformations to compare to conformatins being iterated int the outer loop
        size_t inner_row_number( 0);
        auto itr_rmsd_helpers_p( rmsd_helpers.Begin());
        for
        (
          storage::List< FragmentConformationShared>::const_iterator itr_inner( m_Conformers.Begin());
          itr_inner != itr_outer;
          ++itr_inner, ++inner_row_number, ++itr_rmsd_helpers_p
        )
        {
          // calculate the rmsd between conformations being iterated in outer and inner loop
          double rmsd( itr_rmsd_helpers->SuperimposedRMSD( *itr_rmsd_helpers_p));

          // sum all the rmsd's for the conformation being iterated in the outer loop
          rmsd_sum( row_number) += rmsd;
          rmsd_sum( inner_row_number) += rmsd;
        }
      }

      // find the index with the minimum rmsd sum and store the rmsd and the conformation
      // we would just use min_element here (and skip the find), but doing the find afterwards ensures better consistency
      // across compilers, some of which use <= in std::min_element.
      size_t lowest_index
      (
        rmsd_sum.GetSize() > size_t( 2)
        ? std::find( rmsd_sum.Begin(), rmsd_sum.End(), *std::min_element( rmsd_sum.Begin(), rmsd_sum.End())) - rmsd_sum.Begin()
        : 0
      );
      storage::List< FragmentConformationShared>::iterator itr_cluster_center( m_Conformers.Begin());
      storage::AdvanceIterator( itr_cluster_center, m_Conformers.End(), lowest_index);

      m_Center = FragmentConformationShared( *itr_cluster_center);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RotamerEnsemble::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "Cannot read " + GetClassIdentifier(), -1);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RotamerEnsemble::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      BCL_Exit( "Cannot write " + GetClassIdentifier(), -1);
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
