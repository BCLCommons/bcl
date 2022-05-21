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
#include "chemistry/bcl_chemistry_atom_one_four_interaction_score.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_atom_clash_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_statistics.h"

#define BCL_PROFILE_AtomClashScore
#ifdef BCL_PROFILE_AtomClashScore
#include "util/bcl_util_stopwatch.h"
#endif

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomClashScore::s_Instance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< FragmentComplete, double> >::AddInstance( new AtomClashScore( true))
    );
    const util::SiPtr< const util::ObjectInterface> AtomClashScore::s_HeavyInstance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< FragmentComplete, double> >::AddInstance( new AtomClashScore( false))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor. Can take whether to consider hydrogens or tolerance
    AtomClashScore::AtomClashScore( const bool &CONSIDER_H, const bool &DISPLAY) :
      m_ConsiderH( CONSIDER_H),
      m_Display( DISPLAY)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AtomClashScore
    AtomClashScore *AtomClashScore::Clone() const
    {
      return new AtomClashScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomClashScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AtomClashScore::GetAlias() const
    {
      static const std::string s_name( "AtomClash"), s_heavy_name( "HeavyAtomClash");
      return m_ConsiderH ? s_name : s_heavy_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param MOLECULE molecule that needs to scored
    //! @return clash score for the given molecule
    double AtomClashScore::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      return operator()( static_cast< const ConformationInterface &>( MOLECULE));
    }

    //! @brief evaluate clashes for given atoms
    //! @param MOLECULE molecule that needs to scored
    //! @return clash score for the given atoms
    double AtomClashScore::operator()
    (
      const ConformationInterface &MOLECULE
    ) const
    {
      return GetCachedInfo( MOLECULE).First();
    }

    //! @brief get the clashing pairs
    const storage::Vector< storage::Triplet< size_t, size_t, double> > &AtomClashScore::GetClashingPairs
    (
      const ConformationInterface &MOLECULE
    ) const
    {
      return GetCachedInfo( MOLECULE).Second();
    }

    //! @brief get the score for just a previous set of clashing pairs (saves refinding the clashing pairs, but may miss
    //!        some newly clashing pairs)
    storage::Vector< storage::Triplet< size_t, size_t, double> >::iterator AtomClashScore::UpdateClashingPairs
    (
      const AtomVector< AtomComplete> &MOLECULE,
      storage::Vector< storage::Triplet< size_t, size_t, double> > &CLASHING_PAIRS
    ) const
    {
      auto itr_max( CLASHING_PAIRS.Begin());
      double cs_sum( 0.0);
      for( auto itr( CLASHING_PAIRS.Begin()), itr_end( CLASHING_PAIRS.End()); itr != itr_end; ++itr)
      {
        itr->Third() =
          std::max
          (
            0.0,
            m_VdwSumMaxDistances( itr->First(), itr->Second())
            - linal::Distance( MOLECULE( itr->First()).GetPosition(), MOLECULE( itr->Second()).GetPosition())
          );
        cs_sum += itr->Third();
      }
      double cs_chosen( random::GetGlobalRandom().Random( cs_sum));
      for( auto itr( CLASHING_PAIRS.Begin()), itr_end( CLASHING_PAIRS.End()); itr != itr_end; ++itr)
      {
        cs_chosen -= itr->Third();
        if( cs_chosen < 0.0)
        {
          itr_max = itr;
          break;
        }
      }
      return itr_max;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomClashScore::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Detects clashes between " + std::string( m_ConsiderH ? "all" : "heavy") + " atoms depending on VdW radii "
      );
      return serializer;
    }

    //! @brief Test whether this molecule is the same (constitutionally) as the molecule for which the state of this
    //!        class currently can handle
    bool AtomClashScore::IsMoleculeInfoCached( const ConformationInterface &CONF) const
    {
      if
      (
        m_LastAtomTypes.GetSize() != CONF.GetNumberAtoms()
        || m_LastBondInfo.GetSize() != CONF.GetNumberBonds()
      )
      {
        return false;
      }

      storage::Vector< AtomType>::const_iterator itr_types( m_LastAtomTypes.Begin());
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( CONF.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++itr_types
      )
      {
        if( *itr_types != itr_atoms->GetAtomType())
        {
          return false;
        }
      }

      if( m_LastBondInfo != CONF.GetBondInfo())
      {
        return false;
      }

      return true;
    }

    //! @brief Update molecule change the molecule that this class will compute the clash score for
    //! @param MOL molecule of interest
    void AtomClashScore::UpdateMolecule( const ConformationInterface &CONF) const
    {
      m_LastBondInfo = CONF.GetBondInfo();

      const size_t na( CONF.GetNumberAtoms());
      if( m_VdwSumMaxDistances.GetNumberRows() < na)
      {
        m_VdwSumMaxDistances = linal::Matrix< float>( na, na, float( 0.0));
      }

      const size_t n_atoms( CONF.GetNumberAtoms());
      m_LastAtomTypes.Resize( n_atoms);
      m_MaxVdwRadius = 0.0;
      size_t atom_index( 0);
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( CONF.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++atom_index
      )
      {
        m_LastAtomTypes( atom_index) = itr_atoms->GetAtomType();
      }

      ConformationGraphConverter::t_AtomGraph molgraph
      (
        ConformationGraphConverter::CreateGraphWithAtoms( CONF)
      );

      atom_index = 0;
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( CONF.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++atom_index
      )
      {
        if( !m_ConsiderH && itr_atoms->GetElementType()->GetAtomicNumber() == size_t( 1))
        {
          auto mutable_row( m_VdwSumMaxDistances.GetRow( atom_index));
          mutable_row = float( 0.0);
          continue;
        }
        auto dists_sp( graph::Connectivity::DistancesToOtherVertices( molgraph, atom_index, size_t( 5)));
        const storage::Vector< size_t> &dists( *dists_sp);
        iterate::Generic< const AtomConformationalInterface> itr_atoms_b( CONF.GetAtomsIterator());
        for( size_t atom_index_b( 0); atom_index_b < atom_index; ++atom_index_b, ++itr_atoms_b)
        {
          if( !m_ConsiderH && itr_atoms_b->GetElementType()->GetAtomicNumber() == size_t( 1))
          {
            m_VdwSumMaxDistances( atom_index, atom_index_b) = float( 0.0);
            continue;
          }
          if( dists( atom_index_b) < size_t( 4))
          {
            if( dists( atom_index_b) >= size_t( 2))
            {
              m_VdwSumMaxDistances( atom_index, atom_index_b) = m_VdwSumMaxDistances( atom_index_b, atom_index) =
                itr_atoms->GetElementType()->GetProperty( ElementTypeData::e_CovalentRadius)
                + itr_atoms_b->GetElementType()->GetProperty( ElementTypeData::e_CovalentRadius)
                + 0.1;
            }
            else
            {
              m_VdwSumMaxDistances( atom_index, atom_index_b) = m_VdwSumMaxDistances( atom_index_b, atom_index) = 0.0;
            }
          }
          else
          {
            m_VdwSumMaxDistances( atom_index, atom_index_b) = m_VdwSumMaxDistances( atom_index_b, atom_index) =
              itr_atoms->GetElementType()->GetProperty( ElementTypeData::e_DaltonVdwRadius)
              + itr_atoms_b->GetElementType()->GetProperty( ElementTypeData::e_DaltonVdwRadius)
              - 0.7;
            m_MaxVdwRadius = std::max( m_VdwSumMaxDistances( atom_index, atom_index_b), m_MaxVdwRadius);
          }
        }
      }
    }

    //! @brief prepare the cache entry for a given molecule
    const storage::Pair< double, storage::Vector< storage::Triplet< size_t, size_t, double> > > &AtomClashScore::GetCachedInfo( const ConformationInterface &MOLECULE) const
    {
      if( !m_Display)
      {
        auto itr( m_Cache.Find( size_t( &MOLECULE)));
        if( itr != m_Cache.End())
        {
          return itr->second;
        }
      }
      #ifdef BCL_PROFILE_AtomClashScore
      static util::Stopwatch s_timer( "  AtomClashScore ClashComputation", util::Message::e_Standard, true);
      s_timer.Start();
      #endif

      if( !IsMoleculeInfoCached( MOLECULE))
      {
        UpdateMolecule( MOLECULE);
      }

      storage::Pair< double, storage::Vector< storage::Triplet< size_t, size_t, double> > > &new_val
      (
        m_Cache.Insert
        (
          std::make_pair
          (
            size_t( &MOLECULE),
            storage::Pair< double, storage::Vector< storage::Triplet< size_t, size_t, double> > >()
          )
        ).first->second
      );

      VoxelGridAtom vga( m_MaxVdwRadius);
      vga.SetObjects
      (
        util::SiPtrVector< const AtomConformationalInterface>
        (
          MOLECULE.GetAtomsIterator(),
          MOLECULE.GetAtomsIterator().End()
        )
      );

      double clash_sum( 0.0);
      auto neigh( vga.GetNeighbors( m_MaxVdwRadius));
      for( auto itr_neigh( neigh.Begin()), itr_neigh_end( neigh.End()); itr_neigh != itr_neigh_end; ++itr_neigh)
      {
        size_t index_a( MOLECULE.GetAtomIndex( *itr_neigh->First()));
        size_t index_b( MOLECULE.GetAtomIndex( *itr_neigh->Second()));

        // compute the minimum distance for these atoms to be apart, assuming they meet all the requirements for the
        // vdw to be valid, as was checked earlier, except the requirement that they not have opposite charge
        double vdw_max_distance( m_VdwSumMaxDistances( index_a, index_b));
        if( itr_neigh->Third() < vdw_max_distance)
        {
          double this_score( 0.0);
          this_score = vdw_max_distance - itr_neigh->Third();
          clash_sum += this_score;
          if( index_a > index_b)
          {
            std::swap( index_a, index_b);
          }
          new_val.Second().PushBack( storage::Triplet< size_t, size_t, double>( index_a, index_b, this_score));

          if( ( m_Display && this_score != 0.0) || util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Verbose)
          {
            BCL_MessageStd
            (
              "Type: " + itr_neigh->First()->GetAtomType().GetName() + " " + itr_neigh->Second()->GetAtomType().GetName()
              + " Dist: " + util::Format()( itr_neigh->Third())
              + " Vdw: " + util::Format()( vdw_max_distance)
              + " Score: " + util::Format()( this_score)
              + " A/B " + util::Format()( index_a) + "/" + util::Format()( index_b)
            );
          }
        }
      }
      new_val.First() = clash_sum / double( MOLECULE.GetNumberAtoms());
      new_val.Second().Sort( storage::GreaterThanThird());
      MOLECULE.GetChangeSignal().Connect
      (
        const_cast< AtomClashScore *>( this),
        &AtomClashScore::RemoveResultFromCache
      );
      #ifdef BCL_PROFILE_AtomClashScore
      s_timer.Stop();
      #endif
      return new_val;
    }

    //! @brief remove results for this object from the cache
    //! @param ARGUMENT argument that should be removed
    void AtomClashScore::RemoveResultFromCache( const ConformationInterface &ARGUMENT)
    {
      // find result for that Argument
      auto itr( m_Cache.Find( size_t( &ARGUMENT)));

      // skip if there is no result cached
      if( itr != m_Cache.End())
      {
        m_Cache.RemoveElement( itr);
        ARGUMENT.GetChangeSignal().Disconnect( const_cast< AtomClashScore *>( this));
      }
      else
      {
        BCL_MessageStd( "Caching is messed up for atom clash score");
      }
    }

  } // namespace chemistry
} // namespace bcl
