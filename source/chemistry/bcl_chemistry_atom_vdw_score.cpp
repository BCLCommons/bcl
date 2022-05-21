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
#include "chemistry/bcl_chemistry_atom_vdw_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_statistics.h"

#define BCL_PROFILE_AtomVdwScore
#ifdef BCL_PROFILE_AtomVdwScore
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
    const util::SiPtr< const util::ObjectInterface> AtomVdwScore::s_Instance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< FragmentComplete, double> >::AddInstance( new AtomVdwScore())
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomVdwScore::s_ClashInstance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< FragmentComplete, double> >::AddInstance( new AtomVdwScore( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new AtomVdwScore
    AtomVdwScore *AtomVdwScore::Clone() const
    {
      return new AtomVdwScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomVdwScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AtomVdwScore::GetAlias() const
    {
      static const std::string s_name( "AtomVdw"), s_cname( "AtomClashVdw");
      return m_ClashOnly ? s_cname : s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param MOLECULE molecule that needs to scored
    //! @return clash score for the given molecule
    double AtomVdwScore::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      return operator()( static_cast< const ConformationInterface &>( MOLECULE));
    }

    namespace
    {
      // Score Lennard-Jones potential
      double ScoreLJ( const double &MAX_VDW_SUM, const double &DISTANCE)
      {
        double z( MAX_VDW_SUM / DISTANCE);
        double zsq( z * z);
        double zsix( zsq * zsq * zsq);
        return ( zsix * zsix - 2.0 * zsix);
      }

    }

    //! @brief evaluate clashes for given atoms
    //! @param MOLECULE molecule that needs to scored
    //! @return clash score for the given atoms
    double AtomVdwScore::operator()
    (
      const ConformationInterface &MOLECULE
    ) const
    {
      #ifdef BCL_PROFILE_AtomVdwScore
      static util::Stopwatch s_timer( "vdw score time", util::Message::e_Standard, true);
      s_timer.Start();
      #endif
//      util::ObjectDataLabel label( this->GetCompleteSerializer().GetLabel());
//      auto cache_val( MOLECULE.FindInCache( label));
//      if( cache_val.IsDefined())
//      {
//        #ifdef BCL_PROFILE_AtomVdwScore
//        s_timer.Stop();
//        #endif
//        return cache_val->operator ()( 0);
//      }
      if( !IsMoleculeInfoCached( MOLECULE))
      {
        UpdateMolecule( MOLECULE);
      }
      VoxelGridAtom vga( m_MaxVdwRadius);
      vga.SetObjects
      (
        util::SiPtrVector< const AtomConformationalInterface>
        (
          MOLECULE.GetAtomsIterator(),
          MOLECULE.GetAtomsIterator().End()
        )
      );

      double vdw_sum( 0.0);
      auto neigh( vga.GetNeighbors( m_MaxVdwRadius));
      for( auto itr_neigh( neigh.Begin()), itr_neigh_end( neigh.End()); itr_neigh != itr_neigh_end; ++itr_neigh)
      {
        size_t index_a( MOLECULE.GetAtomIndex( *itr_neigh->First()));
        size_t index_b( MOLECULE.GetAtomIndex( *itr_neigh->Second()));

        // compute the minimum distance for these atoms to be apart, assuming they meet all the requirements for the
        // vdw to be valid, as was checked earlier, except the requirement that they not have opposite charge
        double vdw_max_distance( m_VdwSumMaxDistances( index_a, index_b) + ( m_ClashOnly ? 0.0 : 10.0));
        if( itr_neigh->Third() < vdw_max_distance)
        {
          double this_score( 0.0);
          if( !m_ClashOnly)
          {
            double lj_score( ScoreLJ( m_LJ( index_a, index_b), itr_neigh->Third()));
            lj_score *= m_PVdw( index_a, index_b);
            this_score = -lj_score;
          }
          else
          {
            this_score = vdw_max_distance - itr_neigh->Third();
          }
          vdw_sum += this_score;

          if( m_Display && this_score != 0.0)
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
          else if( !m_ClashOnly)
          {
            BCL_MessageDbg
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

      #ifdef BCL_PROFILE_AtomVdwScore
      s_timer.Stop();
      #endif
      //MOLECULE.Cache( label, linal::Vector< float>( size_t( 1), vdw_sum / MOLECULE.GetNumberAtoms()));
      return vdw_sum / double( MOLECULE.GetNumberAtoms());
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomVdwScore::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Vdw force score between atoms"
      );
      return serializer;
    }

    //! @brief Test whether this molecule is the same (constitutionally) as the molecule for which the state of this
    //!        class currently can handle
    bool AtomVdwScore::IsMoleculeInfoCached( const ConformationInterface &CONF) const
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
    void AtomVdwScore::UpdateMolecule( const ConformationInterface &CONF) const
    {
      m_LastBondInfo = CONF.GetBondInfo();

      const size_t na( CONF.GetNumberAtoms());
      if( m_VdwSumMaxDistances.GetNumberRows() < na)
      {
        m_VdwSumMaxDistances = linal::Matrix< float>( na, na, float( 0.0));
        if( !m_ClashOnly)
        {
          m_PVdw = linal::Matrix< float>( na, na, float( 0.0));
          m_LJ = linal::Matrix< float>( na, na, float( 0.0));
        }
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
        auto dists_sp( graph::Connectivity::DistancesToOtherVertices( molgraph, atom_index, size_t( 5)));
        const storage::Vector< size_t> &dists( *dists_sp);
        iterate::Generic< const AtomConformationalInterface> itr_atoms_b( CONF.GetAtomsIterator());
        for( size_t atom_index_b( 0); atom_index_b < atom_index; ++atom_index_b, ++itr_atoms_b)
        {
          if( dists( atom_index_b) < size_t( 4))
          {
            m_VdwSumMaxDistances( atom_index, atom_index_b) = m_VdwSumMaxDistances( atom_index_b, atom_index) = 0.0;
          }
          else
          {
            m_VdwSumMaxDistances( atom_index, atom_index_b) = m_VdwSumMaxDistances( atom_index_b, atom_index) =
              itr_atoms->GetElementType()->GetProperty( ElementTypeData::e_DaltonVdwRadius)
              + itr_atoms_b->GetElementType()->GetProperty( ElementTypeData::e_DaltonVdwRadius);
            if( m_ClashOnly)
            {
              m_VdwSumMaxDistances( atom_index, atom_index_b) = m_VdwSumMaxDistances( atom_index_b, atom_index) -= 0.7;
              m_MaxVdwRadius = std::max( m_VdwSumMaxDistances( atom_index, atom_index_b), m_MaxVdwRadius);
            }
            else
            {
              m_VdwSumMaxDistances( atom_index, atom_index_b) = m_VdwSumMaxDistances( atom_index_b, atom_index) += 0.7;
              m_PVdw( atom_index, atom_index_b) = m_PVdw( atom_index_b, atom_index) =
                -0.5 *
                (
                  itr_atoms->GetElementType()->GetProperty( ElementTypeData::e_LJEpsilon)
                  + itr_atoms_b->GetElementType()->GetProperty( ElementTypeData::e_LJEpsilon)
                );
              m_LJ( atom_index, atom_index_b) = m_LJ( atom_index_b, atom_index) =
                (
                  itr_atoms->GetElementType()->GetProperty( ElementTypeData::e_LJRadius)
                  + itr_atoms_b->GetElementType()->GetProperty( ElementTypeData::e_LJRadius)
                );
              m_MaxVdwRadius = std::max( m_LJ( atom_index, atom_index_b), m_MaxVdwRadius);
            }
          }
        }
      }
      if( !m_ClashOnly)
      {
        m_MaxVdwRadius += 1.0;
      }
    }

  } // namespace chemistry
} // namespace bcl
