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
#include "chemistry/bcl_chemistry_bond_angle_assignment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_running_average.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BondAngleAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new BondAngleAssignment())
    );

    //! @brief get the pseudocount that is used for all possible rotamers
    //! @return the pseudocount that is used for all possible rotamers
    double &BondAngleAssignment::GetPseudocount()
    {
      static double s_pseudocount( 4.0);
      return s_pseudocount;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BondAngleAssignment::BondAngleAssignment() :
      m_BondAngleCounts(),
      m_BondAngleCountsStorage(), //1
      m_CentralAtomIndex( 0),
      m_AttachedAtomIndices(), //2
      m_AttachedAtomInverseIndices(), //3
      m_AttachedAtomInverseRawAtomIndices(), //4
      m_BondAngleFreeEnergies(), //5
      m_FreeEnergyUnseenRotamer(), //6
      m_ExpectedBondAngleFreeEnergies(), //7
      m_ConsiderAlternateChiralities(), //8
      m_NumberRingBonds(), //9
      m_CanFlip( false),
      m_HasAmide( false),
      m_TrustBondLengths( false)
    {
    }

    //! @brief Constructor from data members
    //! @param FRAGMENT fragment whose rotamer information this class stores
    //! @param FRAGMENT_CSI isomorphism object of fragment and molecule
    BondAngleAssignment::BondAngleAssignment
    (
      const RotamerLibraryInterface::t_BondAnglesWithCounts &BOND_ANGLE_COUNTS,
      const RotamerLibraryInterface::t_BondAngleMapKey &BOND_ANGLE_MAP_KEY,
      const ConformationInterface &MOLECULE,
      const size_t &CENTRAL_ATOM_INDICES,
      const bool &CONSIDER_ALTERNATE_CHIRALITIES,
      const ConformationGraphConverter::AtomComparisonType &ATOM_LEVEL
    ) :
      m_BondAngleCounts( BOND_ANGLE_COUNTS),
      m_CentralAtomIndex( CENTRAL_ATOM_INDICES),
      m_AttachedAtomIndices(),
      m_AttachedAtomInverseIndices( BOND_ANGLE_MAP_KEY.Third().GetSize()),
      m_AttachedAtomInverseRawAtomIndices( BOND_ANGLE_MAP_KEY.Third().GetSize()),
      m_BondAngleFreeEnergies(),
      m_ConsiderAlternateChiralities( CONSIDER_ALTERNATE_CHIRALITIES),
      m_TrustBondLengths( ATOM_LEVEL == ConformationGraphConverter::e_AtomType)
    {
      storage::List< storage::Pair< storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> > > bond_atom_types_to_index;
      auto itr_atom_itr( MOLECULE.GetAtomsIterator());
      itr_atom_itr.GotoPosition( CENTRAL_ATOM_INDICES);
      const storage::Vector< BondConformational> &bonds( itr_atom_itr->GetBonds());
      m_AttachedAtomIndices.Resize( bonds.GetSize());
      m_AttachedAtomIndices.SetAllElements( util::GetUndefinedSize_t());
      m_NumberRingBonds = size_t( 0);
      for( size_t bond( 0), nbonds( bonds.GetSize()); bond != nbonds; ++bond)
      {
        bond_atom_types_to_index.PushBack
        (
          storage::Pair< storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> >
          (
            storage::Pair< size_t, size_t>
            (
              bonds( bond).GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness),
              ConformationGraphConverter::ConvertAtomTypeData( bonds( bond).GetTargetAtom().GetAtomType(), ATOM_LEVEL)
            ),
            storage::Pair< size_t, size_t>
            (
              bond,
              MOLECULE.GetAtomIndex( bonds( bond).GetTargetAtom())
            )
          )
        );
        if( bonds( bond).GetBondType()->IsBondInRing())
        {
          ++m_NumberRingBonds;
        }
      }
      size_t bond_angle_index( 0);
      for
      (
        auto itr( BOND_ANGLE_MAP_KEY.Third().Begin()), itr_end( BOND_ANGLE_MAP_KEY.Third().End());
        itr != itr_end;
        ++itr, ++bond_angle_index
      )
      {
        for
        (
          auto itrlist( bond_atom_types_to_index.Begin()), itrlist_end( bond_atom_types_to_index.End());
          itrlist != itrlist_end;
          ++itrlist
        )
        {
          if( itrlist->First() == *itr)
          {
            // use m_AttachedAtomIndices to go from bond-of-atom index to row-of-matrix index
            m_AttachedAtomIndices( itrlist->Second().First()) = bond_angle_index;
            // use m_AttachedAtomInverseIndices to go from row-of-bond angle matrix to bond-of-atom index
            m_AttachedAtomInverseIndices( bond_angle_index) = itrlist->Second().First();
            m_AttachedAtomInverseRawAtomIndices( bond_angle_index) = itrlist->Second().Second();
            bond_atom_types_to_index.RemoveElement( itrlist);
            break;
          }
        }
      }
      m_CanFlip = false;
      m_HasAmide = itr_atom_itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAmide, size_t( 1));
      int mat_sign_desired( 0);
      if( m_AttachedAtomInverseIndices.GetSize() > size_t( 2) && !m_ConsiderAlternateChiralities)
      {
        util::SiPtrVector< const linal::Vector3D> coord_ptrs( m_AttachedAtomInverseIndices.GetSize());
        for( size_t inv_n( 0), total_n( m_AttachedAtomInverseIndices.GetSize()); inv_n < total_n; ++inv_n)
        {
          coord_ptrs( inv_n) = bonds( m_AttachedAtomInverseIndices( inv_n)).GetTargetAtom().GetPosition();
        }
        auto mat
        (
          GetStandardizedCoordinateMatrix
          (
            itr_atom_itr->GetPosition(),
            coord_ptrs,
            itr_atom_itr->GetAtomType()->GetNumberBonds(),
            m_NumberRingBonds
          )
        );
        int rrow( mat.GetNumberRows() > size_t( 2) ? 2 : 0);
        mat_sign_desired = mat( rrow, 2) > 0.0 ? 1 : mat( rrow, 2) < 0.0 ? -1 : 0;

        if
        (
          itr_atom_itr->GetElementType() == GetElementTypes().e_Nitrogen
          || itr_atom_itr->GetElementType() == GetElementTypes().e_Boron
          || itr_atom_itr->GetElementType() == GetElementTypes().e_Aluminum
        )
        {
          if( itr_atom_itr->GetAtomType()->GetNumberBonds() == size_t( 3))
          {
            if
            (
              itr_atom_itr->GetAtomType()->GetHybridOrbitalType() == GetHybridOrbitalTypes().e_SP3
              || ( mat_sign_desired && itr_atom_itr->GetAtomType()->GetNumberElectronsInBonds() == size_t( 3))
            )
            {
              m_CanFlip = true;
            }
          }
        }

        if( mat_sign_desired && !m_CanFlip)
        {
          for( auto itr( BOND_ANGLE_COUNTS.Begin()), itr_end( BOND_ANGLE_COUNTS.End()); itr != itr_end; ++itr)
          {
            if( itr->First()( rrow, 2) * mat_sign_desired > 0.0)
            {
//              BCL_Debug( itr->First());
              m_BondAngleCountsStorage.PushBack( *itr);
            }
          }
          if( m_BondAngleCountsStorage.IsEmpty())
          {
            // none of the correct chirality was found. Use the incorrect chirality angles, just swapped
            m_BondAngleCountsStorage = BOND_ANGLE_COUNTS;
            for( auto itr( m_BondAngleCountsStorage.Begin()), itr_end( m_BondAngleCountsStorage.End()); itr != itr_end; ++itr)
            {
              itr->First()( rrow, 2) *= -1.0;
              if( itr->First().GetNumberRows() == size_t( rrow + 2))
              {
                itr->First()( rrow + 1, 2) *= -1.0;
              }
            }
          }
          m_BondAngleCounts = util::ToSiPtr( m_BondAngleCountsStorage);
        }
      }

      m_BondAngleFreeEnergies = linal::Vector< float>( m_BondAngleCounts->GetSize(), float( 0.0));
      double total_counts( 0.0);
      bond_angle_index = 0;
      for( auto itr( m_BondAngleCounts->Begin()), itr_end( m_BondAngleCounts->End()); itr != itr_end; ++itr, ++bond_angle_index)
      {
        total_counts += itr->Second().GetWeight() + GetPseudocount();
        m_BondAngleFreeEnergies( bond_angle_index) = itr->Second().GetWeight() + GetPseudocount();
      }
      double ave_counts( total_counts / double( m_BondAngleFreeEnergies.GetSize()));
      for( auto itr( m_BondAngleFreeEnergies.Begin()), itr_end( m_BondAngleFreeEnergies.End()); itr != itr_end; ++itr)
      {
        double original_cnt( *itr);
        *itr = -std::log( original_cnt / ave_counts);
        m_ExpectedBondAngleFreeEnergies.AddWeightedObservation( *itr, original_cnt);
      }
      m_FreeEnergyUnseenRotamer = -std::log( GetPseudocount() / ave_counts);
      m_BondAngleFreeEnergies -= m_ExpectedBondAngleFreeEnergies.GetAverage();
      m_FreeEnergyUnseenRotamer -= m_ExpectedBondAngleFreeEnergies.GetAverage();
    }

    //! @brief Clone function
    //! @return pointer to new BondAngleAssignment
    BondAngleAssignment *BondAngleAssignment::Clone() const
    {
      return new BondAngleAssignment( *this);
    }

    //! @brief get a standardized coordinate matrix for a given atom
    //! @param ATOM the atom of interest
    linal::Matrix< double> BondAngleAssignment::GetStandardizedCoordinateMatrix
    (
      const linal::Vector3D &ATOM_POS,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const size_t &FULL_N_BONDS,
      const size_t &N_RING_BONDS
    )
    {
      linal::Matrix< double> coords( COORDINATES.GetSize() - N_RING_BONDS, size_t( 3), 0.0);
//      BCL_Debug( ATOM_POS);
//      BCL_Debug( COORDINATES);
//      BCL_Debug( FULL_N_BONDS);
      if( COORDINATES.IsEmpty() || COORDINATES.GetSize() == N_RING_BONDS)
      {
        return coords;
      }
      coords( 0, 0) = 1.0;
      if( COORDINATES.GetSize() == size_t( 1))
      {
        return coords;
      }
      else if( COORDINATES.GetSize() == size_t( 2))
      {
        double radians( linal::ProjAngle( ATOM_POS, *COORDINATES( 0), *COORDINATES( 1)));
        coords( 1, 0) = std::cos( radians);
        coords( 1, 1) = std::sin( radians);
        return coords;
      }

      storage::Vector< linal::Vector3D> unitvecs( COORDINATES.GetSize());
      for( size_t i( 0), nb( COORDINATES.GetSize()); i < nb; ++i)
      {
        unitvecs( i) = linal::UnitVector( ATOM_POS, *COORDINATES( i));
      }
      static const linal::Vector3D origin( 0.0);
      if( !N_RING_BONDS)
      {
        math::TransformationMatrix3D tf
        (
          coord::LineSegment3D( origin, linal::Vector3D( 1.0, 0.0, 0.0)),
          coord::LineSegment3D( origin, unitvecs( 0))
        );
        linal::Vector3D vo( unitvecs( 1));
        vo.Transform( tf);
        vo.Normalize();
        double radian( std::atan2( vo.Z(), vo.Y()));
        math::RotationMatrix3D rot( coord::GetAxes().e_X, radian);
        tf( rot);
        for( size_t i( 1), nr( coords.GetNumberRows()); i < nr; ++i)
        {
          vo = unitvecs( i);
          vo.Transform( tf);
          vo.Normalize();
          coords.GetRow( i).CopyValues( vo);
        }
        coords( 1, 2) = 0.0;
        if( math::Absolute( coords( 2, 2)) < 0.05 && FULL_N_BONDS == size_t( 3))
        {
          coords( 2, 2) = 0.0;
          auto row2( coords.GetRow( 2));
          row2.Normalize();
        }
      }
      else
      {
        storage::Vector< linal::Vector3D> identity( N_RING_BONDS + 1, linal::Vector3D( 0.0));
        util::SiPtrVector< const linal::Vector3D> tmp_ref( N_RING_BONDS + 1);
        for( size_t i( 0); i < N_RING_BONDS; ++i)
        {
          unitvecs( i) += ATOM_POS;
          identity( i)( i) = 1.0;
          tmp_ref( i) = unitvecs( i);
        }
        tmp_ref( N_RING_BONDS) = util::ToSiPtr( ATOM_POS);
        math::TransformationMatrix3D tf
        (
          quality::RMSD::SuperimposeCoordinates
          (
            util::SiPtrVector< const linal::Vector3D>( identity.Begin(), identity.End()),
            tmp_ref
          )
        );
//        BCL_Debug( tmp_ref);
//        BCL_Debug( tf);
//        BCL_Debug( ATOM_POS);
        linal::Vector3D vo;
        for( size_t i( N_RING_BONDS), nr( COORDINATES.GetSize()); i < nr; ++i)
        {
          vo = *COORDINATES( i);
          vo.Transform( tf);
          vo.Normalize();
          coords.GetRow( i - N_RING_BONDS).CopyValues( vo);
        }
        coords.ShrinkRows( COORDINATES.GetSize() - N_RING_BONDS);
        if( math::Absolute( coords( 0, 2)) < 0.05 && FULL_N_BONDS == size_t( 3))
        {
          coords( 0, 2) = 0.0;
          auto row2( coords.GetRow( 0));
          row2.Normalize();
        }
//        BCL_Debug( coords);
      }
      return coords;
    }

    //! @brief Transfer standardized coordinate matrix (usually taken from a different molecule) to another
    //! @param MOLECULE the molecule of interest
    //! @param STANDARDIZED_COORDS the coordinates of interest
    void BondAngleAssignment::TransferStandardizedCoordinateMatrixForRingSystem
    (
      const storage::Vector< sdf::AtomInfo> &MOLECULE,
      linal::Matrix< double> &STANDARDIZED_COORDS
    ) const
    {
      const sdf::AtomInfo &atom( MOLECULE( m_CentralAtomIndex));
      const linal::Vector3D &atom_pos( atom.GetCoordinates());
      util::SiPtrVector< const linal::Vector3D> coord_ptrs( m_NumberRingBonds);
      for( size_t inv_n( 0); inv_n < m_NumberRingBonds; ++inv_n)
      {
        coord_ptrs( inv_n) = MOLECULE( m_AttachedAtomInverseRawAtomIndices( inv_n)).GetCoordinates();
      }
      storage::Vector< linal::Vector3D> unitvecs( m_NumberRingBonds);
      for( size_t i( 0); i < m_NumberRingBonds; ++i)
      {
        unitvecs( i) = linal::UnitVector( atom_pos, *coord_ptrs( i));
      }
      storage::Vector< linal::Vector3D> identity( m_NumberRingBonds + 1, linal::Vector3D( 0.0));
      util::SiPtrVector< const linal::Vector3D> tmp_ref( m_NumberRingBonds + 1);
      for( size_t i( 0); i < m_NumberRingBonds; ++i)
      {
        unitvecs( i) += atom_pos;
        identity( i)( i) = 1.0;
        tmp_ref( i) = unitvecs( i);
      }
      tmp_ref( m_NumberRingBonds) = util::ToSiPtr( atom_pos);
      math::TransformationMatrix3D tf
      (
        quality::RMSD::SuperimposeCoordinates
        (
          tmp_ref,
          util::SiPtrVector< const linal::Vector3D>( identity.Begin(), identity.End())
        )
      );
      for( size_t i( 0), nr( STANDARDIZED_COORDS.GetNumberRows()); i < nr; ++i)
      {
        linal::Vector3D ir( STANDARDIZED_COORDS.GetRow( i));
        ir.Transform( tf);
        auto row( STANDARDIZED_COORDS.GetRow( i));
        row.CopyValues( ir);
      }
    }

    //! @brief get the free energy of a particular rotamer
    //! @param ROTAMER rotamer signature of interest
    //! @param CONSIDER_ISOMETRY_CHANGES true - consider changes in isometry
    double BondAngleAssignment::GetFreeEnergy( const FragmentComplete &FRAGMENT) const
    {
      // Key <- atom-type <> Vector< bond type, atom type>
      // Value <- Matrix with unit-vector coordinates of all atoms after the first.
      // The first atom is always moved to 1 0 0, second atom is moved such that it is 0 in
      const auto &atom( FRAGMENT.GetAtomVector()( m_CentralAtomIndex));
      const auto &bonds( atom.GetBonds());

      const size_t heavy_bonds( m_AttachedAtomInverseIndices.GetSize());
      util::SiPtrVector< const linal::Vector3D> coord_ptrs( heavy_bonds);
      for( size_t i( 0); i < heavy_bonds; ++i)
      {
        coord_ptrs( i) = bonds( m_AttachedAtomInverseIndices( i)).GetTargetAtom().GetPosition();
      }
      linal::Matrix< double> coords
      (
        GetStandardizedCoordinateMatrix( atom.GetPosition(), coord_ptrs, atom.GetAtomType()->GetNumberBonds(), m_NumberRingBonds)
      );

      //  static const double s_tolerance( 0.025);
      double closest( 100000.0);
      auto itr_closest( m_BondAngleCounts->Begin());
      size_t index( 0), best_index( 0);
      size_t fpos( m_NumberRingBonds ? 0 : 3);
      // double tolerance( 1000); // s_tolerance * math::Sqrt( coords.GetNumberRows() - ( m_NumberRingBonds ? 0.0 : 1.0)));
      for( auto itr_list( m_BondAngleCounts->Begin()), itr_list_end( m_BondAngleCounts->End()); itr_list != itr_list_end; ++itr_list, ++index)
      {
        double cosine( linal::Distance( itr_list->First().AsVector().Slice( fpos), coords.AsVector().Slice( fpos)));
        if( cosine < closest)
        {
          closest = cosine;
          itr_closest = itr_list;
          best_index = index;
        }
        if( m_NumberRingBonds && itr_list->First()( 0, 2) * coords( 0, 2) < 0.0)
        {
          coords( 0, 2) *= -1.0;
          cosine = linal::Distance( itr_list->First().AsVector().Slice( fpos), coords.AsVector().Slice( fpos));
          coords( 0, 2) *= -1.0;
          if( cosine < closest)
          {
            closest = cosine;
            itr_closest = itr_list;
            best_index = index;
          }
        }
      }

      return m_BondAngleFreeEnergies( best_index);
    }

    //! @brief get a standardized coordinate matrix for a given atom
    //! @param ATOM the atom of interest
    linal::Matrix< double> BondAngleAssignment::GetStandardizedCoordinateMatrix
    (
      const ConformationInterface &MOLECULE
    ) const
    {
      auto itr_atom_itr( MOLECULE.GetAtomsIterator());
      itr_atom_itr.GotoPosition( m_CentralAtomIndex);
      const storage::Vector< BondConformational> &bonds( itr_atom_itr->GetBonds());
      util::SiPtrVector< const linal::Vector3D> coord_ptrs( m_AttachedAtomInverseIndices.GetSize());
      for( size_t inv_n( 0), total_n( m_AttachedAtomInverseIndices.GetSize()); inv_n < total_n; ++inv_n)
      {
        coord_ptrs( inv_n) = bonds( m_AttachedAtomInverseIndices( inv_n)).GetTargetAtom().GetPosition();
      }
      return
        GetStandardizedCoordinateMatrix
        (
          itr_atom_itr->GetPosition(),
          coord_ptrs,
          itr_atom_itr->GetAtomType()->GetNumberBonds(),
          m_NumberRingBonds
        );
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &BondAngleAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BondAngleAssignment::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &BondAngleAssignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl

