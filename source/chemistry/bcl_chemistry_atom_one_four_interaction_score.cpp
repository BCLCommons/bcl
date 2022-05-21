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
#include "chemistry/bcl_chemistry_atom_one_four_interaction_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_statistics.h"

//#define BCL_PROFILE_AtomOneFourInteractionScore
#ifdef BCL_PROFILE_AtomOneFourInteractionScore
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
    const util::SiPtr< const util::ObjectInterface> AtomOneFourInteractionScore::s_Instance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< FragmentComplete, double> >::AddInstance( new AtomOneFourInteractionScore( 1.0, 1.0))
    );

    storage::Map< std::string, storage::VectorND< 3, double> > AtomOneFourInteractionScore::s_InteractionMap =
        storage::Map< std::string, storage::VectorND< 3, double> >();
    namespace
    {
      // initialization for bonds
      const bool map_init( AtomOneFourInteractionScore::InitializeInteractionMap());
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor. Can take whether to consider hydrogens or tolerance
    AtomOneFourInteractionScore::AtomOneFourInteractionScore( double LENGTH_WEIGHT, double INTERACTION_WEIGHT) :
      m_LengthWeight( LENGTH_WEIGHT),
      m_InteractionIdentityWeight( INTERACTION_WEIGHT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AtomOneFourInteractionScore
    AtomOneFourInteractionScore *AtomOneFourInteractionScore::Clone() const
    {
      return new AtomOneFourInteractionScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomOneFourInteractionScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AtomOneFourInteractionScore::GetAlias() const
    {
      static const std::string s_name( "AtomOneFour");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param MOLECULE molecule that needs to scored
    //! @return clash score for the given molecule
    double AtomOneFourInteractionScore::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      return operator()( static_cast< const ConformationInterface &>( MOLECULE));
    }

    //! @brief evaluate clashes for given atoms
    //! @param MOLECULE molecule that needs to scored
    //! @return clash score for the given atoms
    double AtomOneFourInteractionScore::operator()
    (
      const ConformationInterface &MOLECULE
    ) const
    {
      #ifdef BCL_PROFILE_AtomOneFourInteractionScore
      static util::Stopwatch s_timer( "1-4 interaction score", util::Message::e_Standard, true);
      s_timer.Start();
      #endif
      util::ObjectDataLabel label( this->GetCompleteSerializer().GetLabel());
      auto cache_val( MOLECULE.FindInCache( label));
      if( cache_val.IsDefined())
      {
        #ifdef BCL_PROFILE_AtomOneFourInteractionScore
        s_timer.Stop();
        #endif
        return cache_val->operator ()( 0);
      }
      if( !IsMoleculeInfoCached( MOLECULE))
      {
        UpdateMolecule( MOLECULE);
      }
      m_AtomOneFourInteractionScores.SetAllElements( 0.0);

      auto positions( MOLECULE.GetAtomCoordinates());
      size_t n_valid( 0);
      for( size_t i( 0), na( positions.GetSize()); i < na; ++i)
      {
        double closest_distance( 4.6);
        auto itr_closest( m_OneFourNeighbors( i).End());
        double max_contact( 0.0);
        for( auto itr( m_OneFourNeighbors( i).Begin()), itr_end( m_OneFourNeighbors( i).End()); itr != itr_end; ++itr)
        {
          double dist( linal::Distance( *positions( itr->First()), *positions( i)));
          if( dist < closest_distance)
          {
            closest_distance = dist;
            itr_closest = itr;
          }
          max_contact = std::max( max_contact, itr->Second()->First());
        }
        if( itr_closest != m_OneFourNeighbors( i).End())
        {
          ++n_valid;
          m_AtomOneFourInteractionScores( i) = -m_InteractionIdentityWeight * itr_closest->Second()->First() / std::max( max_contact, 0.01);
          if( closest_distance < itr_closest->Second()->Second())
          {
            double length_zscore( ( itr_closest->Second()->Second() - closest_distance) / itr_closest->Second()->Third());
            if( length_zscore > 2.0)
            {
              m_AtomOneFourInteractionScores( i) += m_LengthWeight * ( length_zscore - 2.0);
            }
          }
        }
      }

      double clash_sum( math::Statistics::Sum( m_AtomOneFourInteractionScores.Begin(), m_AtomOneFourInteractionScores.End(), 0.0));
      clash_sum /= std::max( double( n_valid), 1.0);
      #ifdef BCL_PROFILE_AtomOneFourInteractionScore
      s_timer.Stop();
      #endif
      MOLECULE.Cache( label, linal::Vector< float>( size_t( 1), clash_sum));
      return clash_sum;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomOneFourInteractionScore::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Propensity for 1-4 interactions and the lengths between them"
      );
      serializer.AddInitializer
      (
        "length weight",
        "weight assigned to the length z-score term",
        io::Serialization::GetAgent( &m_LengthWeight),
        "1.0"
      );
      serializer.AddInitializer
      (
        "interaction weight",
        "weight for the interaction term",
        io::Serialization::GetAgent( &m_InteractionIdentityWeight),
        "1.0"
      );
      return serializer;
    }

    //! @brief Test whether this molecule is the same (constitutionally) as the molecule for which the state of this
    //!        class currently can handle
    bool AtomOneFourInteractionScore::IsMoleculeInfoCached( const ConformationInterface &CONF) const
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

    namespace
    {
      //! @brief test two (assumed sorted) vectors for any overlapping indices
      bool HaveOverlap( const storage::Vector< size_t> &COMPONENT_A, const storage::Vector< size_t> &COMPONENT_B)
      {
        storage::Vector< size_t>::const_iterator itr_b( COMPONENT_B.Begin()), itr_b_end( COMPONENT_B.End());
        for
        (
          storage::Vector< size_t>::const_iterator itr_a( COMPONENT_A.Begin()), itr_a_end( COMPONENT_A.End());
          itr_a != itr_a_end && itr_b != itr_b_end;
          ++itr_a
        )
        {
          while( itr_b != itr_b_end)
          {
            if( *itr_a > *itr_b)
            {
              ++itr_b;
            }
            else if( *itr_a == *itr_b)
            {
              return true;
            }
            else
            {
              break;
            }
          }
        }
        return false;
      }
    }

    //! @brief Update molecule change the molecule that this class will compute the clash score for
    //! @param MOL molecule of interest
    void AtomOneFourInteractionScore::UpdateMolecule( const ConformationInterface &CONF) const
    {
      m_LastBondInfo = CONF.GetBondInfo();

      const size_t n_atoms( CONF.GetNumberAtoms());
      m_LastAtomTypes.Resize( n_atoms);
      m_AtomOneFourInteractionScores.Resize( n_atoms);
      m_OneFourNeighbors.Resize( n_atoms);
      size_t atom_index( 0);
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( CONF.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++atom_index
      )
      {
        m_AtomOneFourInteractionScores( atom_index) = 0.0;
        m_LastAtomTypes( atom_index) = itr_atoms->GetAtomType();
        m_OneFourNeighbors( atom_index).Reset();
      }

      ConformationGraphConverter::t_AtomGraph molgraph
      (
        ConformationGraphConverter::CreateGraphWithAtoms( CONF)
      );
      FragmentSplitRigid splitter( size_t( 1), false);
      storage::List< storage::Vector< size_t> > components( splitter.GetComponentVertices( CONF, molgraph));
      size_t rigid_id( 0);
      storage::Vector< storage::Vector< size_t> > rigid_components( CONF.GetSize());
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator
          itr_rigid( components.Begin()), itr_rigid_end( components.End());
        itr_rigid != itr_rigid_end;
        ++itr_rigid, ++rigid_id
      )
      {
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_rigid_atom( itr_rigid->Begin()), itr_rigid_atom_end( itr_rigid->End());
          itr_rigid_atom != itr_rigid_atom_end;
          ++itr_rigid_atom
        )
        {
          rigid_components( *itr_rigid_atom).PushBack( rigid_id);
        }
      }

      atom_index = 0;
      auto atom_types( CONF.GetAtomTypesVector());
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( CONF.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++atom_index
      )
      {
        auto dists( graph::Connectivity::DistancesToOtherVertices( molgraph, atom_index));
        for( size_t atom_index_b( 0); atom_index_b < n_atoms; ++atom_index_b)
        {
          if( ( *dists)( atom_index_b) != size_t( 3))
          {
            continue;
          }
          if( HaveOverlap( rigid_components( atom_index), rigid_components( atom_index_b)))
          {
            continue;
          }
          const std::string st
          (
            atom_types( atom_index)->GetTwoLetterCode() < atom_types( atom_index_b)->GetTwoLetterCode()
            ? atom_types( atom_index)->GetTwoLetterCode() +
              atom_types( atom_index_b)->GetTwoLetterCode()
            : atom_types( atom_index_b)->GetTwoLetterCode() +
              atom_types( atom_index)->GetTwoLetterCode()
          );
          auto look( s_InteractionMap.Find( st));

          if( look != s_InteractionMap.End())
          {
            m_OneFourNeighbors( atom_index).PushBack
            (
              storage::Pair< size_t, util::SiPtr< const storage::VectorND< 3, double> > >
              (
                atom_index_b,
                util::SiPtr< const storage::VectorND< 3, double> >( look->second)
              )
            );
            m_OneFourNeighbors( atom_index_b).PushBack
            (
              storage::Pair< size_t, util::SiPtr< const storage::VectorND< 3, double> > >
              (
                atom_index,
                util::SiPtr< const storage::VectorND< 3, double> >( look->second)
              )
            );
          }
          else
          {
            BCL_MessageStd( "Could not find interaction term for: " + st);
          }
        }
      }
    }

    namespace
    {
      double EstimateVdw( const AtomConformationalInterface &A)
      {
        double rad_sum( 0.0);
        if( A.GetAtomType() == GetAtomTypes().H_S && A.GetBonds().GetSize() == size_t( 1))
        {
          const ElementType &bonded_element
          (
            A.GetBonds().FirstElement().GetTargetAtom().GetElementType()
          );
          // New J. Chem. 2005, 29, 371-377, taken as average - std. to give approximate lower bond
          if( bonded_element == GetElementTypes().e_Oxygen)
          {
            rad_sum += 0.45;
          }
          else if( bonded_element == GetElementTypes().e_Nitrogen)
          {
            rad_sum += 0.60;
          }
          else // if( bonded_element == GetElementTypes().e_Carbon)...or any other generic type.
          {
            rad_sum += 0.75;
          }
        }
        else if( A.GetElementType() == GetElementTypes().e_Oxygen)
        {
          rad_sum = 1.3;
        }
        else if( A.GetElementType() == GetElementTypes().e_Fluorine)
        {
          rad_sum = 1.3;
        }
        else if( A.GetElementType() == GetElementTypes().e_Carbon)
        {
          rad_sum = 1.6;
        }
        else if( A.GetElementType() == GetElementTypes().e_Nitrogen && A.GetBonds().GetSize() <= size_t( 3))
        {
          rad_sum = 1.45;
        }
        else
        {
          rad_sum = A.GetElementType()->GetProperty( ElementTypeData::e_VDWaalsRadius);
        }
        return rad_sum;
      }
    }

    //! @brief estimate the minimum distance of two atoms of a given type, bond distance, and presence in the same rigid structure
    double AtomOneFourInteractionScore::EstimateMinimumDistance
    (
      const AtomConformationalInterface &A,
      const AtomConformationalInterface &B,
      const int &BOND_DISTANCE,
      const storage::Vector< size_t> &RIGID_SUBSTRUCTURES_A,
      const storage::Vector< size_t> &RIGID_SUBSTRUCTURES_B
    )
    {
      // ignore bonded atoms
      if( BOND_DISTANCE <= 3)
      {
        return 0.0;
      }
//      double multiplier( 0.81);
//      if( HaveOverlap( RIGID_SUBSTRUCTURES_A, RIGID_SUBSTRUCTURES_B))
//      {
//        multiplier = 0.7;
//      }
      return ( BondLengths::GetAverageCovalentRadius( A) + BondLengths::GetAverageCovalentRadius( B)) * 1.01;
      //return ( EstimateVdw( A) + EstimateVdw( B)) * multiplier;
    }

    //! @brief initialize the interaction map
    //! this is slow and should be performed only once
    bool AtomOneFourInteractionScore::InitializeInteractionMap()
    {
      s_InteractionMap.Reset();

      s_InteractionMap[ "C4N4"] = storage::VectorND< 3, double>( 0.070817, 2.8103,  0.253281);
      s_InteractionMap[ "BrBr"] = storage::VectorND< 3, double>( 0.105263, 3.41984, 0.2);
      s_InteractionMap[ "B3SX"] = storage::VectorND< 3, double>( 0.111111, 4.50248, 0.2);
      s_InteractionMap[ "IXIX"] = storage::VectorND< 3, double>( 0.125,    3.9617,  0.2);
      s_InteractionMap[ "N4S4"] = storage::VectorND< 3, double>( 0.125,    3.17113, 0.2);
      s_InteractionMap[ "N4N4"] = storage::VectorND< 3, double>( 0.131783, 2.7647,  0.2);
      s_InteractionMap[ "SeSe"] = storage::VectorND< 3, double>( 0.133333, 3.4363,  0.2);
      s_InteractionMap[ "N3N4"] = storage::VectorND< 3, double>( 0.135065, 2.57337, 0.250951);
      s_InteractionMap[ "SiSi"] = storage::VectorND< 3, double>( 0.155405, 4.05472, 0.2);
      s_InteractionMap[ "N4SX"] = storage::VectorND< 3, double>( 0.166667, 2.92639, 0.2);
      s_InteractionMap[ "PXPX"] = storage::VectorND< 3, double>( 0.175182, 3.21358, 0.421511);
      s_InteractionMap[ "C2S4"] = storage::VectorND< 3, double>( 0.1875,   3.3622,  0.2);
      s_InteractionMap[ "B4N4"] = storage::VectorND< 3, double>( 0.217391, 2.94818, 0.2);
      s_InteractionMap[ "LiOX"] = storage::VectorND< 3, double>( 0.25,     3.51807, 0.2);
      s_InteractionMap[ "OXS6"] = storage::VectorND< 3, double>( 0.266667, 3.31061, 0.2);
      s_InteractionMap[ "B3NX"] = storage::VectorND< 3, double>( 0.285714, 2.82543, 0.2);
      s_InteractionMap[ "B4SX"] = storage::VectorND< 3, double>( 0.3,      3.52883, 0.2);
      s_InteractionMap[ "BrSi"] = storage::VectorND< 3, double>( 0.3,      3.38557, 0.2);
      s_InteractionMap[ "C4S4"] = storage::VectorND< 3, double>( 0.304628, 2.97592, 0.260738);
      s_InteractionMap[ "B4PX"] = storage::VectorND< 3, double>( 0.327869, 3.28697, 0.355156);
      s_InteractionMap[ "AlSi"] = storage::VectorND< 3, double>( 0.333333, 3.58913, 0.2);
      s_InteractionMap[ "B4FX"] = storage::VectorND< 3, double>( 0.333333, 2.96174, 0.2);
      s_InteractionMap[ "C3Na"] = storage::VectorND< 3, double>( 0.333333, 4.05441, 0.2);
      s_InteractionMap[ "S3S3"] = storage::VectorND< 3, double>( 0.333333, 4.18512, 0.2);
      s_InteractionMap[ "S4Si"] = storage::VectorND< 3, double>( 0.333333, 3.47832, 0.2);
      s_InteractionMap[ "C4IX"] = storage::VectorND< 3, double>( 0.33871,  3.16041, 0.2);
      s_InteractionMap[ "BrC4"] = storage::VectorND< 3, double>( 0.353511, 3.01106, 0.213944);
      s_InteractionMap[ "C4S3"] = storage::VectorND< 3, double>( 0.357394, 2.85919, 0.284947);
      s_InteractionMap[ "SXSX"] = storage::VectorND< 3, double>( 0.361596, 2.88625, 0.310914);
      s_InteractionMap[ "KXOX"] = storage::VectorND< 3, double>( 0.363636, 3.95794, 0.2);
      s_InteractionMap[ "NXNX"] = storage::VectorND< 3, double>( 0.366667, 2.45309, 0.285642);
      s_InteractionMap[ "C3KX"] = storage::VectorND< 3, double>( 0.37037,  3.82281, 0.2);
      s_InteractionMap[ "AlC3"] = storage::VectorND< 3, double>( 0.371429, 3.23738, 0.405485);
      s_InteractionMap[ "BrSX"] = storage::VectorND< 3, double>( 0.375,    3.44836, 0.2);
      s_InteractionMap[ "C4Li"] = storage::VectorND< 3, double>( 0.375,    2.97251, 0.2);
      s_InteractionMap[ "FXS3"] = storage::VectorND< 3, double>( 0.375,    2.95152, 0.2);
      s_InteractionMap[ "N3S3"] = storage::VectorND< 3, double>( 0.387097, 2.93079, 0.2);
      s_InteractionMap[ "C2SX"] = storage::VectorND< 3, double>( 0.392308, 2.93766, 0.299453);
      s_InteractionMap[ "ClCl"] = storage::VectorND< 3, double>( 0.394286, 2.73608, 0.261716);
      s_InteractionMap[ "S4SX"] = storage::VectorND< 3, double>( 0.396226, 3.23484, 0.21667);
      s_InteractionMap[ "BrOX"] = storage::VectorND< 3, double>( 0.398374, 2.91729, 0.2);
      s_InteractionMap[ "B3PX"] = storage::VectorND< 3, double>( 0.4,      3.34562, 0.2);
      s_InteractionMap[ "ClS3"] = storage::VectorND< 3, double>( 0.4,      3.56241, 0.2);
      s_InteractionMap[ "N3Te"] = storage::VectorND< 3, double>( 0.4,      3.82996, 0.2);
      s_InteractionMap[ "C4N3"] = storage::VectorND< 3, double>( 0.402404, 2.52422, 0.289313);
      s_InteractionMap[ "C4C4"] = storage::VectorND< 3, double>( 0.408292, 2.12405, 0.281336);
      s_InteractionMap[ "C4Cl"] = storage::VectorND< 3, double>( 0.415479, 2.87436, 0.243617);
      s_InteractionMap[ "IXN3"] = storage::VectorND< 3, double>( 0.416667, 3.24942, 0.248259);
      s_InteractionMap[ "C4KX"] = storage::VectorND< 3, double>( 0.417722, 3.13618, 0.307729);
      s_InteractionMap[ "B3Si"] = storage::VectorND< 3, double>( 0.41791,  3.09127, 0.218312);
      s_InteractionMap[ "C3N4"] = storage::VectorND< 3, double>( 0.42011,  2.68583, 0.2);
      s_InteractionMap[ "C3S3"] = storage::VectorND< 3, double>( 0.428356, 2.77392, 0.219089);
      s_InteractionMap[ "C4Te"] = storage::VectorND< 3, double>( 0.428571, 3.21121, 0.330505);
      s_InteractionMap[ "HXNa"] = storage::VectorND< 3, double>( 0.428571, 3.17792, 0.2);
      s_InteractionMap[ "C4Na"] = storage::VectorND< 3, double>( 0.433333, 3.31451, 0.2);
      s_InteractionMap[ "C3Li"] = storage::VectorND< 3, double>( 0.434783, 2.80469, 0.2);
      s_InteractionMap[ "HXLi"] = storage::VectorND< 3, double>( 0.434783, 2.70635, 0.2);
      s_InteractionMap[ "AlC4"] = storage::VectorND< 3, double>( 0.439791, 3.22675, 0.215065);
      s_InteractionMap[ "BrC3"] = storage::VectorND< 3, double>( 0.442142, 3.0005,  0.2);
      s_InteractionMap[ "C2N4"] = storage::VectorND< 3, double>( 0.444444, 3.42478, 0.2);
      s_InteractionMap[ "BrN3"] = storage::VectorND< 3, double>( 0.450262, 3.01177, 0.21982);
      s_InteractionMap[ "SeSi"] = storage::VectorND< 3, double>( 0.454545, 3.64919, 0.2);
      s_InteractionMap[ "ClN3"] = storage::VectorND< 3, double>( 0.455652, 2.90909, 0.309346);
      s_InteractionMap[ "N3S4"] = storage::VectorND< 3, double>( 0.457627, 2.94952, 0.285873);
      s_InteractionMap[ "N3Si"] = storage::VectorND< 3, double>( 0.460199, 3.02186, 0.369227);
      s_InteractionMap[ "PXSi"] = storage::VectorND< 3, double>( 0.464789, 3.21903, 0.256917);
      s_InteractionMap[ "N3PX"] = storage::VectorND< 3, double>( 0.470982, 2.92808, 0.377367);
      s_InteractionMap[ "ClSi"] = storage::VectorND< 3, double>( 0.471429, 3.37706, 0.231887);
      s_InteractionMap[ "B4C4"] = storage::VectorND< 3, double>( 0.471487, 2.85739, 0.282397);
      s_InteractionMap[ "N3N3"] = storage::VectorND< 3, double>( 0.472764, 2.57282, 0.239632);
      s_InteractionMap[ "FXS4"] = storage::VectorND< 3, double>( 0.473684, 2.88828, 0.212907);
      s_InteractionMap[ "C3IX"] = storage::VectorND< 3, double>( 0.483221, 3.30496, 0.2);
      s_InteractionMap[ "PXSe"] = storage::VectorND< 3, double>( 0.483871, 3.42283, 0.2);
      s_InteractionMap[ "HXSX"] = storage::VectorND< 3, double>( 0.48866,  2.2436,  0.2);
      s_InteractionMap[ "HXSe"] = storage::VectorND< 3, double>( 0.492105, 2.19564, 0.230706);
      s_InteractionMap[ "BrNX"] = storage::VectorND< 3, double>( 0.492537, 2.98265, 0.2);
      s_InteractionMap[ "B3N4"] = storage::VectorND< 3, double>( 0.5,      3.59316, 0.2);
      s_InteractionMap[ "B3Na"] = storage::VectorND< 3, double>( 0.5,      3.71236, 0.2);
      s_InteractionMap[ "BrS3"] = storage::VectorND< 3, double>( 0.5,      3.66648, 0.2);
      s_InteractionMap[ "C2Cl"] = storage::VectorND< 3, double>( 0.5,      3.00224, 0.2);
      s_InteractionMap[ "C2PX"] = storage::VectorND< 3, double>( 0.5,      3.14231, 0.242434);
      s_InteractionMap[ "ClPX"] = storage::VectorND< 3, double>( 0.5,      3.26375, 0.2);
      s_InteractionMap[ "FXN4"] = storage::VectorND< 3, double>( 0.5,      2.68265, 0.2);
      s_InteractionMap[ "KXN3"] = storage::VectorND< 3, double>( 0.5,      4.00038, 0.2);
      s_InteractionMap[ "NaOX"] = storage::VectorND< 3, double>( 0.5,      3.71812, 0.2);
      s_InteractionMap[ "SXSi"] = storage::VectorND< 3, double>( 0.507463, 3.49989, 0.226763);
      s_InteractionMap[ "C4Si"] = storage::VectorND< 3, double>( 0.507797, 3.01137, 0.314818);
      s_InteractionMap[ "HXKX"] = storage::VectorND< 3, double>( 0.521739, 2.62097, 0.2);
      s_InteractionMap[ "HXN4"] = storage::VectorND< 3, double>( 0.52215,  1.98364, 0.2);
      s_InteractionMap[ "B4Si"] = storage::VectorND< 3, double>( 0.522727, 3.09127, 0.2);
      s_InteractionMap[ "B3C4"] = storage::VectorND< 3, double>( 0.524942, 2.67438, 0.247153);
      s_InteractionMap[ "OXPX"] = storage::VectorND< 3, double>( 0.532151, 2.70601, 0.286731);
      s_InteractionMap[ "ClNX"] = storage::VectorND< 3, double>( 0.536082, 2.79831, 0.226482);
      s_InteractionMap[ "B3N3"] = storage::VectorND< 3, double>( 0.537736, 2.83817, 0.26318);
      s_InteractionMap[ "C3Si"] = storage::VectorND< 3, double>( 0.5384,   3.03272, 0.287411);
      s_InteractionMap[ "PXSX"] = storage::VectorND< 3, double>( 0.545455, 3.11658, 0.343111);
      s_InteractionMap[ "N3SX"] = storage::VectorND< 3, double>( 0.548428, 2.61131, 0.330903);
      s_InteractionMap[ "C3Te"] = storage::VectorND< 3, double>( 0.549223, 2.80085, 0.398038);
      s_InteractionMap[ "C3S4"] = storage::VectorND< 3, double>( 0.549472, 3.00336, 0.2);
      s_InteractionMap[ "C3Cl"] = storage::VectorND< 3, double>( 0.553514, 2.83285, 0.250613);
      s_InteractionMap[ "NXOX"] = storage::VectorND< 3, double>( 0.553781, 2.44625, 0.2);
      s_InteractionMap[ "C3N3"] = storage::VectorND< 3, double>( 0.555796, 2.44181, 0.2);
      s_InteractionMap[ "C4PX"] = storage::VectorND< 3, double>( 0.558119, 2.62054, 0.30986);
      s_InteractionMap[ "BrFX"] = storage::VectorND< 3, double>( 0.565217, 3.03005, 0.2);
      s_InteractionMap[ "IXOX"] = storage::VectorND< 3, double>( 0.568627, 2.95931, 0.2);
      s_InteractionMap[ "AlN3"] = storage::VectorND< 3, double>( 0.571429, 3.52343, 0.2);
      s_InteractionMap[ "B3Cl"] = storage::VectorND< 3, double>( 0.571429, 3.59502, 0.2);
      s_InteractionMap[ "ClS4"] = storage::VectorND< 3, double>( 0.571429, 3.45351, 0.2);
      s_InteractionMap[ "FXS6"] = storage::VectorND< 3, double>( 0.571429, 3.08477, 0.2);
      s_InteractionMap[ "IXSi"] = storage::VectorND< 3, double>( 0.571429, 4.01686, 0.2);
      s_InteractionMap[ "KXNX"] = storage::VectorND< 3, double>( 0.571429, 3.5161,  0.2);
      s_InteractionMap[ "N3Se"] = storage::VectorND< 3, double>( 0.572368, 2.88649, 0.440012);
      s_InteractionMap[ "C4FX"] = storage::VectorND< 3, double>( 0.573686, 2.51087, 0.2);
      s_InteractionMap[ "C4Se"] = storage::VectorND< 3, double>( 0.573957, 2.93906, 0.255662);
      s_InteractionMap[ "HXIX"] = storage::VectorND< 3, double>( 0.573991, 2.80084, 0.2);
      s_InteractionMap[ "C3PX"] = storage::VectorND< 3, double>( 0.576232, 2.91493, 0.242893);
      s_InteractionMap[ "SXSe"] = storage::VectorND< 3, double>( 0.580645, 2.73185, 0.2);
      s_InteractionMap[ "HXTe"] = storage::VectorND< 3, double>( 0.581818, 3.13751, 0.2);
      s_InteractionMap[ "B3C3"] = storage::VectorND< 3, double>( 0.583475, 2.73364, 0.206868);
      s_InteractionMap[ "FXIX"] = storage::VectorND< 3, double>( 0.586207, 3.06177, 0.2);
      s_InteractionMap[ "B4C3"] = storage::VectorND< 3, double>( 0.587531, 2.91816, 0.225399);
      s_InteractionMap[ "NXS3"] = storage::VectorND< 3, double>( 0.588235, 2.84626, 0.2);
      s_InteractionMap[ "C3Se"] = storage::VectorND< 3, double>( 0.593527, 2.9461,  0.251018);
      s_InteractionMap[ "B2C4"] = storage::VectorND< 3, double>( 0.6,      3.33946, 0.2);
      s_InteractionMap[ "NXNa"] = storage::VectorND< 3, double>( 0.6,      3.58826, 0.2);
      s_InteractionMap[ "S4S4"] = storage::VectorND< 3, double>( 0.6,      3.68001, 0.2);
      s_InteractionMap[ "C4OX"] = storage::VectorND< 3, double>( 0.601789, 1.90484, 0.2);
      s_InteractionMap[ "ClHX"] = storage::VectorND< 3, double>( 0.603849, 2.1404,  0.2);
      s_InteractionMap[ "HXOX"] = storage::VectorND< 3, double>( 0.605121, 1.42381, 0.2);
      s_InteractionMap[ "OXOX"] = storage::VectorND< 3, double>( 0.609859, 1.8969,  0.2);
      s_InteractionMap[ "C2Si"] = storage::VectorND< 3, double>( 0.61157,  3.07759, 0.385272);
      s_InteractionMap[ "C4SX"] = storage::VectorND< 3, double>( 0.613546, 2.81937, 0.324675);
      s_InteractionMap[ "C3C4"] = storage::VectorND< 3, double>( 0.614118, 2.19933, 0.277876);
      s_InteractionMap[ "C2S3"] = storage::VectorND< 3, double>( 0.615385, 2.97958, 0.2);
      s_InteractionMap[ "FXPX"] = storage::VectorND< 3, double>( 0.615385, 2.72917, 0.2);
      s_InteractionMap[ "C3C3"] = storage::VectorND< 3, double>( 0.617886, 2.64584, 0.2);
      s_InteractionMap[ "NXSi"] = storage::VectorND< 3, double>( 0.62406,  2.80563, 0.300807);
      s_InteractionMap[ "ClFX"] = storage::VectorND< 3, double>( 0.625,    2.80433, 0.2);
      s_InteractionMap[ "BrHX"] = storage::VectorND< 3, double>( 0.626109, 2.33899, 0.2);
      s_InteractionMap[ "OXS4"] = storage::VectorND< 3, double>( 0.627488, 2.82684, 0.2);
      s_InteractionMap[ "NXPX"] = storage::VectorND< 3, double>( 0.627907, 2.70141, 0.277573);
      s_InteractionMap[ "C2N3"] = storage::VectorND< 3, double>( 0.631884, 2.65011, 0.315205);
      s_InteractionMap[ "OXSe"] = storage::VectorND< 3, double>( 0.633803, 2.87761, 0.244597);
      s_InteractionMap[ "ClOX"] = storage::VectorND< 3, double>( 0.634682, 2.78836, 0.2);
      s_InteractionMap[ "FXOX"] = storage::VectorND< 3, double>( 0.635993, 2.46714, 0.2);
      s_InteractionMap[ "C4HX"] = storage::VectorND< 3, double>( 0.637933, 1.33055, 0.2);
      s_InteractionMap[ "NXS4"] = storage::VectorND< 3, double>( 0.640177, 2.84955, 0.2);
      s_InteractionMap[ "OXS3"] = storage::VectorND< 3, double>( 0.643678, 2.70015, 0.2);
      s_InteractionMap[ "C3OX"] = storage::VectorND< 3, double>( 0.645234, 2.43584, 0.2);
      s_InteractionMap[ "B3OX"] = storage::VectorND< 3, double>( 0.645367, 2.62685, 0.267986);
      s_InteractionMap[ "HXN3"] = storage::VectorND< 3, double>( 0.646248, 1.38664, 0.2);
      s_InteractionMap[ "C3FX"] = storage::VectorND< 3, double>( 0.64778,  2.58316, 0.2);
      s_InteractionMap[ "C3SX"] = storage::VectorND< 3, double>( 0.650968, 2.8139,  0.245687);
      s_InteractionMap[ "B3HX"] = storage::VectorND< 3, double>( 0.653369, 2.14578, 0.26819);
      s_InteractionMap[ "OXTe"] = storage::VectorND< 3, double>( 0.653846, 3.02716, 0.2);
      s_InteractionMap[ "FXNX"] = storage::VectorND< 3, double>( 0.654237, 2.5278,  0.2);
      s_InteractionMap[ "N4OX"] = storage::VectorND< 3, double>( 0.655917, 2.44119, 0.2);
      s_InteractionMap[ "B4N3"] = storage::VectorND< 3, double>( 0.656716, 2.84282, 0.281918);
      s_InteractionMap[ "C2HX"] = storage::VectorND< 3, double>( 0.65675,  2.06082, 0.252472);
      s_InteractionMap[ "FXFX"] = storage::VectorND< 3, double>( 0.657534, 1.81961, 0.2);
      s_InteractionMap[ "B4HX"] = storage::VectorND< 3, double>( 0.658938, 2.40884, 0.257878);
      s_InteractionMap[ "N3OX"] = storage::VectorND< 3, double>( 0.663797, 2.2259,  0.2);
      s_InteractionMap[ "B2C3"] = storage::VectorND< 3, double>( 0.666667, 3.12311, 0.2);
      s_InteractionMap[ "B3Li"] = storage::VectorND< 3, double>( 0.666667, 3.3675,  0.2);
      s_InteractionMap[ "BrPX"] = storage::VectorND< 3, double>( 0.666667, 3.67996, 0.2);
      s_InteractionMap[ "C2FX"] = storage::VectorND< 3, double>( 0.666667, 2.56582, 0.2);
      s_InteractionMap[ "C2IX"] = storage::VectorND< 3, double>( 0.666667, 3.87414, 0.2);
      s_InteractionMap[ "IXS3"] = storage::VectorND< 3, double>( 0.666667, 3.79508, 0.2);
      s_InteractionMap[ "LiN3"] = storage::VectorND< 3, double>( 0.666667, 3.14952, 0.2);
      s_InteractionMap[ "N3Na"] = storage::VectorND< 3, double>( 0.666667, 4.18374, 0.2);
      s_InteractionMap[ "C2C4"] = storage::VectorND< 3, double>( 0.671564, 2.6504,  0.376435);
      s_InteractionMap[ "FXN3"] = storage::VectorND< 3, double>( 0.67426,  2.56685, 0.2);
      s_InteractionMap[ "C3NX"] = storage::VectorND< 3, double>( 0.676417, 2.53908, 0.2);
      s_InteractionMap[ "OXSi"] = storage::VectorND< 3, double>( 0.677985, 2.87194, 0.322436);
      s_InteractionMap[ "HXNX"] = storage::VectorND< 3, double>( 0.678065, 2.06499, 0.2);
      s_InteractionMap[ "NXSe"] = storage::VectorND< 3, double>( 0.679012, 2.73252, 0.583678);
      s_InteractionMap[ "FXHX"] = storage::VectorND< 3, double>( 0.680812, 1.91595, 0.2);
      s_InteractionMap[ "OXSX"] = storage::VectorND< 3, double>( 0.684128, 2.45415, 0.2);
      s_InteractionMap[ "C3HX"] = storage::VectorND< 3, double>( 0.685536, 1.66329, 0.2);
      s_InteractionMap[ "FXSi"] = storage::VectorND< 3, double>( 0.685714, 2.97706, 0.258384);
      s_InteractionMap[ "FXSX"] = storage::VectorND< 3, double>( 0.685714, 2.74225, 0.261758);
      s_InteractionMap[ "B4NX"] = storage::VectorND< 3, double>( 0.696203, 2.82712, 0.248365);
      s_InteractionMap[ "C4NX"] = storage::VectorND< 3, double>( 0.700168, 2.4047,  0.252569);
      s_InteractionMap[ "AlHX"] = storage::VectorND< 3, double>( 0.705882, 2.77006, 0.272129);
      s_InteractionMap[ "ClSX"] = storage::VectorND< 3, double>( 0.705882, 3.03382, 0.2);
      s_InteractionMap[ "HXSi"] = storage::VectorND< 3, double>( 0.713971, 2.49443, 0.256861);
      s_InteractionMap[ "HXHX"] = storage::VectorND< 3, double>( 0.724104, 1.41373, 0.2);
      s_InteractionMap[ "B3B3"] = storage::VectorND< 3, double>( 0.727273, 2.39146, 0.2);
      s_InteractionMap[ "HXS6"] = storage::VectorND< 3, double>( 0.727273, 2.81357, 0.2);
      s_InteractionMap[ "IXNX"] = storage::VectorND< 3, double>( 0.727273, 3.27381, 0.2);
      s_InteractionMap[ "HXS3"] = storage::VectorND< 3, double>( 0.732484, 2.4967,  0.2);
      s_InteractionMap[ "HXPX"] = storage::VectorND< 3, double>( 0.741197, 2.35468, 0.212236);
      s_InteractionMap[ "N3NX"] = storage::VectorND< 3, double>( 0.747584, 2.50898, 0.2);
      s_InteractionMap[ "HXS4"] = storage::VectorND< 3, double>( 0.750133, 2.28853, 0.2);
      s_InteractionMap[ "C2C3"] = storage::VectorND< 3, double>( 0.75759,  2.7023,  0.267591);
      s_InteractionMap[ "N4NX"] = storage::VectorND< 3, double>( 0.758621, 2.56799, 0.227222);
      s_InteractionMap[ "C2OX"] = storage::VectorND< 3, double>( 0.759758, 2.52186, 0.329469);
      s_InteractionMap[ "NXTe"] = storage::VectorND< 3, double>( 0.764706, 2.61505, 0.2);
      s_InteractionMap[ "B4OX"] = storage::VectorND< 3, double>( 0.768421, 2.81293, 0.286012);
      s_InteractionMap[ "C2NX"] = storage::VectorND< 3, double>( 0.776488, 2.56197, 0.317596);
      s_InteractionMap[ "C2C2"] = storage::VectorND< 3, double>( 0.797163, 2.74503, 0.221332);
      s_InteractionMap[ "B4C2"] = storage::VectorND< 3, double>( 0.84,     3.19877, 0.321086);
      s_InteractionMap[ "NXSX"] = storage::VectorND< 3, double>( 0.85038,  2.66054, 0.241621);
      s_InteractionMap[ "BrC2"] = storage::VectorND< 3, double>( 0.875,    3.11339, 0.2);
      s_InteractionMap[ "AlAl"] = storage::VectorND< 3, double>( 0.9,      3.13022, 0.2);
      s_InteractionMap[ "AlC2"] = storage::VectorND< 3, double>( 0.9,      3.33178, 0.2);
      s_InteractionMap[ "AlOX"] = storage::VectorND< 3, double>( 0.9,      3.63278, 0.2);
      s_InteractionMap[ "B2HX"] = storage::VectorND< 3, double>( 0.9,      2.68355, 0.2);
      s_InteractionMap[ "B3C2"] = storage::VectorND< 3, double>( 0.9,      3.19877, 0.2);
      s_InteractionMap[ "B4Te"] = storage::VectorND< 3, double>( 0.9,      3.38697, 0.2);
      s_InteractionMap[ "BrSe"] = storage::VectorND< 3, double>( 0.9,      3.66196, 0.2);
      s_InteractionMap[ "C2Se"] = storage::VectorND< 3, double>( 0.9,      3.26339, 0.2);
      s_InteractionMap[ "ClSe"] = storage::VectorND< 3, double>( 0.9,      3.55229, 0.2);
      s_InteractionMap[ "LiLi"] = storage::VectorND< 3, double>( 0.9,      3.87002, 0.2);

      return true;
    }

  } // namespace chemistry
} // namespace bcl
