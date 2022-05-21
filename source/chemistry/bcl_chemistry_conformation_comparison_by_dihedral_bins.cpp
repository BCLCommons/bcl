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
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_angle.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonByDihedralBins::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonByDihedralBins())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param BIN_SIZE bin size to be used for comparing dihedral bins
    ConformationComparisonByDihedralBins::ConformationComparisonByDihedralBins( const double BIN_SIZE) :
      m_BinSize( BIN_SIZE),
      m_HalfBinSize( BIN_SIZE / 2.0),
      m_MaxBin( 360.0 / BIN_SIZE)
    {
      PriorityDihedralAngles::SetWrappingAngle( m_HalfBinSize);
    }

    //! virtual copy constructor
    ConformationComparisonByDihedralBins *ConformationComparisonByDihedralBins::Clone() const
    {
      return new ConformationComparisonByDihedralBins( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonByDihedralBins::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonByDihedralBins::GetAlias() const
    {
      static const std::string s_Name( "DihedralBins");
      return s_Name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief compare two molecules on the basis of dihedral keys
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_B - second molecule being aligned
    //! @return 0 if molecules are the same conformer
    double ConformationComparisonByDihedralBins::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      // Determine if conformations can be compared
      if( !ConformationComparisonInterface::ConformationsAreComparable( MOLECULE_A, MOLECULE_B))
      {
        // nope, so return an undefined
        return util::GetUndefined< double>();
      }
      if( !m_Mutex.TryLock())
      {
        return ConformationComparisonByDihedralBins( m_BinSize)( MOLECULE_A, MOLECULE_B);
      }
      storage::Vector< double> molecule_a_angles( m_Priority( MOLECULE_A).First());
      storage::Vector< double> molecule_b_angles( m_Priority( MOLECULE_B).First());
      m_Mutex.Unlock();

      // get dihedral keys of two molecule
      linal::Vector< int> keys_molecule_a( DetermineDihedralKeys( molecule_a_angles, true, false));
      linal::Vector< int> keys_molecule_b( DetermineDihedralKeys( molecule_b_angles, true, false));

      linal::Vector< int> difference( keys_molecule_a);
      difference -= keys_molecule_b;
      math::Absolute( difference);

      for
      (
        linal::Vector< int>::iterator itr_diff( difference.Begin()), itr_diff_end( difference.End());
        itr_diff != itr_diff_end;
        ++itr_diff
      )
      {
        *itr_diff = std::min( ( *itr_diff), math::Absolute( int( m_MaxBin) - *itr_diff));
      }

      // if keys differ then we have a different conformers, so return 1
      return math::Statistics::Sum( difference.Begin(), difference.End());
    }

    //! @brief determine the dihedral keys of molecule
    //! @param ANGLE_VECTOR a container containing dihedral bond angle either in radian or degrees
    //! @param ANGLE_IN_DEGREES true if the angle is in degrees
    //! @param IS_WRAPPED_AROUND true if the angle is wrapped around 360 degrees.
    //! @return returns a storage::Vector< double> containing the keys
    storage::Vector< int> ConformationComparisonByDihedralBins::DetermineDihedralKeys
    (
      const linal::Vector< double> &ANGLE_VECTOR,
      bool ANGLE_IN_DEGREES,
      bool IS_WRAPPED_AROUND
    ) const
    {
      linal::Vector< double> dihedral_degrees( DetermineWrappedAroundAngles( ANGLE_VECTOR, !ANGLE_IN_DEGREES));

      storage::Vector< int> dihedral_keys( dihedral_degrees.GetSize(), util::GetUndefined< int>());

      // determine the bin to which dihedral bonds belong
      size_t angle_index( 0);
      for
      (
        linal::Vector< double>::const_iterator itr( dihedral_degrees.Begin()), itr_end( dihedral_degrees.End());
        itr != itr_end;
        ++itr, ++angle_index
      )
      {
        double angle_degree( fmod( *itr + 360.0, 360.0));
        double dihedral_angle_key( angle_degree / m_BinSize);
        dihedral_keys( angle_index) = rint( dihedral_angle_key);
        if( dihedral_keys( angle_index) == 0)
        {
          dihedral_keys( angle_index) = m_MaxBin;
        }
      }
      return dihedral_keys;
    }

    //! @brief determine the dihedral key
    //! @param ANGLE contains a dihedral bond angle either in radian or degrees
    //! @param ANGLE_IN_DEGREES true if the angle is in degrees
    //! @param IS_WRAPPED_AROUND true if the angle is wrapped around 360 degrees.
    //! @return returns a storage::Vector< double> containing the keys
    int ConformationComparisonByDihedralBins::DetermineDihedralKey
    (
      const double &ANGLE,
      bool ANGLE_IN_DEGREES,
      bool IS_WRAPPED_AROUND
    ) const
    {
      double dihedral_degrees( ANGLE_IN_DEGREES ? ANGLE : math::Angle::Degree( ANGLE));
      if( !IS_WRAPPED_AROUND)
      {
        dihedral_degrees = DetermineWrappedAroundAngle( ANGLE, true);
      }

      // determine the bin to which dihedral bonds belong
      double angle_degree( fmod( dihedral_degrees + 360.0, 360.0));
      double dihedral_angle_key( angle_degree / m_BinSize);
      int dihedral_key( rint( dihedral_angle_key));
      if( dihedral_key == 0)
      {
        dihedral_key = m_MaxBin;
      }
      return dihedral_key;
    }

    //! @brief determine the dihedral keys of molecule
    //! @param ANGLE_VECTOR a container of angle values in radians
    //! @param ANGLES_IN_RADIANS true if angles is in radians, otherwise false
    //! @return returns a storage::Vector< double> containing the keys
    const linal::Vector< double> ConformationComparisonByDihedralBins::DetermineWrappedAroundAngles
    (
      const linal::Vector< double> &ANGLE_VECTOR,
      bool ANGLES_IN_RADIANS
    ) const
    {
      // determine the bin to which dihedral bonds belong
      size_t angle_index( 0);
      linal::Vector< double> dihedral_degrees( ANGLE_VECTOR.GetSize());
      for
      (
        linal::Vector< double>::const_iterator itr( ANGLE_VECTOR.Begin()), itr_end( ANGLE_VECTOR.End());
        itr != itr_end;
        ++itr, ++angle_index
      )
      {
        double angle_degree( ANGLES_IN_RADIANS ? math::Angle::Degree( *itr) : *itr);

        if( angle_degree < m_HalfBinSize)
        {
          dihedral_degrees( angle_index) = 360 + ( angle_degree);
        }
        else
        {
          dihedral_degrees( angle_index) = angle_degree;
        }
      }
      return dihedral_degrees;
    }

    //! @brief determine the dihedral key of molecule
    //! @param ANGLE a angle values in radians
    //! @param ANGLES_IN_RADIANS true if angles is in radians, otherwise false
    //! @return returns a storage::Vector< double> containing the keys
    double ConformationComparisonByDihedralBins::DetermineWrappedAroundAngle
    (
      const double &ANGLE,
      bool ANGLES_IN_RADIANS
    ) const
    {
      // determine the bin to which dihedral bonds belong
      double dihedral_degrees( ANGLES_IN_RADIANS ? math::Angle::Degree( ANGLE) : ANGLE);
      if( dihedral_degrees < m_HalfBinSize)
      {
        dihedral_degrees += 360.0;
      }
      return dihedral_degrees;
    }

    //! @brief converts vector of dihedral keys to string
    //! @param DIHEDRAL_KEYS keys which have to be converted to string
    //! @return returns a string containing the keys
    std::string ConformationComparisonByDihedralBins::DihedralKeysString
    (
      const storage::Vector< int> &DIHEDRAL_KEYS
    ) const
    {
      std::stringstream dihedral_keys;
      for
      (
        storage::Vector< int>::const_iterator itr( DIHEDRAL_KEYS.Begin()), itr_end( DIHEDRAL_KEYS.End());
        itr != itr_end;
        ++itr
      )
      {
        dihedral_keys << *itr << ',';
      }
      return dihedral_keys.str();
    }

    //! @brief determine the upper bound and lower bound of bin from angles
    //! @param KEY_VALUE Key value that needs to be converted to range of angles
    //! @return a pair of values, where first corresponds to lower bound and second to upper bound
    const storage::Pair< double, double> ConformationComparisonByDihedralBins::KeyToAngle( size_t KEY_VALUE) const
    {
      BCL_Assert( KEY_VALUE * m_BinSize <= 360, "Key value should correspond to angle less than 360")
      double lower_bound = ( m_BinSize * KEY_VALUE) - m_HalfBinSize;
      if( lower_bound > 180)
      {
        lower_bound = lower_bound - 360;
      }
      double upper_bound = ( m_BinSize * KEY_VALUE) + m_HalfBinSize;
      if( upper_bound > 360)
      {
        upper_bound = upper_bound - 360;
      }
      else if( upper_bound > 180)
      {
        upper_bound = upper_bound - 360;
      }
      return storage::Pair< double, double>( lower_bound, upper_bound);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonByDihedralBins::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Compares dihedral bin strings of two molecules to determine whether they have essentially the same dihedral angles"
      );
      parameters.AddInitializer
      (
        "bin size",
        "size of bins, centered around 0 degrees",
        io::Serialization::GetAgentWithMin( &m_BinSize, 1.0),
        "30.0"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ConformationComparisonByDihedralBins::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_HalfBinSize = m_BinSize / 2.0;
      m_MaxBin = 360.0 / m_BinSize;
      PriorityDihedralAngles::SetWrappingAngle( m_HalfBinSize);
      return true;
    }

  } // namespace chemistry
} // namespace bcl
