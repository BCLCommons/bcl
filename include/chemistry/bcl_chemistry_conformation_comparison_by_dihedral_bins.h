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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_DIHEDRAL_BINS_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_DIHEDRAL_BINS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_priority_dihedral_angles.h"
#include "sched/bcl_sched_mutex.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonByDihedralBins
    //! @brief This class is designed to be used for determining and comparing 3D structures for molecules.
    //!
    //! @see @link example_chemistry_conformation_comparison_by_dihedral_bins.cpp @endlink
    //! @author kothiwsk
    //! @date Nov 21, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonByDihedralBins :
      public ConformationComparisonInterface
    {

    private:

    //////////
    // data //
    //////////

      //!< bin size to be used for comparison of dihedral angles
      double                                      m_BinSize;

      //! half the bin sized, cached for performance
      double                                      m_HalfBinSize;

      //! max bin, cached for performance
      size_t                                      m_MaxBin;

      //!< object to calculate priority dihedral angles
      PriorityDihedralAngles                      m_Priority;

      //! Mutex for changing m_Priority
      mutable sched::Mutex                        m_Mutex;

      //!< single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! @param BIN_SIZE bin size to be used for comparing dihedral bins
      ConformationComparisonByDihedralBins( const double BIN_SIZE = double( 30.0));

      //! virtual copy constructor
      ConformationComparisonByDihedralBins *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief returns the bin size used
      //! @return the bin size used
      const double &GetBinSize() const
      {
        return m_BinSize;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief compare two molecules on the basis of dihedral keys
      //! @param MOLECULE_A - first molecule being aligned
      //! @param MOLECULE_B - second molecule being aligned
      //! @param BINSIZE - binning size for dihedral angle binning
      //! @return 0 if molecules are same conformers
      double operator()
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B
      ) const;

      //! @brief determine the dihedral keys of molecule
      //! @param ANGLE_VECTOR a container containing dihedral bond angle either in radian or degrees
      //! @param ANGLE_IN_DEGREES true if the angle is in degrees
      //! @param IS_WRAPPED_AROUND true if the angle is wrapped around 360 degrees.
      //! @return returns a storage::Vector< double> containing the keys
      storage::Vector< int> DetermineDihedralKeys
      (
        const linal::Vector< double> &ANGLE_VECTOR,
        bool ANGLE_IN_DEGREES = true,
        bool IS_WRAPPED_AROUND = true
      ) const;

      //! @brief determine the dihedral key
      //! @param ANGLE contains a dihedral bond angle either in radian or degrees
      //! @param ANGLE_IN_DEGREES true if the angle is in degrees
      //! @param IS_WRAPPED_AROUND true if the angle is wrapped around 360 degrees.
      //! @return returns a storage::Vector< double> containing the keys
      int DetermineDihedralKey
      (
        const double &ANGLE,
        bool ANGLE_IN_DEGREES = true,
        bool IS_WRAPPED_AROUND = true
      ) const;

      //! @brief determine the dihedral keys of molecule
      //! @param ANGLE_VECTOR a container of angle values in radians
      //! @param ANGLES_IN_RADIANS true if angles is in radians, otherwise false
      //! @return returns a storage::Vector< double> containing the keys
      const linal::Vector< double> DetermineWrappedAroundAngles
      (
        const linal::Vector< double> &ANGLE_VECTOR,
        bool ANGLES_IN_RADIANS = true
      ) const;

      //! @brief determine the dihedral key of molecule
      //! @param ANGLE a angle values in radians
      //! @param ANGLES_IN_RADIANS true if angles is in radians, otherwise false
      //! @return returns a storage::Vector< double> containing the keys
      double DetermineWrappedAroundAngle( const double &ANGLE, bool ANGLES_IN_RADIANS) const;

      //! @brief converts vector of dihedral keys to string
      //! @param DIHEDRAL_KEYS keys which have to be converted to string
      //! @return returns a string containing the keys
      std::string DihedralKeysString
      (
        const storage::Vector< int> &DIHEDRAL_KEYS
      ) const;

      //! @brief determine the upper bound and lower bound of bin from angles
      //! @param KEY_VALUE Key value that needs to be converted to range of angles
      //! @return a pair of values, where first corresponds to lower bound and second to upper bound
      const storage::Pair< double, double> KeyToAngle
      (
        size_t KEY_VALUE
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);
    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_DIHEDRAL_BINS_H_
