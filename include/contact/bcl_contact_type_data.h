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

#ifndef BCL_CONTACT_TYPE_DATA_H_
#define BCL_CONTACT_TYPE_DATA_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TypeData
    //! @brief This is a low level helper class to store contact type properties
    //! @details This class is data class that stores contact type properties, and the class Types enumerates over
    //! over each contact type which is represented with an instance of this class.
    //!
    //! @see @link example_contact_type_data.cpp @endlink
    //! @author karakam
    //! @date 01.29.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TypeData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! This is the radii of windows for first and second windows associated with this contact type
      storage::Pair< size_t, size_t> m_WindowRadii;

      //! This is the length of windows for first and second windows associated with this contact type
      storage::Pair< size_t, size_t> m_WindowLengths;

      //! Whether this type is valid or not
      bool m_IsValid;

      //! Cut off distance for two residues has to be in in order to be in contact according to this type
      double m_ResidueDistanceCutoff;

      //! minimal distance between two contacting sses, so that they are not considered clashing
      double m_MinimalSSEDistance;

      //! distance range for this contact type used for the creating the scoring window from histograms
      math::Range< double> m_DistanceRange;

      //! preferred SSE distance range
      math::Range< double> m_PreferredDistanceRange;

      //! tilt angle range
      math::Range< double> m_TiltAngleRange;

      //! minimal fragment interface length
      double m_MinimalFragmentInterfaceLength;

      //! The types of SSEs which are involved in the Contact
      storage::Set< biol::SSType> m_SSTypes;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined contact type
      TypeData();

      //! @brief construct contact type from provided data
      //! @param WINDOW_RADII radius of window used for representation
      //! @param WINDOW_LENGTH_PAIR window length used for representation
      //! @param IS_VALID whether a valid contact type
      //! @param RESIDUE_DISTANCE_CUTOFF cut off distance for two residues to be considered in contact
      //! @param MINIMAL_SSE_DISTANCE minimal distance of two sses to not clash
      //! @param DISTANCE_RANGE distance range
      //! @param PREFERRED_DISTANCE_RANGE preferred distance range
      //! @param TILT_ANGLE_RANGE tilt angle range
      //! @param MINIMAL_FRAGMENT_INTERFACE_LENGTH minimal fragment interface length
      //! @param SS_TYPES SSTypes associated with this contact type
      TypeData
      (
        const storage::Pair< size_t, size_t> &WINDOW_RADII,
        const storage::Pair< size_t, size_t> &WINDOW_LENGTH_PAIR,
        const bool IS_VALID,
        const double RESIDUE_DISTANCE_CUTOFF,
        const double MINIMAL_SSE_DISTANCE,
        const math::Range< double> &DISTANCE_RANGE,
        const math::Range< double> &PREFERRED_DISTANCE_RANGE,
        const math::Range< double> &TILT_ANGLE_RANGE,
        const double MINIMAL_FRAGMENT_INTERFACE_LENGTH,
        const storage::Set< biol::SSType> &SS_TYPES
      );

      //! @brief virtual copy constructor
      TypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns window length pair
      //! @return window length pair
      const storage::Pair< size_t, size_t> &GetWindowRadii() const
      {
        return m_WindowRadii;
      }

      //! @brief returns window length pair
      //! @return window length pair
      const storage::Pair< size_t, size_t> &GetWindowLengths() const
      {
        return m_WindowLengths;
      }

      //! @brief returns residue distance cutoff
      //! @return residue distance cutoff
      double GetResidueDistanceCutoff() const
      {
        return m_ResidueDistanceCutoff;
      }

      //! @brief returns the minimal distance for the contacting SSEs - anything lower should be considered a clash
      double GetMinimalSSEDistance() const
      {
        return m_MinimalSSEDistance;
      }

      //! @brief returns distance range
      //! @return distance range
      const math::Range< double> &GetDistanceRange() const
      {
        return m_DistanceRange;
      }

      //! @brief returns preferred distance range
      //! @return preferred distance range
      const math::Range< double> &GetPreferredDistanceRange() const
      {
        return m_PreferredDistanceRange;
      }

      //! @brief returns tilt angle range
      //! @return tilt angle range
      const math::Range< double> &GetTiltAngleRange() const
      {
        return m_TiltAngleRange;
      }

      //! @brief returns minimal fragment interface length
      //! @return minimal fragment interface length
      double GetMinimalFragmentInterfaceLength() const
      {
        return m_MinimalFragmentInterfaceLength;
      }

      //! @brief returns vector of involved SSTypes
      //! @return vector of involved SSTypes
      const storage::Set< biol::SSType> &GetSSTypes() const
      {
        return m_SSTypes;
      }

      //! @brief returns whether this contact type is valid or not
      //! @return whether this contact type is valid or not
      bool IsValid() const
      {
        return m_IsValid;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class TypeData

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_TYPE_DATA_H_

