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

#ifndef BCL_DESCRIPTOR_ATOM_VOLUME_H_
#define BCL_DESCRIPTOR_ATOM_VOLUME_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomVolume
    //! @brief is a class which calculates the spherical volume for all atoms in a molecule
    //!
    //! @see @link example_descriptor_atom_volume.cpp @endlink
    //! @author kothiwsk, loweew, mendenjl
    //! @date Feb 18, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomVolume :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      //! Radius to be used to calculate the maximum SA, usually Atom_VDWaalsRadius or Atom_ConvalentRadius
      util::Implementation< Base< chemistry::AtomConformationalInterface, float> > m_MaxRadius;

      //! Radius used to calculate the minimum SA, if undefined, the minimum will be 0 for all atoms
      util::Implementation< Base< chemistry::AtomConformationalInterface, float> > m_MinRadius;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AtomVolume();

      //! @brief default constructor; takes parametric bool for whether to use covalent radius and CSD-derived VDW radii
      AtomVolume
      (
        const util::Implementation< Base< chemistry::AtomConformationalInterface, float> > &MAX_RADIUS,
        const util::Implementation< Base< chemistry::AtomConformationalInterface, float> > &MIN_RADIUS
      );

      //! @brief virtual copy constructor
      AtomVolume *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        // technically the volume depends weakly on the particular conformation chosen since bond lengths can change
        // however, in practice the differences are negligible
        return e_PreferCache;
      }

    }; // class AtomVolume

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ATOM_VOLUME_H_
