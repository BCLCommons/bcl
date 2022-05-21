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

#ifndef BCL_DESCRIPTOR_ATOM_TOPOLOGICAL_POLAR_SURFACE_AREA_H_
#define BCL_DESCRIPTOR_ATOM_TOPOLOGICAL_POLAR_SURFACE_AREA_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "graph/bcl_graph_const_graph.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomTopologicalPolarSurfaceArea
    //! @brief Calculates the spherical surface area for all atoms, convoluted with any other atom property desired
    //!
    //! @see @link example_descriptor_atom_topological_polar_surface_area.cpp @endlink
    //! @author mendenjl, kothiwsk
    //! @date Feb 06, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomTopologicalPolarSurfaceArea :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      graph::ConstGraph< size_t, size_t> m_Graph; //!< Graph of the current molecule

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      AtomTopologicalPolarSurfaceArea *Clone() const;

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

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @return a unique string for the given connectivity information
      static std::string GetTPSATypeString
      (
        const std::string &SYMBOL,
        const size_t &NUMBER_HYDROGEN_BONDS,
        const size_t &NUMBER_SINGLE_BONDS_TO_HEAVY_ATOMS,
        const size_t &NUMBER_DOUBLE_BONDS,
        const size_t &NUMBER_TRIPLE_BONDS,
        const size_t &NUMBER_AROMATIC_BONDS,
        const bool   &IS_IN_3_MEMBERED_RING,
        const short  &CHARGE
      );

      //! @brief create the connectivity string to topological polar surface area
      static storage::Map< std::string, float> MakeTPSATypeStringToPSAMap();

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class AtomTopologicalPolarSurface

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ATOM_TOPOLOGICAL_POLAR_SURFACE_AREA_H_
