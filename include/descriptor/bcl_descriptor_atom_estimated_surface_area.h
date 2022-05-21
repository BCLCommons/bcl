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

#ifndef BCL_DESCRIPTOR_ATOM_ESTIMATED_SURFACE_AREA_H_
#define BCL_DESCRIPTOR_ATOM_ESTIMATED_SURFACE_AREA_H_

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
    //! @class AtomEstimatedSurfaceArea
    //! @brief calculates the effective polarizability for every atom in a molecule
    //!
    //! @see @link example_descriptor_atom_estimated_surface_area.cpp @endlink
    //! @author kothiwsk, mendenjl, geanesar
    //! @date Dec 14, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomEstimatedSurfaceArea :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

        //! bool; whether to compute using covalent radius (true) or vdw radius with
        bool m_UseCovalentRadius;

        //! bool: whether to compute using the CSD-derived van-der waals radii (more exact; but tends to emphasize H more)
        //! or the basic element type statistics
        bool m_UseCSDVdwRadius;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor; takes parametric bool for whether to use covalent radius and CSD-derived VDW radii
      AtomEstimatedSurfaceArea( const bool &COVALENT_RADIUS = true, const bool &USE_CSD_VDW = true);

      //! @brief virtual copy constructor
      AtomEstimatedSurfaceArea *Clone() const;

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
        return e_PreferCache;
      }

    }; // class AtomEstimatedSurfaceArea

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ATOM_ESTIMATED_SURFACE_AREA_H_
