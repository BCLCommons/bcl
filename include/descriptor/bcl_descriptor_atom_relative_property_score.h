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

#ifndef BCL_DESCRIPTOR_ATOM_RELATIVE_PROPERTY_SCORE_H_
#define BCL_DESCRIPTOR_ATOM_RELATIVE_PROPERTY_SCORE_H_

// include the namespace header
#include "bcl_descriptor.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "storage/bcl_storage_map.h"

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
    //! @class AtomRelativePropertyScore
    //! @brief Compute the per-atom contribution to a property relative to the analagous atoms on a reference small molecule
    //!
    //! @see @link example_descriptor_atom_relative_property_score.cpp @endlink
    //! @author brownbp1
    //! @date Dec 24, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomRelativePropertyScore :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      // descriptor to score molecule and reference with
      CheminfoProperty m_Descriptor;

      // reference molecules for relative scoring
      std::string m_ReferenceMolsFilename;
      chemistry::FragmentEnsemble m_ReferenceMols;

      // compute dddG values and map to atoms; do this just once
      std::vector< std::map< size_t, float> > m_FeatureMap;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AtomRelativePropertyScore();

      //! @brief construct with property string
      AtomRelativePropertyScore( const CheminfoProperty &PROPERTY);

      //! @brief construct with reference filename
      AtomRelativePropertyScore( const std::string &REFERENCE_FILENAME);

      //! @brief virtual copy constructor
      AtomRelativePropertyScore *Clone() const;

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

      void SetReferenceMols( chemistry::FragmentEnsemble REFERNCE_MOLS);

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

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

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

#endif // BCL_DESCRIPTOR_ATOM_RELATIVE_PROPERTY_SCORE_H_
