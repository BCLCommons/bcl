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

#ifndef BCL_DESCRIPTOR_ATOM_VCHARGE_H_
#define BCL_DESCRIPTOR_ATOM_VCHARGE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomVcharge
    //! @brief calculates the Vcharge for every atom in a molecule
    //!
    //! @see http://pubs.acs.org/doi/full/10.1021/ci034148o
    //! @see https://structbio.vanderbilt.edu/twiki/bin/view/MeilerLab/UnderstandingHowV-ChargeIsCalculated
    //! @see https://structbio.vanderbilt.edu/twiki/bin/view/MeilerLab/InformationAboutTheNewVchargeImplementation
    //! @see @link example_descriptor_atom_vcharge.cpp @endlink
    //!
    //! @author mendenjl, kothiwsk
    //! @date Feb 07, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomVcharge :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {
    private:

      //! Cached charges, created whenever Calculate is called after a new molecule is set
      storage::Vector< float> m_Vcharge;

      //! extrapolation version - 0 = old version (prior to March, 2015), 1 = revised version (March 2015 onwards)
      size_t m_ExtrapolationVersion;

    public:

      //! the instance of this class is created in chemistry::AtomProperties, so no s_Instance is necessary

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, defaults to old extrapolation version
      AtomVcharge() :
        m_ExtrapolationVersion( 0)
      {
      }

      //! @brief constructor from extrapolation version
      explicit AtomVcharge( const size_t &VERSION) :
        m_ExtrapolationVersion( VERSION)
      {
      }

      //! @brief virtual copy constructor
      AtomVcharge *Clone() const;

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

      //! @brief Recalculates the charges for the current molecule
      void RecalculateCharges();

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

      //! @return a unique string for the given connectivity information
      //! @param SPECIAL_FEATURE One of the following letters: A: in aromatic ring, P: in planar ring, (space): otherwise
      //! @note: Planar ring feature is only applicable for N with 2-3 single bonds, no float  bonds
      //! @param EXTRAPOLATION_VERSION version of extrapolated parameters to use
      static std::string GetVchargeTypeString
      (
        const size_t &NUMBER_SINGLE_BONDS,
        const size_t &NUMBER_DOUBLE_BONDS,
        const size_t &NUMBER_TRIPLE_BONDS,
        const char &SPECIAL_FEATURE,
        const size_t &EXTRAPOLATION_VERSION
      );

      //! @brief create the map from element type & bond counts to Vcharge Electronegativity and Hardness
      //! @return the actual Map with the Vcharge Electronegativity and Hardness
      static storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> >
        MakeVchargeTypeStringToVchargeParametersMap();

    }; // class AtomVcharge

  } // namespace descriptor
} // namespace bcl

#endif //BCL_DESCRIPTOR_ATOM_VCHARGE_H_
