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

#ifndef BCL_DESCRIPTOR_AA_BLAST_WEIGHTED_PROPERTY_H_
#define BCL_DESCRIPTOR_AA_BLAST_WEIGHTED_PROPERTY_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "biol/bcl_biol_aa_base.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_reference.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AABlastWeightedProperty
    //! @brief Average of a given AA property, weighted by a particular AA's blast profile
    //!
    //! @see @link example_descriptor_aa_blast_weighted_property.cpp @endlink
    //! @author mendenjl
    //! @date Mar 07, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AABlastWeightedProperty :
      public BaseElement< biol::AABase, float>
    {

    public:

      // all blast-weighting methods
      enum Method
      {
        e_Probability,         //!< Weight by probabilities
        e_LogPOffset,          //!< Weight by log probability, with an offset to make the final weight == 1
        e_LogPSign,            //!< Average all values where the blast profile gave >= 0 (x+)
        e_LogPSignDiff,        //!< Average all values where the blast profile gave >= 0 + (x+ - x-)*conservation
        e_LogPSignTValueSign,  //!< signed t-value;  between all values >= 0 and those < 0 (population)
        e_PearsonCorrelation,  //!< Pearson correlation value (signed) between log blast profile and descriptor
        e_SpearmanCorrelation, //!< Spearman correlation value (signed) between log blast profile and descriptor
        s_NumberMethods
      };

      //! @brief get the string for the method
      //! @param METHOD the method to retrieve the name for
      static const std::string &GetMethodName( const Method &METHOD);

      //! Typedef for the method enum wrapper
      typedef util::WrapperEnum< Method, &GetMethodName, s_NumberMethods> MethodEnum;

    private:

    //////////
    // data //
    //////////

      //! set of PropertyType to be used
      biol::AATypeData::PropertyTypeEnum m_Property;

      //! Method to use
      MethodEnum m_Method;

      //! descriptor used to calculate AA blast profile; must return 20 values
      util::Implementation< Base< biol::AABase, float> > m_TypeCalculator;

      //! Simple pointer to vector with 31 floats, goes between blast profile values of [-15 15], gives relative weight
      //! to be applied for each member of the blast profile, for use with e_PearsonCorrelation only
      linal::Vector< float> m_BlastProfileWeighting;

      //! Offset from actual blast profile value to value in *m_BlastProfileWeighting
      static const int s_BlastProfileWeightingOffset = 15;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param METHOD weighting method to use
      AABlastWeightedProperty( const Method &METHOD = e_LogPOffset);

      //! @brief constructor from a property to use
      //! @param PROPERTY property to be used
      //! @param METHOD weighting method to use
      AABlastWeightedProperty( const biol::AATypeData::PropertyType &PROPERTY, const Method &METHOD = e_LogPOffset);

      //! @brief Clone function
      //! @return pointer to new AABlastProfile
      AABlastWeightedProperty *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to forward calls from SetObject on to all internal implementations
      //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
      //! that function is called
      iterate::Generic< Base< biol::AABase, float> > GetInternalDescriptors();

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

    }; // class AABlastWeightedProperty

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_BLAST_WEIGHTED_PROPERTY_H_
