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

#ifndef BCL_DESCRIPTOR_BAND_PASS_FILTER_H_
#define BCL_DESCRIPTOR_BAND_PASS_FILTER_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_replace_undefined_values.h"
#include "bcl_descriptor_window.h"
#include "bcl_descriptor_window_alignment_type.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BandPassFilter
    //! @brief computes the periodogram at specified frequency intervals, see http://en.wikipedia.org/wiki/BandPassFilter
    //! @details This class allows detection of frequency amplitude for noisy data regardless of phase
    //! @details The windowing function is currently hard-coded to being the wWlch window
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_band_pass_filter.cpp @endlink
    //! @author mendenjl
    //! @date Oct 16, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class BandPassFilter :
      public BaseElement< t_DataType, float>
    {

    private:

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, float> > m_Descriptor;    //!< Actual descriptor implementation
      size_t                                           m_Size;          //!< Window radius
      float                                            m_HighPassPeriod; //!< Minimal period desired; <= 1  for a high-pass filter
      float                                            m_LowPassPeriod;  //!< Maximal period desired; >= window size for a low-pass filter

      //! Alignment of the window relative to the current element
      WindowAlignmentEnum                              m_Alignment;

      //! matrix containing the wavelets used
      linal::Vector< float>                            m_Filter;

      //! bool; whether to use a reflecting window to minimize edge effects
      bool                                             m_Reflect;

      //! handles replacement of undefined values
      ReplaceUndefinedValues< t_DataType>              m_UndefinedReplacer;

      //! object used to create the window weights
      util::Implementation< WindowWeightingInterface> m_WindowWeightsCreator;

    public:

    //////////
    // data //
    //////////

      //! instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_ReflectingInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor from whether to be springy at the end
      BandPassFilter( const bool &REFLECTING);

      //! @brief Clone function
      //! @return pointer to new BandPassFilter
      BandPassFilter< t_DataType> *Clone() const;

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

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors();

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

    private:

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const t_DataType> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    }; // class BandPassFilter

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API BandPassFilter< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API BandPassFilter< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_BAND_PASS_FILTER_H_
