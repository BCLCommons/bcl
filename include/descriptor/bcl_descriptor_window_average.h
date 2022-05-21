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

#ifndef BCL_DESCRIPTOR_WINDOW_AVERAGE_H_
#define BCL_DESCRIPTOR_WINDOW_AVERAGE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_window.h"
#include "bcl_descriptor_window_weighting_interface.h"
#include "util/bcl_util_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WindowAverage
    //! @brief generates descriptions for a window around an element, weigthed by a desired weighting function
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_window_average.cpp @endlink
    //! @author mendenjl
    //! @date Mar 15, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class WindowAverage :
      public BaseElement< t_DataType, float>
    {

    private:

    //////////
    // data //
    //////////

      //! object used to create the window weights
      util::Implementation< WindowWeightingInterface> m_WindowWeightsCreator;

      //! actual weights for the window
      linal::Vector< float> m_WindowWeights;

      //! object used to retrieve the window
      Window< t_DataType, float> m_Window;

      //! True - window return a window containing the cumulative average out to N positions
      bool m_Cumulative;

      //! number of features returned ( = # of features per position in the window)
      size_t m_InternalDescriptorSize;

      //! bins to place given values in, used only for cumulative averages
      storage::Vector< size_t> m_Bins;

      //! bin averages, used only for cumulative averages
      storage::Vector< math::RunningAverage< float> > m_BinAverages;

    public:

    //////////
    // data //
    //////////

      //! instances of that class
      static const util::SiPtr< const util::ObjectInterface> s_ReflectiveInstance;
      static const util::SiPtr< const util::ObjectInterface> s_OpaqueInstance;
      static const util::SiPtr< const util::ObjectInterface> s_ReflectiveCumulativeInstance;
      static const util::SiPtr< const util::ObjectInterface> s_OpaqueCumulativeInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor from whether to be springy at the end and whether to do a cumulative average
      WindowAverage( const bool &REFLECTIVE, const bool &CUMULATIVE);

      //! @brief Clone function
      //! @return pointer to new WindowAverage
      WindowAverage< t_DataType> *Clone() const;

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

    private:

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

    }; // class WindowAverage

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API WindowAverage< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WindowAverage< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WINDOW_AVERAGE_H_
