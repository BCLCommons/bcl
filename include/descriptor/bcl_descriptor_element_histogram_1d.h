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

#ifndef BCL_DESCRIPTOR_ELEMENT_HISTOGRAM_1D_H_
#define BCL_DESCRIPTOR_ELEMENT_HISTOGRAM_1D_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ElementHistogram1D
    //! @brief Computes histogram in two dimensions over two properties
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_element_histogram_1d.cpp @endlink
    //! @author mendenjl
    //! @date May 31, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class ElementHistogram1D :
      public BaseElement< t_DataType, float>
    {

    private:

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, float> > m_Descriptor; //!<  descriptor
      float m_Min;         //!< Min value
      float m_Max;         //!< Max value
      float m_BinSize;     //!< Bin size
      size_t m_NumberBins; //!< # of bins

      //! if non-zero, return a grid with points smoothed with a gaussian kernel, e.g.
      //! Ae^(-BinDistance/SmoothingDistance), where BinDistance is the euclidean distance to center of the given bin
      //! A is adjusted to ensure that the sum of the histogram's returned values == 1
      float m_SmoothingDistance;

      //! Option to catch all instances that fall outside the histogram's boundaries and placing them in the nearest bin
      bool m_Catchall;

    public:

    //////////
    // data //
    //////////

      //! instances of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ElementHistogram1D();

      //! @brief Clone function
      //! @return pointer to new ElementHistogram1D
      ElementHistogram1D< t_DataType> *Clone() const;

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    protected:

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

    }; // class ElementHistogram1D

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API ElementHistogram1D< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ElementHistogram1D< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ElementHistogram1D< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ElementHistogram1D< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ELEMENT_HISTOGRAM_1D_H_
