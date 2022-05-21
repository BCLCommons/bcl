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

#ifndef BCL_DESCRIPTOR_WITHIN_RANGE_H_
#define BCL_DESCRIPTOR_WITHIN_RANGE_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "math/bcl_math_trigonometric_transition.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WithinRange
    //! @brief determines whether a descriptor's mean value lies between a certain range
    //!
    //! @see @link example_descriptor_within_range.cpp @endlink
    //! @author geanesar
    //! @date Oct 23, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class WithinRange :
      public Base< t_DataType, float>
    {
    private:

    //////////
    // data //
    //////////

      //! the property to calculate a result for
      util::Implementation< Base< t_DataType, float> > m_Descriptor;

      //! range beginning
      float m_RangeBegin;

      //! range end
      float m_RangeEnd;

      //! whether to include end points
      bool m_InclusiveRange;

    public:

      //! single ascending instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      WithinRange();

      //! @brief constructor
      //! @param RANGE_BEGIN where return values of 1 should start
      //! @param RANGE_END where return values of 1 should end
      //! @param INCLUSIVE whether the range is inclusive on the ends
      WithinRange
      (
        const float &RANGE_BEGIN,
        const float &RANGE_END,
        const bool &INCLUSIVE = true
      );

      //! @brief virtual copy constructor
      WithinRange *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief set descriptor to be queried
      //! @param DESCRIPTOR the descriptor to use
      void SetDescriptor( const util::Implementation< Base< t_DataType, float> > &DESCRIPTOR);

      //! @brief get the descriptor that is queried
      //! @return an implementation of the descriptor
      const util::Implementation< Base< t_DataType, float> > &GetDescriptor() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void RecalculateImpl( const Iterator< t_DataType> &ITR, linal::VectorReference< float> &STORAGE);

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const;

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const;

    }; // class WithinRange

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API WithinRange< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WithinRange< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WithinRange< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WithinRange< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WITHIN_RANGE_H_
