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

#ifndef BCL_DESCRIPTOR_MIN_MAX_INDEX_H_
#define BCL_DESCRIPTOR_MIN_MAX_INDEX_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MinMaxIndex
    //! @details Computes the min or max index of a multi-valued descriptor
    //!
    //! @see @link example_descriptor_min_max_index.cpp @endlink
    //! @author mendenjl
    //! @date Mar 06, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MinMaxIndex :
      public Base< t_DataType, float>
    {
    private:

    //////////
    // data //
    //////////

      //! the internal descriptor to find the min/max index of
      util::Implementation< Base< t_DataType, float> > m_Property;

      bool m_ComputeMax; //!< whether to compute max index (true) or min index (false)

    public:

      //! single ascending instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_InstanceMax;

      //! single descending instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_InstanceMin;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, accepts bool of whether to compute max
      MinMaxIndex( const bool &COMPUTE_MAX);

      //! @brief virtual copy constructor
      MinMaxIndex *Clone() const;

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

      //! @brief set properties to be logarithmized
      //! @param PROPERTY SmallMoleculeProperty
      void SetProperty( const util::Implementation< Base< t_DataType, float> > &PROPERTY);

      //! @brief get atom property of code
      //! @return atom property mapped in 2da code
      const util::Implementation< Base< t_DataType, float> > &GetProperty() const;

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

    }; // class MinMaxIndex

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API MinMaxIndex< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MinMaxIndex< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MinMaxIndex< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MinMaxIndex< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MIN_MAX_INDEX_H_
