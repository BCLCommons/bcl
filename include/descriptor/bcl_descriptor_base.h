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

#ifndef BCL_DESCRIPTOR_BASE_H_
#define BCL_DESCRIPTOR_BASE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_cache_map.h"
#include "bcl_descriptor_cache_preference.h"
#include "bcl_descriptor_sequence_interface.h"
#include "bcl_descriptor_type.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Base
    //! @brief A generic base class for descriptors of arbitrary number of tuples of objects from a given sequence or
    //!        of the sequence class itself
    //!
    //! @tparam t_DataType type of items inside the sequences that this descriptor is for
    //! @tparam t_ReturnType elemental type that this descriptor returns
    //! Normally t_ReturnType is either float, for numeric descriptors, or char, for string descriptors
    //! @note this class makes use of the Non-Virtual-Interface idiom, see http://www.gotw.ca/publications/mill18.htm
    //! @note the term dimension in this class refers to the # of sequence-elements that go into calculating one feature
    //! @note for the descriptor.  Every descriptor has a native dimension (what the descriptor would like to be) and
    //! @note a set dimension (how we are actually treating it).  Dimension conversion is currently limited to these cases:
    //! @note 0 -> x (using a sequence descriptor in an element/pairwise/etc context)
    //! @note 1 -> x (using an elementwise descriptor in a pairwise or higher context)
    //! @note x -> 0 (summing an element/pairwise or triplet descriptor to make a sequence descriptor)
    //! @note All other conversions are forbidden
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Dec 11, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Base :
      public util::SerializableInterface
    {

    /////////////
    // friends //
    /////////////

      //! a class in the hpp that is used to handle caching if t_ReturnType is the same as the type held in the cache
      friend class DescriptorHelper;

    private:

    //////////
    // data //
    //////////

      //! The dimension that *this class will act like it is for via operator()
      //! For example, a sequence descriptor may be driven like an elementwise descriptor, in which case the type setting
      //! is elementwise (dimension = 1)
      size_t m_DimensionSetting;

      //! storage for the descriptor's last calculated values.  This vector is usually unused if the preference is cached
      linal::Vector< t_ReturnType> m_Storage;

      //! si-ptr to the object this base is currently attached to
      util::SiPtr< const SequenceInterface< t_DataType> > m_SequencePtr;

      //! cache descriptor reference
      linal::MatrixConstReference< t_ReturnType> m_CachedDescriptor;

      //! label in cache
      util::ObjectDataLabel m_CacheLabel;

      static bool s_HaveUpdatedHelp;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Base();

      //! @brief constructor from dimension
      Base( const size_t &DIMENSION);

      //! @brief copy constructor
      Base( const Base &BASE);

      //! @brief Clone function
      //! @return pointer to new Base
      virtual Base< t_DataType, t_ReturnType> *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the type of this descriptor
      //! @return the type of this descriptor (should ignore dimension setting)
      Type GetType() const;

      //! @brief get the type of this descriptor, considering dimension setting
      //! @return the type of this descriptor, considering dimension setting
      Type GetEffectiveType() const;

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      virtual Type::Symmetry GetSymmetry() const
      {
        return Type::e_Symmetric;
      }

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      virtual bool ConsiderRepeatedElements() const
      {
        return false;
      }

      //! @brief get the size of the features for this descriptor
      //! @return the size of features for this descriptor
      size_t GetSizeOfFeatures() const
      {
        // the storage already has the correct size, if it has been initialized
        return m_Storage.GetSize() ? m_Storage.GetSize() : GetNormalSizeOfFeatures();
      }

      //! @brief get the cache preference under the current dimension setting (e.g. m_DimensionSetting)
      //! @return the cache preference for the descriptor
      CachePreference GetCachePreference() const;

      //! @brief get the label of the current descriptor if it were placed in the cache
      //! @return the label of the current descriptor when it is placed in the cache
      const util::ObjectDataLabel &GetCacheLabel() const
      {
        return m_CacheLabel;
      }

      //! @brief override the dimension associated with this descriptor
      //! @param NEW_DIMENSION the new dimension to use
      void SetDimension( const size_t &NEW_DIMENSION);

      //! @brief set the sequence object and reset the storage
      //! @param SEQUENCE the sequence object of interest
      void SetObject( const SequenceInterface< t_DataType> &SEQUENCE);

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      virtual size_t GetNormalSizeOfFeatures() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest
      //! @return generated description for given argument
      linal::VectorConstReference< t_ReturnType> Recalculate( const Iterator< t_DataType> &ITR);

    ////////////////
    // operators //
    ////////////////

      //! @brief operator to calculate the descriptors for a given sequence iterator position, considering cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest
      //! @return generated description for given argument
      linal::VectorConstReference< t_ReturnType> operator()( const Iterator< t_DataType> &ITR);

      //! @brief operator to calculate the descriptors for a given sequence; all inner descriptors will be appended
      //! @param SEQ sequence of interest
      //! @return generated description for given argument
      //! returns GetNormalSizeOfFeatures() * SEQ.GetSize() values (GetNormalSizeOfFeatures for each element)
      linal::Vector< t_ReturnType> CollectValuesOnEachElementOfObject( const SequenceInterface< t_DataType> &SEQ);

      //! @brief operator to calculate the descriptors for a given sequence; all inner descriptors will be appended
      //! @param SEQ sequence of interest
      //! @return generated description for given argument
      //! returns GetNormalSizeOfFeatures() * SEQ.GetSize() values (GetNormalSizeOfFeatures for each element)
      linal::Vector< t_ReturnType> CollectValuesOnEachElementOfObject( const SequenceInterface< t_DataType> &SEQ) const;

      //! @brief operator to calculate the descriptors for a given sequence; all inner descriptors will be summed
      //! @param SEQ sequence of interest
      //! @return generated description for given argument
      //! returns GetNormalSizeOfFeatures(); if the descriptor is multi-dimensional; results are summed
      linal::Vector< t_ReturnType> SumOverObject( const SequenceInterface< t_DataType> &SEQ);

      //! @brief operator to calculate the descriptors for a given sequence; all inner descriptors will be summed
      //! @param SEQ sequence of interest
      //! @return generated description for given argument
      //! returns GetNormalSizeOfFeatures(); if the descriptor is multi-dimensional; results are summed
      linal::Vector< t_ReturnType> SumOverObject( const SequenceInterface< t_DataType> &SEQ) const;

      //! @brief helper function to write help for all instances of this class
      //! @param STREAM stream to write the help to
      //! @param INSTANCES all instances of the class
      static io::FixedLineWidthWriter &WriteInstancesHelp
      (
        io::FixedLineWidthWriter &STREAM,
        const storage::Map< std::string, util::OwnPtr< Base> > &INSTANCES,
        bool FULL_HELP
      );

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      virtual bool DimensionIsWellDefined() const
      {
        return !InjectDimensions();
      }

      //! @brief True if a derived class expects all calls to SetDimension to also call SetDimension on internally-held classes
      //! @note users do not normally need to override this function
      virtual bool InjectDimensions() const
      {
        return true;
      }

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @return the dimension (# of sequence elements) this descriptor is currently configured to accept
      size_t GetDimensionSetting() const
      {
        return m_DimensionSetting;
      }

      //! @brief get the object that descriptors are currently being calculated for
      //! @return simple pointer to the object that descriptors are currently being calculated for
      util::SiPtr< const SequenceInterface< t_DataType> > GetCurrentObject() const
      {
        return m_SequencePtr;
      }

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to forward calls from SetObject on to all internal implementations
      //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
      //! that function is called
      virtual iterate::Generic< Base< t_DataType, t_ReturnType> > GetInternalDescriptors();

      //! @brief get the logical name for the complete sequence; evaluates to e.g. sequence, molecule, or string
      static const std::string &GetObjectName();

      //! @brief get the logical name for each element of the sequence; evaluates to e.g. amino acid, atom, or character
      static const std::string &GetElementName();

    private:

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      virtual size_t GetNormalDimension() const = 0;

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      virtual void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< t_ReturnType> &STORAGE
      ) = 0;

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook()
      {
      }

      //! @brief hook that derived classes can override to add behavior after every time SetDimension is called
      virtual void SetDimensionHook()
      {
      }

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      virtual CachePreference GetNormalCachePreference() const;

      //! @brief True if the intended dimension can only be determined if this descriptor is used with other descriptors
      //! @note users do not normally need to override this function
      virtual bool DimensionIsContextual() const
      {
        return false;
      }

    }; // template class Base

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Base< chemistry::AtomConformationalInterface, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Base< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Base< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Base< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Base< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Base< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Base< char, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Base< char, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_BASE_H_
