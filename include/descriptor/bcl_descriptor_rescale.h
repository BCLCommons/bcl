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

#ifndef BCL_DESCRIPTOR_RESCALE_H_
#define BCL_DESCRIPTOR_RESCALE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_average.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Rescale
    //! @brief Rescales the value for the sequence; equivalent to writing ReplaceUndefinedValuesWithSequenceMean((X-SequenceMean(X)) / SequenceStd(X))
    //!
    //! @see @link example_descriptor_rescale.cpp @endlink
    //! @author mendenjl
    //! @date Oct 02, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Rescale :
      public Base< t_DataType, float>
    {
    private:

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, float> > m_Descriptor; //!< primary property of interest

      //! Property used internally to calculate actual rescaled value
      util::Implementation< Base< t_DataType, float> > m_Rescaled;

    public:

    //////////
    // data //
    //////////

      //! instances of the class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Rescale() :
        m_Descriptor(),
        m_Rescaled()
      {
      }

      //! @brief constructor from SOURCE descriptor
      explicit Rescale( const util::Implementation< Base< t_DataType, float> > &SOURCE) :
        m_Descriptor( SOURCE),
        m_Rescaled()
      {
        ReadInitializerSuccessHook( SOURCE.GetLabel(), util::GetLogger());
      }

      //! @brief copy constructor
      Rescale *Clone() const
      {
        return new Rescale( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const
      {
        static std::string s_alias( "Rescale");
        return s_alias;
      }

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const
      {
        return m_Descriptor->GetSymmetry();
      }

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const
      {
        return m_Descriptor->ConsiderRepeatedElements();
      }

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to forward calls from SetObject on to all internal implementations
      //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
      //! that function is called
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors()
      {
        // if auto-replacing with the mean of the non-nan values, then do not return the replacement property
        return
          iterate::Generic< Base< t_DataType, float> >( &m_Rescaled, &m_Rescaled + 1);
      }

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_Descriptor->GetSizeOfFeatures();
      }

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const
      {
        return m_Descriptor->GetType().GetDimension();
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer parameters;

        parameters.AddInitializer
        (
          "",
          "descriptor to rescale",
          io::Serialization::GetAgent( &m_Descriptor)
        );
        parameters.SetClassDescription
        (
          "Rescales values relative to the " + this->GetObjectName() + " values, specifically, "
          "computes: (X-" + this->GetObjectName() + "Mean(X))/" + this->GetObjectName() + "Std(X). "
          "In the event that std is 0, returns 0"
        );

        return parameters;
      }

    private:

      //! @brief True if a derived class expects all calls to SetDimension to also call SetDimension on internally-held classes
      //! @note users do not normally need to override this function
      bool InjectDimensions() const
      {
        return false;
      }

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      bool DimensionIsWellDefined() const
      {
        return false;
      }

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< float> &STORAGE
      )
      {
        STORAGE.CopyValues( m_Rescaled->operator()( ITR));
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
      {
        // create an object data label to initialize the rescaled implementation with
        const std::string internal_label_string( m_Descriptor.GetString());
        const util::ObjectDataLabel rescaled_label
        (
          "DefineNaN("
          "  replacement=Constant(0),"
          "  Divide("
          "    lhs=Subtract(lhs=" + internal_label_string
          +               ",rhs=" + Base< t_DataType, float>::GetObjectName() + "Mean(" + internal_label_string + ")),"
          "    rhs="  + Base< t_DataType, float>::GetObjectName() + "StandardDeviation(" + internal_label_string + ")"
          "  )"
          ")"
        );
        const bool success( m_Rescaled.TryRead( rescaled_label, ERR_STREAM));
        return success;
      }

    }; // class Rescale

    BCL_EXPIMP_TEMPLATE template class BCL_API Rescale< char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Rescale< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Rescale< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Rescale< chemistry::AtomConformationalInterface>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_RESCALE_H_
