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

#ifndef BCL_CHEMISTRY_DESCRIPTOR_TO_SCORE_ADAPTOR_H_
#define BCL_CHEMISTRY_DESCRIPTOR_TO_SCORE_ADAPTOR_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DescriptorToScoreAdaptor
    //! @brief Adapts descriptor::Base< AtomConformationalInterface> to math::FunctionInterfaceSerializable< FragmentComplete, double>
    //!        for use in scoring
    //!
    //! @see @link example_chemistry_descriptor_to_score_adaptor.cpp @endlink
    //! @author mendenjl
    //! @date Mar 18, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorToScoreAdaptor :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>
    {
    private:

      //! descriptor to compute
      mutable descriptor::CheminfoProperty m_Descriptor;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      DescriptorToScoreAdaptor()
      {
      }

      //! @brief constructor
      //! @param PROPERTY property to compute
      explicit DescriptorToScoreAdaptor( const descriptor::CheminfoProperty &PROPERTY);

      //! @brief Clone function
      //! @return pointer to new DescriptorToScoreAdaptor
      DescriptorToScoreAdaptor *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetAlias() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief test whether or not the implementation is defined
      //! @return true if the implementation is defined
      bool IsDefined() const;

      //! @brief get the current property
      //! @return the current property
      const descriptor::CheminfoProperty &GetProperty() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate clashes for given atom pair
      //! @param MOLECULE molecule that needs to scored
      //! @return propensity score for observing rotamers that exist in conformation
      double operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class DescriptorToScoreAdaptor

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_DESCRIPTOR_TO_SCORE_ADAPTOR_H_
