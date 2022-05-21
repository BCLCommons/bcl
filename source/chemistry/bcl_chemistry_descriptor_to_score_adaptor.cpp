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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_descriptor_to_score_adaptor.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DescriptorToScoreAdaptor::s_Instance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< FragmentComplete, double> >::AddInstance( new DescriptorToScoreAdaptor())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param PROPERTY property to compute
    DescriptorToScoreAdaptor::DescriptorToScoreAdaptor( const descriptor::CheminfoProperty &PROPERTY) :
      m_Descriptor( PROPERTY)
    {
      m_Descriptor->SetDimension( 0);
    }

    //! @brief Clone function
    //! @return pointer to new DescriptorToScoreAdaptor
    DescriptorToScoreAdaptor *DescriptorToScoreAdaptor::Clone() const
    {
      return new DescriptorToScoreAdaptor( *this);
    }

    /////////////////fragment_probabilities
    // data access //
    /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DescriptorToScoreAdaptor::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DescriptorToScoreAdaptor::GetAlias() const
    {
      static const std::string s_name( "DescriptorAdaptor");
      return s_name;
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DescriptorToScoreAdaptor::GetScheme() const
    {
      return m_Descriptor.IsDefined() ? m_Descriptor->GetAlias() : GetAlias();
    }

    //! @brief test whether or not the implementation is defined
    //! @return true if the implementation is defined
    bool DescriptorToScoreAdaptor::IsDefined() const
    {
      return m_Descriptor.IsDefined();
    }

    //! @brief get the current property
    //! @return the current property
    const descriptor::CheminfoProperty &DescriptorToScoreAdaptor::GetProperty() const
    {
      return m_Descriptor;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param MOLECULE molecule that needs to scored
    //! @return propensity score for observing rotamers that exist in conformation
    double DescriptorToScoreAdaptor::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      m_Descriptor->SetObject( MOLECULE);
      descriptor::Iterator< AtomConformationalInterface> itr( descriptor::Type(), MOLECULE);
      return m_Descriptor->operator ()( itr)( 0);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool DescriptorToScoreAdaptor::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      if( m_Descriptor.IsDefined())
      {
        m_Descriptor->SetDimension( 0);
        if( m_Descriptor->GetSizeOfFeatures() != size_t( 1))
        {
          ERR_STREAM << GetAlias() << " expects a single-valued descriptor, but a descriptor returning "
                     << m_Descriptor->GetSizeOfFeatures() << " values was given";
          return false;
        }
      }
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DescriptorToScoreAdaptor::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Allows use of descriptors in the scoring framework");
      serializer.AddInitializer
      (
        "",
        "descriptor to use; must return only a single value",
        io::Serialization::GetAgent( &m_Descriptor)
      );
      return serializer;
    }

  } // namespace chemistry
} // namespace bcl
