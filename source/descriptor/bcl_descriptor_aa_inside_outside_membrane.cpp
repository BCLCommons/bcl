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
#include "descriptor/bcl_descriptor_aa_inside_outside_membrane.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AAInsideOutsideMembrane::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AAInsideOutsideMembrane()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAInsideOutsideMembrane::AAInsideOutsideMembrane() :
      m_Method( sspred::GetMethods().e_OCTOPUS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AAInsideOutsideMembrane *AAInsideOutsideMembrane::Clone() const
    {
      return new AAInsideOutsideMembrane( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAInsideOutsideMembrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAInsideOutsideMembrane::GetAlias() const
    {
      static const std::string s_name( "InsideOutsideMembrane");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAInsideOutsideMembrane::GetNormalSizeOfFeatures() const
    {
      return 3;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A the element of interest
    //! @param STORAGE storage for the descriptor
    void AAInsideOutsideMembrane::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // get the prediction map for this AA
      util::SiPtr< const sspred::MethodInterface> prediction_method_result( ELEMENT->GetSSPrediction( m_Method));

      // handle the case where the prediction is unavailable
      if( !prediction_method_result.IsDefined())
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }

      biol::EnvironmentType environment_type( prediction_method_result->GetOneStateTMPrediction());
      STORAGE = float( 0.0);
      if( environment_type == biol::GetEnvironmentTypes().e_MembraneCore)
      {
        STORAGE( 2) = 1.0;
      }
      else if( environment_type == biol::GetEnvironmentTypes().e_SolutionInside)
      {
        STORAGE( 0) = 1.0;
      }
      else if( environment_type == biol::GetEnvironmentTypes().e_SolutionOutside)
      {
        STORAGE( 1) = 1.0;
      }
      else if( environment_type == biol::GetEnvironmentTypes().e_Solution)
      {
        STORAGE( 0) = STORAGE( 1) = 0.5;
      }
      else
      {
        // undefined
        STORAGE = util::GetUndefined< float>();
      }
      return;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAInsideOutsideMembrane::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Returns whether a secondary structure prediction for an AA is inside, outside, or within the membrane. "
        "This method works only for methods that consider the inside and outside of the membrane separately, such as "
        "OCTOPUS and BOCTOPUS"
      );
      parameters.AddInitializer
      (
        "method",
        "secondary structure prediction or analysis method to use to obtain result",
        io::Serialization::GetAgent( &m_Method),
        "OCTOPUS"
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
