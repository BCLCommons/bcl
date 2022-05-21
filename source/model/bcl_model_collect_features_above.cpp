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

// include forward header of this class
#include "model/bcl_model_collect_features_above.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> CollectFeaturesAbove::s_Instance
    (
      util::Enumerated< CollectFeaturesInterface>::AddInstance( new CollectFeaturesAbove( 0.0))
    );

    //! @brief constructor from parameter
    //! @param THRESHOLD threshold above which feature values will be kept
    CollectFeaturesAbove::CollectFeaturesAbove( const float &THRESHOLD) :
      m_Threshold( THRESHOLD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectFeaturesAbove
    CollectFeaturesAbove *CollectFeaturesAbove::Clone() const
    {
      return new CollectFeaturesAbove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CollectFeaturesAbove::GetClassIdentifier() const
    {
      return GetStaticClassName< CollectFeaturesAbove>();
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &CollectFeaturesAbove::GetAlias() const
    {
      static const std::string s_Name( "Above");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief collect the top m_NumberToKeep value indices by score
    //! @param SCORES any values
    //! @return top m_NumberToKeep value indices
    storage::Vector< size_t> CollectFeaturesAbove::Collect( const linal::Vector< float> &SCORES) const
    {
      // find the indices of values above the given threshold
      storage::Vector< size_t> to_keep;
      to_keep.AllocateMemory( SCORES.GetSize());
      for( size_t i( 0), i_max( SCORES.GetSize()); i < i_max; ++i)
      {
        if( SCORES( i) > m_Threshold)
        {
          to_keep.PushBack( i);
        }
      }
      return to_keep;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectFeaturesAbove::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Selects the indices of the values above a threshold");

      member_data.AddInitializer
      (
        "",
        "threshold value",
        io::Serialization::GetAgent( &m_Threshold),
        "0.0"
      );

      return member_data;
    }

  } // namespace model
} // namespace bcl
