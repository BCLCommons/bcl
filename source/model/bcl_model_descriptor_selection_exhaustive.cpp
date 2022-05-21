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
#include "model/bcl_model_descriptor_selection_exhaustive.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math.h"
#include "util/bcl_util_message.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DescriptorSelectionExhaustive::s_Instance
    (
      util::Enumerated< DescriptorSelectionInterface>::AddInstance
      (
        new DescriptorSelectionExhaustive()
      )
    );

    //! @brief Clone function
    //! @return pointer to new DescriptorSelectionExhaustive
    DescriptorSelectionExhaustive *DescriptorSelectionExhaustive::Clone() const
    {
      return new DescriptorSelectionExhaustive();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &DescriptorSelectionExhaustive::GetAlias() const
    {
      static const std::string s_Name( "ExhaustiveSelection");
      return s_Name;
    }

    namespace
    {
      //! @brief update chosen with the next combination n-choose-k corresponds to N-Choose-chosen.getsize()
      bool NextCombination( storage::Vector< size_t> &CHOSEN, const size_t &N)
      {
        const size_t k( CHOSEN.GetSize());
        size_t mx_value( N - 1);
        for( size_t i( 0); i < k; ++i, --mx_value)
        {
          if( CHOSEN( i) < mx_value)
          {
            ++CHOSEN( i);
            for( size_t j( 0); j < i; ++j)
            {
              CHOSEN( j) = CHOSEN( i) + i - j;
            }
            return true;
            break;
          }
        }
        return false;
      }
    }

    //! @brief assemble all descriptor combinations based on an initial and a total descriptor set as an object label
    //! @param INITIAL initial descriptor set as an object label
    //! @param TOTAL all available descriptor groups in descriptor selection process
    //! @return container with all possible descriptor combinations based on initial descriptor set
    const storage::Vector< util::ObjectDataLabel> DescriptorSelectionExhaustive::operator()
    (
      const util::ObjectDataLabel &INITIAL,
      const util::ObjectDataLabel &TOTAL
    ) const
    {
      storage::Vector< util::ObjectDataLabel> labels;
      const size_t n( TOTAL.GetNumberArguments());
      for( size_t k( 1), max_k( std::min( m_MaxFeatures, TOTAL.GetNumberArguments() - 1)); k <= max_k; ++k)
      {
        storage::Vector< size_t> indices( k);
        for( size_t i( 0); i < k; ++i)
        {
          indices( i) = k - i - 1;
        }

        labels.AllocateMemory( labels.GetSize() + math::BinomialCoefficient( n, k));
        do
        {
          storage::Vector< util::ObjectDataLabel> sublabels;
          sublabels.AllocateMemory( k);
          for( size_t i( 0); i < k; ++i)
          {
            sublabels.PushBack( TOTAL.GetArgument( indices( i)));
          }
          labels.PushBack( util::ObjectDataLabel( "Combine", sublabels));
        } while( NextCombination( indices, n));
      }

      // return list of descriptor combinations
      return labels;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DescriptorSelectionExhaustive::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetCommandLineIdentifier( "Exhaustive Forward feature selection");
      serial.AddInitializer( "max features", "max # of features to use", io::Serialization::GetAgent( &m_MaxFeatures), "2");
      return serial;
    }

  } // namespace model
} // namespace bcl
