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

#ifndef BCL_DESCRIPTOR_MUTATION_NEAREST_IN_SEQUENCE_H_
#define BCL_DESCRIPTOR_MUTATION_NEAREST_IN_SEQUENCE_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "bcl_descriptor_iterator.h"
#include "bcl_descriptor_type.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_mutation.h"
#include "iterate/bcl_iterate_generic.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutationNearestInSequence
    //! @brief MutationNearestInSequence gives the normalized position of an amino acid within the membrane the values are from
    //! (0-1.0) starting at inner membrane going to outer membrane. Based on Octopus TM predictions i = 0 o = 1.0
    //!
    //! @see @link example_descriptor_mutation_nearest_in_sequence.cpp @endlink
    //! @author mendenjl
    //! @date Jan 30, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutationNearestInSequence :
      public Base< biol::Mutation, float>
    {
    public:

    //////////
    // data //
    //////////

    private:

      //! Iterators sorted by position along the sequence
      storage::Map< biol::Mutation, Iterator< biol::Mutation> > m_Mutations;

      //! bool of whether to include the distance to each mutant
      bool m_IncludeDistances;

      //! bool of whether to include the direction of each mutant (-1 - N-terminal to this mutation, 0 at same residue, 1 - C-terminal)
      bool m_IncludeDirection;

      //! Maximum # of different Mutations to find, e.g. if set to 1, only find the first Mutations with that type, for 2 find the
      //! first two of that type, etc.
      size_t m_MaxToFind;

      //! whether to ignore/skip mutants for which the descriptor is undefined
      bool m_SkipUndefined;

      //! Actual descriptor to compute. Can be undefined
      util::Implementation< Base< biol::Mutation, float> > m_Descriptor;

      //! instances of the class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutationNearestInSequence();

      //! @brief virtual copy constructor
      MutationNearestInSequence *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< biol::Mutation, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void RecalculateImpl( const Iterator< biol::Mutation> &ITR, linal::VectorReference< float> &STORAGE);

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const;

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const;

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

      //! @brief find mutants
      void FindMutants( const Iterator< biol::Mutation> &ITR);

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    }; // class MutationNearestInSequence

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MUTATION_NEAREST_IN_SEQUENCE_H_
