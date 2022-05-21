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

#ifndef BCL_DESCRIPTOR_MUTATION_AA_DENSITY_H_
#define BCL_DESCRIPTOR_MUTATION_AA_DENSITY_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "bcl_descriptor_iterator.h"
#include "bcl_descriptor_type.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_mutation.h"
#include "iterate/bcl_iterate_generic.h"
#include "math/bcl_math_trigonometric_transition.h"
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
    //! @class MutationAADensity
    //! @brief MutationAADensity gives the normalized position of an amino acid within the membrane the values are from
    //! (0-1.0) starting at inner membrane going to outer membrane. Based on Octopus TM predictions i = 0 o = 1.0
    //!
    //! @see @link example_descriptor_mutation_aa_density.cpp @endlink
    //! @author mendenjl
    //! @date Feb 01, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutationAADensity :
      public Base< biol::Mutation, float>
    {
    public:

    //////////
    // data //
    //////////

    private:

      //! Iterators sorted by position along the sequence
      storage::Map< storage::Pair< int, char>, Iterator< biol::AABase> > m_Mutations;

      //! voxel grid for fast lookup of AAs in the given distance
      assemble::VoxelGridAA m_VoxelGrid;

      //! bool of whether to include the weight of the output (number of mutants, if none are in the sigmoidal region)
      bool m_IncludeWeight;

      //! width of the cutoff
      float m_TransitionWidth;

      //! cutoff distance. Values beyond this will be subjected to the transition width
      float m_CutoffDistance;

      //! trigonometric transition
      math::TrigonometricTransition m_Transition;

      //! Actual descriptor to compute. Can be undefined
      util::Implementation< Base< biol::AABase, float> > m_Descriptor;

      //! instances of the class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutationAADensity();

      //! @brief virtual copy constructor
      MutationAADensity *Clone() const;

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

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief find mutants
      void FindMutants( const Iterator< biol::Mutation> &ITR);

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    }; // class MutationAADensity

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MUTATION_AA_DENSITY_H_
