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
#include "chemistry/bcl_chemistry_descriptor_dimension.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_type.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief DescriptorDimension as string
    //! @param PREFERENCE the descriptor dimension desired
    //! @return the dimension as string
    const std::string &GetDescriptorDimensionName( const DescriptorDimension &PREFERENCE)
    {
      static const std::string s_names[ s_NumberDescriptorDimensions + 1] =
      {
        "Molecule",                   // Whole molecule descriptor generation; one row per molecule
        "Atom",                       // Per-atom descriptor generation; one row per atom
        "UniquePairsOfDistinctAtoms", // All unique pairs of atoms that are distinct, e.g. A-B, A-C, B-C (not A-A)
        "UniquePairsOfAtoms",         // All unique pairs of atoms, e.g. A-A, A-B, A-C, B-B, B-C, C-C
        "PairsOfDistinctAtoms",       // All pairs of atoms that are distinct, e.g. A-B, A-C, B-A, B-C, C-A, C-B
        "PairsOfAtoms",               // All pairs of atoms, A-A, A-B, A-C, B-A, B-B, B-C, C-A, C-B, C-C
        "GivenByResultIdDescriptors", // Use the highest-dimension of the result/id descriptors to determine type
        "GivenByDescriptors",         // Use the highest-dimension of the feature/result/id descriptors to determine type
      };
      return s_names[ PREFERENCE];
    }

    //! @brief determine the descriptor type given the dimension preference stated here, and the types of the
    //!        feature, result, and id descriptors
    //! @param PREFERENCE the descriptor dimension/type setting provided by the user
    //! @param FEATURE_TYPE the native dimension/type of the feature descriptor
    //! @param RESULT_TYPE the native dimension/type of the result descriptor
    //! @param ID_TYPE the native dimension/type of the id descriptor
    descriptor::Type DetermineDescriptorDimension
    (
      const DescriptorDimension &PREFERENCE,
      const descriptor::Type &FEATURE_TYPE,
      const descriptor::Type &RESULT_TYPE,
      const descriptor::Type &ID_TYPE
    )
    {
      switch( PREFERENCE)
      {
        case e_Molecule:                   return descriptor::Type( 0, false, descriptor::Type::e_Symmetric);  break;
        case e_Atom:                       return descriptor::Type( 1, false, descriptor::Type::e_Symmetric);  break;
        case e_UniquePairsOfDistinctAtoms: return descriptor::Type( 2, false, descriptor::Type::e_Symmetric);  break;
        case e_UniquePairsOfAtoms:         return descriptor::Type( 2, true,  descriptor::Type::e_Symmetric);  break;
        case e_PairsOfDistinctAtoms:       return descriptor::Type( 2, false, descriptor::Type::e_Asymmetric); break;
        case e_PairsOfAtoms:               return descriptor::Type( 2, true,  descriptor::Type::e_Asymmetric); break;
        case e_GivenByResultIdDescriptors: return descriptor::Type( RESULT_TYPE).GeneralizeToHandle( ID_TYPE); break;
        case s_NumberDescriptorDimensions:
        case e_GivenByDescriptors:
        default:
          return descriptor::Type( RESULT_TYPE).GeneralizeToHandle( ID_TYPE).GeneralizeToHandle( FEATURE_TYPE); break;
      }
      return descriptor::Type();
    }

  } // namespace chemistry
} // namespace bcl
