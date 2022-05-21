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

#ifndef BCL_CHEMISTRY_DESCRIPTOR_DIMENSION_H_
#define BCL_CHEMISTRY_DESCRIPTOR_DIMENSION_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_chemistry_descriptor_dimension.h
  //! @brief provides names for different dimensions of descriptors
  //!
  //! @remarks example unnecessary
  //! @author mendenjl
  //! @date Feb 24, 2014
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace chemistry
  {

    enum DescriptorDimension
    {
      e_Molecule,                   //!< Whole molecule descriptor generation; one row per molecule
      e_Atom,                       //!< Per-atom descriptor generation; one row per atom
      e_UniquePairsOfDistinctAtoms, //!< All unique pairs of atoms that are distinct, e.g. A-B, A-C, B-C (not A-A)
      e_UniquePairsOfAtoms,         //!< All unique pairs of atoms, e.g. A-A, A-B, A-C, B-B, B-C, C-C
      e_PairsOfDistinctAtoms,       //!< All pairs of atoms that are distinct, e.g. A-B, A-C, B-A, B-C, C-A, C-B
      e_PairsOfAtoms,               //!< All pairs of atoms, A-A, A-B, A-C, B-A, B-B, B-C, C-A, C-B, C-C
      e_GivenByResultIdDescriptors, //!< Use the highest-dimension of the result/id descriptors to determine type
      e_GivenByDescriptors,         //!< Use the highest-dimension of the feature/result/id descriptors to determine type
      s_NumberDescriptorDimensions
    };

    //! @brief DescriptorDimension as string
    //! @param PREFERENCE the descriptor dimension desired
    //! @return the dimension as string
    BCL_API const std::string &GetDescriptorDimensionName( const DescriptorDimension &PREFERENCE);

    //! DescriptorDimensionEnum simplifies the usage of the DescriptorDimension enum of this class
    typedef util::WrapperEnum< DescriptorDimension, &GetDescriptorDimensionName, s_NumberDescriptorDimensions>
      DescriptorDimensionEnum;

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
    );

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_DESCRIPTOR_DIMENSION_H_
