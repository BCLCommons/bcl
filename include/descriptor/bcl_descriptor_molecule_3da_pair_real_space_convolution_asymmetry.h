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

#ifndef BCL_DESCRIPTOR_MOLECULE_3DA_PAIR_REAL_SPACE_CONVOLUTION_ASYMMETRY_H_
#define BCL_DESCRIPTOR_MOLECULE_3DA_PAIR_REAL_SPACE_CONVOLUTION_ASYMMETRY_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "bcl_descriptor_molecule_3da_pair_convolution_asymmetry.h"
#include "bcl_descriptor_molecule_3da_smooth_sign_code.h"
#include "bcl_descriptor_window.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Molecule3DAPairRealSpaceConvolutionAsymmetry
    //! @brief code object for 3DARealSpace
    //! @details Calculates the 3DASign of the contacting residues between a protein and a ligand
    //!  such that the 3DA distance bins are partitioned across contact distance.
    //!
    //! @see @link example_descriptor_molecule_3DA_pair_convolution.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Aug 28, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule3DAPairRealSpaceConvolutionAsymmetry :
      public Molecule3DAPairConvolutionAsymmetry
    {

    //////////
    // data //
    //////////

      //! Distance cutoff for detecting protein atom neighbors
      float m_Cutoff;

      //! Protein voxel grid
      chemistry::VoxelGridAtom m_VoxelGridPocket;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Molecule3DAPairRealSpaceConvolutionAsymmetry();

      //! @brief constructor from number of steps, and mapped atom property
      Molecule3DAPairRealSpaceConvolutionAsymmetry
      (
        const CheminfoProperty &ATOM_PROPERTY_A,
        const CheminfoProperty &ATOM_PROPERTY_B,
        const size_t NUMBER_STEPS = 20,
        const float STEP_SIZE = 0.50,
        const size_t WINDOW_SIZE = 4,
        const float CUTOFF = 7.0
      );

      //! @brief virtual copy constructor
      Molecule3DAPairRealSpaceConvolutionAsymmetry *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      //! @brief reduce pocket to atoms contacting the ligand
      virtual void ReduceProteinPocket( chemistry::FragmentComplete &POCKET);

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    }; // class Molecule3DAPairRealSpaceConvolutionAsymmetry

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_3DA_PAIR_REAL_SPACE_CONVOLUTION_ASYMMETRY_H_
