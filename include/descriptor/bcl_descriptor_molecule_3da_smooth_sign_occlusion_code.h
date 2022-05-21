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

#ifndef BCL_DESCRIPTOR_MOLECULE_3DA_SMOOTH_SIGN_OCCLUSION_CODE_H_
#define BCL_DESCRIPTOR_MOLECULE_3DA_SMOOTH_SIGN_OCCLUSION_CODE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element_or_sequence.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "coord/bcl_coord_line_segment_3d.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Molecule3DASmoothSignOcclusionCode
    //! @brief code object for 3DA/RDF Smooth Sign Code
    //! @details This class provides methods for calculating 3DA related input code
    //!          This version of the 3DA returns 3 values for every bin:
    //!          Bin 0: Sum of RDF kernel values for when both atom's property is < 0,
    //!          Bin 1: Sum of RDF kernel values for when both atom's property is > 0,
    //!          Bin 2: Sum of RDF kernel values for when the atom properties have opposite sign
    //!
    //! @see @link example_descriptor_molecule_3da_smooth_sign_occlusion_code.cpp @endlink
    //! @author mendenjl
    //! @date Apr 22, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule3DASmoothSignOcclusionCode :
      public BaseElementOrSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      CheminfoProperty     m_AtomProperty; //!< the atom property encode in the 3D autocorrelation function
      size_t               m_NumberSteps;  //!< number of steps in 2D autocorrelation function
      float                m_StepSize;     //!< step size for 3d autocorrelation function
      size_t               m_MaxOcclusions;
      float                m_ClashDistance;
      bool                 m_Sqrt;

      //! Vector of bonds; one line per bond
      storage::Vector< coord::LineSegment3D> m_BondLines;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Molecule3DASmoothSignOcclusionCode();

      //! @brief constructor from property
      Molecule3DASmoothSignOcclusionCode
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const size_t NUMBER_STEPS,
        const float STEP_SIZE,
        const bool &SQRT = false
      );

      //! @brief virtual copy constructor
      Molecule3DASmoothSignOcclusionCode *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_NumberSteps * 3 * ( m_MaxOcclusions + 1);
      }

      //! @brief get number of steps of code
      //! @return number of steps in 2da code
      size_t GetNumberSteps() const;

      //! @brief get step size of code
      //! @return step size of 3DA code
      float GetStepSize() const;

      //! @brief get atom property of code
      //! @return atom property mapped in 2da code
      const CheminfoProperty &GetAtomProperty() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param ELEMENT_A, ELEMENT_B the elements of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_A,
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_B,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

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

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      //! In this case, updates the line segments
      void SetObjectHook();

    protected:

      //! @brief add an observed distance/property value to m_DiscreteCode
      //! @param where to accumulate the data at
      //! @param DISTANCE actual distance of the two atoms
      //! @param PROP_A property from atom A
      //! @param PROP_B property from atom B
      //! @param N_CLASHES number of clashes
      void Accumulate
      (
        linal::VectorReference< float> &STORAGE,
        const float &DISTANCE,
        const float &PROP_A,
        const float &PROP_B,
        const size_t &N_CLASHES
      );

      //! @brief determine the # of clashes between two atom pairs
      size_t CountClashes
      (
        const chemistry::AtomConformationalInterface &A,
        const chemistry::AtomConformationalInterface &B
      ) const;

    }; // class Molecule3DASmoothSignOcclusionCode

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_3DA_SMOOTH_SIGN_OCCLUSION_CODE_H_
