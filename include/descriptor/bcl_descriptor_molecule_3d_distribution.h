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

#ifndef BCL_DESCRIPTOR_MOLECULE_3D_DISTRIBUTION_H_
#define BCL_DESCRIPTOR_MOLECULE_3D_DISTRIBUTION_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_base_element_or_sequence.h"
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Molecule3DDistribution
    //! @brief code object for 3DA Soft max c
    //! @details Calculates normal 3da, except that rather than sum, the max of each bin is found
    //!          After all 3da values have been computed, the maxes are smoothed over the remaining bins
    //!
    //! @see @link example_descriptor_molecule_3d_distribution.cpp @endlink
    //! @author mendenjl
    //! @date Aug 09, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule3DDistribution :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      CheminfoProperty     m_AtomProperty; //!< the atom property encode in the 3D autocorrelation function
      CheminfoProperty     m_CenterProperty; //!< the atom property encode in the 3D autocorrelation function
      size_t               m_NumberSteps;  //!< number of steps in 2D autocorrelation function
      float                m_StepSize;     //!< step size for 3d autocorrelation function

      //! Temporary vector used in calculations; stored here to avoid reallocating it during every call to Calculate
      linal::Vector< float> m_DiscreteCode;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Molecule3DDistribution();

      //! @brief constructor from number of steps, and mapped atom property
      Molecule3DDistribution
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const CheminfoProperty &CENTER_ATOM_PROPERTY,
        const size_t NUMBER_STEPS = 128,
        const float STEP_SIZE = 0.1
      );

      //! @brief virtual copy constructor
      Molecule3DDistribution *Clone() const;

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
        return m_NumberSteps;
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

      //! @brief Get a simple pointer to a smoothing coefficient vectors
      //!        Because the same parameters are often reused across many instances of this class, it is helpful
      //!        to have only a single smoothing coefficients vector
      //! @param NUMBER_STEPS number of steps needed
      //! @param TEMPERATURE gaussian-smoothing temperature
      //! @param STEP_SIZE size between each step
      static linal::VectorConstReference< float> GetSmoothingCoefficientVector
      (
        const size_t &NUMBER_STEPS,
        const float &TEMPERATURE,
        const float &STEP_SIZE
      );

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

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

      //! @brief create a gaussian-smoothed signal from m_DiscreteCode and store it in STORAGE
      //! @param STORAGE storage for the gaussian-smoothed signal
      void Smooth( linal::VectorReference< float> &STORAGE) const;

      //! @brief add an observed distance/property value to m_DiscreteCode
      //! @param DISTANCE actual distance of the two atoms
      //! @param PROP_A property from atom A
      void Accumulate
      (
        const float &DISTANCE,
        const float &PROP_A
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    }; // class Molecule3DDistribution

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_3D_DISTRIBUTION_H_
