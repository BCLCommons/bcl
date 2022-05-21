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

#ifndef BCL_DESCRIPTOR_MOLECULE_3DA_SOFT_MAX_H_
#define BCL_DESCRIPTOR_MOLECULE_3DA_SOFT_MAX_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Molecule3DASoftMax
    //! @brief code object for 3DA Soft max c
    //! @details Calculates normal 3da, except that rather than sum, the max of each bin is found
    //!          After all 3da values have been computed, the maxes are smoothed over the remaining bins
    //!
    //! @see @link example_descriptor_molecule_3da_soft_max.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Feb 12, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule3DASoftMax :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      CheminfoProperty     m_AtomProperty; //!< the atom property encode in the 3D autocorrelation function
      size_t               m_NumberSteps;  //!< number of steps in 2D autocorrelation function
      float                m_StepSize;     //!< step size for 3d autocorrelation function
      float                m_Temperature;  //!< Exponent of the gaussian function
      bool                 m_Smooth;       //!< Whether to perform smoothing; otherwise, linear kernel is used

      //! Temporary vector used in calculations; stored here to avoid reallocating it during every call to Calculate
      linal::Vector< float> m_DiscreteCode;

      //! Cached vector of smoothing coefficients; used to avoid excessive computations of the gaussian kernel
      linal::VectorConstReference< float> m_SmoothingCoefficients;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Molecule3DASoftMax();

      //! @brief constructor from number of steps, and mapped atom property
      Molecule3DASoftMax
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const size_t NUMBER_STEPS = 128,
        const float STEP_SIZE = 0.1,
        const float &TEMPERATURE = 100.0
      );

      //! @brief virtual copy constructor
      Molecule3DASoftMax *Clone() const;

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

      //! @brief get temperature of code
      //! @return const float  temperature of 3DA code
      const float &GetTemperature() const;

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

    }; // class Molecule3DASoftMax

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_3DA_SOFT_MAX_H_
